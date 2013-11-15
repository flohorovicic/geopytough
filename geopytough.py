"""Create a mulgrid mesh, apply to a TOUGH2 template file and update rocktypes from GeoModeller """

from matplotlib import use
use("Agg")

from t2grids import *
from t2data import *
from mulgrids import *
import copy
from os import path
import numpy
import os
import shutil
import subprocess
import time

class GeoT2Data(t2data):
    """Base class for extension of PyTOUGH t2data module with Geomodeller links
    
    overwrites functions (and other relevant changes): 
    self.write() : contains optional kwds to only save ELEME and CONNE blocks (for iTOUGH)
    self.write_blocks() : add keyword: xyz = True/ False to include point centre into 
                          Element list (pos. 51 - 80)
        Caution: this required a change in t2data_parser: the type blocks_xyz was added
                 to the specification dictionary!!
    """    
    def __init__(self, filename = '', meshfilename = '', **kwds):
        """Initialisation of GeoT2Data object for extension of PyTOUGH with Geomodeller functionality
        
        **Optional Keywords**:
            - *print_debug* = True/False : set flag to print debug messages (default = False)
        """ 
        if kwds.has_key('print_debug'):
            self.print_debug = kwds['print_debug']
        else:
            self.print_debug = False
        
        super(GeoT2Data, self).__init__(filename, meshfilename)
       
    
    def write_bak(self, filename='', **kwds):
        """Writes data to file
        overwrites standard t2data.write() method with additional keywords to
        adjust output for iTOUGH2
        optional keywords:
        iTOUGH2 = True/ False: only save geometry definition, i.e. ELEME and CONNE blocks
        """
        if filename: self.filename=filename
        if self.filename=='': self.filename='t2data.dat'
        outfile=t2data_parser(self.filename,'w')
        if kwds.has_key('iTOUGH2') and kwds['iTOUGH2']:
            # save in format optimal for iTOUGH2
            write_fn=[self.write_blocks_xyz, self.write_connections]
        else:
            write_fn=[self.write_title,self.write_simulator,self.write_rocktypes,self.write_parameters,self.write_start,
                      self.write_rpcap,self.write_lineq,self.write_multi,self.write_times,self.write_blocks,
                      self.write_connections,self.write_generators,self.write_short_output,self.write_incons]
        for fn in write_fn: fn(outfile)
        if kwds.has_key('iTOUGH2') and kwds['iTOUGH2']:
            pass
        else:
            outfile.write('ENDCY\n')
        outfile.close()  
    

    def update_model_from_geomodeller_xml_file(self, geomodeller_xml_file, **kwds):
        """use gemodeller direct dll access to directly read geology data into shemat object
        idea: also update poro, perm, and all other parameters directly?
        Either: evaluate from shemat-file or read from (external) parameter file?
        optional kwds:
        compute = True/False: (re-) compute geological model before processing
        lower_left_x = float : x-position of lower-left corner (default: model range)
        lower_left_y = float : y-position of lower-left corner (default: model range)
        lower_left_z = float : z-position of lower-left corner (default: model range)
        save_grid = True/ False : save (mul)grid to file
        """
        import geomodeller_api as g_api
        g_api.initialise_geomodeller_api(r'C:\Geomodeller\GeoModeller_1502\bin')
        g_api.load_model(geomodeller_xml_file)
#        if kwds.has_key('compute') and kwds['compute']==True:
#            g_api.compute()
#        (geo_dx,geo_dy,geo_dz) = g_api.get_model_extent()
#        (xmin,ymin,zmin,xmax,ymax,zmax) = g_api.get_model_bounds()
        if kwds.has_key('lower_left_x') and kwds['lower_left_x'] != '':
            xmin = kwds['lower_left_x']
        if kwds.has_key('lower_left_y') and kwds['lower_left_y'] != '':
            ymin = kwds['lower_left_y']
        if kwds.has_key('lower_left_z') and kwds['lower_left_z'] != '':
            zmin = kwds['lower_left_z']
        
        i = 0
        total_i = len(self.grid.blocklist)-1
        for blk in self.grid.blocklist:
            print "%07d of %07d blocks, %02.f Percent" % (i,total_i,i/float(total_i)*100)
            if blk.name == 'ATM 0' or blk.name == 'at  0': continue
            rock_id = g_api.get_lithology(blk.centre[0], blk.centre[1], blk.centre[2])
            rock_tough_name = self.geomodel.tough_formation_names[self.geomodel.stratigraphy_list[rock_id-1]]
            try:
                blk.rocktype = self.grid.rocktype[rock_tough_name]
            except KeyError:
                pass
#                if rock_tough_name == '  OUT': # maybe atmosphere not defined in this model!
#                    continue
#                else:
#                    raise KeyError
            i += 1
#        geo_new = [] # new geology array
#        n_cells = len(self.centre_z) * len(self.centre_y) * len(self.centre_x)
#        i = 0
#        for z in self.centre_z:
#            for y in self.centre_y:
#                for x in self.centre_x:
#                    # geo_new.append(g_api.get_lithology(x,y,-geo_dz+z)) # old version
#                    print "process cell %8d of %d, %4.1f Percent" % (i, n_cells, (float(i)/float(n_cells)*100))
#                    geo_new.append(g_api.get_lithology(xmin + x, ymin + y, zmin + z))
#                    i += 1
#        self.set_array("GEOLOGY", geo_new)

        # remove unused rocktypes
        self.grid.clean_rocktypes()
        if kwds.has_key('save_grid') and kwds['save_grid']:
            print "Save grid to file"
            self.geo.write("grid.dat")
        if kwds.has_key("save_rocktypes") and kwds['save_rocktypes']:
            print("Save rocktype indices array to text file")
            if kwds.has_key("rock_filename"):
                f = open(kwds['rock_filename'],'w')
            else:
                f = open("rocktype_indices.csv",'w')
            for ind in self.grid.rocktype_indices:
                f.write("%d," % ind)
            f.close()
        print self.grid.rocktypelist
        
    def update_model_from_exported_grid(self, grid_file="exported_grid.txt", **kwds):
        """Update the self.grid.rocktypelist according to exported geomodel grid
        
        **Arguments**:
            - *grid_file*: filename of exported grid file (default="exported_grid.txt")
        
        **Optional keywords**:
            - *save_grid* = True/False : save grid to file in mulgrid format (grid.dat)
            - *save_rocktypes* = True/False : save rocktype list (for visualisation/ later analysis!)
            - *rock_filename* = filename of rocklist (default = "rocktype_indices.csv")
        """
        from numpy import mod
        total_i = len(self.grid.blocklist)-1
        i = total_i
        f = open(grid_file, "r")
        lines = f.readlines()
        f.close()
        if self.print_debug:
            print("Rocktypes:")
            print(self.grid.rocktype)
            print("TOUGH formation names:")
            print self.tough_formation_names
        for line in lines:
            l = line.split(',')
            l = l[:-1]
            l = l[::-1]
            for l1 in l:
                if self.print_debug:
                    if mod(i/float(total_i)*100,10) == 0:
                        print "%07d of %07d blocks, %02.f Percent" % (i,total_i,i/float(total_i)*100)
                try:
                    rock_tough_name = self.geomodel.tough_formation_names[self.geomodel.stratigraphy_list[int(l1)-1]]
                except AttributeError: # for the case that only loaded from grid (and not from a Geomodeller model)
                    rock_tough_name = self.tough_formation_names[int(l1)]
                self.grid.blocklist[i].rocktype = self.grid.rocktype[rock_tough_name]
                i -= 1

        # remove unused rocktypes
        self.grid.clean_rocktypes()
        if kwds.has_key('save_grid') and kwds['save_grid']:
            print "Save grid to file"
            self.geo.write("grid.dat")
        if kwds.has_key("save_rocktypes") and kwds['save_rocktypes']:
            print("Save rocktype indices array to text file")
            if kwds.has_key("rock_filename"):
                # f = open(kwds['rock_filename'],'w')
                rock_filename = kwds['rock_filename']
            else:
                # f = open("rocktype_indices.csv",'w')
                rock_filename = "rocktype_indices.txt"
            numpy.savetxt(rock_filename,self.grid.rocktype_indices,"%d")
            # for ind in self.grid.rocktype_indices:
            #     f.write("%d," % ind)
            # f.close()
        # print self.grid.rocktypelist
        
    
    def create_regular_mesh_from_geomodel(self, nx, ny, nz, **kwds):
        """create a regular mesh based on dimensions of Geomodeller Model
        => t2data.geomodel GeoModeller-Object has to be defined!
        nx, ny, nz = integer : Nodes in x, y, z-direction
        **Optional Keywords**:
            - *convention* = int : mesh convention (default: 2)
            - *atmos_type* = int : athmosphere type (default: 0)
            - *keep_rocktpyes* = True/False : keep current rocktypes (usually cleaned when grid created with mulgrids)
            - *save_mesh* = True/ False: save mesh to file (mesh.geo)
        """
        try:
            self.geomodel
        except AttributeError:
            print "\n\n" + 80 * "*" + "\n"
            print "\tGeomodel Object is not defined, please load Geomodel and try again"
            print "\tt2data.load_geomodel(geomodeller_xml_file)\n"
            print "\n" + 80 * "*" + "\n"
            raise
        # check atmosphere convention
        if kwds.has_key('atmos_type'):
            atmos_type = kwds['atmos_type']
        else:
            atmos_type = 0
        if kwds.has_key('convention'):
            convention = kwds['convention']
        else:
            convention = 2
        # determine extend and origin of geomodel and create mesh accordingly
        (range_x, range_y, range_z) = self.geomodel.get_model_range()
        (x_min, x_max, y_min, y_max, z_min, z_max) = self.geomodel.get_model_extent()
        # determine origin
        # Attention: origin in t2data defined form z top downwards?
        self.geomodel.origin = (x_min, y_min, z_max)
        if self.print_debug:
            print("Geomodel Origin:")
            print(self.geomodel.origin)
            print("Geomodel Range:")
            print(range_x, range_y, range_z)
        # now create mulgrid mesh
        spacing_x = range_x / nx
        spacing_y = range_y / ny
        spacing_z = range_z / nz
        if kwds.has_key('keep_rocktypes') and kwds['keep_rocktypes']:
            rtypes = [copy.deepcopy(self.grid.rocktype[s.name]) for s in self.grid.rocktypelist]
            if self.print_debug:
                print("Keep rocktypes")
        self.geo = mulgrid().rectangular([spacing_x] * nx, 
                                    [spacing_y] * ny, 
                                    [spacing_z] * nz, 
                                    origin=self.geomodel.origin, 
                                    convention=convention,
                                    atmos_type = atmos_type)
        if kwds.has_key('keep_rocktypes') and kwds['keep_rocktypes']:
            for rtype in rtypes:
                self.grid.add_rocktype(rtype)
                
        if kwds.has_key('save_mesh') and kwds['save_mesh']:
            self.geo.write("mesh.geo")
        
        # update mesh in tough2 object
        self.grid.fromgeo(self.geo)
        
    def create_rectilinear_mesh_from_geomodel(self, dx, dy, dz, **kwds):
        """Create a rectilinear mesh based on dimensions of Geomodeller Model
        
        Model is exported *relative* to passed array lists dx, dy, dz; i.e.: if dz is linearly
        increasing from 1 to 10 in 10 steps, then the mesh dimensions in z-direction
        will be linearly increasing in the same relative step sizes - but the step value
        itself is calculated from the extent of the geomodel.
        Example:
        dx = 10 * [1.]
        dy = 20 * [1.]
        dz = linspace(1,10,10)
        
        **Arguments**:
            - dx, dy, dz = float(3) : relative increemnts in each direction (will be scaled to model extent!)
            
        ..note:: absolute spacing can be forced with keyword absoulte = True

        **Optional Keywords**:
            - *origin* = float(x,y,z) : origin of grid (for subgrid of whole model range)
            - *extent* = float(extent_x,extent_y,extent_z) : extent of model in each direction (for submodel)
            - *convention* = int : mesh convention (default: 2)
            - *atmos_type* = int : athmosphere type (default: 0)
            - *keep_rocktpyes* = True/False : keep current rocktypes (usually cleaned when grid created with mulgrids)
            - *save_mesh* = True/ False: save mesh to file (mesh.geo)
            - *absolute* = True/ False: use absolute instead of relative spacing (default: False)

        ..note:: naming convention for tough block names (keyword convention):
                    0 : max. 17575 nodes in layer, max. 99 layers (total 1,739,925 nodes)
                    1 : max. 99 nodes in layer, max. 17575 layers (total 1,739,925 nodes)
                    2 : max. 999 nodes in layer, max. 675 layers (total 674,325 nodes) 


        """
        try:
            self.geomodel
        except AttributeError:
            print "\n\n" + 80 * "*" + "\n"
            print "\tGeomodel Object is not defined, please load Geomodel and try again"
            print "\tt2data.load_geomodel(geomodeller_xml_file)\n"
            print "\n" + 80 * "*" + "\n"
            raise
        if kwds.has_key('extent'):
            (extent_x, extent_y, extent_z) = kwds['extent']
        else:
            # determine extend and origin of geomodel and create mesh accordingly
            (extent_x, extent_y, extent_z) = self.geomodel.get_model_range()

        # determine origin
        # Attention: origin in t2data defined form z top downwards?
        if kwds.has_key('origin'):
            self.geomodel.origin = kwds['origin']
        else:
            (x_min, x_max, y_min, y_max, z_min, z_max) = self.geomodel.get_model_extent()
            self.geomodel.origin = (x_min, y_min, z_max)

        if kwds.has_key('convention'):
            convention = kwds['convention']
        else:
            convention = 0
        if kwds.has_key('atmos_type'):
            atmos_type = kwds['atmos_type']
        else:
            atmos_type = 0

        # now create mulgrid mesh
        # make sure that variables are arrays
        from numpy import array, sum
        # first: scale sum of array to 1, then multiply by range
        dx = array(dx)/sum(dx) * extent_x
        dy = array(dy)/sum(dy) * extent_y
        dz = array(dz)/sum(dz) * extent_z
        print "Create grid with %d cells\n" % (len(dx)*len(dy)*len(dz))
        if kwds.has_key('keep_rocktypes') and kwds['keep_rocktypes']:
            rtypes = [copy.deepcopy(self.grid.rocktype[s.name]) for s in self.grid.rocktypelist]
            if self.print_debug:
                print("Keep rocktypes")
        self.geo = mulgrid().rectangular(dx, dy, dz, 
                                    origin = self.geomodel.origin, 
                                    convention = convention,
                                    atmos_type = atmos_type)
        # update mesh in tough2 object
        if kwds.has_key('keep_rocktypes') and kwds['keep_rocktypes']:
            for rtype in rtypes:
                self.grid.add_rocktype(rtype)
                
                
        if kwds.has_key('save_mesh') and kwds['save_mesh']:
            self.geo.write("mesh.geo")
        

#        self.geo.write("grid.geo")
        self.grid.fromgeo(self.geo)

    def create_rectilinear_mesh_from_delxyz_file(self, delxyz_file, **kwds):
        """create a rectilinear mesh based on dimensions defined in an external (delxyz.txt) file
        
        convention = int : naming convention for tough block names (and generators?); default: 0
                    0 : max. 17575 nodes in layer, max. 99 layers (total 1,739,925 nodes)
                    1 : max. 99 nodes in layer, max. 17575 layers (total 1,739,925 nodes)
                    2 : max. 999 nodes in layer, max. 675 layers (total 674,325 nodes) 
        atmos_type = int : atmosphere type (default: 2)
        """

        if kwds.has_key('convention'):
            convention = kwds['convention']
        else:
            convention = 0
        if kwds.has_key('atmos_type'):
            atmos_type = kwds['atmos_type']
        else:
            atmos_type = 0

        # load dx,dy and dz from delxyz file
        lines = open(delxyz_file,"r").readlines()
        dels = []
        for line in lines:
            d_tmp = []
            l = line.split(",")
            for l2 in l:
                if "*" in l2:
                    l3 = l2.split("*")
                    for l4 in range(int(l3[0])):
                        d_tmp.append(float(l3[1]))
                else:
                    d_tmp.append(float(l2))
            print d_tmp
            dels.append(d_tmp)

        print "Create grid with %d cells\n" % (len(dels[0])*len(dels[1])*len(dels[2]))
        self.geo = mulgrid().rectangular(dels[0], dels[1], dels[2], 
                                    convention = convention,
                                    atmos_type = atmos_type)
        # update mesh in tough2 object
        # self.geo.write("grid.geo")
        self.grid.fromgeo(self.geo)


    def load_geomodel(self, geomodeller_xml_file):
        """Load geomodel and store object as attribute"""
        import Geomodeller_Obj_4 as GeoModel
        self.geomodel = GeoModel.GeomodellerClass()
        self.geomodel.load_geomodeller_file(geomodeller_xml_file)
        self.geomodel.get_stratigraphy_list(out='upper')
        self.geomodel.create_TOUGH_formation_names(out='upper')
        
    def update_properties_from_csv_list(self, csv_file_name, **kwds):
        """update the rock properties from csv list
        use same keywords as for SHEMAT files, i.e.:
        compulsary:
        NAME = str : Original formation name (as defined in GeoModel!)
        GEOLOGY = int : Geology identifier (not used if NAME is defined)
        PERM = float : Permeability in z-direction
        POR = float : porosity
        WLFM0 = float : thermal conductivity of rock
        
        optional:
        ANISOI = float : anisotropy factor for permeability in x-direction (perm_x = perm_z / ANISOI)
        ANISOJ = float : anisotropy factor for permeability in y-direction (perm_y = perm_z / ANISOJ)
        
        **Optional keywords**:
            - *generate_tough_names* = True/False : generate TOUGH2 compliant names from "NAME" column in csv file
        
        """
        
        """ read csv file """
        
        csv = open(csv_file_name,'r')
        csv_lines = csv.readlines()
        csv.close()
        csv_header = csv_lines[0].split(',')
        # remove end of line character from last entry
        csv_header[-1] = csv_header[-1].rstrip()
        # test if file contains another column (e.g. a 'Name' column) before GEOLOGY column
        if csv_header[0] == 'GEOLOGY':
            start = 1
        else:
            start = 2 # in this case, GEOLOGY column HAS TO BE in the second column!!
        csv_datalines = []
        for i_csv,l in enumerate(csv_lines[1:]):
#            print "line %d: %s" % (i_csv,l)
            """ check implementation of permeability functions for TOUGH!
            if "PERM_FUNC" in l:
                print "permeability function defined"
                # test correct definition
                if csv_lines[i_csv+2].split(',')[0] != "Type":
                    print "please define type (kwd Type) in input line %d" % i_csv+2
                    break
                if csv_lines[i_csv+3].split(',')[0] != "Lines":
                    print "please define number of lines (kwd Lines) in input line %d" % i_csv+3
                    break
                # deconstruct input file and create permeability_function_dict with formation id as key
                type = csv_lines[i_csv+2].split(',')[1]
                lines = csv_lines[i_csv+3].split(',')[1]
                perm_func = {}
                perm_func['type'] = type
                if type == 'linear_global':
                    for n_line in range(int(lines)):
                        line = csv_lines[i_csv+5+n_line].split(',')
                        func_dict = {'k_min' : float(line[start]),
                                     'k_min_depth' : float(line[start+1]),
                                     'k_max' : float(line[start+2]),
                                     'k_max_depth' : float(line[start+3])}
                        perm_func[line[start-1]] = func_dict
            if "PERM_RANDOM" in l:
                print "permeability random distribution defined"
                # test correct definition
                if csv_lines[i_csv+2].split(',')[0] != "Type":
                    print "please define type (kwd Type) in input line %d" % i_csv+2
                    break
                if csv_lines[i_csv+3].split(',')[0] != "Lines":
                    print "please define number of lines (kwd Lines) in input line %d" % i_csv+3
                    break
                # deconstruct input file and create permeability_random_dict with formation id as key
                type = csv_lines[i_csv+2].split(',')[1]
                lines = csv_lines[i_csv+3].split(',')[1]
                perm_rand = {}
                perm_rand['type'] = type
                if type == 'lognormal':
                    for n_line in range(int(lines)):
                        line = csv_lines[i_csv+5+n_line].split(',')
                        rand_dict = {'log_stdev' : line[start]}
                        perm_rand[line[start-1]] = rand_dict
            else:
            """
            if l[-1] == '\n':
                l = l[:-1]
            csv_datalines.append(l.split(','))
        if csv_header[0] != "GEOLOGY" and csv_header[1] != "GEOLOGY":
            print "first or second column should contain GEOLOGY data!"
            return None
        
        """ update properties in tough grid file """
        
        for b in self.grid.blocklist:
            # print b
            self.grid.delete_block(b)
        
        # to store tough names (in case they are generated)
        tough_name_list = []
        self.tough_formation_names = {}

        for line in csv_datalines:
            # csv input file is parsed based on SHEMAT variable names!
            
            # option: create tough formation names directly from name in csv file, if "NAME" is defined
            if kwds.has_key("generate_tough_names") and kwds['generate_tough_names']:
                # create tough2 compiant names from 'NAME' in csv_header
                formation_name = line[csv_header.index('NAME')] 
                cut_name = formation_name[0:5] 
                if cut_name in tough_name_list:
                    for j in range(100):
                        if "%4s%d" % (cut_name[0:4],j) in tough_name_list:
                            continue
                        else: 
                            cut_name = "%4s%1d" % (cut_name[0:4],j)
                            tough_name_list.append(cut_name)
                            break
                else:
                    tough_name_list.append(cut_name)        
                rockname = "%5s" % str.upper(cut_name)
                self.tough_formation_oames[int(line[csv_header.index('GEOLOGY')])] = rockname
                
                
                print rockname
            
            elif 'NAME' in csv_header: 
                try:       
                    rockname = self.geomodel.tough_formation_names[line[csv_header.index('NAME')]]
                except KeyError: # might be an empty line in csv file
                    if line[csv_header.index('NAME')]== '': continue
                    else: 
                        print csv_header.index('NAME')
                        print line
                        raise KeyError
                    
            else:
                print "\n\n\tNo column with title NAME found in csv-file! Please check and try again!\n\n"
                raise NameError
            if 'PERM' in csv_header:        
                perm_z = float(line[csv_header.index('PERM')])
            else:
                print "\n\n\tNo column with title PERM found in csv-file! Please check and try again!\n\n"
                raise NameError
            if 'ANISOI' in csv_header:
                anisoi = float(line[csv_header.index('ANISOI')])
                perm_x = perm_z * anisoi
            else:
                perm_x = perm_z
            if 'ANISOJ' in csv_header:
                anisoj = float(line[csv_header.index('ANISOJ')])
                perm_y = perm_z * anisoj
            else:
                perm_y = perm_z
            if 'POR' in csv_header:        
                porosity = float(line[csv_header.index('POR')])
            else:
                print "\n\n\tNo column with title POR found in csv-file! Please check and try again!\n\n"
                raise NameError
            if 'HPR' in csv_header:        
                hpr = float(line[csv_header.index('HPR')])
            else:
                if self.print_debug:
                    print "\n\n\tNo column with title HPR found in csv-file."
            if 'WLFM0' in csv_header:        
                th_cond = float(line[csv_header.index('WLFM0')])
            else:
                print "\n\n\tNo column with title WLFM0 found in csv-file! Please check and try again!\n\n"
                raise NameError
            new_rock = rocktype(rockname,
                                permeability = [perm_x,perm_y,perm_z],
                                conductivity = th_cond,
                                porosity = porosity)
            self.grid.add_rocktype(new_rock)
            
            
    def basal_heat_flux(self,q):
        """ assign generators to all cells in bottom of model to define
        a uniform basal heat flux of q W/m^2
        """
#        self.clear_generators()
        layer = self.geo.layerlist[-1] # get bottom layer
        cols = [col for col in self.geo.columnlist] #  if 10.e3 <= col.centre[0] <= 20.e3]
        i = 0
        for col in cols:
            blockname = self.geo.block_name(layer.name, col.name)
            # Attention with naming convention for generators: 
            # first three characters are arbitrary, the last two have to be numbers!!
            gen = t2generator(name='q%4d' % i, block=blockname, type='HEAT',gx=q*col.area)
#            gen = t2generator(name=' q'+col.name, block=blockname, type='HEAT',gx=q*col.area)
#            gen = t2generator(name=('q%04s' % col.name), block=blockname, type='HEAT', gx=q*col.area)
            self.add_generator(gen)
            i += 1
    
    def assign_basal_layer(self):
        """assing another rocktype to baselayer for constant temperature BC at bottom
        (setting density to a very high value)
        """
        # add base rock type
        new_rock = rocktype('BASE ',
                            permeability = [1E-17,1E-17,1E-17],
                            conductivity = 3.,
                            porosity = 0.1,
                            density = 1E50)
        self.grid.add_rocktype(new_rock)

        layer = self.geo.layerlist[-1] # get bottom layer
        cols = [col for col in self.geo.columnlist] #  if 10.e3 <= col.centre[0] <= 20.e3]
        i = 0
        for col in cols:
            blockname = self.geo.block_name(layer.name, col.name)
#            self.grid.blocklist[blockname].rocktype
            self.grid.blocklist[self.grid.block_index(blockname)].rocktype = 'BASE '
            i += 1
        
    def add_well_block(self, well_position):
        """add a well block at this position (not a generator, but different rock type!)
        well_position = (float,float,float) : coordinates of block
        attention: works only with mulgrid grid structure!
        """
#        layer = self.geo.layer_containing_elevation(well_position[2])
        # define new rocktype
        new_rock = rocktype('WELL ',
                            permeability = [1E-18,1E-18,1E-18],
                            conductivity = 3.,
                            porosity = 0.01,
                            density = 2.6E3)
        self.grid.add_rocktype(new_rock)
        # now change block at well position to this rock type
        blockname = self.geo.block_name_containing_point(well_position)
        self.grid.blocklist[self.grid.block_index(blockname)].rocktype = 'WELL '
        
    def set_atmosphere_block_name(self,atmosphere_rocktype):
        """set atmopshere block to this name (e.g. to define another rock type for it"""
        try:
            self.grid.blocklist[self.grid.block_index('at  0')].rocktype = atmosphere_rocktype
        except TypeError:
            self.grid.blocklist[self.grid.block_index('ATM 0')].rocktype = atmosphere_rocktype

#    def write_blocks_xyz(self,outfile):
#        """Write blocks with block centre values, e.g. for use in iTOUGH2"""
#        if self.grid.num_blocks>0:
#            outfile.write('ELEME\n')
#            from copy import copy
#            for blk in self.grid.blocklist:
#                blkw=copy(blk.__dict__)
#                blkw['name']=unfix_blockname(blkw['name'])
#                if blk.centre==None: outfile.write_value_line(blkw,'blocks')
#                else:
#                    vals=[blkw['name'],blk.nseq,blk.nadd,blk.rocktype.name,blk.volume,
#                          blk.ahtx,blk.pmx]+list(blk.centre)
#                    outfile.write_values(vals,'blocks')
#            outfile.write('\n')


class GeoMulgrid(mulgrid):
    """Overwrite some functionality of mulgrids class for own requirements"""
    
    def slice_plot(self, line=None, variable=None, variable_name=None, unit=None, block_names=None, colourmap=None, linewidth=0.2,
                   linecolour='black', aspect='auto', plt=None, subplot=111, title=None, xlabel=None, ylabel='elevation (m)',
                   contours=False, contour_label_format='%3.0f', contour_grid_divisions=(100,100), colourbar_limits=None,
                   plot_limits=None, column_axis = False, layer_axis = False, wells = None, well_names = True,
                   hide_wells_outside = False, wellcolour = 'blue', welllinewidth = 1.0, wellname_bottom = False,
                   rocktypes = None, allrocks = False, rockgroup = None, flow = None, grid = None, flux_matrix = None,
                   flow_variable_name = None, flow_unit = None, flow_scale = None, flow_scale_pos = (0.5, 0.02),
                   flow_arrow_width = None, **kwds):
        """Produces a vertical slice plot of a Mulgraph grid, shaded by the specified variable (an array of values for each block).
       
       ..Note: this version is slightly adapted for geopytough but mainly the original mulgrids.slice_plot()
       
       A unit string can be specified for annotation.  Block names can be optionally superimposed, and the colour 
       map, linewidth, aspect ratio, colour-bar limits and plot limits specified.
       
       If no variable is specified, only the grid is drawn, without shading.  If no line is specified, a slice
       through the grid bounds is made (bottom left to top right).
       
       If a string 'x' or 'y' is passed in instead of a line, a plot is made through the centre of the grid along
       the x- or y-axes, and the coordinate along the slice represents the actual x- or y- coordinate.  
       
       If a northing
       (float, in degrees) is passed instead of a line, a plot is made through the centre along the specified northing direction.
       
       Adaptations for geopytough all controlled through additional keywords:
       
       **Optional keywords (for geopytough)**:
           - *cbar_orientation* = 'vertical'/ 'horizontal' : orientation of colourbar
           - *show* = True/ False: show figure (set to false for non-interactive sessions)
           - *rocknames* = [string] : custom rocknames instead of TOUGH2 names 
       
       """
        if kwds.has_key('cbar_orientation'):
            cbar_orientation = kwds['cbar_orientation']
        else:
            cbar_orientation = 'vertical'
       
        if line is None:
            l=self.bounds
            default_title='vertical slice across grid bounds'
        elif isinstance(line,str):
            axislines={'x':[np.array([self.bounds[0][0],self.centre[1]]),np.array([self.bounds[1][0],self.centre[1]])],
                       'y':[np.array([self.centre[0],self.bounds[0][1]]),np.array([self.centre[0],self.bounds[1][1]])]}
            if line in axislines:
                l=axislines[line]
                default_title='vertical slice along '+line+' axis'
            else:
                l=self.bounds
                default_title='vertical slice across grid bounds'
        elif isinstance(line,(float,int)):
            r=0.5*norm(self.bounds[1]-self.bounds[0])
            from math import radians,cos,sin
            theta=radians(line)
            d=r*np.array([sin(theta),cos(theta)])
            l=[self.centre-d,self.centre+d]
            default_title='vertical slice '+("%3.0f"%float(line)).strip()+'$^\circ$N'
        else:
            l=line
            default_title='vertical slice from ('+("%7.0f"%l[0][0]).strip()+','+("%7.0f"%l[0][1]).strip()+') to ('+("%7.0f"%l[1][0]).strip()+','+("%7.0f"%l[1][1]).strip()+')'
        if norm(l[1]-l[0])>0.0:
            # matplotlib.rcParams.update({'mathtext.default': 'regular','figure.figsize':(12,9)}) 
            if kwds.has_key('ax'):
                ax = kwds['ax']
            else:
                if plt is None:
                    import matplotlib
                    import matplotlib.pyplot as plt
                    loneplot=True
                else: loneplot=False
                ax = plt.subplot(subplot,aspect=aspect)
            if variable is not None:
                if len(variable)==self.num_columns: variable=self.column_values_to_block(variable)
            if variable_name: varname=variable_name
            else: varname='Value'
            if rocktypes: variable, varname = rocktypes.rocktype_indices, 'Rock type'
            if block_names:
                if not isinstance(block_names,list): block_names=self.block_name_list
            else: block_names=[]

            track=self.column_track(l)
            if track:

                if xlabel is None:
                    if line=='x': xlabel='x (m)'
                    elif line=='y': xlabel='y (m)'
                    else: xlabel='distance (m)'
                ax.set_xlabel(xlabel)
                if column_axis: colnames, colcentres = [],[]
                verts,vals=[],[]
                if not isinstance(contours,bool): contours=list(contours)
                if contours != False or flow is not None: xc,yc=[],[]
                if flow is not None:
                    if flow_variable_name is None: flow_variable_name = 'Flow'
                    if flow_unit is None: flow_unit = 'units'
                    if grid is None:
                        from t2grids import t2grid
                        grid = t2grid().fromgeo(self)
                    if flux_matrix is None: flux_matrix = grid.flux_matrix(self)
                    blkflow = flux_matrix * flow
                    blkflow = blkflow.reshape((self.num_underground_blocks,3))
                    natm = self.num_atmosphere_blocks
                    U,V = [],[]
                    slice_dirn = (l[1] - l[0]).T
                    slice_dirn /= np.linalg.norm(slice_dirn) # normal vector in slice direction

                for trackitem in track:
                    col,points=trackitem[0],trackitem[1:]
                    inpoint=points[0]
                    if len(points)>1: outpoint=points[1]
                    else: outpoint=inpoint
                    if line=='x': din,dout=inpoint[0],outpoint[0]
                    elif line=='y': din,dout=inpoint[1],outpoint[1]
                    else: din,dout=norm(inpoint-l[0]),norm(outpoint-l[0])
                    dcol = 0.5*(din+dout)
                    if column_axis: colnames.append(col.name); colcentres.append(dcol)
                    for lay in self.layerlist[1:]:
                        if col.surface>lay.bottom:
                            blkname=self.block_name(lay.name,col.name)
                            if variable is not None: val=variable[self.block_name_index[blkname]]
                            else: val=0
                            vals.append(val)
                            top = self.block_surface(lay,col)
                            centre = self.block_centre(lay,col)
                            verts.append(((din,lay.bottom),(din,top),(dout,top),(dout,lay.bottom)))
                            if blkname in block_names:
                                ax.text(dcol, centre[2], blkname, clip_on = True, horizontalalignment = 'center')
                            if contours != False or flow is not None:
                                xc.append(dcol); yc.append(centre[2])
                            if flow is not None:
                                blkindex = self.block_name_index[blkname] - natm
                                q = blkflow[blkindex]
                                qslice = np.dot(slice_dirn, q[:2])
                                U.append(qslice)
                                V.append(q[2])

                import matplotlib.collections as collections
                if variable is not None: facecolors=None
                else: facecolors=[]
                if rocktypes: 
                    vals, rocknames, colourmap, colourbar_limits = \
                        self.setup_rocktype_plot(rocktypes, vals, colourmap, allrocks, rockgroup)
                    if kwds.has_key('rocknames'):
                        rocknames = kwds['rocknames']
                else: rocknames, rocktypes = None, None
                col=collections.PolyCollection(verts,cmap=colourmap,linewidth=linewidth,facecolors=facecolors,edgecolors=linecolour)
                if variable is not None: col.set_array(np.array(vals))
                if colourbar_limits is not None: col.norm.vmin,col.norm.vmax=tuple(colourbar_limits)
                ax.add_collection(col)
                if plot_limits is not None:
                    ax.set_xlim(plot_limits[0]); ax.set_ylim(plot_limits[1])
                else: ax.autoscale_view()
                if contours<>False:
                    from matplotlib.mlab import griddata
                    xc,yc=np.array(xc),np.array(yc)
                    valc=np.array(vals)
                    bds=((np.min(xc),np.min(yc)),(np.max(xc),np.max(yc)))
                    xgrid=np.linspace(bds[0][0],bds[1][0],contour_grid_divisions[0])
                    ygrid=np.linspace(bds[0][1],bds[1][1],contour_grid_divisions[1])
                    valgrid=griddata(xc,yc,valc,xgrid,ygrid)
                    if isinstance(contours,list): cvals=contours
                    else: cvals=False
                    CS=plt.contour(xgrid,ygrid,valgrid,cvals,colors='k')
                    if contour_label_format is not None: plt.clabel(CS, inline=1,fmt=contour_label_format)
                ax.set_ylabel(ylabel)
                scalelabel=varname
                if unit: scalelabel+=' ('+unit+')'
                if variable is not None:
                    self.plot_colourbar(plt, col, scalelabel, rocktypes, rocknames, cbar_orientation = cbar_orientation)
                    default_title=varname+' in '+default_title
                if column_axis:
                    ax.set_xticks(colcentres)
                    ax.set_xticklabels(colnames)
                    ax.set_xlabel('column')
                if layer_axis:
                    ax.set_yticks([lay.centre for lay in self.layerlist])
                    ax.set_yticklabels([lay.name for lay in self.layerlist])
                    ax.set_ylabel('layer')
                self.slice_plot_wells(plt, ax, line, l, wells, well_names, hide_wells_outside, wellcolour, welllinewidth, wellname_bottom)
                if flow is not None: self.plot_flows(plt, xc, yc, U, V, flow_variable_name, flow_unit, flow_scale,
                                                     flow_scale_pos, flow_arrow_width)
                if title is None: title=default_title
                if kwds.has_key("plot_points"):
                    for i,point in enumerate(kwds['plot_points']):
                        ax.plot(point[0],point[1],'rs')
                        if kwds.has_key('point_labels'):
                            ax.text(point[0]+30,point[1]+105,kwds['point_labels'][i], 
                                    va = 'center',
                                    fontsize = 16, color='blue', weight = 'bold', 
                                    bbox={'color' : 'blue','facecolor' : 'white', 'pad' : 10})
                plt.title(title)
                if kwds.has_key("show") and (kwds['show'] == False):
                    pass
                else:
                    if loneplot: plt.show()
            else: print 'Slice',str(line),'does not intersect the grid.'
#        return plt

    def plot_colourbar(self, plt, col, scalelabel, rocktypes, rocknames, **kwds):
        """Draws colour bar on a layer or slice plot."""
        if kwds.has_key('cbar_orientation'):
            cbar_orient = kwds['cbar_orientation']
        else:
            cbar_orient = 'vertical'
        cbar = plt.colorbar(col, orientation = cbar_orient)
        cbar.set_label(scalelabel)
        if rocktypes:
            cbar.set_ticks([i+0.5 for i in range(len(rocknames))])
            cbar.set_ticklabels(rocknames)
            cbar.ax.invert_yaxis() # to get in same top-down order as in the data file



def create_tough_sim_dir(**kwds):
    """create directory to store created tough2 and mulgrid files, 
    default: .\tough_sim_n
    n is automatically asigned as continuous numbers, starting with "1"
    but can also be manually set with kwds: dir = str
    optional kwds:
    basename = string : basename of directory, number added
    n_format = int : number of integers (for %0nd format description)
    """
    from os import chdir, path, mkdir, getcwd
    if kwds.has_key(dir):
        mkdir(dir)
    else:
        if kwds.has_key('basename'):
            basename = kwds['basename']
        else:
            basename = 'tough_sim'
        for n in range(1,100000): # sets a maximum of 100 result directories...
            if not path.isdir(path.join(".","%s_%d" % (basename,n))):
                print "create results directory %s\\%s_%02d" % (getcwd(), basename, n)
                resultsdir = path.join(".","%s_%d" % (basename, n))
                mkdir(resultsdir)
                break
    return resultsdir



def copy_files_export_model(geomodel_dir, geomodeller_file, new_dir, nx, ny, nz, **kwds):
    """Copy all generated models into a temporary directory, compute geomodel and export (regular) grid
    
    **Arguments**:
        - *geomodel_dir* = string : path to directory with all (original) geomodel files (xml + others)
        - *geomodeller_file* = filename : filename of Geomodeller xml file 
        - *nx, ny, nz* = integers : number of cells in each coordinate direction (grid is regular per default)
        - *new_dir* = dir-path : new directory for model export
        
    **Optional Keywords**:
        - *nz_list* = list of int : list of nz-values for discretization study
    """
    
    
    
    os.chdir(geomodel_dir)
    
    
    # first step: determine all Geomodeller files in original directory
    geomodel_files = []
    for l in os.listdir("."):
        if os.path.splitext(l)[1] == ".sec":
            geomodel_files.append(l)
        if os.path.splitext(l)[1] == ".s3d":
            geomodel_files.append(l)
        if os.path.splitext(l)[1] == ".md5":
            geomodel_files.append(l)
        if os.path.splitext(l)[1] == ".wsp":
            geomodel_files.append(l)
        if os.path.splitext(l)[1] == ".xml":
            # This is the original file - keep for later
            geomodel_xml_ori = l
    
    print geomodel_xml_ori
    os.chdir(new_dir)

    for f in geomodel_files:
        shutil.copyfile(os.path.join(geomodel_dir, f), os.path.join(".",f))

    for l2 in os.listdir("."):
        if l2 == geomodeller_file:
            print l2
            shutil.copyfile(l2,os.path.join(".", geomodel_xml_ori))
            if "original" in l2:
                new_grid = "exported_grid_original.txt"
                new_delxyz = "delxyz_original.txt"
            else:
                new_grid = "exported_grid_" + os.path.splitext(l2)[0] + ".txt"
                new_delxyz = "delxyz_" + os.path.splitext(l2)[0] + ".txt"
            
            if kwds.has_key('nz_list'):  # cell discretization study
                try:
                    i = int(l2[-7:-4])
                except ValueError:
                    os.chdir('..') 
                    continue # original model
                subprocess.call("/home/flo/bin/export_model %s %d %d %d" % (geomodel_xml_ori, nx, ny, kwds['nz_list'][i]), shell=True)
            else:

                while True:
                    proc = subprocess.Popen("/home/flo/bin/export_model %s %d %d %d" % (geomodel_xml_ori, nx, ny, nz), shell=True)
                    for i in range(10):
                        time.sleep(5)
                        if proc.poll() == 0:
                            print("Process finished successfully")
                            break
                    if proc.poll() == None:
                        proc.kill()
                        print("\n\n\tProcess killed, restart!\n\n")
                        continue
                    elif proc.poll() == 0:
                        print("Process finished successfully")
                        break
                    print("\n\n\tProcess not finished correctly!\n\n\n")
                    break
                # subprocess.call("/home/flo/bin/export_model %s %d %d %d" % (geomodel_xml_ori, nx, ny, nz), shell=True)
            shutil.copyfile("exported_grid.txt", os.path.join(".", new_grid))
            shutil.copyfile("delxyz.txt", os.path.join(".", new_delxyz))





