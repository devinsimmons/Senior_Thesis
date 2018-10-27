#this script creates the class modis_et. This class can take the shapefile containing flux tower locations,
#select a specific flux tower, create a buffer around it at the desired distance, and then clip the MOD16
#raster to that buffer. It can also determine the mean of the raster pixels in this buffer area. This average is 
#compared to the flux tower's measurements over the same 8 day period

import arcpy 
from os import listdir

arcpy.env.overwriteOutput = True

class modis_et:
    '''
    initalization of the class requires a shapefile with flux towers, a tower name,
    and an output folder
    the tower in question is identified with the tower_name string.
    tower_name should equal the tower's fluxnetid
    
    '''
    def __init__(self, flux_towers, tower_name, output_env):
        #flux towers is a shapefile
        self.flux_towers = flux_towers
        #tower name is a string that matches the "fluxnetid" attribute
        self.tower_name = tower_name
        #folder path
        self.output_env = output_env

    '''
    buffer inputs should be floats and the units should be meters
    this function buffers around the tower using the first_buffer, clips the 
    raster to the extent of the buffer
    '''
    def make_buffer(self, first_buffer):
        self.first_buffer = first_buffer
        #first this turns the tower shapefile into a layer
        #then it selects the specified tower and makes a layer out of it
        tower_layer = arcpy.MakeFeatureLayer_management(self.flux_towers, "tower_layer")
        tower = arcpy.SelectLayerByAttribute_management(tower_layer, "NEW_SELECTION",
                                                        ''' "fluxnetid" = '{}' '''.format(self.tower_name))
        
        #makes a layer out of the chosen tower, buffers it, creates a shapefile of the buffer
        arcpy.CopyFeatures_management(tower, self.output_env + r"/chosen_tower.shp")
        tower = self.output_env + r"/chosen_tower.shp"
        
        self.buffer_fc = self.output_env + r"/buffered_tower.shp"
        arcpy.Buffer_analysis(tower, self.buffer_fc,
                              str(self.first_buffer) + " Meters")

        
    '''
    clips raster to the specified buffer for the flux tower
    '''
    def clip_raster(self, mod16_raster):
        self.mod16_raster = mod16_raster
        #this extracts the julian day the raster starts from the file name
        #it looks better this way
        self.mod16_ET = self.output_env + r"/{}".format(self.mod16_raster[9:16] + "ET.tif")

        #extracts just the ET layer from the MOD16 hdf. this is located at index 0
        arcpy.ExtractSubDataset_management(self.mod16_raster, self.mod16_ET, "0")

        self.clipped_et = self.mod16_ET[0:-4] + "{}_clip.tif"
        arcpy.Clip_management(self.mod16_ET,
                              '0 0 0 0', self.clipped_et,
                              self.buffer_fc, '#', "ClippingGeometry")

    '''
    returns mean et in millimeters for the clipped raster
    '''
    @property
    def mean_et(self):
        mean = arcpy.GetRasterProperties_management(self.clipped_et, "MEAN")
        return float(mean.getOutput(0))/10

    

in_towers = r'C:\Users\Devin Simmons\Desktop\GEOL393\GIS\MOD16_2014_ET_Annual\flux_towers\fluxnet_sites\fluxnet_sites.shp'
in_rasters = r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\mod16_tiles'
out_folder = r"C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\clip_features"

test1 = modis_et(in_towers, "CA-Qcu", out_folder)
test1.make_buffer(1692)

#without these lines the program had trouble reading the HDF files
arcpy.env.workspace = in_rasters
rasters = arcpy.ListRasters("*", "HDF")

for raster in rasters:
    test1.clip_raster(raster)
    print test1.mod16_raster[9:16]
    print test1.mean_et
    




        
        
        
        
                                                
                                        
        
        

    
        


