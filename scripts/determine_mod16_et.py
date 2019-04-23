#This script downloads MOD16 files for a specific site, then clips them to the buffer zone around that site
#it outputs MOD16 values that I then copied and pasted into another script (analysis.py)
from flux_module import *

#this shapefile can be found in the repository
in_towers = r'C:\Users\Devin Simmons\Desktop\GEOL393\GIS\MOD16_2014_ET_Annual\flux_towers\fluxnet_sites\fluxnet_sites.shp'
in_rasters = r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\mod16_tiles\h13v04'
out_folder = r"C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\clip_features"

#downloading rasters for h13v04 from 2004 to 2010 at tower ca_na1
ca_na1_years = [i for i in range(2004, 2010)]
#creates a list of the julian days that correspond to MOD16 8day periods
ca_na1_days = [str(i).zfill(3) for i in range(1, 365, 8)]

#sets parameters needed to download the right files
#filepath is the folder where the rasters will be downloaded
#this folder needs a subfolder with the name of the tile being downloaded, in this case, h13v04
ca_na1_rasters = download_modis(ca_na1_years, ca_na1_days, "h13v04",
                                r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\mod16_tiles')
#actually downloads these files
ca_na1_rasters.download_mod16()

in_rasters = ca_na1_rasters.mod16_folder

#makes a buffer around the CA-Na1 tower
test1 = modis_et(in_towers, "CA-Na1", out_folder)
test1.make_buffer(1692)


#without these lines the program had trouble reading the HDF files
arcpy.env.workspace = in_rasters
rasters = arcpy.ListRasters("*", "HDF")


mod16_values = {}
#clips the raster to the tower, determines its ET value, adds the julian day
#that the 8 day period starts to a dictionary as a key with its ET value as the value
for raster in rasters:
    test1.clip_raster(raster)
    print("clipped " + raster)
    mod16_values[test1.julian_day] = test1.mean_et

print(mod16_values)
