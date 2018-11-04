import arcpy 

arcpy.env.overwriteOutput = True

class modis_et:
    '''
    initalization of the class requires a shapefile with flux towers, a tower name,
    and an output folder
    the tower in question is identified with the tower_name string.
    tower_name should equal the tower's fluxnetid
    
    '''
    def __init__(self, flux_towers, tower_name, output_env):
        self.flux_towers = flux_towers
        self.tower_name = tower_name
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
        self.julian_day = self.mod16_raster[9:16]
        self.mod16_ET = self.output_env + r"/{}".format(self.julian_day + "ET.tif")

        #extracts just the ET layer from the MOD16 hdf. this is located at index 0
        arcpy.ExtractSubDataset_management(self.mod16_raster, self.mod16_ET, "0")

        self.clipped_et = self.mod16_ET[0:-4] + "_clip.tif"
        arcpy.Clip_management(self.mod16_ET,
                              '0 0 0 0', self.clipped_et,
                              self.buffer_fc, '#', "ClippingGeometry")

    '''
    returns mean et in millimeters
    '''
    @property
    def mean_et(self):
        mean = arcpy.GetRasterProperties_management(self.clipped_et, "MEAN")
        return float(mean.getOutput(0))/10

    




import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

class flux_reader:
    
    '''
    infile is a filepath string
    read_file opens the input csv file, gets content into a list,
    then closes the file
    '''
    def __init__(self, infile):
        self.infile = infile
        tower_name = infile.split("\\")[-1]
        self.tower = tower_name[4:10]
        fileobj = open(self.infile, 'r')
        self.content = fileobj.readlines()
        fileobj.close()
    
    '''
    this function reteurns a dictionary where each year is a key with a dict value
    the dict value contains error based on energy balance closure discrepancy
    '''
    def yearly_error(self):
        self.year_dict = {}
        for i in range(5, len(self.content)):
            entries = self.content[i].split(',')    
            #slice year from datetime start column
            year = entries[0][0:4]
            
            if year not in self.year_dict:
                self.year_dict[year] = {}
                self.year_dict[year]['left error'] = 0
                self.year_dict[year]['right error'] = 0
            
            #variables that we care about for determining error
            h = float(entries[9]) #sensible heat turbulent flux, w/m^2
            g = float(entries[13]) #soil heat flux, w/m^2
            rn = float(entries[24]) #net radiation, w/m^2
            le = float(entries[11]) #latent heat turbulent flux
            #this may be useful for other sites
            #s = float(entries[8]) #heat storage flux in biomass
            #checks if these variables were measured and determines the error
            if h != -9999 and g != -9999 and rn != -9999 and le != -9999:
                #according to the energy balance closure error equation, Rn - G should equal LE + H
                left = rn - g
                right = le + h
                self.year_dict[year]['left error'] += left
                self.year_dict[year]['right error'] += right
        
        #after the year dict is filled out this determines the error based on energy closure
        #balance discrepancy
        for i in self.year_dict:
            left = self.year_dict[i]['left error']
            right = self.year_dict[i]['right error']
            if left != 0 and right != 0:
                #this was my method of determining energy balance closure error. it probably is bad
                #year_dict[i]['error'] = (abs(left - right)/right + abs(left - right)/left)/2
                #method that follows twine et al. (2002)
                self.year_dict[i]['error'] = abs(1 - right/left)
        
        return self.year_dict
    
    '''
    goes through all lines in the file, gets ET for each day based on half hour measurements,
    fills gaps based on values from adjacent days
    '''
    def et_by_day(self):
        self.et_dict = {}
        
        for i in range(100, len(self.content) - 100):
            entries = self.content[i].split(',')
    
            #slices timestamp to just the date  
            date = entries[0][0:8]
            
            #adds each date to the et dictionary
            if date not in self.et_dict:
                self.et_dict[date] = {}
                self.et_dict[date]['et'] = 0
                self.et_dict[date]['measurements'] = 0
                    
            #variables needed to determine evapotranspiration
            le = float(entries[11]) #latent heat turbulent flux
            ta = float(entries[3]) #air temp
    
            #determines if there are missing values for LE or temperature. if not the
            #evapotranspiration is calculated
            
            #commented out. for now I'm looking into using a latent heat
            #of vaporization constant
            if le != -9999 and ta != -9999:
            #if float(entries[11]) != -9999:
                #evapotranspiration = float(entries[11]) * 60 * 30 / ((2.501 - 0.002361 * float(entries[3])) * 10 ** 6)
                evapotranspiration = le * 60 * 30 / ((2.501) * 10 ** 6)
                #cumulates the et measurements for each day and adds a value
                #to show how many measurements were taken
                self.et_dict[date]['et'] += evapotranspiration
                self.et_dict[date]['measurements']+= 1
       
        
            #if there are missing values ET is determined based on values in adjacent rows
            else:
                adjacent_et = []
                
                #two days prior and previous day at the same time, respectively
                entries_neg_2 = self.content[i - 96].split(',')
                entries_neg_1 = self.content[i - 48].split(',')
                
                #next two days
                entries1 = self.content[i + 48].split(',')
                entries2 = self.content[i + 96].split(',')
                
                adjacent_days = [entries_neg_2, entries_neg_1, entries1, entries2]
                
                #determines if there are et values in adjacent days
                #if there are it adds them to a list
                for day in adjacent_days:
                     le = float(day[11])
                     ta = float(day[3])
                     #commented out. for now I'm looking into using a latent heat
                     #of vaporization constant
                     if le != -9999 and ta != -9999:
                     #if float(day[11]) != -9999:
                         #adjacent_et.append(float(i[11]) * 60 * 30 / ((2.501 - 0.002361 * float(i[3])) * 10 ** 6))
                         adjacent_et.append(le * 60 * 30 / ((2.501) * 10 ** 6))
                #we can only use this cocrrection if we get two values
                if len(adjacent_et) >= 2:
                    evapotranspiration = sum(adjacent_et) / len(adjacent_et)
                    #cumulates the et measurements for each day and adds a value
                    #to show how many measurements were taken
                    self.et_dict[date]['et'] += evapotranspiration
                    self.et_dict[date]['measurements']+= 1
        return self.et_dict
    
    '''
    this function can only be called after et_by_day has been called
    it returns a dictionary that displays all the valid 8 day ranges recorded
    by the flux tower for a given year. year can be an int or a string
    '''
    def valid_8day_by_year(self, year):
        from datetime import datetime
        
        self.year = year
        #will contain the julian day an 8-day period starts, averaged cumulative
        #ET over that period
        self.flux_et = {}
        
        #julian day counter, step value of 8 corresponds to 8 day period
        for i in range (1, 365, 8):
            
            et_counter = 0
            measurement_counter = 0
            
            i = str(i)
            #add leading zeroes to julian days less than 100
            if len(i) < 3:
                if len(i) < 2:
                    i = "00" + i
                else:
                    i = "0" + i
            #julian day
            i = str(year) + i
            
            for day in range(int(i), int(i) + 8):
                #jd361 corresponds to a period shorter than 8 days
                if day < 2008361:
                    #convert 
                    jd = datetime.strptime(str(day), '%Y%j')
                    dt = datetime.strftime(jd, '%Y%m%d')
                    #need all 48 measurements filled to be 
                    if self.et_dict[dt]['measurements'] == 48:
                        et_counter += self.et_dict[dt]['et']
                        measurement_counter += 1
            #need at least five days with all 48 measurements filled in
            if measurement_counter > 4:
                et_averaged = (et_counter * 8) / measurement_counter
                self.flux_et[i] = et_averaged
        
        return self.flux_et
    
    '''
    makes plots that compare et data to another dataset (ie MOD16 data)
    date_list is a list of ints that represent the start julian day for the
    8 day comparison period. can only be called after et_by_day and yearly_error.
    both flux_et and mod16_et should be dictionaries where the key is
    the julian day (YYYYDDD) and the value is ET. 
    '''
    def linreg_et_plot(self, flux_et, mod16_et):
        self.flux_et = flux_et
        self.mod16_et = mod16_et

        flux_keys = [key for key in self.flux_et]
        flux_keys.sort()
        
        flux_et_sorted = []
        modis_et_sorted = []

        for key in flux_keys:
            flux_et_sorted.append(self.flux_et[key])
            modis_et_sorted.append(self.mod16_et[key])

        slope, intercept, r_value, p_value, std_err = stats.linregress(flux_et_sorted, modis_et_sorted)
        modis_et_sorted = np.array(modis_et_sorted)
        flux_et_sorted = np.array(flux_et_sorted)

        plt.style.use('seaborn-darkgrid')

        #scatter plot of points
        plt.scatter(flux_et_sorted, modis_et_sorted, color = '#31a354')
        #error bars
        plt.errorbar(flux_et_sorted, modis_et_sorted, 
                     #error equation: error is a function of 2008's error, flux tower measurement
                     xerr = flux_et_sorted * self.year_dict['2008']['error'], 
                     linewidth = 0, elinewidth = 1, ecolor = 'black',
                     capsize = 2, zorder = -100)

        #extend linreg line
        flux_et_sorted = np.append(flux_et_sorted, 0)
        flux_et_sorted = np.append(flux_et_sorted, 40)

        plt.xlim(-1, 30)
        plt.ylim(0, 28)
        
        #lin regression line
        line = slope * flux_et_sorted + intercept
        plt.plot(flux_et_sorted, line, 'r', label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(slope,intercept, r_value), color = '#377eb8')

        plt.xlabel('Flux tower ET measurement in mm')
        plt.ylabel('MOD16 estimated ET in mm')
        plt.title('Comparison of 8-Day ET at CA-Qcu Flux Tower')

        plt.legend()
        plt.show()
        plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\all_of_2008.png', dpi = 350)
        
    '''
    this function takes a string that represents a year as an input
    it returns a linear regression graph of LE + H vs. Rn - H
    '''
    def yearly_uncertainty_plot(self, year):
        self.year = year
        
        #x values, Rn - G
        rn_g = []
        #y values, LE + H
        le_h = []
        
        for i in range(5, len(self.content)):
            entries = self.content[i].split(',')
            if self.year == entries[0][0:4]:
                #variables that we care about for determining error
                h = float(entries[9]) #sensible heat turbulent flux, w/m^2
                g = float(entries[13]) #soil heat flux, w/m^2
                rn = float(entries[24]) #net radiation, w/m^2
                le = float(entries[11]) #latent heat turbulent flux
                #this may be useful for other sites
                #s = float(entries[8]) #heat storage flux in biomass
                #checks if these variables were measured and determines the error
                if h != -9999 and g != -9999 and rn != -9999 and le != -9999:
                    #according to the energy balance closure error equation, Rn - G should equal LE + H
                    left = rn - g
                    right = le + h
                    rn_g.append(left)
                    le_h.append(right)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(rn_g, le_h)
        
        rn_g = np.array(rn_g)
        le_h = np.array(le_h)
        
        #scatter plot of data
        plt.scatter(rn_g, le_h)
        line = slope * rn_g + intercept
        
        plt.style.use('seaborn-pastel')
        #plots linear regression
        plt.plot(rn_g, line, 'r', label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(slope,intercept, r_value), color = '#377eb8')
                
        plt.xlabel('Rn - G, W/m$^2$')
        plt.ylabel('LE + H, W/m$^2$')
        plt.title('Energy balance closure error at {} for {}'.format(self.tower, self.year))
        
        
        plt.legend()
        plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\{}_closure_error_{}.png'.format(self.tower, self.year), dpi = 350)

                
        
        
test_file = r"C:\Users\Devin Simmons\Desktop\GEOL393\flux_towers\figure_for_1st_pres\AMF_CA-Qcu_BASE_HH_1-1.csv"

test_flux = flux_reader(test_file)
test_flux.et_by_day()
in_flux = test_flux.valid_8day_by_year(2008)


in_towers = r'C:\Users\Devin Simmons\Desktop\GEOL393\GIS\MOD16_2014_ET_Annual\flux_towers\fluxnet_sites\fluxnet_sites.shp'
in_rasters = r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\mod16_tiles\ca_qcu'
out_folder = r"C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\clip_features"

test1 = modis_et(in_towers, "CA-Qcu", out_folder)
test1.make_buffer(1692)

#without these lines the program had trouble reading the HDF files
arcpy.env.workspace = in_rasters
rasters = arcpy.ListRasters("*", "HDF")

mod16_values = {}
for raster in rasters:
    test1.clip_raster(raster)
    mod16_values[test1.julian_day] = test1.mean_et


test_flux.yearly_error()
test_flux.linreg_et_plot(in_flux, mod16_values)


        
        
        
                                                
                                        
        
        

    
        


