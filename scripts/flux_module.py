from bs4 import BeautifulSoup
import requests
import arcpy 
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from datetime import datetime

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
    ebr closure is evaluated on a yearly scale
    '''
    def ebr_error(self):
        self.year_dict = {}
        for i in range(5, len(self.content)):
            entries = self.content[i].split(',')    
            #slice year from datetime start column
            year = entries[0][0:4]
            
            if year not in self.year_dict:
                self.year_dict[year] = {}
                self.year_dict[year]['left error'] = 0
                self.year_dict[year]['right error'] = 0
                #the number of half hour periods used to determine error
                self.year_dict[year]['n'] = 0
            
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
                self.year_dict[year]['n'] += 1
        
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
    this function also determines error based on gap-filling and random error in LE measurements
    for each day
    '''
    def et_by_day(self):
        self.et_dict = {}
        
        for i in range(100, len(self.content) - 100):
            entries = self.content[i].split(',')
    
            #slices timestamp to just the date  
            date = entries[0][0:8]
            
            #turns date into a datetime object so that it can be converted to julian
            dt = datetime.strptime(date, '%Y%m%d')
            #converts to julian day which is a lot better than yyyymmdd
            date = datetime.strftime(dt, '%Y%j')
            
            #adds each date to the et dictionary
            if date not in self.et_dict:
                self.et_dict[date] = {}
                self.et_dict[date]['et'] = 0
                #measurements that weren't gap filled
                self.et_dict[date]['real measurements'] = []
                #filled and real measurements combined
                self.et_dict[date]['measurements'] = 0
                #a list that will contain the difference between et measurements
                #on the given day and et measurements taken at the same time on the
                #previous and/or preceding day
                self.et_dict[date]['paired et'] = []
                #will contain random error, in mm, for et measurements
                self.et_dict[date]['random error'] = 0
                
                    
            #variables needed to determine evapotranspiration
            le = float(entries[11]) #latent heat turbulent flux, w/m^2
            
            ta = float(entries[3]) #air temp, degrees c
            ws = float(entries[6]) #wind speed, m/s
            ppfd = float(entries[25]) - float(entries[28]) #difference between incoming and outgoing ppfd, µmolPhoton m-2 s-1
            
            #this handles the issue at some towers where every single midnight value
            #is missing
            if le == -9999 and str(entries[0][-4:]) == '0000':
                #set et equal to zero
                et = 0
                #cumulates the et measurements for each day and adds a value
                #to show how many measurements were taken
                self.et_dict[date]['et'] += et
                self.et_dict[date]['real measurements'].append(et)
                self.et_dict[date]['measurements']+= 1
            
            #determines if there are missing values for LE or temperature. if not, the
            #evapotranspiration is calculated

            elif le != -9999:
                #evapotranspiration
                et = le * 60 * 30 / ((2.501) * 10 ** 6)
                #cumulates the et measurements for each day and adds a value
                #to show how many measurements were taken
                self.et_dict[date]['et'] += et
                self.et_dict[date]['real measurements'].append(et)
                self.et_dict[date]['measurements']+= 1
                
                #these lines will determine the difference between le and le on the two adjdacent days
                #this is used to determine random error
                
                #previous day
                entries_neg_1 = self.content[i - 48].split(',')
                #the next day
                entries1 = self.content[i + 48].split(',')
                
                adjacent_days = [entries_neg_1, entries1]
                
                for day in adjacent_days:
                    #these values will be compared to le, ta etc. 
                    le2 = float(day[11])
                    ta2 = float(day[3]) 
                    ws2 = float(day[6]) 
                    ppfd2 = float(day[25]) - float(day[28])
                    #check if comparison of le can be made
                    #the environmental conditions need to be nearly identical
                    if le2 != -9999 and (abs(ta - ta2) < 3) and (abs(ppfd - ppfd2) < 75) and (abs(ws - ws2) < 1):
                        et2 = le2 * 60 * 30 / ((2.501) * 10 ** 6)
                        #adds the difference between et, et on adjacent day to this list
                        self.et_dict[date]['paired et'].append(et - et2)
        
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
                     #ta = float(day[3])
                     #commented out. for now I'm looking into using a latent heat
                     #of vaporization constant
                     #if le != -9999 and ta != -9999:
                     if le != -9999:
                         #adjacent_et.append(float(i[11]) * 60 * 30 / ((2.501 - 0.002361 * float(i[3])) * 10 ** 6))
                         adjacent_et.append(le * 60 * 30 / ((2.501) * 10 ** 6))
                #we can only use this cocrrection if we get two values
                if len(adjacent_et) >= 2:
                    evapotranspiration = sum(adjacent_et) / len(adjacent_et)
                    #cumulates the et measurements for each day and adds a value
                    #to show how many measurements were taken
                    self.et_dict[date]['et'] += evapotranspiration
                    self.et_dict[date]['measurements']+= 1
                
                #adding this keeps the real measurements list the same length for all days
                #which makes it easier to calculate gap filling error
                self.et_dict[date]['real measurements'].append('null')
                
        #this goes through each day, determines the standard deviation of
        #paired et measurements, and uses this value to determine random error
        #for that day
        for day in self.et_dict:
            if len(self.et_dict[day]['paired et']) > 1:
                #calculate random error
                self.et_dict[day]['random error'] = 2 * float(statistics.stdev(self.et_dict[day]['paired et']) / (2**0.5))
                
        
        #all the days in the study period, sorted
        study_days = [day for day in self.et_dict]
        study_days.sort()
        #this dictionary will contain keys that correspond to the percentage
        #of gaps filled, while the values will be lists of the corresponding error
        #for each day that was tested
        self.gap_error = {12.5: [], 25: [], 50:[], 75:[], 87.5:[]}
        
        #goes through days with all real measurements, creates artifical, gap-fillled 
        #versions, determines a relationship between % gaps filled and % deviation from
        #true value
        for i in range(0, len(study_days)):
            #list of real measurements for day in question
            measurements = self.et_dict[study_days[i]]['real measurements']
            #does not do this on days that are at the edge of the study period, as they don't have 
            #a four day window to compare with
            if "null" not in measurements and i not in [0, 1, len(study_days), len(study_days)-1, len(study_days) - 2]:
                #these will contains artifical data sets where 12.5%, 25% etc.
                #of the gaps are filled in from the surrounding days
                
                i125 = []
                i25 = []
                i50 = []
                i75 = []
                i875 = []
                
                #represents the measurements taken in a 4-day window surrounding the 
                #artificial dataset
                measurements_neg_2 = self.et_dict[study_days[i - 2]]['real measurements']
                measurements_neg_1 = self.et_dict[study_days[i - 1]]['real measurements']
                measurements_1 = self.et_dict[study_days[i + 1]]['real measurements']
                measurements_2 = self.et_dict[study_days[i + 2]]['real measurements']
                
                window = [measurements_neg_2, measurements_neg_1, measurements_1, 
                          measurements_2]
                
                #there is a pesky day in one of the files that is somehow not getting an array
                #of length 48 for 'real measurements'
                if sum([len(i) for i in window]) != (48*4):
                    pass
                
                else:
                    #goes through each half hour measurement in order
                    for x in range(0, len(measurements)):
                        #the four values of ET measured at the same time on the adjacent days
                        #only adds values that aren't "null"
                        
                        adjacent_days = [day[x] for day in window if type(day[x]) != str]
                        
                        if len(adjacent_days) < 2:
                            #the script will not use any days that have a null fill value to 
                            #determine error
                            window_avg = 'null'
                        else:
                            #the value that will be used to fill artificial gaps
                            window_avg = sum(adjacent_days)/len(adjacent_days)
                        
                        #fills the lists with a certain percentage of evenly spaced gap-fills                  
                        if x % 8 == 0:
                            i125.append(window_avg)
                            i875.append(measurements[x])
                        else:
                            i875.append(window_avg)
                            i125.append(measurements[x])
                        
                        if x % 4 == 0:
                            i75.append(measurements[x])
                            i25.append(window_avg)
                        else:
                            i25.append(measurements[x])
                            i75.append(window_avg)
                            
                        if x % 2 == 0:
                            i50.append(measurements[x])
                        else:
                            i50.append(window_avg)
                    
                    #determine what percent error there is on the artificial data 
                    #when compared to the real data. this is a real mess and should be fixed
                    if "null" not in i125:
                        #if abs(1 - sum(i125)/sum(measurements)) < 10: 
                            self.gap_error[12.5].append(abs(1 - sum(i125)/sum(measurements)))
                    if "null" not in i25:
                        #if abs(1 - sum(i25)/sum(measurements)) < 20: 
                            self.gap_error[25].append(abs(1 - sum(i25)/sum(measurements)))
                    if "null" not in i50:
                        #if abs(1 - sum(i50)/sum(measurements)) < 40: 
                            self.gap_error[50].append(abs(1 - sum(i50)/sum(measurements)))
                    if "null" not in i75:
                        #if abs(1 - sum(i75)/sum(measurements)) < 50: 
                            self.gap_error[75].append(abs(1 - sum(i75)/sum(measurements)))
                    if "null" not in i875:
                        #if abs(1 - sum(i875)/sum(measurements)) < 50: 
                            self.gap_error[87.5].append(abs(1 - sum(i875)/sum(measurements)))
        
        #will contain the x/y data for percent gaps, percent error
        percent_gaps = []
        errors = []
        for i in [.125, .25, .50, .75, .875]:
            for x in self.gap_error[i*100]:
                values = np.array(self.gap_error[i*100])
                #excludes the 5% most extreme values from each list
                if x < np.percentile(values, 95):
                    percent_gaps.append(i)
                    errors.append(x)
                
        self.gap_slope, self.gap_intercept, r_value, p_value, std_err = scipy.stats.linregress(percent_gaps, errors)
            
        percent_gaps = np.array(percent_gaps)
        errors = np.array(errors)
        
        percent_gaps *= 100
        errors *= 100
        
        #commented out because this plot isn't terribly useful
        
        #scatter plot of points
        # plt.scatter(percent_gaps, errors, color = '#e41a1c',
        #             zorder = 4, edgecolors = 'black')
        # #linreg relationship between percent gaps, percent error
        # line = self.gap_slope * percent_gaps + self.gap_intercept
        
        # plt.title("Relationship between % Gaps Filled, % Error at {}".format(self.tower))
        # plt.xlabel("Percentage of gaps created and filled in artificial datasets")
        # plt.ylabel("Percent error")
        # plt.plot(percent_gaps, line, 'r', zorder = 5, 
        #          label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(self.gap_slope, self.gap_intercept, r_value ** 2), 
        #          color = '#fb9a99', linewidth = 2.5)
        # plt.legend()
        
        # plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\gap_error\{}_gap_error.png'.format(self.tower), dpi = 350)
        # plt.show()
        
        
        #goes through et_dict and applies equation to determine percent error on
        #gap-filled days
        for day in self.et_dict:
            #number of real measurements
            real_measurements = len([i for i in self.et_dict[day]['real measurements'] if i != 'null'])
            #gap error in mm represents the range of uncertainty
            if real_measurements < 48:
                self.et_dict[day]['gap error in mm'] = (self.gap_slope * (1 - real_measurements/48) + self.gap_intercept) * self.et_dict[day]['et']
            else:
                self.et_dict[day]['gap error in mm'] = 0
            
        return self.et_dict
        
    
    '''
    this function can only be called after et_by_day has been called
    it returns a dictionary that displays all the valid 8 day ranges recorded
    by the flux tower for a given time range. year is a list containing strings 
    or ints
    '''
    def valid_8day_by_year(self, years):
             
        self.years = years
        #will contain the julian day an 8-day period starts, averaged cumulative
        #ET over that period, gap-filling error, and random error. Vqlues are tuples. ET is the first
        #tuple value
        self.flux_et = {}
        
        #contains julian day 8-day period starts, number of days the 8-day average 
        #represents
        self.measurement_dict = {}
        
        for year in years:
            #julian day counter, step value of 8 corresponds to 8 day period
            for i in range (1, 365, 8):
                
                et_counter = 0
                days_counter = 0
                error_counter = 0
                
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
                    if int(str(day)[-3:]) < 361:
                        #convert to julian day
                        jd = datetime.strptime(str(day), '%Y%j')
                        jd = datetime.strftime(jd, '%Y%j')
                        
                        if jd in self.et_dict:
                            #a day needs all 48 measurements filled to be valid
                            if self.et_dict[jd]['measurements'] == 48:
                                et_counter += self.et_dict[jd]['et']
                                days_counter += 1
                                #dependent errors, propagated in quadrature
                                error_counter += (self.et_dict[jd]['gap error in mm']**2 + self.et_dict[jd]['random error']**2)**0.5
                #need at least five days with all 48 measurements filled in
                if days_counter > 4:
                    et_averaged = (et_counter * 8) / days_counter
                    error_percent = error_counter/et_averaged
                    self.flux_et[i] = et_averaged, error_percent
            
                self.measurement_dict[i] = days_counter
        return self.flux_et
    
    '''
    makes plots that compare et data to another dataset (ie MOD16 data)
    date_list is a list of ints that represent the start julian day for the
    8 day comparison period. can only be called after et_by_day.
    both flux_et and mod16_et should be dictionaries where the key is
    the julian day (YYYYDDD) and the value is ET. 
    years should be a list of years
    '''
    def linreg_et_plot(self, mod16_et):

        self.mod16_et = mod16_et

        flux_keys = [key for key in self.flux_et]
        flux_keys.sort()
        
        
        flux_et_sorted = []
        modis_et_sorted = []
        
        #error is calculated using gap-filled error percentages 
        error_sorted = [abs(self.flux_et[i][1]) for i in flux_keys]

        for key in flux_keys:
            flux_et_sorted.append(self.flux_et[key][0])
            modis_et_sorted.append(self.mod16_et[key])

        slope, intercept, r_value, self.p_value_8day, std_err = scipy.stats.linregress(flux_et_sorted, modis_et_sorted)
        
        #making these accessible outside of the function so that
        #points can be aggregated for east and west
        self.modis_et_sorted = np.array(modis_et_sorted)
        self.flux_et_sorted = np.array(flux_et_sorted)
        self.error_sorted = np.array(error_sorted)

        plt.style.use('seaborn-darkgrid')

        #scatter plot of points
        plt.scatter(self.flux_et_sorted, self.modis_et_sorted, color = '#008837',
                    zorder = 4, edgecolors = 'black')
                    
        #error bars
        plt.errorbar(self.flux_et_sorted, self.modis_et_sorted, 
                     #error equation: error is a function gap filling error, random error, flux et
                     xerr = flux_et_sorted * self.error_sorted,
                     linewidth = 0, elinewidth = 1, ecolor = 'black',
                     capsize = 2, zorder = 3)
        
        #lin regression line
        line = slope * self.flux_et_sorted + intercept
        plt.plot(self.flux_et_sorted, line, 'r', zorder = 5, 
                 label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(slope, intercept, r_value ** 2), 
                 color = '#7b3294', linewidth = 2.5)

        
        plt.xlabel('Flux tower ET measurement in mm')
        plt.ylabel('MOD16 estimated ET in mm')
        plt.title('Comparison of 8-Day ET at {} Flux Tower'.format(self.tower))

        plt.legend()
        plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\linreg\{}_linreg.png'.format(self.tower), dpi = 350)
        plt.show()
    '''
    groups the assessesd 8-day periods by year and then does linear regression.
    can only be called after et_by_day.
    both flux_et and mod16_et should be dictionaries where the key is
    the julian day (YYYYDDD) and the value is ET. 
    years should be a list of years.
    makes a linreg of 8 day modis et compared to flux tower et for a single year
    '''
    def annual_et_linreg(self, mod16_et, year, make_plot = False):
        self.year = year
        self.mod16_et = mod16_et
        
        
        flux_keys = [key for key in self.flux_et]
        flux_keys.sort()
        
        flux_et = []
        mod_et = []
        
        error = []
        #check 8 day periods for given year
        for key in flux_keys:
            if key[0:4] == str(self.year):

                flux_et.append(self.flux_et[key][0])
                mod_et.append(self.mod16_et[key])
                #error == et * percent error
                error.append(self.flux_et[key][0] * abs(self.flux_et[key][1]))

        #stats
        #year_slope and year_intercept are declared so that they can be accessed and 
        #plotted later
        year_slope, year_intercept, r_value, p_value, std_err = scipy.stats.linregress(flux_et, mod_et)
        
        #turn lists to np arrays
        flux_et = np.array(flux_et)
        mod_et = np.array(mod_et)
        error = np.array(error)
        
        if make_plot == True:
            plt.style.use('seaborn-darkgrid')
            
            #scatter plot of points
            plt.scatter(flux_et, mod_et, color = '#008837',
                        zorder = 4, edgecolors = 'black')
    
            #error bars
            plt.errorbar(flux_et, mod_et, 
                         xerr = error,
                         linewidth = 0, elinewidth = 1, ecolor = 'black',
                         capsize = 2, zorder = 3)
        
        
            
            #lin regression line
            line = year_slope * flux_et + year_intercept
            plt.plot(flux_et, line, 'r', zorder = 5, 
                     label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(year_slope, year_intercept, r_value ** 2), 
                     color = '#7b3294', linewidth = 2.5)
                     
            plt.xlabel('Flux tower ET measurement in mm')
            plt.ylabel('MOD16 estimated ET in mm')
            plt.title('Comparison of {} ET at {} Flux Tower'.format(self.year, self.tower))
    
            plt.legend()
            plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\{}_linreg_{}.png'.format(self.tower, self.year), dpi = 350)
            plt.show()
        return year_slope, year_intercept
        
    '''
    make a time series plot of flux data vs. modis data
    year should just be one year (may change this)
    '''
    def time_series_plot(self, flux_et, mod16_et, ppt_data, year):
        self.flux_et = flux_et
        self.mod16_et = mod16_et
        self.ppt_data = ppt_data
        
        #sorts the data in chronological order
        flux_keys = [key for key in self.flux_et]
        flux_keys.sort()
        mod_keys = [key for key in self.mod16_et]
        mod_keys.sort()
        
        if ppt_data is not None:
            ppt_keys = [key for key in self.ppt_data]
        
        flux_et_sorted = np.array([])
        modis_et_sorted = np.array([])
        
        if ppt_data is not None:
            ppt_sorted = np.array([])
        
        for key in mod_keys:
            if key[0:4] == str(year):
                modis_et_sorted = np.append(modis_et_sorted, self.mod16_et[key])
                if key in flux_keys:
                    flux_et_sorted = np.append(flux_et_sorted, self.flux_et[key][0])
                else:
                    flux_et_sorted = np.append(flux_et_sorted, np.nan)
                if key in ppt_keys:
                    ppt_sorted = np.append(ppt_sorted, self.ppt_data[key]['ppt in mm'])
                else:
                    ppt_sorted = np.append(ppt_sorted, np.nan)

        
        #gets rid of jd 361, which has no flux value associated
        #modis_et_sorted = modis_et_sorted[0:-1]
        plt.style.use('seaborn-darkgrid')
        
        julian_days = [i for i in range(1, 365, 8)]
        
        if self.ppt_data is not None:
            #ppt data graph
            plt.bar(julian_days, ppt_sorted, color = '#984ea3', zorder = 1, width = 1.6)
        
        #time series plots
        plt.plot(julian_days, flux_et_sorted, color = "#74a9cf", zorder = 2)
        plt.plot(julian_days, modis_et_sorted, color = "#045a8d", zorder = 2)
        plt.scatter(julian_days, flux_et_sorted, color = "#74a9cf", s = 10, zorder = 3)
        plt.scatter(julian_days, modis_et_sorted, color = "#045a8d", s = 10, zorder = 3)
                  
        #error is calculated as gap-filling error
        error_sorted = np.array([])
        for i in julian_days:
            jd_year = str(year) + str(i).zfill(3)
            if jd_year in flux_keys:
                #error is a function of gap filling error, random error
                error_sorted = np.append(error_sorted, abs(self.flux_et[jd_year][1]))
            else:
                #nan means that nothing will plot
                error_sorted = np.append(error_sorted, np.nan)
        
        #error bars
        plt.errorbar([i for i in range(1, 365, 8)], flux_et_sorted, 
                     #error equation: error is a function of 2008's error, flux tower measurement
                     yerr = flux_et_sorted * error_sorted, 
                     linewidth = 0, elinewidth = 1, ecolor = "black", capsize = 2,
                     zorder = 2)
        
        blue = matplotlib.patches.Patch(color = '#74a9cf', label = 'Flux tower ET')
        dark_blue = matplotlib.patches.Patch(color = '#045a8d', label = 'MOD16 estimations')
        purp = matplotlib.patches.Patch(color = '#984ea3', label = 'Precipitation')
        
        plt.xticks([i for i in range (1, 356, 24)], [i for i in range (1, 356, 24)])
        plt.xlabel("Julian day that 8-day period began")
        plt.ylabel("ET in mm for 8-day period")
        plt.title("{} Time Series of MOD16 and Flux Tower ET at {}".format(year, self.tower))
        
        plt.legend(handles = [blue, dark_blue, purp])
        
        plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\time_series\{}_{}_time_Series.png'.format(self.tower, year), dpi = 350)
    
    '''
    function that compares ppt to evapotranspiration measures via linear regression
    '''
    def ppt_et_linreg(self, mod16_et, ppt_data):
        
        #we're only considering 8 day periods that aren't missing precipitation data
        self.ppt_keys = [i for i in ppt_data if ppt_data[i]['missing days'] == 0]
        self.ppt_keys.sort()
        
        ppt_sorted = []
        flux_et_sorted = []
        modis_et_sorted = []
        
        #error is calculated using gap-filled error percentages 
        flux_error_sorted = [abs(self.flux_et[i][1]) for i in self.ppt_keys if i in self.flux_et]

        for key in self.ppt_keys:
            #8 day period needs valid data from 
            if key in self.flux_et:
                flux_et_sorted.append(self.flux_et[key][0])
                modis_et_sorted.append(mod16_et[key])
                ppt_sorted.append(ppt_data[key]['ppt in mm'])
        
        #added ppt as an addendum bc these arrays have less data since they are constricted
        #to 8 day periods during which ppt isn't missing data
        self.modis_et_sorted_ppt = np.array(modis_et_sorted)
        self.flux_et_sorted_ppt = np.array(flux_et_sorted)
        self.flux_error_sorted_ppt = np.array(flux_error_sorted)
        self.ppt_sorted = np.array(ppt_sorted)
        
        #data for linreg that compares flux et to ppt
        slope, intercept, r_value, p_value_8day, std_err = scipy.stats.linregress(flux_et_sorted, ppt_sorted)
        
        #making these accessible outside of the function so that
        #points can be aggregated for east and west
        

        plt.style.use('seaborn-darkgrid')

        #scatter plot of points
        plt.scatter(self.flux_et_sorted_ppt, self.ppt_sorted, color = '#74a9cf',
                    zorder = 4, edgecolors = 'black')
                    
        #error bars
        plt.errorbar(self.flux_et_sorted_ppt, self.ppt_sorted, 
                     #error equation: error is a function gap filling error, random error, flux et
                     xerr = self.flux_et_sorted_ppt * self.flux_error_sorted_ppt,
                     linewidth = 0, elinewidth = 1, ecolor = 'black',
                     capsize = 2, zorder = 3)
        
        #lin regression line
        line = slope * self.flux_et_sorted_ppt + intercept
        plt.plot(self.flux_et_sorted_ppt, line, 'r', zorder = 5, 
                 label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(slope, intercept, r_value ** 2), 
                 color = '#984ea3', linewidth = 2.5)

        
        plt.xlabel('Flux tower ET measurement in mm')
        plt.ylabel('Precipitation in mm')
        plt.title('Comparison of Flux Tower ET and Precipitation at {} Flux Tower'.format(self.tower))

        plt.legend()
        plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\linreg\{}_flux_ppt_linreg.png'.format(self.tower), dpi = 350)
        plt.show()
        
        #data for linreg that compares modis et to ppt
        slope, intercept, r_value, p_value_8day, std_err = scipy.stats.linregress(modis_et_sorted, ppt_sorted)
        #scatter plot of points
        plt.scatter(self.modis_et_sorted_ppt, self.ppt_sorted, color = '#045a8d',
                    zorder = 4, edgecolors = 'black')

        #lin regression line
        line = slope * self.modis_et_sorted_ppt + intercept
        plt.plot(self.modis_et_sorted_ppt, line, 'r', zorder = 5, 
                 label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(slope, intercept, r_value ** 2), 
                 color = '#984ea3', linewidth = 2.5)

        
        plt.xlabel('MOD16 ET estimation in mm')
        plt.ylabel('Precipitation in mm')
        plt.title('Comparison of MOD16 ET and Precipitation at {} Flux Tower'.format(self.tower))

        plt.legend()
        plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\linreg\{}_mod16_ppt_linreg.png'.format(self.tower), dpi = 350)
        plt.show()
                
    
    '''
    this function takes a string that represents a year as an input
    it returns a linear regression graph of LE + H vs. Rn - H
    plot is a boolean and should equal True or False
    '''
    def yearly_uncertainty_plot(self, year, plot):
        self.year = str(year)
        
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
        
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(rn_g, le_h)
        
        rn_g = np.array(rn_g)
        le_h = np.array(le_h)
        
        if plot == True:
            #scatter plot of data
            plt.scatter(rn_g, le_h)
            line = slope * rn_g + intercept
            
            plt.style.use('seaborn-pastel')
            #plots linear regression
            plt.plot(rn_g, line, 'r', label='y = {:.2f}x + {:.2f}, R$^2$ = {:.2f}'.format(slope,intercept, r_value ** 2), color = '#377eb8')
                    
            plt.xlabel('Rn - G, W/m$^2$')
            plt.ylabel('LE + H, W/m$^2$')
            plt.title('Energy balance closure error at {} for {}'.format(self.tower, self.year))
            
            
            plt.legend()
            plt.savefig(r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\figures\{}_closure_error_{}.png'.format(self.tower, self.year), dpi = 350)

        #returns the slope and intercept of OLS analysis for EBR error given a particular year,
        #this is used in the ebr_table function
        return [slope, intercept]
    '''
    returns the root mean square error of the data
    annual should be a boolean value.
    false, the default, calculates rmse for 8-day periods
    true calculates it on a yearly scale
    '''
    def rmse(self, mod16_et, seasonal = False):
        
        difference = 0
        counter = 0
        
        if seasonal == False:
            for key in self.flux_et:
                difference += (self.flux_et[key][0] - mod16_et[key])**2
                counter += 1
            self.root_mse = (difference/counter)**0.5
            return self.root_mse
        
        else: 
            
            #difference between mod16, flux et
            spring_dif = 0
            summer_dif = 0
            fall_dif = 0
            wint_dif = 0
            
            #counters
            spring_ct = 0
            summer_ct = 0
            fall_ct = 0
            wint_ct = 0
            
            for key in self.flux_et:
                #determines season based on julian day. special cases for leap years
                #seasons are meteorological seasons
                if (int(key[-3:]) in range(61, 153) and int(key[0:4]) % 4 == 0) or (int(key[-3:]) in range(60, 152) and int(key[0:4]) % 4 != 0):
                    spring_dif += (mod16_et[key] - self.flux_et[key][0]) ** 2
                    spring_ct += 1
                
                elif (int(key[-3:]) in range(153, 245) and int(key[0:4]) % 4 == 0) or (int(key[-3:]) in range(152, 244) and int(key[0:4]) % 4 != 0):
                    summer_dif += (mod16_et[key] - self.flux_et[key][0]) ** 2
                    summer_ct += 1
                
                elif (int(key[-3:]) in range(245, 336) and int(key[0:4]) % 4 == 0) or (int(key[-3:]) in range(244, 335) and int(key[0:4]) % 4 != 0):
                    fall_dif += (mod16_et[key] - self.flux_et[key][0]) ** 2
                    fall_ct += 1
                
                else: 
                    wint_dif += (mod16_et[key] - self.flux_et[key][0]) **2 
                    wint_ct += 1
            
            self.spring_rmse = (spring_dif/spring_ct) ** 0.5
            self.summer_rmse = (summer_dif/summer_ct) ** 0.5
            self.fall_rmse = (fall_dif/fall_ct) ** 0.5
            self.wint_rmse = (wint_dif/wint_ct) ** 0.5
            
            return self.spring_rmse, self.summer_rmse, self.fall_rmse, self.wint_rmse
        
    
    '''
    determines bias and percent bias of MOD16 relative to flux tower data
    seasonal should be a boolean value.
    false, the default, calculates bias for all 8-day periods
    true calculates it for all four seasons
    '''
    def bias(self, mod16_et, seasonal = False):

        #flux tower et, in mm
        fet = 0
        #difference between MOD16 and flux et 
        difference = 0
        #number of comparisons, or N 
        counter = 0
        
        if seasonal == False:
            for key in self.flux_et:
                #eq for bias
                difference += (mod16_et[key] - self.flux_et[key][0])
                counter += 1
                fet += self.flux_et[key][0]
            self.bias_calc = difference/counter
            #eq for pbias
            self.percent_bias = 100 * self.bias_calc/((1/counter)*fet)
            
            return self.bias_calc, self.percent_bias
            
        else:
            #flux tower et in mm
            spring_fet = 0
            summer_fet = 0
            fall_fet = 0
            wint_fet = 0
            
            #difference between mod16, flux et
            spring_dif = 0
            summer_dif = 0
            fall_dif = 0
            wint_dif = 0
            
            #counters
            spring_ct = 0
            summer_ct = 0
            fall_ct = 0
            wint_ct = 0
            
            
            for key in self.flux_et:
                #determines season based on julian day. special cases for leap years
                #seasons are meteorological seasons
                if (int(key[-3:]) in range(61, 153) and int(key[0:4]) % 4 == 0) or (int(key[-3:]) in range(60, 152) and int(key[0:4]) % 4 != 0):
                    spring_dif += (mod16_et[key] - self.flux_et[key][0])
                    spring_ct += 1
                    spring_fet += self.flux_et[key][0]
                
                elif (int(key[-3:]) in range(153, 245) and int(key[0:4]) % 4 == 0) or (int(key[-3:]) in range(152, 244) and int(key[0:4]) % 4 != 0):
                    summer_dif += (mod16_et[key] - self.flux_et[key][0])
                    summer_ct += 1
                    summer_fet += self.flux_et[key][0]
                
                elif (int(key[-3:]) in range(245, 336) and int(key[0:4]) % 4 == 0) or (int(key[-3:]) in range(244, 335) and int(key[0:4]) % 4 != 0):
                    fall_dif += (mod16_et[key] - self.flux_et[key][0])
                    fall_ct += 1
                    fall_fet += self.flux_et[key][0]
                
                else: 
                    wint_dif += (mod16_et[key] - self.flux_et[key][0])
                    wint_ct += 1
                    wint_fet += self.flux_et[key][0]
                
            self.spring_bias = spring_dif/spring_ct
            self.spring_pbias = 100 * self.spring_bias/((1/spring_ct) * spring_fet)
                
            self.summer_bias = summer_dif/summer_ct
            self.summer_pbias = 100 * self.summer_bias/((1/summer_ct) * summer_fet)
                
            self.wint_bias = wint_dif/wint_ct
            self.wint_pbias = 100 * self.wint_bias/((1/wint_ct) * wint_fet)
                
            self.fall_bias = fall_dif/fall_ct
            self.fall_pbias = 100 * self.fall_bias/((1/fall_ct) * fall_fet)

                
            return self.spring_bias, self.spring_pbias, self.summer_bias, self.summer_pbias, self.fall_bias, self.fall_pbias, self.wint_bias, self.wint_pbias
    
    '''
    writes a csv file that contains a summary table of MOD16 data, flux tower data.
    out_file should be a file path that ends in .csv
    i should edit this to automatically adjust the number of sig figs
    '''
    def summary_table(self, mod16_et, out_file):
        
        #all days that are compared
        mod_keys = [key for key in mod16_et]
        mod_keys.sort()
        
        fileobj = open(out_file, "w")
        fileobj.write("8-day period start date (Julian),MOD16 estimated ET in mm,Flux Tower Measured ET in mm,Measurement Uncertainty as a Percentage,Days of Flux Tower Measurement in 8 day period \n")
        for i in mod_keys:
            fileobj.write("{}, {},".format(i, mod16_et[i]))
            #writes nothing for days where (measurements == 48) < 5
            if i in self.flux_et:
                error = abs(self.flux_et[i][1])
                fileobj.write("{}, {},".format(self.flux_et[i][0], error))
            else:
                fileobj.write(",,")
            if i in self.measurement_dict:
                fileobj.write('{}\n'.format(self.measurement_dict[i]))
            else:
                fileobj.write('\n')
        fileobj.close()
    
    '''
    this makes a summary table of the ebr error for each year in the years arg 
    (which should be an array)
    ols_slope and ols_intercept can be acquired using the yearly_uncertainty_plot function
    they should be arrays
    '''
    def ebr_table(self, years, out_file, ols_slopes, ols_intercepts):
        
        fileobj = open(out_file, 'w')
        fileobj.write("Year, n, OLS Intercept, OLS SLope, Annual EBR\n")
        
        counter = 0
        for i in years:
            i = str(i)
            fileobj.write('{},{},{},{},{}\n'.format(i, self.year_dict[i]['n'],
                                                    ols_intercepts[counter], ols_slopes[counter],
                                                    self.year_dict[i]['right error']/self.year_dict[i]['left error']))
            counter += 1
    
        fileobj.close()
        
class ppt_reader:
    
    def __init__(self, infile):
        self.infile = infile
        fileobj = open(self.infile, 'r')
        self.content = fileobj.readlines()
        self.content = self.content[1:]
        fileobj.close()
    
    #this function creates a dictionary where ppt is summed up over 8 day 
    #periods
    @property
    def eight_day_ppt(self):
        
        counter = 0
        eight_day_ppt_dict = {}
        

        for i in self.content:     
            
            date = i.split(",")[0]
            
            #when ppt = '' there is a missing value
            ppt = i.split(",")[19]
            if ppt != '':
                ppt = float(ppt)

            
            #turns date into a datetime object so that it can be converted to julian
            dt = datetime.strptime(date, '%m/%d/%Y')
            #converts to julian day which is a lot better than yyyymmdd
            date = str(datetime.strftime(dt, '%Y%j'))

            #create new dictionary key for each 8 day period
            if counter % 8 == 0:
                
                eight_day_ppt_dict[date] = {'ppt in mm': 0, 'missing days': 0}
                
                if ppt != '':
                    eight_day_ppt_dict[date]['ppt in mm'] = ppt
                elif ppt == '':
                    eight_day_ppt_dict[date]['missing days'] += 1
            else:    
                
                if ppt != '':
                    #matching this ppt value to the correct key
                    eight_day_ppt_dict[str(int(date) - (counter % 8))]['ppt in mm'] += ppt
                elif ppt == '':
                    eight_day_ppt_dict[str(int(date) - (counter % 8))]['missing days'] += 1
             
            #the counter needs to be reset each year because 365 % 8 != 0
            #i am indexing the days as Jan 1 == 0 on the counter
            if counter >= 364:
                #checks for leap year
                if counter == 364 and int(date[0:4]) % 4 == 0:
                    counter += 1
                else: 
                    counter = 0
            elif counter < 364:
                counter += 1
        
        self.eight_day_ppt_dict = eight_day_ppt_dict
        return self.eight_day_ppt_dict
    
    @property
    def ordered_ppt(self):
        
        eight_day_periods = [i for i in self.eight_day_ppt_dict]
        eight_day_periods.sort()
        
        ppt_values = [self.eight_day_ppt_dict[i]['ppt in mm'] for i in eight_day_periods]
        
        return eight_day_periods, ppt_values

class download_modis:
    '''
    this class downloads modis data
    the __init__ function takes a list of years and a list of julian days
    julian days should be strings and should have three digits, with leading 
    zeros if jd < 100
    years can be int or str
    the function can be applied to one tile at a time
    it returns the tile folder where all the new rasters were saved
    '''
    def __init__(self, years, days, tile, out_folder):
        self.years = years
        self.days = days
        self.tile = tile
        self.out_folder = out_folder
    def download_mod16(self):
        #file type we want to download
        ext = 'hdf'
        base_url = r'http://files.ntsg.umt.edu'
        
        for year in self.years:
            for day in self.days:
                url = r'http://files.ntsg.umt.edu/data/NTSG_Products/MOD16/MOD16A2.105_MERRAGMAO/Y{}/D{}/'.format(str(year), day)
                page = requests.get(url).text
                #parses the url's html file to find links to the rasters
                soup = BeautifulSoup(page, 'html.parser')
                #makes a list of the urls for the hdf files
                files = [base_url + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]
                
                for file in files:
                    if file.split(".")[-4] == self.tile:
                        r = requests.get(file, allow_redirects = True)
                        #save file to outpath
                        tile_filepath = self.out_folder + r'\{}\{}'.format(self.tile, file.split("/")[-1])
                        open(tile_filepath, 'wb').write(r.content)
                        print("wrote to {}".format(tile_filepath))
    @property
    def mod16_folder(self):
        return self.out_folder + r'\{}'.format(self.tile)
    
    


#Here are lines that show how I execute the code


in_towers = r'C:\Users\Devin Simmons\Desktop\GEOL393\GIS\MOD16_2014_ET_Annual\flux_towers\fluxnet_sites\fluxnet_sites.shp'
in_rasters = r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\mod16_tiles\h12v04'
out_folder = r"C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\clip_features"

#downloading rasters for h13v04 from 2004 to 2010
ca_na1_years = [2003]
#creates a list of the julian days that correspond to MOD16 8day periods
ca_na1_days = [str(i).zfill(3) for i in range(1, 365, 8)]

#sets parameters needed to download the right files
ca_na1_rasters = download_modis(ca_na1_years, ca_na1_days, "h13v04",
                                r'C:\Users\Devin Simmons\Desktop\GEOL393\figure_for_rd\gis_things\mod16_tiles')
#actually downloads these files
ca_na1_rasters.download_mod16()

in_rasters = ca_na1_rasters.mod16_folder

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

