
file = r"figure_for_1st_pres/AMF_CA-Qcu_BASE_HH_1-1.csv"

fileobj = open(file, 'r')

content = fileobj.readlines()

fileobj.close()
et_dict = {}
#will contain each year and its associated error
year_dict = {}

for i in range(100, len(content) - 100):
    entries = content[i].split(',')
    
    year = entries[0][0:4]
    #slices timestamp to just the date  
    date = entries[0][0:8]
    
    #adds each date to the et dictionary
    if date not in et_dict:
        et_dict[date] = {}
        et_dict[date]['et'] = 0
        et_dict[date]['measurements'] = 0
        et_dict[date]['left error'] = 0
        et_dict[date]['right error'] = 0
    
    #adds each year to the year dictionary
    if year not in year_dict:
        year_dict[year] = {}
        year_dict[year]['left error'] = 0
        year_dict[year]['right error'] = 0
    
    h = float(entries[9]) #sensible heat turbulent flux, w/m^2
    g = float(entries[13]) #soil heat flux, w/m^2
    rn = float(entries[24]) #net radiation, w/m^2
    le = float(entries[11]) #latent heat turbulent flux
    #checks if these variables were measured and determines the error
    if h != -9999 and g != -9999 and rn != -9999 and le != -9999:
        left = rn - g
        right = le + h
        year_dict[year]['left error'] += left
        year_dict[year]['right error'] += right
            
    
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
        et_dict[date]['et'] += evapotranspiration
        et_dict[date]['measurements']+= 1
       
        
    #if there are missing values ET is determined based on values in adjacent rows
    else:
        adjacent_et = []
        
        #two days prior and previous day at the same time, respectively
        entries_neg_2 = content[i - 96].split(',')
        entries_neg_1 = content[i - 48].split(',')
        
        #next two days
        entries1 = content[i + 48].split(',')
        entries2 = content[i + 96].split(',')
        
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
            et_dict[date]['et'] += evapotranspiration
            et_dict[date]['measurements']+= 1
    
for i in year_dict:
    left = year_dict[i]['left error']
    right = year_dict[i]['right error']
    if left != 0 and right != 0:
        year_dict[i]['error'] = (abs(left - right)/right + abs(left - right)/left)/2
    
'''shows all days with 48 measurements
counter = 0
for i in et_dict:
    if et_dict[i]['measurements'] == 48:
        print("{} {}".format(i, et_dict[i]['et']))
        counter += 1
'''
days = [str(i) for i in range(20080601, 20080609)]
days1 = [str(i) for i in range(20080609, 20080617)]
days2 = [str(i) for i in range(20080703, 20080711)]
days3 = [str(i) for i in range(20080711, 20080719)]
days4 = [str(i) for i in range(20080524, 20080532)]
days5 = [str(i) for i in range(20080516, 20080524)]
days6 = [str(i) for i in range(20080804, 20080812)]
days7 = [str(i) for i in range(20080508, 20080516)]
days8 = [str(i) for i in range(20080820, 20080828)]
days9 = [str(i) for i in range(20080422, 20080430)]


all_days = [days, days1, days2, days3, days4, days5, days6, days7, days8, days9]


from datetime import datetime

for day in all_days:
    et_counter = 0
    measurement_counter = 0
    #these represent the left and right sides of the energy closure balance error equation
    
    #takes first day from the list of days, converts it to
    #a julian day. this matches the format 
    dt = datetime.strptime(day[0], '%Y%m%d')
    jd = datetime.strftime(dt, '%Y%j')
    print(day[0])
    print(jd)
    for i in day:
        if et_dict[i]['measurements'] == 48:
            et_counter += et_dict[i]['et']
            measurement_counter += 1
            
    et_averaged = (et_counter * 8) / measurement_counter
    print(et_averaged)
    print('\n')
print(year_dict)

#mod = [197, 228, 235, 237, 239, 240, 242, 242, 246, 253, 271, 275]
#print(sum(mod)/len(mod)/10)
    
'''tells how much et in a given year
counter = 0
for i in et_dict:
    if str(i[0:4]) == '2010':
        counter += et_dict[i]['et']
        print(et_dict[i]['et'])
print(counter)
'''