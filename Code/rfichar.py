####
##  (1) rfi characteristics stored directly in a dataframe. 
##  (2) for each rfi type, get fr

import pandas as pd

"""
Structure of the RFI source dictionary
--------------------------------------
types : Type of RFI signal
freqrange : Frequency ranges (units:GHz) 
freqres : Channel width of the signal  (units:kHz)
timeres : Time resolution of the signals. 
timegap : Time space between samples. 0.0 for 'continuous'
timefrac : Fraction of time (0-1) it is expected to be on. ( A list of size N )
arrayfrac : Fraction of the array that sees this RFI, or range of RFI signal ( A list of size N )
             - outlier:10-30km, few outlier ants 
             - core: 50km, 90 antennas
             - full: 500km range (all except the outliers)
"""


rfi = {
    'People': {'types':['wifi','bluetooth','wireless usb','wireless video','car radar'], 
               'freqrange':[[0.9,0.93],[2.4,2.5],[3.6,3.7],[4.9,5.0],[5.1,5.9],[61,61.5],[3.0,11.0],[5.7,5.9],[76.0,81.0]] , ## GHz
               'freqres': 20000.0,  # kHz
               'timeres': 1e-03 , ## millisecond signal packets
               'timegap': 0.01, ## Gaps of 100msec (Nearly continuous when on.)
#               'timefrac':[0.4] , ## Always, someone radiating near some antennas.
#               'arrayfrac' : ['core'] ## 'outlier', 'core', 'full'   
               'timefrac':[0.9,0.4] , ## Always, someone radiating near some antennas.
               'arrayfrac' : ['outlier','core'] ## 'outlier', 'core', 'full'   
           },

    'UWB': {'types':['wifi','bluetooth','wireless usb','wireless video','car radar','surveillance','handheld comm'], 
               'freqrange':[[1.99,10.6],[7.5,8.5],[22.0,29.0]] , ## GHz
               'freqres': 900.0*1e+3,  # 900 MHz  ( also need 8 GHz at 24 GHz for vehichles )
               'timeres': 1e-09, # 
               'timegap': 1e-07, ## 10 MHz pulse rate
               'timefrac':[0.9,0.2] , ## Always, someone radiating near some antennas.
               'arrayfrac' : ['outlier','core'] ## 'outlier', 'core', 'full'   
           },

    'UWBcar': {'types':['car radar'], 
               'freqrange':[[22.0,29.0]] , ## GHz
               'freqres': 8000.0*1e+3,  # 8 GHz at 24 GHz for vehicles
               'timeres': 1e-09, # 
               'timegap': 1e-07, ## 10 MHz pulse rate
               'timefrac':[0.2,0.03] , ## Cars near outlier antennas, and passing on Hwy 60.
               'arrayfrac' : ['outlier','core'] ## 'outlier', 'core', 'full'   
           },

    'Cell 5G': {'types':['cell phones','cars on hwy60', 'aircraft 5G'],
                'freqrange':[[1.427,1.518],[3.3,3.8],[5.15,5.925],[24.25,27.5],[31.8,33.4],[37.0,43.5],[45.5,50.2],[50.4,52.6],[66.0,76.0],[81.0,86.0]] , ## GHz
                'freqres': 200.0,  # kHz
                'timeres': 1e-04, #millisec signal packets
                'timegap': 0.0, ## continuous when on
                'timefrac':[1.0, 0.2, 0.2] , ## Fraction of time per day, from 0 to 1
                'arrayfrac' : ['outlier','core', 'full'] ## people, cars on hwy 60, 5G illuminated aircraft
            },
    
    'LEO Sat': {'types':['spacex','oneweb','boeing','Iridium Data', 'sirius XM'],
                       'freqrange':[[1.61,1.63],[2.2,2.33],[10.7,12.7],[13.85,14.5],[17.8,18.6],[18.8,19.3],[27.5,29.1], [29.5,30.0],[17.8,19.3],[27.5,29.1],[29.5,30.0],[37.5,42.0],[47.2,50.2],[50.4,51.4]] , ## GHz
                       'freqres': 200.0,  # kHz
                       'timefrac':[0.9] , ## Fraction of time per day, from 0 to 1
                       'timeres': 1e-04 , ## Time res at which coherence is preserved## Continuous when on
                       'timegap': 0.0 ,  ## continuous when on
                       'arrayfrac' : ['full'] ## 'outlier', 'core', 'full'   
                   },
    
    'Aircraft Comm': {'types':['ground radar','aircraft comm','aircraft data','aircraft nav','aircraft radar','drones'],
                      'freqrange':[[1.24,1.37],[2.7,3.0],[5.6,5.65],[9.0,9.2],[15.7,16.2],[0.978,1.213],[4.2,4.4],[5.3,5.5],[8.7,8.9],[9.3,9.5],[8.0,12.0],[12.0,18.0],[13.2,13.4]] , ## GHz
                      'freqres': 100.0,  # kHz
                      'timefrac':[0.5,0.9] , ## Fraction of time per day, from 0 to 1 that several of these may be present
                      'timeres': 2e-05 , ## 10s of micro-seconds
                      'timegap': 5e-05 , ## one and off protocols
                      'arrayfrac' : ['core','outlier'] ## 'outlier', 'core', 'full'   
                  },
    
    'Sat Comm': {'types':['GNSS','SAR','Iridium Comm','Military','CloudSat'],
                       'freqrange':[[1.1,1.8],[2.0,2.3],[2.5,2.7],[3.4,4.2],[7.25,7.75],[8.0,9.0],[10.7,12.8],[18.0,20.0], [23.0,27.0], [94.01,94.09] ] , ## GHz
                       'freqres': 100.0,  # kHz
                       'timefrac':[0.4] , ## Fraction of time per day, from 0 to 1
                       'timeres': 2e-05 , ## 
                       'timegap': 1e-04 , 
                       'arrayfrac' : ['full'] ## 'outlier', 'core', 'full'   
                   }
    
}


### LEO sats
###  One hemisphere is 180 * 180 = 32400 square degrees.
### If there are 2000 satellites covering one hemisphere,
### Spacing between satellites is  32400/2000 = 16.2 square degrees.
### => every 4 x 4 degrees there is one satellite.
### ====> NEED a hole above the ngVLA.
###  Assume that even with a hole, need 20deg to not saturate receivers.
###  ==> assume max FOV at which RFI will enter for LEOs is 20 deg.

### About 40 satellites visible at any given time. 

###  USA : 6% of world area.
###  LEO footprint of one satellite is between 2-4% of world area.
### 100 min orbits => 3.6arcmin/sec speed.
###                          => Move through beam in a few seconds. 
### Multiple satellites visible at once, at 4-deg separation.
###
###  Area of USA : 3,600,000 sq miles
###  Area of planet : 4*3.14*3,950^2  
### Fraction of area covered by USA :  = 1.84 %  
### Assume 2% coverage of the ngVLA. 
###  ==> Satellites will see the entire array at once. 
####




def get_rfi_df(rfi=None):
    rfistr = {}
    for akey in rfi.keys():
        rfistr[akey] = rfi[akey]
        rfistr[akey]['name'] = akey
        for bkey in rfistr[akey].keys():
            rfistr[akey][bkey] = str( rfistr[akey][bkey] )
    return (pd.DataFrame(data=rfistr)).transpose()

def get_rfi_dict(indf=None):
    tdict = indf.transpose().to_dict()
    newdict = {}
    for akey in tdict.keys():
        keyname = tdict[akey]['name']
        newdict[keyname] = tdict[akey]
        for bkey in newdict[keyname].keys():
            if bkey != 'name':
                newdict[keyname][bkey] = eval(newdict[keyname][bkey] )
    return newdict

def old_get_rfi_df():
    newrfi = {}
    for akey in rfi.keys():
        onesource = rfi[akey]
        nent = len( onesource['arrayfrac'] )
        asrc = onesource.copy()
        for ii in range(0,nent):
            asrc['arrayfrac'] = [ onesource['arrayfrac'][ii] ] 
            asrc['timefrac'] = [ onesource['timefrac'][ii] ] 
            newrfi[akey+str(ii)] = asrc.copy()
            newrfi[akey+str(ii)].pop('freqrange')
            newrfi[akey+str(ii)].pop('types')
            newrfi[akey+str(ii)]['name'] = akey
            newrfi[akey+str(ii)]['calc'] = "On"

    df = (pd.DataFrame(data=newrfi)).transpose()

    qqs = df.to_html(sparsify=False)
    af = open('out_rfichar.html','w')
    af.writelines(qqs)
    af.close()

    return df
    

def printrfi(outname='rfilist.tex'):
    
    fp = open(outname,'w')

    fp.write("\\begin{tabular}{|p{3cm}|p{4cm}|c|c|p{1.8cm}|p{1.8cm}|}\n")
    fp.write("\\hline\n")
    fp.write("Type & Freq Range & Freq Res & Time Res & Time Frac & Array Frac\\\\\n")
    fp.write("\\hline\n")

    for atom in rfi.keys():
        
        fp.write(atom + "\\newline\\newline " + str(rfi[atom]['types']) + " & " + str(rfi[atom]['freqrange']) + " & " + str(rfi[atom]['freqres']) + " kHz & " + str(rfi[atom]['timeres']) + " s & " + splitlines(rfi[atom]['timefrac']) + " & " + splitlines(rfi[atom]['arrayfrac']) + "\\\\\n")
        fp.write("\\hline\n")

    #fp.writeline(asrc)
    fp.write("\\end{tabular}\n")
    fp.close()


def splitlines(astr):
    newstr = ''
    for aw in astr:
        newstr = newstr + str(aw) + " \\newline "
    return newstr


