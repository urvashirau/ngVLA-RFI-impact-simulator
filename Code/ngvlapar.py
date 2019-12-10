import numpy as np
import pylab as pl
import matplotlib
import pandas as pd

ngvlabands = [[1.2,3.9],[3.9,12.6],[12.6,21.0],[21.0,35.0],[30.5,50.5],[70.0,116.0]]

arrtype = 'ngvla'

if arrtype=='basic':
    ## Basic Array
    arrayfracinfo={
        'core': {'nant':100.0, 'maxbase':2.0},  # maxbase in km
        'outlier':{'nant':10, 'maxbase':50.0},
        'full':{'nant':214, 'maxbase':500.0}
    }
    totalants = 214.0
    totalvis = totalants * ( totalants - 1)/2.0
    maxbase = 500.0 # km
    ## Baseline count profile
    base_lengths = np.arange(0.02,maxbase,0.1)   # 1-500 km baselines
    base_counts = base_lengths.copy()
    
    ## Relative numbers of baselines for each baseline length range.
    base_counts.fill(0.3)
    base_counts[ base_lengths < 10 ] = 1.0
    base_counts[ base_lengths >=250  ] = 0.1
    
if arrtype=='ngvla' or arrtype=='approx_ngvla':
    
    ## ngVLA Definition
    ##   94  x  18m  :  Main : Core  ( 1.25 km baselines )
    ##   74 x 18m  : Main : Plains Spiral  ( 70 km baselines )
    ##   46 x 18m : Main : Mid Baselines (upto 1000km )
    ##   19 x 6m  : Compact Core ( 0.1 km baselines )
    ##   233 in total so far.
    ##   30 x 18m  : Long Baselines  ( 1000km baselines, in 10 clusters/stations )
    arrayfracinfo={
        'core': {'nant':94+19+30, 'maxbase':5.0},   ## 143 Core and compact core and 20 plains ants
        'outlier':{'nant':30, 'maxbase':1000.0},    ## baselines that all see the RFI.
        'full': {'nant':263, 'maxbase':1000.0} ## Main and Compact arrays
    }
    totalants = 263
    totalvis = totalants * ( totalants - 1)/2.0
    maxbase = 8000.0 # km
    ## Baseline count profile
    base_lengths = np.arange(0.02,maxbase,0.1)   # 1-500 km baselines
    base_counts = base_lengths.copy()
    
    ## Relative numbers of baselines for each baseline length range.
    if arrtype=='approx_ngvla':
        ### 119 within 1km baselines
        ### 114 in spirals out to 100k baselines
        ### 30 outliers upto 1000k baselines
        nbase_upto_1k = 119*118/2.0
        nbase_1k_to_100k = 114*113/2.0 + 114*119  
        nbase_beyond_100k = 30*29/2.0 + 30*233

        base_counts.fill( nbase_1k_to_100k )
        base_counts[ base_lengths < 1.0 ] = nbase_upto_1k
        base_counts[ base_lengths >=100  ] = nbase_beyond_100k

    else:   ## main ngvla ( arrtype = 'ngvla' )

        nbase_upto_5k = 143*142/2.0   # core : 143 ants
        nbase_5k_to_1000k = 90*89.0/2.0 + 90*143  # 233 - 143 = 90 (rest of plains, and midbase)
        nbase_beyond_1000k = 30*29/2.0 + 30*233  # all long baselines

        base_counts.fill( nbase_5k_to_1000k )
        base_counts[ base_lengths < 5.0 ] = nbase_upto_5k
        base_counts[ base_lengths >=1000  ] = nbase_beyond_1000k

    
## Normalize such that the integrated total is 'totalvis'
base_counts = (base_counts / np.sum(base_counts)) * totalvis


#
pl.ion()
pl.plot(base_lengths, base_counts)

