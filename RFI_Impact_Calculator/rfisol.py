import numpy as np
import pylab as pl
import matplotlib
import pandas as pd
pl.ion()

#execfile('rfichar.py')
#execfile('ngvlapar.py')

from rfichar import *
from ngvlapar import *


"""
RFI solutions
"""
solutions = {
'A':{'desc':"Antenna based high resolution flagging", 
     "Effort":"Medium",  "Risk" : "low",  "Gains" : "High", 
     'datares':np.array([1e-10, 1e-09, 1e-08]),## Data windows are 100 times higher.
     'chanres':1e+6},  ## (in kHz) nanosec pulses, 1GHz bw.
'B':{'desc':"Baseline based high resolution flagging", 
     "Effort" : "Medium",  "Risk" : "low",  "Gains" : "Medium",
#    'datares' : np.array([1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02]),
     'datares' : np.array([1e-05, 1e-04, 1e-03, 1e-02]),
     'chanres' : 100.0},
'C':{'desc':"High resolution modeling and subtraction", 
     "Effort" :"High", "Risk" : "high", "Gains" : "High",
     'datares' : np.array([1e-06, 1e-05, 1e-04, 1e-03]),
     'chanres' : 10.0 },
'D': {'desc':"Post-processing flagging",
      "Effort" : "Low", "Risk" : "low", "Gains" : "Low",
      'datares' : np.array([0.1, 1.0, 10.0]),
      'chanres' : 1e+3}
}

#Type E : smart scheduling in time/freq/direction ( Effort : Medium, Risk : Low, Gains : High )

def cost_of_solutions():
    
    ## Base cost is type D, defined in terms of numbers of visibilities, with 0.1 sec timesteps and 2k channels per band. . 
    timeress = pl.exp(pl.arange(-10,3))
    
    pl.figure(3)
    pl.clf()

    cols = {'A':'red', 'B':'blue', 'C': 'orange', 'D':'green'}
    labs = {'A':'Antenna-based Flagging', 'B':'In-correlator Flagging', 'C': 'Modeling and subtraction', 'D':'Post Processing Flagging'}

    for sol in ['A','B','C','D']:
    ## ntimesteps per 1 sec and nchan per 1MHz (1e+3 kHz)

        tres = 1.0      ## ntimes in 1 second.
        fres = 1e+5   ## nchans in 100 MHz ( 1e+5 kHz ) : 1 spw. 

        ntimes = tres / solutions[sol]['datares'] 
        nchans = fres / solutions[sol]['chanres'] 
        npol = 4

        volperchunk = 1
        parwidth = 1  ## Parallelization by N times

        totalvol = totalants * ntimes * nchans * npol

        if sol=='A':
            parwidth = 'antenna'
            parquant = totalants
        if sol=='B':
            parwidth = 'baseline'
            parquant = totalvis
        if sol=='C':
            parwidth = 'time,freq'
            parquant = ntimes*nchans
        if sol=='D':
            parwidth = 'baseline'
            parquant = totalvis

        volperchunk = totalvol / (parquant)


#        pl.subplot(211)
        pl.plot(solutions[sol]['datares'], (totalvol), '.-', color=cols[sol], label=labs[sol])
        pl.plot(solutions[sol]['datares'], volperchunk, '.--', color=cols[sol] , label=labs[sol] + ' [ parallelize on ' + parwidth  + ']')
        
#        pl.subplot(212)
#        pl.plot(solutions[sol]['datares'], parquant, '.-', color=cols[sol] , label=sol+', partitioned by '+parwidth)
        

#    pl.subplot(211)
    pl.xscale('log')
    pl.yscale('log')
    pl.legend()
#    pl.xlabel('Data time resolution (sec)')
    pl.ylabel('N samples [ per '+str(tres)+ ' sec, '+str(fres/1e+3) + ' MHz ]')
    pl.title('Data rates for RFI mitigation')
    pl.ylim(0.1,10e+18)
    pl.xlabel('Time resolution (sec) of RFI mitigation solutions')

#    pl.subplot(212)
#    pl.xscale('log')
#    pl.yscale('log')
#    pl.legend()
#    pl.xlabel('Data time resolution (sec)')
#    pl.ylabel('Parallelization width')
    
    pl.savefig('fig_rfi_mitigation_cost.png')
