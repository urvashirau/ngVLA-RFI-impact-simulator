import numpy as np
import pylab as pl
import matplotlib
import pandas as pd
pl.ion()

#execfile('rfichar.py')
#execfile('rfisol.py')
#execfile('ngvlapar.py')

from rfichar import *
from rfisol import *
from ngvlapar import *


def calc_base_loss_frac(arrayfrac='core', asol='D',siglen=1e-05,datares=1e-05, decorr='no_decorr', dynr='1e-04'):
    """
    Calculate fraction of baselines lost to RFI
          -- Fraction of baselines over which RFI is seen and remains correlated. 
          -- Account for antennas that see RFI versus those that do not
          -- Account for decorrelation due to baseline length, timestep and channel width. 
    """
    ## Number of antennas that see the RFI
    nant_see_rfi = arrayfracinfo[arrayfrac]['nant']

    ## Nbaselines over which RFI is seen in both antennas per pair
    nbase_see_rfi = ( nant_see_rfi*(nant_see_rfi - 1)/2.0 )  

    # Fraction of the data remains correlated (given baseline length and time/freq res and signal length) -- upto maxbase. 
    ## range : 0 to 1.0.  For a small maxbase and mostly short baselines, this is close to 1.0
    if decorr=='no_decorr':
        frac_corr_rfi = 1.0
    else:
        ## Max angular offset from phase reference center. This is half of the image field of view. 
        fov=90.0 ## Full Sky
        if decorr=='prac_decorr':
            fov=20.0 ###  average expected dist to a satellite or ground RFI
                            ### no point considering imaging decorrelation (RFI varies)
        frac_corr_rfi = base_corr_frac(fov=fov, 
                                       dnu=solutions[asol]['chanres'], 
                                       maxb=arrayfracinfo[arrayfrac]['maxbase'],
                                       siglen=siglen,
                                       datares=datares,
                                       dynr=dynr)

    ## Nbaselines over which RFI remains correlated.
    nbase_corr_rfi = nbase_see_rfi * frac_corr_rfi

    ## Nbaselines over which RFI gets decorrelated.
    nbase_decorr_rfi = nbase_see_rfi * ( 1.0 - frac_corr_rfi )

    ## Return the fractions of total baselines that see correlated RFI and that see decorrelated RFI (higher noise)
    return nbase_corr_rfi / totalvis , nbase_decorr_rfi / totalvis


def base_corr_frac(fov=90.0, 
                   dnu=50.0, 
                   maxb=1.0, 
                   siglen=1e-05,
                   datares=1e-05,
                   ret2=False, 
                   dynr='1e-04'):
    """
    Fraction of baselines with correlated RFI is less than 1.0 only if :  
      -- We are working with correlated visibilities : datares >= 1e-06.  Smaller timescales are assumed to be antenna-based voltages. 
      -- Data time resolution is smaller than the signal time resolution :   TODO : Implement from Rick's memo...

    For a given baseline length, calculate the attenuation of the signal due to decorrelation. 

    fov = Field of View in degrees
    dnu = Del nu, in kHz
    base = longest baseline, in km

    Implements eqn 12 (and Fig 3) from NGVLA Memo #38
    """

    ######  If siglen is smaller than datares, attenuation can happen across siglen.
    ######  fov 90 implies all baselines see the RFI end-on. This is not the case. Pick an average effective 'fov'.
    ######  light travel time between baselines is 3.5 msec for 1000km baselines. => How to even get higher time-res visibilities for those baselines. 

    ## By default, assume that the RFI remains correlated across all baselines.
    fracb = 1.0
    attenuation = base_lengths.copy()
    attenuation.fill(1.0)

    ## Calculate RFI decorrelation fraction 
#    if datares <= siglen and datares >= 1e-06:
    if datares >= 1e-06:

        atten_len = min(datares,siglen)

        ### QUESTION :  SIN or COS ?!!!!

        ## Light travel time difference between pair of antennas
        tdel = (1e+3 * base_lengths) * np.sin(fov*np.pi/180.0) / 3.0e+8
        ## Phase change across the channel bandwidth when phase referenced to the middle (? factor of 2 here ? )
        freqphi = 2*np.pi* tdel * (dnu*1e+3)
        ## Attenuation due to decorrelation (from channel bandwidth)
        freqdecorr = np.fabs(np.sinc( freqphi / (2*np.pi) ))

        ## Fringe frequency
        fdel = 3.0e+8 / ( (1e+3 * base_lengths) * np.sin(fov*np.pi/180.0) )
        ## Phase change across timestep when phase-references to the middle of the timestep. ( ? factor of 2 here ? )
        timephi = 2*np.pi* atten_len * fdel
        ## Attentuation due to decorrelation from fringe rotation in time (per baseline)
        timedecorr = np.fabs(np.sinc( timephi / (2*np.pi) ))
        
        #timedecorr=1.0

        ## Total attenuation from channel and time averaging to produce a single single visibility
        attenuation = freqdecorr * timedecorr

        ## Amount of attenuation at which point RFI 'disappears' is the same as target imaging dynamic range.
        decor_lim = eval(dynr)

        attenuation[ attenuation < decor_lim ] = 0.0  ## RFI effectively goes away, but noise increases.
        attenuation[ attenuation >= decor_lim ] = 1.0   ## RFI is attenuated, but is still going to be there. 

        ## Fraction of baselines over which RFI is correlated
        fracb = np.sum(attenuation[base_lengths<maxb] * base_counts[base_lengths<maxb]   )  / np.sum(base_counts[base_lengths<maxb])  

#        print('Baseline fraction over which RFI is correlated : ' + str(fov) + ' deg ' + str(dnu)+'kHz : max b : ' + str(maxb) + " : " + str(fracb) + " --- " + str(np.mean(attenuation)) )

    if ret2==True:
        return fracb, base_lengths, base_counts, attenuation
    else:
        return fracb #, attenuation



### Calculate RFI duty cycle, when viewed at the data resolutoin
def calcfilling(siglen=1e-05,siggap=5e-05, datares=1e-05,asol='D',arrayfrac='core', rfitype=''):
    """
    Sol = C :  modeling and subtraction
            Assume that the timerange is fully-filled with RFI.
            But, if the datares is smaller than the signal length over which coherence is preserved (siglen),
            assume a 60% RFI excision efficiency. i.e.  60% of the timesteps are recovered 'perfectly' and 40% are lost.

            TODO : Alternate criteria :  All RFI is only partially subtracted, making it similar to 'decorrelation' and equivalent to 'higher noise'

    """
    ## If modeling and subtraction
    ##     and data resolution is finer than signal 
    ##     and data res is what is available after correlation
    if asol=='C' and datares <= siglen and datares >= 1e-06 : 
        frac = 0.5  # 50% efficiency in removal, for baselines where correlation is maintained
        if rfitype.count('LEO')>0:
            frac = 0.9
    else:
        ## Filling factor for the chosen datares. 
        frac = max(siglen,datares)/(siglen+siggap)

        ## Completely filled at the datares time resolution. Usually the case for low time resolution (big datares)
        if frac>1.0:   
            frac=1.0

        ## For very small datares, there needs to be enough RFI on/off across the data window available at that datares.
        ## If the visible window itself is smaller than 3 * siglen, then assume that there is insufficient information to
        ## distinguish between RFI and non-RFI and so all the data within datares will be lost to RFI. 
        datawin = datares * 100.0
        if datawin < 3*siglen:
            frac=1.0
   
    return frac


def calc_noise_equiv_data_loss(dnoise=0.2):
    ### Imnoise = Tsys / sqrt(N)
    ### N * dataloss where dataloss is between 0 and 1
    ### Imnoise = Tsys / sqrt(  (N * dataloss) * 1/dataloss )   --- need 1/dloss more time to get to same sensitivity.
    ###
    ### Inmoise = Tsys * 1.x / sqrt( N * 1.x * 1.x )  --- need  (1.x)^2 more time to get to same sensitivity
    ### This is the same as an equivalent dataloss of 1/(1.x)^2

    ### Noise goes up by dnoise (between 0 and 1)
    noise = 1.0+dnoise
    
    ### Fraction of data loss ( or extra time needed to compensate for it )
    nfrac = (noise*noise) - 1.0
    return nfrac
    
    


################################################

def get_data_loss_fraction(rfichar={}, frange=[1.0,120.0],sol=['D'],verbose=False, decorr='no_decorr', dynr='1e-04'):

    totalx = pl.exp(pl.arange(pl.log(frange[0]),pl.log(frange[1]),0.001))

    maxtotal={}
    thresh = [0.05,0.25,0.5]
    for th in thresh:
        maxtotal[str(th)] = pl.zeros(len(totalx)) + 1.0

    totals={}
    rfinames=[]
    procpars={}
    ## For each type of RFI
    for atom in rfichar.keys():
        rfinames.append(atom)
        totals[atom] = pl.zeros(len(totalx))
        procpars[atom] = {'datares':pl.zeros(len(totalx)), 'sol':['D']*len(totalx)}

        ## Calculate fractions of baselines where RFI is correlated vs uncorrelated
        timebasefrac={}
        timebasefrac_uncorr={}
           ## For correlated RFI, what's the best we can do ? 
        bestsol='D'
        bestres=1.0 #sec
        bestlossfrac = 1.0
        
        ## Loop over chosen RFI mitigation options and pick the one that gives the smallest data loss. 
        for asol in sol:
            timebasefrac=0.0
            ## Add together contributions 
            for ii in range(0,len(rfichar[atom]['timefrac'])):
                loss=1.0
                ## Try all viable time resolutions for the chosen mitigation option, and pick the one that gives the smallest data loss. 
                for dress in solutions[asol]['datares']:
                    tbfrac, tbfrac_uncorr =calc_base_loss_frac(
                        arrayfrac=rfichar[atom]['arrayfrac'][ii],
                        asol=asol,
                        siglen=rfichar[atom]['timeres'], 
                        datares=dress,
                        decorr=decorr,
                        dynr=dynr)

                    tbfrac = tbfrac * rfichar[atom]['timefrac'][ii]
                    tbfrac_uncorr = tbfrac_uncorr * rfichar[atom]['timefrac'][ii]

                    aloss =  (calcfilling(siglen=rfichar[atom]['timeres'], siggap=rfichar[atom]['timegap'], datares=dress, asol=asol,arrayfrac=rfichar[atom]['arrayfrac'][ii], rfitype = atom)  )

                    ### If RFI filling is low, then assume that 'A' can take it out completely, so no decorr noise.
                    aloss_uncorr = 0.1 if ( asol=='A' and aloss < 0.5 ) else aloss

                    tlosspart = tbfrac*aloss + tbfrac_uncorr*aloss_uncorr * calc_noise_equiv_data_loss(dnoise=0.2)
                    
#                    if verbose==True:
#                        print (atom + ' : Signal res/gap = ' + rfichar[atom]['timeres'] , "/", rfichar[atom]['timegap'] ,' -->  Sol', asol, ' at datares ',dress, 's has RFI corr frac ', tbfrac, ' uncorr', tbfrac_uncorr*calc_noise_equiv_data_loss(0.2), ' at loss ', aloss

                    ## Pick the datares that gives the best RFI mitigation
                    if tlosspart < loss:
                        loss = tlosspart
                        bestres = dress
                
                ## Add the contributions from different RFI sources
                timebasefrac = timebasefrac + loss

            if timebasefrac < bestlossfrac:
                bestlossfrac = timebasefrac
                bestsol = asol

        for ran in rfichar[atom]['freqrange']:
            for tt in range(0,len(totalx)):
                if totalx[tt]>=ran[0] and totalx[tt]<=ran[1]:
                    totals[atom][tt] = bestlossfrac
                    procpars[atom]['datares'][tt] = dress
                    procpars[atom]['sol'][tt] = bestsol
                    
                    ## Special case for geo-stat satellites.... (iridium and sirius XM)
                    if atom.count('LEO')>0 and ran[1] < 8.0 and 'C' in sol:
                        totals[atom][tt] = bestlossfrac*(0.5/0.9)   ## To match 0.5 for non-LEO (line 146)

    #print('Total RFI loss with '+ str(sol)+ ' across freqs '+ str(frange) + ' = ('+ str(pl.mean(maxfrac)) + ' to ' + str(pl.mean(totalfrac)) )

    return totalx, totals, procpars


def get_observing_time_ratio(totalx, totals):

    totalfrac = np.zeros( len(totalx) )
    maxfrac = pl.zeros( len(totalx) )

    for ii in range(0,len(totalfrac)):
        for atom in totals.keys():
            if totals[atom][ii] > maxfrac[ii]:
                maxfrac[ii] = totals[atom][ii] 
            totalfrac[ii] = totalfrac[ii] + totals[atom][ii]

    totalfrac[ totalfrac > 1.0 ] = 1.0

    ttime_min = pl.zeros(len(totalfrac))
    ttime_max = pl.zeros(len(totalfrac))
    for ii in range(0,30):
       ttime_min = ttime_min + maxfrac**ii
       ttime_max = ttime_max + totalfrac**ii

#    print ('Required total observing time to compensate for data loss and reach the same continuum imaging sensitivity: ' + str( pl.mean(ttime) ) )

    return ttime_min, ttime_max, np.mean(maxfrac), np.mean(totalfrac)

def get_compute_intensity(totalx, procpars):
    
    compute_intensity = {}
    for atom in procpars.keys():
        compute_intensity[atom] = np.zeros( len(totalx) )

    tres = 1.0    ## ntimes in 1 second
#    fres = 1e+5  ## nchans in 100 MHz

    for ii in range(0,len(totalx)):        ## for each freq
        for atom in procpars.keys():    ## for each RFI type
            datares = procpars[atom]['datares'][ii]
            sol = procpars[atom]['sol'][ii]

            if datares>0.0:
                ntimes = tres / float(datares)
                #            nchans = fres/solutions[sol]['chanres']
                nchans=1
                npol=4
            
                ## Number of visibilities per 1 second (and per channel)
                nvis = totalvis * ntimes * nchans * npol

                if sol=='A':
                    nvis = totalants * ntimes * nchans + npol/2
                
                ## Number of operations per datum ( relative ! )
                nops = 1   # if sol is 'D' it's just a one-step average calculation
                if sol=='B':
                    nops = 20    # 20 times more operations : several averages, etc.
                if sol=='C':
                    nops = 200  # 200 times more operations : matrix decompositions
                    
                comp_int_atom = nvis * nops

                compute_intensity[atom][ii] = comp_int_atom

    return compute_intensity



