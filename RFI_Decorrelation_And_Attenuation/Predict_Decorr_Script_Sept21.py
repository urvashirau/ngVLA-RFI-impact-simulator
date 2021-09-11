import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import log10,pi
import math
import random
from random import seed
from random import gauss
from matplotlib import pyplot as plt
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord
from itertools import combinations_with_replacement
from IPython.display import display, HTML
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from matplotlib.gridspec import GridSpec

def decorr_predict_tool(
        ha_src=0.0,  dec_src=34.0, 
        ha_rfi=90.0, dec_rfi=34.0,
        speed_rfi=0.0,
        rec_power_spfd=1.0, 
        freq_obs=2.345e9,
        obs_bw=62.5e3,
        tau=1.0, 
        cfg='cfgfiles/vlab.cfg',
        high_thresh=0.9,
        low_thresh=0.1,
        ylog=False):
    """
    Input Parameters :

    ha_src, dec_src :  HourAngle, Declination of the PhaseCenter (degrees)

    ha_rfi, dec_rfi : HourAngle, Declination of the RFI source (degrees)
                              - HourAngle=0 is the N-S local Meridian.
                              - HourAngle=-90, +90 are on the horizon
                              - Declination=latitude is the Zenith at HA=0
                              - Declination=90 is the North Celestial Pole
                              - Dec=0 is the Ecliptic.

    rec_power_spfd : Received SPFD ( Jy )
    
    speed_rfi : Speed of the source (degrees/sec)
                      - For a stationary object (geo-stat, terrestrial) 
                        this is the Earth Rotation rate (0.004 deg/sec)
                      - For a moving satellite, this is its speed (0.06 deg/sec)
                      - For an astronomical source, this is zero.

    freq_obs : Central observing frequency (Hz)
    obs_bw : channel bandwidth (Hz)

    tau : integration time (or UV cell crossing time - for imaging). (sec)
    
    cfg : Array config file
    
    high_thresh : Vis Amp above this implies no decorrelation
    low_thresh : Vis Amp below this implies perfect decorrelation 
    ylog : False, plot amps on a linear scale (power).  True, plot amps on a log scale. None : No plots.
    """
    
    ha_base,dec_base,base_len, ant1, ant2 = baseline_ha_dec_calc(cfg)

    print("\n====================================")
    print("Phasecenter :\t  Hour angle : %3.4f deg \t  Declination : %3.4f deg"%(ha_src,dec_src))
    print("RFI location:\t  Hour angle : %3.4f deg \t  Declination : %3.4f deg"%(ha_rfi,dec_rfi))

    rad=np.max( [np.max(np.abs(ant1)) ,np.max(np.abs(ant2)) ] ) * 1.2
    x_src, y_src = convert_ha_dec_to_xy(ha_src,dec_src,rad)
    x_rfi, y_rfi = convert_ha_dec_to_xy(ha_rfi,dec_rfi,rad)

    print("RFI speed : %3.4f deg/sec"%(speed_rfi))
    print("\nObserving freq : %3.4f GHz \t Channel width : %3.4f kHz"%(freq_obs/1e+9,obs_bw/1e+3))
    print("Array Config : %s\n"%(cfg))
    #plot_ant_base(base_len, antx, anty)

    totphase, theta, phi = total_phase_delay(ha_src,ha_base,ha_rfi,
                                             dec_src,dec_base,dec_rfi,
                                             freq_obs,base_len)

    #vis_amp_config = sinc_theory_rfi(base_len,theta,obs_bw,freq_obs,phi)

    vis_amp_config,uvdist = calc_attenuation(ha_src=ha_src, dec_src=dec_src,
                                             ha_rfi=ha_rfi, dec_rfi=dec_rfi,
                                             ha_base=ha_base,
                                             dec_base=dec_base,
                                             b = base_len,
                                             rfi_speed = speed_rfi, ## 0.00418019,  # Earth rot vel (rad/sec)
                                             tau_int = tau,  #sec
                                             beta = obs_bw)

    ## Multiply by received power level 
    print("Received SPFD at each antenna : %3.2e Jy"%(rec_power_spfd))
    vis_amp_config = vis_amp_config * rec_power_spfd

    base_high, ant1_high, ant2_high, amp_high = get_baseline_subset( base_len, ant1, ant2, 
                                                           vis_amp_config, high_thresh, checktype='high')
    base_low, ant1_low, ant2_low,amp_low = get_baseline_subset( base_len, ant1, ant2, 
                                                           vis_amp_config, low_thresh, checktype='low')

    print("\nBaseline Count :  Above %3.1e Jy : %d \t  Below %3.1e Jy : %d \t  In-between %d "%(high_thresh, len(base_high), low_thresh, len(base_low), len(base_len)-len(base_high)-len(base_low) ))

    median_vis_amp = np.median(vis_amp_config).value
    max_vis_amp = np.max(vis_amp_config).value
    min_vis_amp = np.min(vis_amp_config).value
    print("Median Vis Amplitude : %3.1e Jy\t [Max:%3.1e Jy, Min:%3.1e Jy]"%(median_vis_amp, max_vis_amp, min_vis_amp))

    ret = {'Visibility Amplitude':{}, 'Baseline Counts':{} }
    ret['Visibility Amplitude']['median_vis_amp'] = median_vis_amp
    ret['Visibility Amplitude']['max_vis_amp'] = max_vis_amp
    ret['Visibility Amplitude']['min_vis_amp'] = min_vis_amp
    ret['Baseline Counts']['nbase_high'] = len(base_high)/len(base_len)
    ret['Baseline Counts']['nbase_low'] = len(base_low)/len(base_len)
    ret['Baseline Counts']['high_thresh'] = high_thresh
    ret['Baseline Counts']['low_thresh'] = low_thresh

    if ylog==None:
        return ret

    extraplots=False

    plt.close('all')
    if extraplots==True:
        fig = plt.figure(figsize=(8,8),constrained_layout=True);
        gs = GridSpec(2,2, figure=fig);
        ax1 = fig.add_subplot( gs[0,:] )
        ax2 = fig.add_subplot( gs[1,0] )
        ax3 = fig.add_subplot( gs[1,1] )
    else:
        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot()

    ## Plot Vis Amp
    ax1.plot(base_len, np.abs(vis_amp_config), '.-', color='lightgrey',markeredgecolor='grey',);
    ax1.plot(base_high, np.abs(amp_high), 'b.');
    ax1.plot(base_low, np.abs(amp_low), 'r.');  
    ax1.set_title('Predicted Visibility Amplitude : Jy  [1e-26 W/m2/Hz]')
    ax1.set_xlabel('Baseline Length')
    #ax1.set_ylim(np.min(np.abs(vis_amp_config)), np.max(np.abs(vis_amp_config)))
    if ylog==True:
        ax1.set_yscale('log',nonposy='clip')

    if extraplots==False:
        return ret
   
    for i in range(0,len(ant1)):
        ax2.plot( [ant1[i][0]], [ant1[i][1]] , 'ko' )

    ax2.set_title('Antennas and baselines') 
    ax2.plot( [x_src], [y_src], 'co',markersize=10)
    ax2.plot( [x_rfi], [y_rfi], 'mo',markersize=7)
    
    circ_x = []
    circ_y = []
    for th in range(0,361,20):
        circ_x.append( rad * np.cos(np.radians(th)) )
        circ_y.append( rad * np.sin(np.radians(th)) )
    ax2.plot( circ_x, circ_y, 'g-' )


    ## Plot baseline histogram
    bins = np.arange(0, np.max(base_len), np.max(base_len)/15)
    ax3.hist(x=base_len,bins=bins,color='lightgrey')
    ax3.hist(x=base_low,bins=bins,color='r',alpha=0.5)
    ax3.hist(x=base_high,bins=bins,color='b',alpha=0.5)
    ax3.set_title('Baseline Length Distribution')

    return ret

def baseline_ha_dec_calc(filename):
    
    x,y,z= np.loadtxt(filename,dtype= 'float',usecols=(0,1,2),unpack = True)
    m = np.prod(x.shape) #number of antennas in the file
    
    blist=[]

    HA  = []
    Dec =[]
    base_list=[]
    ant1 = []
    ant2 = []

    #Averaging each column to get a center
    cx = np.mean(x)
    cy = np.mean(y)
    cz = np.mean(z)
    #print('Averages',cx,cy,cz)

    xm =np.array(x- cx)
    ym = np.array(y - cy)
    zm = np.array(z- cz)
 
    for i in range(0,m):
        for j in range(i+1,m):    
         
            dx = x[i]-x[j]   
            dy = y[i]-y[j]
            dz = z[i]-z[j]  

            dxm = xm[i]-xm[j]   
            dym = ym[i]-ym[j]   
            dzm = zm[i]-zm[j]   

             #calculating baselines
            base = np.sqrt(dxm**2+dym**2+dzm**2)
            # print('Baseline',base)

            #calculating ra and dec

            dec = np.arcsin(dz/base)
            ra = np.arccos(dx/(base*np.cos(np.radians(dec))))

            HA.append(np.degrees(ra))
            Dec.append(np.degrees(dec))
            base_list.append(base)  ## Baseline length.
            ant1.append( (xm[i], ym[i]) )
            ant2.append( (xm[j], ym[j]) )

   
    HA = (np.asarray(HA))
    Dec = (np.asarray(Dec))
    #base_list = np.asarray(base_list)

    return sort_on_first(base_list, HA, Dec, ant1, ant2)


def sort_on_first( base_list, HA, Dec, ant1, ant2 ):
    #Sort in increasing baseline order
    from operator import itemgetter
    indices, base_list_sorted = zip(*sorted(enumerate(base_list), key=itemgetter(1)))
    HA_sorted = []
    Dec_sorted = []
    ant1_sorted = []
    ant2_sorted = []
    for i in indices:
        HA_sorted.append( HA[i] )
        Dec_sorted.append( Dec[i] )
        ant1_sorted.append( ant1[i] )
        ant2_sorted.append( ant2[i] )

    return np.array(HA_sorted), np.array(Dec_sorted), np.array(base_list_sorted), np.array(ant1_sorted), np.array(ant2_sorted)


def get_angles(ra,dec):
    radec = SkyCoord(ra, dec, frame='icrs')
    ra = radec.ra.value
    dec = radec.dec.value
    return ra,dec
    
def total_phase_delay(H_s,H_b,H_rfi,d_s,d_b,d_rfi,F_obs,b):
    #H's are the hourangle values in degrees and d's are the declination values in degrees
   
    A = np.sin(np.radians(d_s)) * np.sin(np.radians(d_b)) + np.cos(np.radians(d_s)) * np.cos(np.radians(d_b)) * np.cos(np.radians(H_s-H_b))
    B = np.sin(np.radians(d_rfi)) * np.sin(np.radians(d_b)) + np.cos(np.radians(d_rfi)) * np.cos(np.radians(d_b)) * np.cos(np.radians(H_rfi-H_b))
  
    theta_rfi = np.degrees(np.arccos(B))
    theta_source = np.degrees(np.arccos(A))
    #Calulating total phase delay
    
    c = 3e8 # in m/s
#    D_t = np.degrees(((b *np.abs(A-B) * F_obs * 2*np.pi) / c))  #in degrees
    D_t = ((b *(theta_rfi-theta_source) * F_obs * 2*np.pi) / c)  #in degrees

    tau_delay = b * (A-B) / c

    return D_t,theta_source,theta_rfi


def sinc_theory_rfi(b,theta,F_del,F_obs,phi ): #for the case of an RFI referenced at phase center

    c = 3e8    
    sinc_it = np.sinc(((b *np.sin(np.radians(phi))* F_del * np.pi) / c))  
    cos_factor = np.cos(((b *np.sin(np.radians(phi)) *np.cos(np.radians(theta))* F_obs * 2 * np.pi) / c))
    theory = (cos_factor) * (sinc_it)
 #   y =np.abs(theory) 
 
    return theory

def calc_attenuation(ha_src=0.0, dec_src=34.0,
                      ha_rfi=0.0, dec_rfi=-20.0,
                      ha_base=np.array([-90.0,+90.0]),
                      ###
                      dec_base=np.array([+30.0,-20.0]),
                      b = np.array([200.0, 500.0]),
                      ###
#                      lx = np.array([50,100]),
#                      ly = np.array([150,-300]),
                      ###
                      rfi_speed = 0.00418019,  # Earth rot vel (rad/sec)
                      tau_int = 1.0,  #sec
                      beta = 50e+3):
    """
    ha_src : Hour Angle of the phasecenter  
    dec_src : Declination of the phasecenter 
    ha_rfi : Hour Angle of the phasecenter  
    dec_rfi : Declination of the phasecenter 

    ha_base : [V] Hour angle vector for all baselines
    dec_base : [V] Dec vector for all baselines
    b : [V] baseline lengths (this is a vector over all baselines)

###    lx, ly : [V] Components each baseline in the (H=0,d=0) , (H=-90,d=0) directions. Vectors
    rfi_speed : deg/sec.  
                 Angular velocity : 7.2921159e-5 radians/second == 0.00418019 deg/sec

    tau_int : integration time. ( Or, cell-crossing time.... as in baseline-based averaging )
    
    beta : channel bandwidth
    

    Calculate attenuation from both the fringe averaging and the bandwidth-decorrelation effects.
    att1 : Fringe washing effect
    att2 : Bandwidth decorrelation
    """
    ### Should this be for the phasecenter or the rfi source, or a diff ? 
    #Eqn 2
    lx = b * np.cos(np.radians(ha_base))*np.cos(np.radians(dec_base))
    ly = -1.0 * np.array(b) * np.sin(np.radians(ha_base))*np.cos(np.radians(dec_base))
    lz = b *np.cos(np.radians(ha_base))+np.sin(np.radians(dec_base))
        
    u = lx*np.sin(np.radians(ha_src)) + ly*np.cos(np.radians(ha_src))
    v = lz * np.cos(np.radians(dec_src)) - lx*np.sin(np.radians(dec_src))*np.cos(np.radians(ha_src)) + ly*np.sin(np.radians(dec_src))*np.sin(np.radians(ha_src))
        
    u = lx*np.sin(np.radians(ha_src)) + ly*np.cos(np.radians(ha_src))
    # Eqn 4
    fringe_freq = np.radians(rfi_speed) * u * np.cos(dec_src)
    # Eqn 3
    att1 = np.sinc( np.pi * fringe_freq * tau_int )

    ## Eqn 17 to get the time delay
    A = np.sin(np.radians(dec_src)) * np.sin(np.radians(dec_base)) + np.cos(np.radians(dec_src)) * np.cos(np.radians(dec_base)) * np.cos(np.radians(ha_src-ha_base))
    B = np.sin(np.radians(dec_rfi)) * np.sin(np.radians(dec_base)) + np.cos(np.radians(dec_rfi)) * np.cos(np.radians(dec_base)) * np.cos(np.radians(ha_rfi-ha_base))

    # Eqn 18
    tau_delay = b * np.abs(A-B) / 3e8
    att2 = np.sinc(np.pi * beta * tau_delay)

    # Eqn 19 (combining both effects).
    attenuation = np.sqrt(att1**2 * att2**2)
    uvdist = np.sqrt(u**2 + v**2)

    return attenuation,uvdist

def set_ha_dec_phasecenter_Cconfigdata():
    """
    Set parameters for the C-config data
    """
    ## Source Location (changes from 2.36 hrs to 2.373 hrs for C config data) 
    ha_s,dec_s = get_angles(ra='2h20m0s', dec='16d22m')
    #print("Phasecenter :  Hour angle : %3.4f deg \t  Declination : %3.4f deg"%(ha_s,dec_s))
    return ha_s, dec_s

def set_ha_dec_satellite_Cconfigdata():
    """
    Set parameters for the C-config data
    """
    #VLA Location : 34.0784 deg N, 107.6184 deg W
    #Satellite Locator Website:  https://www.groundcontrol.com/Satellite_Look_Angle_Calculator.htm
    
    #Address: VLA
    lat= 34.0784 #deg
    long= -107.6184 #deg
    
    E= 5.4 #deg #elevation
    Azi= 22#186.2 #deg,Azimuth
    
    declination_sat = (np.degrees(np.arcsin(np.sin(np.radians(E))*np.sin(np.radians(lat))+ np.cos(np.radians(E))*np.cos(np.radians(lat))*np.cos(np.radians(Azi)))))
    RA_sat = np.degrees(np.arcsin(-np.sin(np.radians(Azi))*np.cos(np.radians(E))/np.cos(np.radians(declination_sat))))
    
    ha_rfi = RA_sat
    dec_rfi = declination_sat
    
    #print("RFI :  Hour angle : %3.4f deg \t  Declination : %3.4f deg"%(ha_rfi,dec_rfi))
    return ha_rfi, dec_rfi


                                                        
def get_baseline_subset( base_len, ant1, ant2, 
                         vis_amp_config, 
                         thresh, 
                         checktype='low'):

    base_sub = []
    ant1_sub = []
    ant2_sub = []
    amp_sub = []

    for i in range(0,len(vis_amp_config)):
        if checktype=='low':
            if np.abs(vis_amp_config[i]) < thresh:
                base_sub.append(base_len[i])
                ant1_sub.append(ant1[i])
                ant2_sub.append(ant2[i])
                amp_sub.append(vis_amp_config[i])
        if checktype=='high':
            if np.abs(vis_amp_config[i]) > thresh:
                base_sub.append(base_len[i])
                ant1_sub.append(ant1[i])
                ant2_sub.append(ant2[i])
                amp_sub.append(vis_amp_config[i])

    return base_sub, ant1_sub, ant2_sub, amp_sub


def convert_ha_dec_to_xy(ha,dec,rad):
    """
    convert from HA DEC to an approximate location on the array config plot, for plotting in a vertical projection. 
     - 0 deg HA is the N-S line going through (x.y) = (0,0) 
     - 90 deg HA should map to the max radius of the horizon circle.
     - Similarly, max radius represents 90 deg of declination, shifted by latitude.... since the NCP is at 34deg Alt.
    """
    lat = dec - 34.0
    lon = ha

    circum = 2 * np.pi * 6369 #Earth circum, in km
    #circum = 2 * np.pi * rad/1000.0   # Circumference of our horizon circle.... (km)

    x = lon * circum * np.cos(lat*np.pi/360.0)/360.0
    y = lat * circum/360.0

    r = np.sqrt(x**2 + y**2)
    th = np.angle(np.complex(x,y))

    if r > rad:
        x = rad * np.cos(th)
        y = rad * np.sin(th)

    return x,y




import matplotlib.pyplot as plt
import numpy as np
from pycraf import conversions as cnv
from astropy import units as u
from astropy import constants
import math

def calc_vla_antenna_gain(ha_src, dec_src, ha_rfi, dec_rfi, obs_freq):
    """
    #Calculate the angular distance between the phasecenter and the rfi source. 
    ## Angular distance : https://en.wikipedia.org/wiki/Angular_distance
    """
    ang_dist = np.degrees(np.arccos( np.sin(np.radians(dec_src)) * np.sin(np.radians(dec_rfi))  + 
                          np.cos(np.radians(dec_src)) * np.cos(np.radians(dec_rfi)) * np.cos(np.radians(ha_src-ha_rfi))  ))
    ## Primary Beam HPBW for a 25m dish (VLA)
    hpbw = np.degrees((3e+8/obs_freq)/25.0 )

    ### Check if the VLA sees the RFI through a mainlobe or a sidelobe..... 
    if ang_dist > hpbw * 2:  # Assume x2 to get to first sidelobe.
        vla_lobe='SideLobe'
        vla_gain = 0 
    else:    # MainLobe ( https://www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain)
        vla_lobe='MainLobe'
        d=25.0 # diameter in m
        lam = 3e8/obs_freq # lamda in m
        vla_gain = 10*log10( 0.75 * (pi*d/lam)**2 )  

    print("Phasecenter and RFI are %3.2f deg apart. PB main lobe at %3.2f GHz : %3.2f deg."%(ang_dist,obs_freq/1e9,hpbw))
    return vla_lobe, vla_gain,ang_dist


def calc_received_power(ha_src, dec_src, 
                        ha_ut,dec_ut, dist_ut, freq_ut, 
                        ha_sat,dec_sat, freq_sat, obs_bw):   
    """
    Calculate the received power, given the source power and the gains of the transmitting and receiving antennas.

    ha_src, dec_src : HA/DEC of the observation phasecenter (deg)
    ha_ut, dec_ut : HA/DEC of the UT (deg)
    dist_ut : Distance betwen the UT and a VLA antenna
    freq_ut : Observing frequency for the UT uplink test (Hz)
    ha_sat, dec_sat : HA/DEC of the SAT (deg)
    freq_sat : Observing frequency for the SAT downlink test (Hz)
    obs_bw : Channel bw (Hz)

    """

    ################################################
    #### VLA antenna gain : MainLobe or SideLobe ?  Gain ?
    ################################################

    ################################################
    #### UT SideLobe Gain, and SPFD at the VLA antenna, from the UT at a given distance.
    ################################################
    gtx= 0 * cnv.dBi  ## Always a sidelobe
    ptx = 3.2 * u.W   ## Transmitted power (confirmed by SpaceX as their maximum EIRP)
    dist = dist_ut * u.km  ## Distance
    rec_pfd = cnv.powerflux_from_ptx(ptx, dist, gtx)   ## W/m2   in a 60MHz transmission band
    rec_ut_spfd = rec_pfd /(60e6 * u.Hz)   ## W/m2/Hz
    rec_ut_spfd_db = 10*log10(rec_ut_spfd/(1.0*u.W/(u.m*u.m*u.Hz)))
    #print("SPFD incident at the VLA antenna (from UT SideLobe at %3.2f km distance): %3.2f [dB W/m2/Hz] "%(dist_ut,rec_ut_spfd_db))

    #vla_lobe_ut, vla_gain_ut = calc_vla_antenna_gain(ha_src, dec_src, ha_ut, dec_ut, freq_ut)
    #print("UT Uplink : VLA %s gain : %3.2f dB"%(vla_lobe_ut,vla_gain_ut))

    ################################################
    #### SAT MainLobe or SideLobe ?  Calculate SPFD at the VLA antenna.
    ################################################
    ## PFD from the SAT (at the Earth's surface). MainLobe : -146.0 dB W/m2/4kHz  ( from SpaceX )
    ## Satellite forward gain = 38.3 dB.  Power in a sidelobe = Power in a mainlobe (dB) - 38.3 (dB)
    sat_pfd={'MainLobe':-146,  'SideLobe':-146 - 38.3}  ## dB : W/(m2(4kHz))

    ## Do we see a SAT Mainlobe or sidelobe ? 
    if dist_ut<22.0: ##  SpaceX Sats expect to have a 22km footprint on the ground.
        sat_lobe='MainLobe'
    else:
        sat_lobe='SideLobe'

    pfd_db = sat_pfd[sat_lobe]
    pfd_w = 10**(pfd_db/10.0)
    spfd_w_hz = pfd_w/4e+3   ## W/(m2 Hz)
    spfd_db_hz = 10*log10(spfd_w_hz)
    #print("SPFD incident at the VLA antenna (from SAT %s): %3.2f [dB W/m2/Hz]"%(sat_lobe,spfd_db_hz))
    sat_spfd_db=spfd_db_hz

    #vla_lobe_sat, vla_gain_sat = calc_vla_antenna_gain(ha_src, dec_src, ha_sat, dec_sat, freq_sat)
    #print("SAT Downlink : VLA %s gain : %3.2f dB"%(vla_lobe_sat,vla_gain_sat))

    ######################################################
    ### Choose what VLA and RFI mainlobe sidelobe combination is relevant here. 
    sigtypes = {'UT':{'rfi_dist':dist_ut,'rfi_spfd':rec_ut_spfd_db,'freq':freq_ut,'rfi_lobe':'SideLobe'}, 
                       'SAT':{'rfi_dist':570.0,'rfi_spfd':sat_spfd_db,'freq':freq_sat, 'rfi_lobe':sat_lobe}}

    for asig in ['UT','SAT']:
        spfd_db = sigtypes[asig]['rfi_spfd'] ## dB
        freq = sigtypes[asig]['freq'] * u.Hz


        if asig=='UT':
            vla_lobe, vla_gain,ang_dist = calc_vla_antenna_gain(ha_src, dec_src, ha_ut, dec_ut, freq_ut)
        else:
            vla_lobe, vla_gain,ang_dist = calc_vla_antenna_gain(ha_src, dec_src, ha_sat, dec_sat, freq_sat)
        sigtypes[asig]['vla_lobe'] = vla_lobe
        sigtypes[asig]['vla_gain'] = vla_gain
        sigtypes[asig]['ang_dist'] = ang_dist


        grx = sigtypes[asig]['vla_gain'] * cnv.dBi
        rfi_lobe = sigtypes[asig]['rfi_lobe'] 
        dist_rfi = sigtypes[asig]['rfi_dist']

        area = 0.75 * 3.14*(12.5**2)  * u.m*u.m   ## Collecting area of a VLA dish :  m2
        chanbw = obs_bw * u.Hz    ## Chan bw  : Hz

        print("\n%s %s --> Beam Interaction : %s %s and VLA %s "%(asig, ("Up-Link" if asig=='UT' else "Down-Link"),asig,rfi_lobe, vla_lobe ) )
        print("================================================================")
        #print("Beam interaction : %s %s and VLA %s "%(asig,rfi_lobe, vla_lobe ))
        
        print("SPFD incident at the VLA antenna from %3.2f km distance : %3.2f [dB W/m2/Hz] "%(dist_rfi,spfd_db))
        print("VLA %s gain : %3.2f dB"%(vla_lobe,sigtypes[asig]['vla_gain']))

        spfd_w = 10**(spfd_db/10.0) * u.W / (u.m*u.m*u.Hz)  ##  W/m2/Hz
        pfd_rec_chan = spfd_w * chanbw   ## Power Flux Density : W/m2 in one 125kHz channel.... 

        rec_pwr = cnv.prx_from_powerflux(powerflux=pfd_rec_chan, grx=grx, freq=freq)   # W  in a 125kHz channel
        rec_pwr_W = rec_pwr / (1.0*u.W)   # W
        rec_pwr_WdB = 10*log10(rec_pwr_W) # dB W
        print("Recd power in one 125kHz channel :\t %3.2f dB W \t\t[%3.2e W]"%(rec_pwr_WdB, rec_pwr_W))

        rec_spfd = rec_pwr / (area*chanbw)     # W / m2 /Hz    --> to get to Units of Jy...
        rec_spfd_W = rec_spfd / (1.0*u.W/(u.m*u.m*u.Hz))
        rec_spfd_dB = 10*log10(rec_spfd_W)
        print("Recd spectral flux density :\t\t %3.2f dB W/m2/Hz \t[%3.2e W/m2/Hz] [%3.2e Jy]"%(rec_spfd_dB,rec_spfd_W,rec_spfd_W/1e-26))

        sigtypes[asig]['recd_spfd_db'] = rec_spfd_dB
        sigtypes[asig]['recd_spfd_jy'] = rec_spfd_W/1e-26

    return sigtypes


def get_lst_range(obs_date = '21-Sep-09 18:50:25'):
    from casatools import measures,quanta
    me = measures()
    qa =quanta()
    
    t1 = me.epoch('utc', obs_date)
    me.doframe(me.observatory('VLA'))
    me.doframe(t1)
    t2 = me.measure(t1,'LAST')
    
    d1=me.riseset(me.direction('sun'))
    t_rise = d1['rise']['last']
    t_set = d1['set']['last']
    
    print('LST now: ' + str(qa.time(qa.sub(t2['m0'],qa.floor(t2['m0'])))))
    print('LST sunrise: '+str(qa.time(qa.sub(t_rise['m0'],qa.floor(t_rise['m0'])))))
    print('LST sunset: '+str(qa.time(qa.sub(t_set['m0'],qa.floor(t_set['m0'])))))

    lst_sunrise = qa.time(qa.sub(t_rise['m0'],qa.floor(t_rise['m0'])))
    lst_sunset = qa.time(qa.sub(t_set['m0'],qa.floor(t_set['m0'])))

    return lst_sunrise, lst_sunset

def print_summary(summary={},attenuation_ut={},attenuation_sat={}):
    ## Accumulate summary
    summary['UT']['Visibility Amplitude'] = attenuation_ut['Visibility Amplitude']
    summary['UT']['Baseline Counts'] = attenuation_ut['Baseline Counts']
    summary['SAT']['Visibility Amplitude'] = attenuation_sat['Visibility Amplitude']
    summary['SAT']['Baseline Counts'] = attenuation_sat['Baseline Counts']

    ## Print version
    psum={}
    for asig in ['UT','SAT']:
        psum[asig]={}
        psum[asig]['Beam Interaction'] = "VLA %s + %s %s"%(summary[asig]['vla_lobe'], asig, summary[asig]['rfi_lobe'])
        psum[asig]['Distance : VLA to RFI [km]'] = "%3.2f km"%(summary[asig]['rfi_dist'])
        psum[asig]['Angular Separation : Obs to RFI [deg]'] = "%3.2f deg"%(summary[asig]['ang_dist'])
        psum[asig]['SPFD Measured per antenna [Jy,dB]'] = "%3.1e Jy [%3.1f dB]"%(summary[asig]['recd_spfd_jy'],summary[asig]['recd_spfd_db'])
        psum[asig]['Median Visibility Amp [Jy, dB]'] = "%3.1e Jy [%3.1f dB]"%(
            summary[asig]['Visibility Amplitude']['median_vis_amp'],
            10*log10(summary[asig]['Visibility Amplitude']['median_vis_amp']*1e-26))
        psum[asig]['Max Visibility Amp [Jy],dB '] = "%3.1e Jy [%3.1f dB]"%(
            summary[asig]['Visibility Amplitude']['max_vis_amp'],
            10*log10(summary[asig]['Visibility Amplitude']['max_vis_amp']*1e-26))
        psum[asig]['Min Visibility Amp [Jy],dB '] = "%3.1e Jy [%3.1f dB]"%(
            summary[asig]['Visibility Amplitude']['min_vis_amp'],
            10*log10(summary[asig]['Visibility Amplitude']['min_vis_amp']*1e-26))

        psum[asig]['Baseline Count : > %s Jy (%3.1f dB)'%(summary[asig]['Baseline Counts']['high_thresh'],10*log10(summary[asig]['Baseline Counts']['high_thresh']*1e-26))] = "%3.2f %%"%(summary[asig]['Baseline Counts']['nbase_high']*100)
        psum[asig]['Baseline Count : < %s Jy (%3.1f dB)'%(summary[asig]['Baseline Counts']['low_thresh'],10*log10(summary[asig]['Baseline Counts']['low_thresh']*1e-26))] = "%3.2f %%"%(summary[asig]['Baseline Counts']['nbase_low']*100)

    pdf = pd.DataFrame(psum)
    display(HTML(pdf.to_html()))

