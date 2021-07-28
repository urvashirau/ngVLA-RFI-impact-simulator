import matplotlib.pyplot as plt
import numpy as np
import math
import random
from random import seed
from random import gauss
from matplotlib import pyplot as plt
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord
from itertools import combinations_with_replacement

from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from matplotlib.gridspec import GridSpec

def decorr_predict_tool(
        ha_src=0.0,  dec_src=34.0, 
        ha_rfi=90.0, dec_rfi=34.0,
        speed_rfi=0.0,
        bw=62.5e3,
        tau=1.0, 
        cfg='cfgfiles/vlab.cfg',
        high_thresh=0.9,
        low_thresh=0.1):
    """
    Input Parameters :

    ha_src, dec_src :  HourAngle, Declination of the PhaseCenter (degrees)
    ha_rfi, dec_rfi : HourAngle, Declination of the RFI source (degrees)
                              - HourAngle=0 is the N-S local Meridian.
                              - HourAngle=-90, +90 are on the horizon
                              - Declination=latitude is the Zenith at HA=0
                              - Declination=90 is the North Celestial Pole, Dec=0 is the Ecliptic.
    
    speed_rfi : Speed of the source (degrees/sec)
                      - For a geo-stationary or terrestrial object, this Earth rotation (0.004 deg/s)
                      - For a moving satellite, this is its speed (LEO : 0.06 deg/s)
                      - For an astronomical source, this is zero.
    
    bw : channel bandwidth

    tau : integration time (or UV cell crossing time - for imaging). 
    
    cfg : Array config file
    
    high_thresh : Vis Amp above this implies no decorrelation
    low_thresh : Vis Amp below this implies perfect decorrelation 
    """
    
    ha_base,dec_base,base_len, ant1, ant2 = baseline_ha_dec_calc(cfg)

    print("Phasecenter :\t  Hour angle : %3.4f deg \t  Declination : %3.4f deg"%(ha_src,dec_src))
    print("RFI location:\t  Hour angle : %3.4f deg \t  Declination : %3.4f deg"%(ha_rfi,dec_rfi))

    rad=np.max( [np.max(np.abs(ant1)) ,np.max(np.abs(ant2)) ] ) * 1.2
    x_src, y_src = convert_ha_dec_to_xy(ha_src,dec_src,rad)
    x_rfi, y_rfi = convert_ha_dec_to_xy(ha_rfi,dec_rfi,rad)

    print("RFI speed : %3.4f deg/sec"%(speed_rfi))
    print("\nChannel width : %3.4f kHz"%(bw/1e+3))
    print("Array Config : %s"%(cfg))

    vis_amp_config,uvdist = calc_attenuation(ha_src=ha_src, dec_src=dec_src,
                                             ha_rfi=ha_rfi, dec_rfi=dec_rfi,
                                             ha_base=ha_base,
                                             dec_base=dec_base,
                                             b = base_len,
                                             rfi_speed = speed_rfi, ## 0.00418019,  # Earth rot vel (rad/sec)
                                             tau_int = tau,  #sec
                                             beta = bw)

    #base_len = uvdist

    base_high, ant1_high, ant2_high, amp_high = get_baseline_subset( base_len, ant1, ant2, 
                                                           vis_amp_config, high_thresh, checktype='high')
    base_low, ant1_low, ant2_low,amp_low = get_baseline_subset( base_len, ant1, ant2, 
                                                           vis_amp_config, low_thresh, checktype='low')

    print("\nBaseline Count (decorrelation) :  None : %d \t  Full : %d \t  Partial %d "%(len(base_high), len(base_low), len(base_len)-len(base_high)-len(base_low) ))

    plt.close('all')
    fig = plt.figure(figsize=(8,8),constrained_layout=True);
    gs = GridSpec(2,2, figure=fig);
    ax1 = fig.add_subplot( gs[0,:] )
    ax2 = fig.add_subplot( gs[1,0] )
    ax3 = fig.add_subplot( gs[1,1] )

    ## Plot Vis Amp
    ax1.plot(base_len, np.abs(vis_amp_config), '.-', color='lightgrey',markeredgecolor='grey',);
    ax1.plot(base_high, np.abs(amp_high), 'b.');
    ax1.plot(base_low, np.abs(amp_low), 'r.');  
    ax1.set_title('Predicted Visibility Amplitude')

    ## Plot antennas and baselines
    for i in range(0,len(ant1)):
        ax2.plot( [ant1[i][0], ant2[i][0]], [ant1[i][1], ant2[i][1]] , 'o-' ,alpha=0.1,color='grey')
    for i in range(0,len(ant1_low)):
        ax2.plot( [ant1[i][0], ant2[i][0]], [ant1[i][1], ant2[i][1]] , 'r-' ,alpha=0.3)
    for i in range(0,len(ant1_high)):
        ax2.plot( [ant1[i][0], ant2[i][0]], [ant1[i][1], ant2[i][1]] , 'b-' ,alpha=0.3)
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
    #Eqn 2
    lx = b * np.cos(np.radians(ha_base))*np.cos(np.radians(dec_base))
    ly = -1.0 * np.array(b) * np.sin(np.radians(ha_base))*np.cos(np.radians(dec_base))
    lz = b *np.cos(np.radians(ha_base))+np.sin(np.radians(dec_base))
        
    u = lx*np.sin(np.radians(ha_src)) + ly*np.cos(np.radians(ha_src))
    v = lz * np.cos(np.radians(dec_src)) - lx*np.sin(np.radians(dec_src))*np.cos(np.radians(ha_src)) + ly*np.sin(np.radians(dec_src))*np.sin(np.radians(ha_src))
        
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
    lon= -107.6184 #deg
    
    E= 5.4 #deg #elevation
    Azi= 22#186.2 #deg,Azimuth
    
    declination_sat = (np.degrees(np.arcsin(np.sin(np.radians(E))*np.sin(np.radians(lat))+ np.cos(np.radians(E))*np.cos(np.radians(lat))*np.cos(np.radians(Azi)))))
    RA_sat = np.degrees(np.arcsin(-np.sin(np.radians(Azi))*np.cos(np.radians(E))/np.cos(np.radians(declination_sat))))
    
    ha_rfi = RA_sat  ## Nope !
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
