import matplotlib.pyplot as plt
import numpy as np
import math
import random
from random import seed
from random import gauss
from pandas import Series  
from pandas.plotting  import autocorrelation_plot
from matplotlib import pyplot as plt


def sample_calc(Bandwidth, F_obs, b,o_fact):
    
    #Calculating Samplesize
    
    #b - baseline length (m)
    #c - speed of light (m/s)
    #F_obs = Center Frequency Observed (Hz)
    #Fs = Sampling rate 
    c = 3e8
    samplesize = round(1000* F_obs*(b/c)) #unitless ,later only pick shorter baseline len, b = 100 m as maximum
    #1.2- arbitrary
    
    actual_samplesize = round(200 * samplesize) #5 was the prev value
    print('actual samplesize', actual_samplesize)
    
    #sampling rate
    #o_fact = 250
    dt = 1/(o_fact*F_obs)  # in seconds, 1000 is the oversampling factor, 400
    print('dt',dt,'s')
    return samplesize,dt,actual_samplesize
    

def sig_gen(Bandwidth, F_obs, b, dt, samplesize):
    
    
    trange = np.zeros(samplesize, 'float')

    Time = np.arange(0.0,samplesize) * dt #in seconds

    Signal = np.cos(2*np.pi* Time *  F_obs) # in radians
    
    # to randomize phase in the sine signal add 'random.uniform(0,np.pi*dt*freq*360) )'
    
    return Signal,Time

def mul_sig_gen(Bandwidth, F_obs,F_end,F_del, b,o_fact):
    
    samplesize,dt,actual = sample_calc(Bandwidth, F_obs, b,o_fact)
    
    SigB =  np.zeros(actual, 'float')
    nfreq = 0
    for ff in np.arange(F_obs,F_end,F_del):
       sig3 , time3 = sig_gen(Bandwidth, ff, b, dt, actual)
       SigB = SigB + sig3 #*np.random.uniform(len(sig3))
       nfreq = nfreq+1
    SigB = SigB/nfreq
    
    return SigB,time3,samplesize,dt,actual


def tot_phase(F_obs, b, Bandwidth, theta):
    #Calulating total phase delay
    c = 3e8 # in m/s
    D_t = np.degrees(((b *np.sin(np.radians(theta))*np.cos(np.radians(phi)) * F_obs * 2*np.pi) / c))  #in degrees

        
    return D_t

def tot_phase_2(F_obs, b, Bandwidth, theta, phi,v_rfi):
    #Calulating total phase delay
    c = 3e8 # in m/s
    D_t = np.degrees(((b *np.sin(np.radians(theta))*np.sin(np.radians(phi)) * F_obs *v_rfi* 2*np.pi) / c))  #in degrees

        
    return D_t
    
    
def tshift_calc(Signal, F_obs, b, Bandwidth, theta ,dt):
      
    
    #Calculating the shifted samples version 2
    c = 3e8 # in m/s
    lam = c/F_obs # Wavelength in m
    
    D_t = tot_phase(F_obs, b, Bandwidth, theta) #change phase to calling it theta
        
    lambda_frac = (D_t)/360  #unlitless

    tshift = int( (lambda_frac *lam)/(c*dt)) #unitless
    
    #print('tshift',tshift)
    
    #Calculating samples to shift for the imaginary magnitude 
    
    imtshift =  int( (((D_t+90)/360) *lam/c)/dt)#+ int( ((90/360) *lam/c)/dt) #unitless
    
    #print('imtshift', imtshift)
    

    return tshift, imtshift

def shift_sigs(Signal,samplesize, tshift,imtshift,actualsamplesize):
    if (imtshift)>0:
       # sig_orig = Signal[0:samplesize]
       # sig_real= Signal[tshift:samplesize+tshift] 
       # sig_im = Signal[imtshift:samplesize+imtshift] 
        sig_orig = Signal[int(actualsamplesize/2)+50:int(actualsamplesize/2)+ samplesize] 
        sig_real= Signal[int(actualsamplesize/2)+tshift:int(actualsamplesize/2)+ samplesize+tshift] 
        sig_im= Signal[int(actualsamplesize/2)+imtshift:int(actualsamplesize/2)+ samplesize+imtshift]
 
    elif (imtshift)< 0 :  
        sig_orig = Signal[int(actualsamplesize/2):int(actualsamplesize/2)+ samplesize] 
        sig_real= Signal[int(actualsamplesize/2)+tshift:int(actualsamplesize/2)+ samplesize+tshift] 
        sig_im= Signal[int(actualsamplesize/2)+imtshift:int(actualsamplesize/2)+ samplesize+imtshift]
    else : 
        sig_orig = Signal[int(actualsamplesize/2):int(actualsamplesize/2)+ samplesize] 
        sig_real= Signal[int(actualsamplesize/2)+tshift:int(actualsamplesize/2)+ samplesize+tshift] 
        sig_im= Signal[int(actualsamplesize/2)+imtshift:int(actualsamplesize/2)+ samplesize+imtshift]
    
    if (tshift)>0:
       # sig_orig = Signal[0:samplesize] 
       # sig_real= Signal[tshift:samplesize+tshift] 
       # sig_im = Signal[imtshift:samplesize+imtshift] 
        sig_orig = Signal[int(actualsamplesize/2):int(actualsamplesize/2)+ samplesize] 
        sig_real= Signal[int(actualsamplesize/2)+tshift:int(actualsamplesize/2)+ samplesize+tshift] 
        sig_im= Signal[int(actualsamplesize/2)+imtshift:int(actualsamplesize/2)+ samplesize+imtshift]
 
    elif (tshift)< 0 :  
        sig_orig = Signal[int(actualsamplesize/2):int(actualsamplesize/2)+ samplesize] 
        sig_real= Signal[int(actualsamplesize/2)+tshift:int(actualsamplesize/2)+ samplesize+tshift] 
        sig_im= Signal[int(actualsamplesize/2)+imtshift:int(actualsamplesize/2)+ samplesize+imtshift]
    else : 
        sig_orig = Signal[int(actualsamplesize/2):int(actualsamplesize/2)+ samplesize] 
        sig_real= Signal[int(actualsamplesize/2)+tshift:int(actualsamplesize/2)+ samplesize+tshift] 
        sig_im= Signal[int(actualsamplesize/2)+imtshift:int(actualsamplesize/2)+ samplesize+imtshift]
 
    return sig_orig,sig_real,sig_im 

def mac(a1,a2, step, samplesize):
    
    corrs = np.zeros(samplesize, 'float')
    stddev = np.zeros(samplesize, 'float')
    corrs = []
    stddev= []
    
    for n in range(0,samplesize,step):

        multiply =  (a1[n:n+step]) * (a2[n:n+step]) #multiplying the phase shifted signals
        
        avg = np.nanmean(multiply)

        corrs.append(avg) 

        std = np.nanstd(multiply)
        
        stddev.append(std)
    
     
    return corrs, stddev

def correlate(a1,a2,a3,step,chanfreq,bandwidth,samplesize):
    
    #real-cosine
    corr_r,sd_r = mac(a1,a2, step, samplesize)
    meanr = np.nanmean(np.array(corr_r))
    
    #imaginary sine
    corr_i,sd_i = mac(a1,a3, step, samplesize) 
    meani = np.nanmean(np.array(corr_i))
   
    complexvis = np.complex(meanr,meani)
    
    #plt.plot(corr_r)
    #plt.plot(corr_i)
    #plt.show()
    
    return meanr, meani, complexvis


def loop(samplesize, Bandwidth, F_obs,step, phinterval,startphase,endphase,b,Sig,dt,actual_samplesize):

    #phase recovery loop using Sinc Definition
    phasestored = []
    amplitude = []
    D_t_store=[]
    theta_input = []
    print('Integration Time(s)',dt*step)
    #print('Coherence Time(s)',(1/F_del))



    for i in np.arange(startphase,endphase,phinterval):
        theta= i
        
        #calculating total phase for each theta
        D_t = tot_phase(F_obs, b, Bandwidth, i)
        
        #calculating total angle and the shifted samples        
        tshift,imtshift= tshift_calc(Sig, F_obs, b, Bandwidth, i,dt)

        #Getting the signals- original, shifted with given phase, shfited with given phase + 90deg

        orig,real,im = shift_sigs(Sig,samplesize, tshift,imtshift,actual_samplesize)
        
        
        #correlation
        realvis,imagvis,compvis = correlate(orig,real,im,step,F_obs,Bandwidth,samplesize) 
     
        #phase and amplitude output 
        visamp = np.abs(compvis)
        outphase = np.angle(realvis-imagvis*1j,deg=True) +(360 * np.round( D_t/360))

        phasestored.append(outphase)
        amplitude.append(visamp)
        D_t_store.append(D_t)
        theta_input.append(i)
    return theta_input,phasestored, amplitude,D_t_store

def sinc_theory(b,theta,F_del,F_obs): 

    c = 3e8    
    sinc_it = np.sinc(((b *np.sin(np.radians(theta))* F_del * np.pi) / c))  
    cos_factor = np.cos(((b *np.sin(np.radians(theta)) * F_obs * 2 * np.pi) / c))
    theory = (cos_factor) * (sinc_it)
    y =np.abs(theory) 
 
    return y

def loop_base(Bandwidth, F_obs,step, start,end,interval,theta,F_end,F_del,o_fact):

    #phase recovery loop using Sinc Definition
    phasestored = []
    amplitude = [] 
    D_t_store=[]
    baselines= []
   
    
    Sig,time3,samplesize,dt,actual=  mul_sig_gen(Bandwidth, F_obs,F_end,F_del, end,o_fact)

    print('Integration Time(s)',dt*step)
    print('Coherence Time',(1/F_del))

    
    for i in np.arange(start,end,interval):
        #Sig,time3,samplesize,dt,actual=  mul_sig_gen(Bandwidth, F_obs,F_end,F_del, i)

        #calculating total phase for each theta only to be recorded in later step
        D_t = tot_phase(F_obs, i, Bandwidth, theta)

        #calculating total angle and the shifted samples        
        tshift,imtshift= tshift_calc(Sig, F_obs, i, Bandwidth, theta,dt)
        
        #Getting the signals- original, shifted with given phase, shfited with given phase + 90deg

        orig,real,im = shift_sigs(Sig,samplesize, tshift,imtshift,actual)
        
        
        #correlation
        realvis,imagvis,compvis = correlate(orig,real,im,step,F_obs,Bandwidth,samplesize) 
     
        #phase and amplitude output 
        visamp = np.abs(compvis)
        outphase = np.angle(realvis-imagvis*1j,deg=True) +(360 * np.round( D_t/360))

        phasestored.append(outphase)
        amplitude.append(visamp)
        D_t_store.append(D_t)
        b = i
        baselines.append(b)
    
    baselines_array = np.asarray(baselines)
    phasestored_array = np.asarray(phasestored)
    amplitude_array = np.asarray(amplitude)
    D_t_store_array = np.array(D_t_store)
    
    return phasestored_array, amplitude_array,D_t_store_array,baselines_array


def dofft(asig):
    return np.fft.ifftshift(np.fft.fft(np.fft.fftshift(asig)))

def doifft(asig):
    return np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(asig)))

def chan_filter(a1, atime, Bandwidth, F_obs, chanwidth,doplot=False):
    ## Take FFT
    fa1 = dofft(a1)
    nrt = len(a1)
    
    ## Get frequency units 
    dnu = 1.0/(2*np.max(atime))  # in hz
    print('dnu',dnu)
    nu = np.zeros(nrt)
    mid = int(nrt/2)
    nu = np.arange(-1*dnu*mid, dnu*mid, dnu)
    
    off = int(F_obs/dnu + 0.5)
    wid =  int(0.5*chanwidth/dnu + 0.5)#int(chanwidth)
    
    
    print('Filter Width',np.abs(mid-off-wid+1-(mid-off+wid+1)))
    print('Coherence Time',np.abs(1/(mid-off-wid+1-(mid-off+wid+1))))


    ffa1 = np.zeros(nrt,'complex')
     
    ffa1[mid-off-wid+1:mid-off+wid+1]=fa1[mid-off-wid+1:mid-off+wid+1]
    ffa1[mid+off-wid:mid+off+wid]=fa1[mid+off-wid:mid+off+wid]

    if doplot==True:
        pl.figure(figsize=(20,15))
        pl.subplot(211)
        #print('Displaying only frequency ranges upto +/- 3e+5')
        #pl.plot(nu[mid-disprange: mid+disprange], np.abs(ffa1[mid-disprange :mid+disprange]), label='Filtered FT of signal')
        pl.plot(nu, np.abs(ffa1), label='Filtered FT of signal')
        pl.legend(fontsize=20)

    # Take iFFT
    iffa1 = doifft(ffa1)

    if doplot==True:
        pl.subplot(212)
        pl.plot(atime, np.real(iffa1), label='Filtered signal - real')
        pl.plot(atime, np.imag(iffa1), label='Filtered signal - imag')
        pl.legend(fontsize=20)

    filt_a1 = np.real(iffa1)
   
    return filt_a1

def make_signal(actual_samplesize, bandwidth,dt, noise_sfd, sig_sfd, F_obs,chanwidth):

    atime = np.arange(0,actual_samplesize) * dt  
    
    #seed random number generator
    seed(1)
    ## Noise 
    series = [gauss(0.0,noise_sfd) for i in range(actual_samplesize)]
    noise1 = Series(series)
    series2 = [gauss(0.0,noise_sfd) for i in range(actual_samplesize)]
    noise2 = Series(series2)
   
    ## Signal
    series3 = [gauss(0.0,sig_sfd) for i in range(actual_samplesize)]
    sig = Series(series3)


    ## Pass the noise and signal through a baseband filter.
    bb_sig = chan_filter(sig, atime, bandwidth, F_obs, chanwidth,doplot=False)
    bb_noise1 = chan_filter(noise1, atime, bandwidth,F_obs, chanwidth,doplot=False)
    bb_noise2 = chan_filter(noise2, atime, bandwidth, F_obs, chanwidth,doplot=False)
    
    ## Scale the noise 
    bb_sig = bb_sig * sig_sfd * bandwidth
    bb_noise1 = bb_noise1 * noise_sfd * bandwidth
    bb_noise2 = bb_noise2 * noise_sfd * bandwidth

    
    return bb_sig, bb_noise1,bb_noise2 

def loop_base_white_noise(Bandwidth, F_obs,step, start,end,interval,theta,F_del,sig_sfd, noise_sfd,o_fact):

    #phase recovery loop using Sinc Definition
    phasestored = []
    amplitude = [] 
    D_t_store=[]
    baselines= []
    
    samplesize,dt,actual = sample_calc(Bandwidth, F_obs, end,o_fact)


    Sig,noise1,noise2=make_signal(actual, Bandwidth,dt, noise_sfd, sig_sfd, F_obs,F_del)

    print('Integration Time',step*dt)


    for i in np.arange(start,end,interval):
        #Sig,time3,samplesize,dt,actual=  mul_sig_gen(Bandwidth, F_obs,F_end,F_del, i)


        b= i
        #calculating total phase for each theta only to be recorded in later step
        D_t = tot_phase(F_obs, i, Bandwidth, theta)

        #calculating total angle and the shifted samples        
        tshift,imtshift= tshift_calc(Sig, F_obs, i, Bandwidth, theta,dt)
        


        #Getting the signals- original, shifted with given phase, shfited with given phase + 90deg

        orig,real,im = shift_sigs(Sig,samplesize, tshift,imtshift,actual)
        
        
        #correlation
        realvis,imagvis,compvis = correlate(orig,real,im,step,F_obs,Bandwidth,samplesize) 
     
        #phase and amplitude output 
        visamp = np.abs(compvis)
        outphase = np.angle(realvis-imagvis*1j,deg=True) +(360 * np.round( D_t/360))

        phasestored.append(outphase)
        amplitude.append(visamp)
        D_t_store.append(D_t)
        baselines.append(b)
        
    baselines_array = np.asarray(baselines)
    phasestored_array = np.asarray(phasestored)
    amplitude_array = np.asarray(amplitude)
    D_t_store_array = np.array(D_t_store)
    
    return phasestored_array, amplitude_array,D_t_store_array,baselines_array

#Baseline Calculator - use text files
def baseline_calc(filename):

    x,y,z= np.loadtxt(filename,dtype= 'float',usecols=(0,1,2),unpack = True)
    #Averaging each column to get a center
    cx = np.mean(x)
    cy = np.mean(y)
    cz = np.mean(z)
    #print('Averages',cx,cy,cz)

    x = x- cx
    y = y - cy
    z = z- cz
    m = np.prod(x.shape) #number of antennas in the file
    base_list = []

    for i in np.arange(0,m):
        dx = x[i]-x[:]   
        dy = y[i]-y[:]   
        dz = z[i]-z[:]   

        #calculating baselines
        base = np.sqrt(dx**2+dy**2+dz**2)
        # print('Baseline',base)
        base_list.append(base)
    
    base_list = np.asarray(base_list)
    w = np.prod(base_list.shape)
    reshaped_baselines = np.reshape(base_list,(w,))
     
    return reshaped_baselines

def loop_base_white_noise_version2(Bandwidth, F_obs,step, start,end,interval,theta,F_del,sig_sfd, noise_sfd,baseline_set,o_fact):

    #phase recovery loop using Sinc Definition
    phasestored = []
    amplitude = [] 
    D_t_store=[]
    baselines= []
    
    samplesize,dt,actual = sample_calc(Bandwidth, F_obs, end,o_fact)


    Sig,noise1,noise2=make_signal(actual, Bandwidth,dt, noise_sfd, sig_sfd, F_obs,F_del)

    print('Integration Time',step*dt)

    for i in np.arange(start,end,interval):
        #Sig,time3,samplesize,dt,actual=  mul_sig_gen(Bandwidth, F_obs,F_end,F_del, i)


        b= baseline_set[i]
        #calculating total phase for each theta only to be recorded in later step
        D_t = tot_phase(F_obs, b, Bandwidth, theta)

        #calculating total angle and the shifted samples        
        tshift,imtshift= tshift_calc(Sig, F_obs, b, Bandwidth, theta,dt)
        


        #Getting the signals- original, shifted with given phase, shfited with given phase + 90deg

        orig,real,im = shift_sigs(Sig,samplesize, tshift,imtshift,actual)
        
        
        #correlation
        realvis,imagvis,compvis = correlate(orig,real,im,step,F_obs,Bandwidth,samplesize) 
     
        #phase and amplitude output 
        visamp = np.abs(compvis)
        outphase = np.angle(realvis-imagvis*1j,deg=True) +(360 * np.round( D_t/360))

        phasestored.append(outphase)
        amplitude.append(visamp)
        D_t_store.append(D_t)
        baselines.append(b)
        
    baselines_array = np.asarray(baselines)
    phasestored_array = np.asarray(phasestored)
    amplitude_array = np.asarray(amplitude)
    D_t_store_array = np.array(D_t_store)
    
    return phasestored_array, amplitude_array,D_t_store_array,baselines_array

def smoothing(inputarray,step):
    
    samplesize = len(inputarray)
    newarr = np.zeros(samplesize, 'float')
    newarr=[]
    for n in range(0,samplesize,step):
        
        avg = np.nanmean(inputarray[n:n+step])

        newarr.append(avg) 
        
        
    return newarr

