# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 08:10:42 2021

@author: G. Hovey, 10 Dec 2021;Onsala Space Observatory, Sweden
@revised: G. Hovey, 14 Apr 2022; changed to use rfiMonFileIO to get data
@revised: G. Hovey, 09 Nov 2022; Modified to plot CXA rfi spectra
@modified: F. DiVruno, 21/04/2023
@modified: H. Qiu, 05/2024 modified to loading script

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import argparse
import os
from datetime import datetime
import cysgp4
import astropy.units as u
from astropy import time

import glob
from astropy.table import vstack
from astropy.io import ascii
from astropy.table import vstack,Table


from rfiMonFileRO import RfiMonFileIO
import argparse

font = {'family' : 'DejaVu Sans','weight' : 'normal','size'   : 10}
import matplotlib
matplotlib.rc('font', **font)



bytes2Str = lambda xBytes : xBytes.decode('UTF-8')

bytes2Datetime = lambda xBytes : datetime.strptime(bytes2Str(xBytes), 
                            "%Y-%m-%dT%H:%M:%S.%f")

bytes2mdate = lambda xBytes : mdates.date2num(bytes2Datetime(xBytes))

bytes2mdates = np.vectorize(bytes2mdate)



# main data loading script
import argparse

def main():
    parser = argparse.ArgumentParser(description='load file location')
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                        help='a list of files to process')

    args = parser.parse_args()

    for filename in args.files:
        print(f'Processing file: {filename}')


def bytes2astrotime(utctime):
    return  [time.Time(bytes2Datetime(utctime[i])) for i in range(len(utctime))]
    




def logLinLogMean(data_dBm, axis=0):
    dataLin = 10.**(data_dBm/10.)
    aveLin = dataLin.mean(axis=axis)
    return 10.*np.log10(aveLin)


def getMonData(filespec):
    global fileId
    
    fileId = RfiMonFileIO(None)
    fileId.filespec = filespec
    f_GHz = fileId.rdField('f_GHz')[0]
    ave_dBm = fileId.rdField('ave_dBm')
    max_dBm = fileId.rdField('max_dBm')
    min_dBm = fileId.rdField('min_dBm')
    tUtcBytes = fileId.rdField('datetimeUTC')
    az_deg = fileId.rdField('az_deg')
    el_deg = fileId.rdField('el_deg')
    mjd = fileId.rdField('mjd')
    satTles = fileId.rdField('satTles')
    
    hdrD = fileId.rdHdr()

    return ave_dBm[:], max_dBm[:], min_dBm, f_GHz, tUtcBytes[:], az_deg, el_deg, mjd, satTles, hdrD

def getCorrData(filespec):
    global fileId
    
    fileId = RfiMonFileIO(None)
    fileId.filespec = filespec
    aveX = fileId.rdField('antX/ave')
    az_degX = fileId.rdField('antX/az_deg')
    el_degX = fileId.rdField('antX/el_deg')
    mjdX = fileId.rdField('antX/mjd')
    aveY = fileId.rdField('antY/ave')
    az_degY = fileId.rdField('antY/az_deg')
    el_degY = fileId.rdField('antY/el_deg')
    mjdY = fileId.rdField('antY/mjd')
    XYave = fileId.rdField('xy/ave')
    
    tUtcBytes = fileId.rdField('datetimeUTC')
    satTles = fileId.rdField('satTles')
    f_GHz = fileId.rdField('f_GHz')[0]
    
    hdrD = fileId.rdHdr()

    return aveX[:], f_GHz, tUtcBytes[:], az_degX, el_degX, mjdX,\
             satTles, aveY[:], az_degY, el_degY, mjdY,\
             XYave, hdrD


def closeFile():
    fileId.closeFile()

def downSample(arry, func, byKrows):
    shape = arry.shape
    nRows = shape[0]
    remainder = nRows%byKrows

    nD = len(shape)
    nCols = 1
    if 1 < nD:
        nCols = shape[1]
        newArry = func(arry[:-remainder].reshape(-1,byKrows,nCols), axis=1)
        return np.append(newArry, [func(arry[-remainder:], axis=0)], axis=0)

    nCols = 1
    newArry = func(arry[:-remainder].reshape(-1,byKrows), axis=1)
    return np.append(newArry, [func(arry[-remainder:], axis=0)], axis=0)




def slicer(arr, t, vmin, vmax, ax = 0 ):
    '''
    Slice the 2D matrix arr between Vmin and vmax in the ax axis

    Parameters
    ----------
    arr : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.
    vmin : TYPE
        DESCRIPTION.
    vmax : TYPE
        DESCRIPTION.
    ax : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    arr_new : TYPE
        DESCRIPTION.
    t_new : TYPE
        DESCRIPTION.

    '''
    print('Slicing data...')
    idx = (vmin<t)*(t<vmax)
    idx = np.where(idx)
    t_new = t[idx]
    if ax==0:
        arr_new = arr[idx[0],:]
    else:
        arr_new = arr[:,idx[0]]
    print('Sliced data from %.2f to %.2f'%(vmin,vmax))
    print('Original size: ' + str(arr.shape))
    print('New size: ' + str(arr_new.shape))
    
    return arr_new, t_new


def binner(arr, t, t_interval, ax = 0, bin_method = 'max'):
    '''
    bin a 2 dimensional matrix in one dimmension 
    in the time_inteval 
    

    Parameters
    ----------
    arr : TYPE
        DESCRIPTION.
    t : TYPE
        DESCRIPTION.
    t_interval : TYPE
        DESCRIPTION.
    bin_method : string, 'max' or 'avg'
        DESCRIPTION. The default is 'max'.

    Returns
    -------
    None.

    '''
    if ax==0:
        if (arr.shape[0] != len(t)):
            print('dimmensions 0 are not equal!')
            return -1,-1
    if ax==1:
        if (arr.shape[1] != len(t)):
            print('dimmensions 1 are not equal!')
            return -1,-1
    
    print('Starting binning...')
    ind0 = 0 #original array index
    to = t[0]
    arr_bin = [] #np.array([t//t_interval,arr.shape[1]]) #new binned array
    t_new = []
    
    if bin_method == 'max':
        method = np.nanmax
    if bin_method == 'min':        
        method = np.nanmin
    if bin_method == 'mean':
        method = np.nanmean
    
    if ax == 0:
        print('Axis 0')
        for i in range(t.shape[0]):
            if (t[i] >= to+t_interval):
                to = t[i]
                t_new.append(t[i]-t_interval/2)
                arr_bin.append(method(arr[ind0:i,:],axis=0))
                ind0 = i #update the starting index for the original array
                
                
        arr_bin = np.array(arr_bin)
        
    if ax == 1:
        print('Axis 1')
        for i in range(t.shape[0]):
            if (t[i] >= to+t_interval):
                to = t[i]
                t_new.append(t[i]-t_interval/2)
                arr_bin.append(method(arr[:,ind0:i],axis=1).T)
                ind0 = i #update the starting index for the original array
        arr_bin = np.array(arr_bin).T
    
    t_new = np.squeeze(np.array(t_new))
    
    return arr_bin, t_new

def pow_per_channel(data, freq, f0, delta_f,N,axis):
    Plist = []
    flist = []
    for i in range(N):
        P1, f1 = slicer(data,freq,f0+i*delta_f,\
                        f0+(i+1)*delta_f,ax=axis)
        P1_tot = 10*np.log10(np.nansum(10**(P1/10),axis=axis))
        Plist.append(P1_tot)
        flist.append(np.mean(f1))
    return Plist, flist


def propagate_satellites(const_name,TLE_type,obs_mjds):
    '''

    Parameters
    ----------
    const_name : TYPE
        'starlink' or 'oneweb' if TLE_type is 'SUP'
        'starlink','oneweb','active' or 'geo' if TLE_type is 'GP'
    TLE_type : String
        'SUP' or 'GP'
    obs_mjds: numpy array
        time steps of the observation in mjd
        

    Returns
    -------
    result: output of the cysgp4 function "propagate_many"
    noradID:
    satname:

    '''
    tle_folder = r'C:\Users\f.divruno\Dropbox (SKAO)\TLEs\TLEs'
    
    if TLE_type=='SUP':
        if const_name == 'starlink':
            const = '\SUP\starlink'
        elif const_name == 'oneweb':
            const = '\SUP\oneweb'
        else:
            print('Error in constellation name')
            return 0,0

        tle_files = sorted(glob.glob(tle_folder +const+ r'\*.npz'))
        
        #find which file to load using the mjd1 of the capture
        files_datetime = []
        for f in tle_files:
            DT = str.split(str.split(f,'\\')[-1],'.')[0]
            files_datetime.append(datetime.strptime(DT,"%Y%m%d_%H%M%S"))
        files_datetime = np.array(files_datetime)    
     
        # use chunks of 4 hours to find TLEs closer to the time:
        meas_datetime = time.Time(obs_mjds[0],format='mjd')
        
        # index in the files array
        id_tle = np.where(meas_datetime>=files_datetime)[0][-1] 
                  
        
        TLEs = np.load(tle_files[id_tle], allow_pickle=True)['arr_0'].item().text
        

    if TLE_type == 'GP': #load all TLEs
        tle_files = sorted(glob.glob(tle_folder +'\GP'+ r'\*.npz'))
        
        #find which file to load using the mjd1 of the capture
        files_datetime = []
        for f in tle_files:
            DT = str.split(str.split(f,'\\')[-1],'.')[0]
            files_datetime.append(datetime.strptime(DT,"%Y%m%d"))
        files_datetime = np.array(files_datetime)    

        # use chunks of 4 hours to find TLEs closer to the time:
        meas_datetime = time.Time(obs_mjds[0],format='mjd')
        
        # index in the files array
        id_tle = np.where(meas_datetime>=files_datetime)[0][-1] 
                  
        
        TLEs = np.load(tle_files[id_tle],\
                        allow_pickle=True)['arr_0'].item().text
    
        meas_datetime = time.Time(obs_mjds[0],format='mjd')
        id_tle = np.where(meas_datetime>=files_datetime)[0][-1] # index in the files array
        
        aux = np.load(tle_files[id_tle], allow_pickle=True)
        
        TLEs_GP_SL = aux['arr_0'].item().text
        TLEs_GP_OW = aux['arr_1'].item().text
        TLEs_GP_AC = aux['arr_2'].item().text
        TLEs_GP_GEO = aux['arr_3'].item().text
        
        if const_name == 'starlink':
            TLEs = TLEs_GP_SL
        elif const_name == 'oneweb':
            TLEs = TLEs_GP_OW
        elif const_name == 'active':
            TLEs = TLEs_GP_AC
        elif const_name == 'geo':
            TLEs = TLEs_GP_GEO
        else:
            print('Error in GP constellation name')
            return 0,0,0
    
    #generate tles with cysgp4
    tles = np.array(cysgp4.tles_from_text(TLEs))
                
    
    # check all tles files have the same number of satellites
    
    
    #satellite noradIDs
    noradID = np.array([tles[k].catalog_number for k in range(len(tles))])
    satname = str.split(TLEs,'\r\n')[0::3]
    
    # ONSALA observatory as observer:
    obs_name = 'Onsala'
    observer = cysgp4.PyObserver( 11.917778,57.393056, 0.02) # long,lat,alt, [deg, deg, km]

    print('Propagating '+const_name+' '+TLE_type+
          ' tles for %d timesteps starting at MJD= %.3f'%
          (len(obs_mjds),obs_mjds[0]))
    #propagate satellites in the obs_mjds time steps
    # use the chunks of the TLEs
    
    result = cysgp4.propagate_many(
        obs_mjds[:, np.newaxis],
        tles[np.newaxis,:],
        observer,
        do_sat_azel=False,
        on_error = 'coerce_to_nan')
    
    # topo_pos_az, topo_pos_el, topo_pos_dist, _ = [result['topo'][..., j] for j in range(4)]
    # topo_pos_az = (topo_pos_az + 180.) % 360. - 180.    
    
#    vis_mask = topo_pos_el>0
    
    print('\nDone!')
    print('%d satellites'%result['topo'][...,0].shape[1])
    print('%d timesteps'%result['topo'][...,0].shape[0])
    result['noradID'] = noradID
    result['satname'] = satname
    result['epochs'] = obs_mjds
    return result

def plot_waterfall(data, t, freq, title='', savefolder='' ,savename='', savefig=False):
    '''
    function to plot a waterfall plot    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    time : TYPE
        DESCRIPTION.
    freq : TYPE
        DESCRIPTION.
    title : TYPE
        DESCRIPTION.
    savefolder : TYPE
        DESCRIPTION.
    savename : TYPE
        DESCRIPTION.
    savefig : boolean
        flag to save the figure

    Returns
    -------
    None.

    '''                    
    vmax = np.nanpercentile(data,99.9)
    vmin = np.nanpercentile(data,10)
    
    fig, ax = plt.subplots(2,2,figsize=(20, 16),
                            gridspec_kw={'height_ratios': [1, 0.25], 'width_ratios':[1,0.25]})
    
    ax[1,0].sharex(ax[0,0])
    ax[0,1].sharey(ax[0,0])
    
    img = ax[0,0].imshow(data, aspect='auto', vmin=vmin, vmax=vmax, 
                          origin='lower',extent=[freq[0], freq[-1],t[0], t[data.shape[0]-1]])
      
    ax[0,0].set_ylabel("time")
                      
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("Power (dB)", rotation=270)
    

    Pmax_f = 10*np.log10(np.nanmax(10**(data/10),axis=0))
    Pmean_f = 10*np.log10(np.nanpercentile(10**(data/10),50,axis=0))
    P99_f = 10*np.log10(np.nanpercentile(10**(data/10),99.9,axis=0))
    ax[1,0].plot(freq,Pmax_f,'-')
    ax[1,0].plot(freq,Pmean_f,'-')
    # ax[1,0].plot(freq,P99_f,'-')
    ax[1,0].set_ylabel('dBm')
    ax[1,0].set_xlabel("Frequency (GHz)")
    ax[1,0].grid('both',alpha=0.5)
    ax[1,0].legend(['Max','Mean','99.9 percentile'],fontsize=10)
    
    Pm_t = 10*np.log10(np.nanmean(10**(data/10),axis=1))
    Pmax_t = 10*np.log10(np.nanmax(10**(data/10),axis=1))
    ax[0,1].plot(Pm_t,t,'-') #np.arange(0,data.shape[0]),'-')
    ax[0,1].plot(Pmax_t,t,'-') #np.arange(0,data.shape[0]),'-')
    ax[0,1].set_xlabel("dBm")
    
    ax[0,1].grid('both',alpha=0.5)
    
    ax[0,0].set_title(title)
    
    ax[1,1].axis('off')
    
    if savefig: #if the name is provided, save the figure
        plt.savefig(savefolder+'\\'+savename, dpi=200)    
        # plt.close()

    return fig


def calibrate(data,mjds,Tel_az, Tel_el,filename='',mode='',
              plot_figs=False):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    mjds : TYPE
        DESCRIPTION.
    Tel_az : TYPE
        DESCRIPTION.
    Tel_el : TYPE
        DESCRIPTION.
    filename : TYPE, optional
        DESCRIPTION. The default is ''.
    mode : TYPE, optional
        DESCRIPTION. The default is ''.
    plot_figs : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    corr_factor : TYPE
        DESCRIPTION.

    '''
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.time import Time
    # see if any observation was near Cas A
    casA = SkyCoord.from_name('casA')

    OSO = EarthLocation(lon= 11.917778*u.deg,lat=57.393056*u.deg, height=0.02*u.km)


    epochs = Time(np.array(mjds),format='mjd')
    time_human = epochs.iso

    casAaltaz = casA.transform_to(AltAz(obstime=epochs,location=OSO))

    # find times where pointing was near Cas A:

    ang_dist = (abs(casAaltaz.az-Tel_az)**2+abs(casAaltaz.alt-Tel_el)**2)**0.5
    match =  ang_dist < 0.1*u.deg

    ind = np.where(match)
    start = [ind[0][0]]
    end = []
    
    for i in range(len(ind[0])-1):
        if (ind[0][i+1]-ind[0][i])>10:
            end.append(ind[0][i])
            start.append(ind[0][i+1])
    end.append(ind[0][-1])
    
    #start every hour
    # mjd_hour = 3600/60/60/24
    # start = np.where(((np.array(mjd)-mjd_hour) - (np.array(mjd)-mjd_hour)//mjd_hour*mjd_hour) < 2/60/60/24)[0]
    
    #start after every hour since the start of the obs
    start = np.where( ((t-3600) - (t-3600)//3600*3600) < 2)[0]
        
    print('Found %d calibration sequences'%(len(start)))
    L = 200
    for i in range(len(start)):
        plt.figure()
        t1 = t[start[i]:start[i]+L] - t[start[i]]
        plt.plot(t1, casAaltaz.az[start[i]:start[i]+L],'.')
        plt.plot(t1, casAaltaz.alt[start[i]:start[i]+L],'.')
        plt.plot(t1, Tel_az[start[i]:start[i]+L])
        plt.plot(t1, Tel_el[start[i]:start[i]+L])
        plt.title(epochs[start[i]].isot)
        plot_waterfall(data[start[i]:start[i]+L], t[start[i]:start[i]+L]-t[start[i]], freq)
        
    #the calibration in CAS A is 120 sec on source and 120 sec off source
    casA = []
    off_source_ind = int(100//(t[1]-t[0])) #
    for i  in range(len(start)):
        ind_start = max(start[i],0)
        ind_end = min(end[i]+off_source_ind,len(match))

        #if its the SouthWest, has a lot of RFI in the 10 GHz range.
        # use a smaller freq range.
        if 'South' in filename:
            i_f = (f_GHz>10.8) & (f_GHz<12) #due to RFI
        elif 'Corr' in filename:
            i_f = (f_GHz>11.96) & (f_GHz<12.2) #due to RFI
        else:
            i_f = f_GHz>0


        if plot_figs:
            fig,axs = plt.subplots(2,1,figsize = [20,15])
            axs[0].sharex(axs[1])
            axs[0].plot(epochs[ind_start:ind_end].mjd,az_deg[ind_start:ind_end],'.',label='tel pointing')
            axs[0].plot(epochs[ind_start:ind_end].mjd,casAaltaz.az[ind_start:ind_end],'.',label='Cas A')
        
            axs[1].plot(epochs[ind_start:ind_end].mjd,el_deg[ind_start:ind_end],'.',label='tel pointing')
            axs[1].plot(epochs[ind_start:ind_end].mjd,casAaltaz.alt[ind_start:ind_end],'.',label='Cas A')   
            axs[0].grid()
            axs[1].grid()
            axs[0].set_title('Azimuth')
            axs[1].set_title('Elevation')
            axs[0].set_xlabel('Samples')
            axs[1].set_xlabel('Samples')
            axs[0].set_ylabel('deg')
            axs[1].set_ylabel('deg')
            fig.legend()
            plt.title(filename)
            # plt.savefig(fig_folder+'Cal pointing - start:'+ str(start[i]) + ' -' +filename[:-5]+'.png')
            
            # normalizing with the mean in frequency domain to remove the frequency structure
            data_norm = data[ind_start:ind_end,i_f] - np.mean(data[ind_start:ind_end,i_f],axis=0)
            
            fig = plot_waterfall(data_norm,t,f_GHz[i_f],
                            title = 'Calibration on and off source',
                            savefolder = fig_folder,
                            savename = 'Cal_waterfall - start:'+ str(start[i]) + ' -'+  filename[:-5],
                            savefig = False)
    
    
        D1 = 10**(data[ind_start:ind_end-off_source_ind,i_f]/10) #dat on source
        D2 = 10**(data[ind_end-off_source_ind:ind_end,i_f]/10) #data off source
        
        
        #remove RFI - finding the mask
        D1_n = D1/np.mean(D1,axis=0)-1
        D2_n = D2/np.mean(D2,axis=0)-1
        
        mask1 = D1_n>np.std(D1_n)*3
        mask2 = D2_n>np.std(D2_n)*3
    
        # apply the mask to the original D1 and D2
        D1[mask1] = np.nan
        D2[mask2] = np.nan
        
 
        if (D2.shape[0]>0) * (D1.shape[0]>0): #check that both sequences have anything inside
            if plot_figs:
                fig = plot_waterfall(D1,
                                t,f_GHz[i_f],
                                title = 'Cal sequence nr %d No RFI'%i,
                                savefolder = fig_folder,
                                savename = '',
                                savefig = False)
    
                fig = plot_waterfall(D2,
                                t,f_GHz[i_f],
                                title = 'Cal sequence nr %d No RFI'%i,
                                savefolder = fig_folder,
                                savename = '',
                                savefig = False)
            
            
            casA.append((np.nanmean(D1) - np.nanmean(D2))) #subtracting the system noise
            # ignoring the atmosphere


    corr_factor = (400*u.Jy / (casA*u.mW/1e6/u.Hz)).to(1/u.m**2)
    
    # print('Correction Factor: %.4f %s'%(corr_factor.value,corr_factor.unit))
    return corr_factor

def rebin(x,y,bins,func=np.nanmean):
    '''
    

    Parameters
    ----------
    x : TYPE
        x coordinate of the existing points, ordered.
    y : TYPE
        y coordinate of the existing points
    bins : TYPE
        subrange of x where to apply the function.
    func : TYPE, optional
        Function to be applied. The default is np.nanmean.

    Returns
    -------
    x1 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.

    '''
    y1 = []
    x1 = []
    N = len(bins)
    
    i=0
    while i < (N-1):
        j = 1
        ind = np.where((x>=bins[i])*(x<=bins[i+j]))
        while (len(ind[0])==0)*((i+j)<N):
            ind = np.where((x>=bins[i])*(x<=bins[i+j]))
            j+=1
        if (i+j)<N:
            x1.append((bins[i]+bins[i+j])/2)
            y1.append(func(y[ind]))
            i+=j
        else:
            i=N
    return x1,y1
 
def rebin_90perc(x,y,bins,func=np.nanmean):
    '''
    

    Parameters
    ----------
    x : TYPE
        x coordinate of the existing points, ordered.
    y : TYPE
        y coordinate of the existing points
    bins : TYPE
        subrange of x where to apply the function.
    func : TYPE, optional
        Function to be applied. The default is np.nanmean.

    Returns
    -------
    x1 : TYPE
        DESCRIPTION.
    y1 : TYPE
        DESCRIPTION.

    '''
    y1 = []
    x1 = []
    N = len(bins)
    
    i=0
    while i < (N-1):
        j = 1
        ind = np.where((x>=bins[i])*(x<=bins[i+j]))
        while (len(ind[0])==0)*((i+j)<N):
            ind = np.where((x>=bins[i])*(x<=bins[i+j]))
            j+=1
        if (i+j)<N:
            x1.append((bins[i]+bins[i+j])/2)
            y1.append(np.nanpercentile(y[ind],90))
            i+=j
        else:
            i=N
    return x1,y1




# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()