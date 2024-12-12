# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 16:27:39 2023

@author: f.divruno
"""

import numpy as np
import matplotlib.pyplot as plt
    
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

