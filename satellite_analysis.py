import numpy as np
import matplotlib.pyplot as plt
### import the example plot function


def readnpz(data,index=0):
    npz=np.load(data)
    files=npz.files
#     print(f"current files available{index}")   
    data=npz[files[index]]
    print(f"extracting file:{files[index]}")
    return(data)

def convertdbm(dbm):
    """convert dBm/dB units to mW/W
    """
    p=1*10**(dbm/10)
    return p ### mW
def changetodbm(p):
    """convert mw/W units to dBm/dB
    """
    dbm=10*np.log10(p)
    return dbm

class telescope:
    def __init__(self,file):
        """initiate class for single dish telescope data analysis
        Parameters
        ----------
        file : str
            File name of npz file containing the following:
            data : array/float
                Telescope data
            tel_az : array/float
                Telescope azimuth
            tel_el : array/float
                Telescope elevation
            mjd : array/float
                1-d array to record MJD time
            f_GHz : array/float
                1-d array to record frequency in GHz
        """
        data=readnpz(file',0)
        tel_az=readnpz(file',1)
        tel_el=readnpz(file',2)
        f_GHz=readnpz(file',3)
        mjd=readnpz(file',4)
        self.data=data
        self.tel_az=tel_az
        self.tel_el=tel_el
        self.mjd=mjd
        self.f_GHz=f_GHz
        self.nchan=len(f_GHz)
        self.max_epoch=len(mjd)
        self.tsamp=(mjd[1]-mjd[0])*24*3600
        print("Dataset info:")
        totaltime=(mjd[-1]-mjd[0])*24
        print(f"frequency range {np.min(f_GHz)} - {np.max(f_GHz)} GHz, total channels {self.nchan}, total time {totaltime}")
        print(f"nsamp = {self.max_epoch}, tsamp = {self.tsamp} seconds")

    def catalogue(self,file):
        """read in catalogue of sources
        Parameters
        ----------
        file : str
            File name of npz file containing the following:
            sat_az : array/float
                Source azimuth
            sat_el : array/float
                Source elevation
            sat_range : array/float
                1-d array to record distance
        """
        self.sat_az=readnpz(file,0)
        self.sat_el=readnpz(file,1)
        self.sat_range=readnpz(file,2)
        print(f"number of sources {self.sat_az.shape[1]}")
        
    def bandpass(self,len=150,exampleinterval=0,show=False):
        i=exampleinterval
        bandpass=np.mean(self.data[0:100],axis=0)
        self.bandpass=bandpass
        exampledata=self.data[i]
        stdnoise=np.std(convertdbm(exampledata)-convertdbm(bandpass))
        self.stdnoise=stdnoise
        print(self.stdnoise)
        if show:
            gs = GridSpec(2,2,hspace=0.5)
            plt.figure(figsize=(10,6))
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[1,0])
            ax3= plt.subplot(gs[:,1],projection='polar')
            ax1.set_title('averaged zenith observation data')
            # ax1.plot(f_GHz,exampledata)
            ax1.plot(self.f_GHz,bandpass)
            ax1.set_ylabel("dBm")
            ax1.set_xlabel('Frequency (GHz)')
            # ax1.set_ylim(-53,-36)
            ax2.set_title('example calibrated bandpass')

            ax2.plot(self.f_GHz,(convertdbm(exampledata)-convertdbm(bandpass))/stdnoise)
            ax2.set_ylabel("S/N")
            ax2.set_xlabel('Frequency (GHz)')
            # ax2.set_ylim(-5,20)
            ax3.set_ylim(90,20)
            ax3.scatter(self.sat_az[i],self.sat_el[i],color='seagreen') 
            ax3.scatter(self.tel_az[i],self.tel_el[i],marker='^',color='C1',s=100,alpha=0.5,label='telescope')
            ax3.set_xlabel('Azimuth (Deg)')
            # plt.ylabel('Elevation (Deg)')
            ax3.legend(loc=1)
            plt.savefig("baselinecalibration.png")
            plt.show()
    def playobs(self,length=30000,step=100,fps=10):
        livetrack(self.data,self.bandpass,self.stdnoise,self.f_GHz,self.sat_az,self.sat_el,length,step,fps)
    

def livetrack(data,bandpass,stdnoise,f_GHz,sat_az,sat_el,length=30000,step=100,fps=10)
    fig = plt.figure(figsize=(8, 6),dpi=80)
    gs = GridSpec(2,2,hspace=0.5)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])
    ax3= plt.subplot(gs[:,1],projection='polar')

    ax1.plot(f_GHz,data[0])
    ax1.set_ylim(-53,-32)
    ax1.set_ylabel("dBm")
    ax1.set_xlabel('Frequency (GHz)')

    ax2.set_title('calibrated bandpass')
    # stdnoise=np.std(10**(data[0]-bandpass))
    ax2.plot(f_GHz,(data[0]-bandpass)/stdnoise)
    ax2.set_ylabel("Log S/N")
    ax2.set_xlabel('Frequency (GHz)')

    ax3.set_ylim(90,20)
    ax3.scatter(sat_az[0], sat_el[0], c="seagreen")
    ax3.scatter(tel_az[0], tel_el[0],marker='^',color='C1',s=100,alpha=0.5,label='telescope')
    x=np.arange(1,length,step)
    ax3.legend(loc=1)

    def animate(i):
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax1.set_title('observation data')

        ax1.plot(f_GHz,data[i])
        ax1.set_ylabel("dBm")
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylim(-53,-32)
        ax2.set_title('calibrated bandpass')
        ax2.plot(f_GHz,(convertdbm(data[i])-convertdbm(bandpass))/stdnoise)
        ax2.set_ylabel("S/N")
        ax2.set_xlabel('Frequency (GHz)')
        ax2.set_ylim(-5,30)
        ax3.set_title(f"MJD={mjd[i]}")
        ax3.set_ylim(90,20)
        ax3.scatter(sat_az[i],sat_el[i],color='seagreen',label='Satellites') 
        ax3.scatter(tel_az[i],tel_el[i],marker='^',color='C1',s=100,alpha=0.5,label='Telescope pointing')
        ax3.set_xlabel('Azimuth (Deg)')
        ax3.legend(loc=1)
    ani = animation.FuncAnimation(fig,animate,repeat=False,frames=x)
    ani.save("tracking.gif", writer='Pillow',fps=fps)

    plt.show()