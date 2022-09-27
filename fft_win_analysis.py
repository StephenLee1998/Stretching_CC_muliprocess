# Program to do fft analysis of windowed ncfs

from fft_analysis import sacs,get_std,smooth,get_fft_std
import struct
import pandas as pd
from scipy.fftpack import fft,ifft
import numpy as np
import os,glob
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter 
from obspy import read
import matplotlib.gridspec as gridspec

def cal_correlation(data1,data2):
    D1 = pd.Series(data1)
    D2 = pd.Series(data2)
    cc = D1.corr(D2)
    return(cc)

def select_fft(standard_fft,fft):
    for i in range(len(fft)):
        cc = cal_correlation(standard_fft[0:10],fft[i,0:10])
        print(cc)
        #print(max(standard_fft) , max(fft[i,:]))
        if cc < 0.98 or max(standard_fft) > 1.2*max(fft[i,:]) or max(standard_fft) < 0.2*max(fft[i,:]):
            fft[i,:]=0
    return standard_fft

def get_fft(data,flag,trange,T_band):
    sac = sacs()
    sac.readsac(data, T_band)
    sac.window(trange)
    [frq,Y] = sac.fft('k',False,flag)
    return frq,Y

def get_fft_std_win(dayrange,yearrange,flag,trange,T_band):
    fft = []
    Daynum = 0
    for year in range(yearrange[0],yearrange[1]):
        print(year)
        for iday in range(0,365):
            if year == 2015 and iday < 78:
                continue
            file = glob.glob('*day'+str(iday)+'_*_' + str(year) + '*SAC')
            if len(file) == 0:
                continue
            Daynum += 1
            print(file,Daynum)
            sac = sacs()
            sac.readsac(file[0],T_band)
            sac.window(trange)
            [frq,Y] = sac.fft('k',False,flag)
            if Daynum == 1:
                fft = Y
            else:
                fft = np.vstack((fft,Y))
    std = get_std(fft)
    std_smooth = smooth(std, 5)
    return std_smooth,frq,fft

def Days():
    Station_pairs = ['AXCC1_AXBA1','AXEC2_AXBA1','AXEC1_AXBA1','AXAS1_AXBA1','AXID1_AXBA1','AXAS2_AXBA1','AXEC3_AXBA1']
    Station_pairs = ['AXCC1_AXBA1','AXAS2_AXBA1','AXAS1_AXBA1','AXID1_AXBA1']
    #Station_pairs = ['AXEC2_AXBA1','AXEC1_AXBA1','AXEC3_AXBA1']
    fig = plt.figure(1,figsize=(20,4))
    ifig = 0


    for station_pair in Station_pairs:
        path = '/work/li_chao/work/Axial_Seamount/NCFs/Distant_Station/'+station_pair+'/Z-Z/ascii'
        print(path)

        ifig+=1
        ax = plt.subplot(2,2,ifig)
        dayrange = [1,365]
        T_band = [2,8]
        yearrange = [2015,2021]
        #win = [15,35]
        win = [15,35]
        os.chdir(path)

        #standard_file = glob.glob('/work/li_chao/work/Axial_Seamount/NCFs/Standard/DS/Z-Z/ascii/*'+station_pair.split('_')[0]+'-'+station_pair.split('_')[1]+'*.SAC')[0]
        #[frq,standard_fft] = get_fft(standard_file,'Nega_win',win,T_band)
        
    
        [std_positive,frq,fft] = get_fft_std_win(dayrange,yearrange, 'Nega_win',win,T_band)
        #select_fft(standard_fft, fft)
        
        plt.plot([-80,80],[114,114],'red',dashes=[6, 2],linewidth='2')
        plt.plot([0,0],[0,365],'k',linewidth='1')
        plt.ylim(78,150)
        plt.xlim(-50,50)
        #plt.contourf(np.arange(78,150),frq,fft.T)
        levels = np.arange(0,np.max(fft),0.001)
        
        cs = ax.contourf(np.arange(0,len(fft)),frq,fft.T,levels,cmap=plt.get_cmap('seismic'))
        

        #plt.title('Positive')
        #plt.xlabel('Days')
        plt.title('%s' % (station_pair))
        if ifig == 1 or ifig == 3 or ifig == 5 or ifig == 7 :
            plt.ylabel('Frequency (Hz)')
        if ifig == 2 or ifig == 4 or ifig == 6:
            plt.yticks([])
        if ifig < 3:
            plt.xticks([])

        #plt.plot([114,114],[-10,10],'white')
        #plt.plot([99,99],[-10,10],'white',dashes=[6,2])
        plt.plot([287,287],[-10,10],'white',dashes=[6,2])
        plt.plot([653,653],[-10,10],'white',dashes=[6,2])
        plt.plot([1018,1018],[-10,10],'white',dashes=[6,2])
        plt.plot([1383,1383],[-10,10],'white',dashes=[6,2])
        plt.xlim(0,len(fft))
        plt.ylim(0.1,0.6)

    #cbar = fig.colorbar(cs)
    
    
    print(len(fft))
    path = '/work/li_chao/work/Axial_Seamount/process_code/Distant_Station/code'
    plt.savefig(path+'/fft_nega_normalized.png',dpi=1000)
    


if __name__ == "__main__":
    Days()
