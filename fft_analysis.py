# Program to do fft analysis of ncfs

import struct
from scipy.fftpack import fft,ifft
import numpy as np
import os,glob
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter 
from obspy import read

plt.switch_backend('agg')

def get_std(data):
    std = []
    for iday in range(len(data[1,:])):
        std.append(np.std(data[:,iday]))
    return std

def smooth(data,N):
    n = len(data)
    data_smooth = np.zeros(n)
    for i in range(n):
        if i-N < 0:
            data_smooth[i] = sum(data[0:i+N])/len(data[0:i+N])
        elif i+N > n-1:
            data_smooth[i] = sum(data[i-N:n])/len(data[0:i+N])
        else:
            data_smooth[i] = sum(data[i-N:i+N])/len(data[i-N:i+N])
    return data_smooth


class sacs:
    def window(self,t_win):
        t = np.linspace(-self.half_len*self.delta,self.half_len*self.delta,self.npts)
        iwp = np.where((t>t_win[0])&(t<t_win[1]))
        iwn = np.where((t<-t_win[0])&(t>-t_win[1]))
        self.winp = self.data[iwp]
        self.winn = self.data[iwn]
        
        #plt.plot(t,self.data/max(abs(0.7*self.data))+day,'k',linewidth='1')
        
        #plt.plot(t[iwp],self.winp/max(0.7*abs(self.data))+day,'r',linewidth='1')
        #plt.plot(t[iwn],self.winn/max(0.7*abs(self.data))+day,'r',linewidth='1')

    def sides(self):
        #print(t_win[0],t_win[1])
        t = np.linspace(-self.half_len*self.delta,self.half_len*self.delta,self.npts)
        iP = np.where(t>=0)
        iN = np.where(t<=0)
        dP = self.data[iP]
        dN = self.data[iN]

        self.DataP = dP
        self.DataN = dN
        #plt.plot(t[iP],dP,'k',linewidth='1')
        #plt.plot(t[iN],dN,'r',linewidth='1')

        

    def readsac(self,sacfile,T_band):
        tr = read(sacfile)
        st = tr[0]

        ### Convert cf data to egf data ###
        self.half_len = (st.stats.sac.npts-1)/2
        self.npts = st.stats.sac.npts
        self.delta = st.stats.sac.delta

        cf_positve = st.data[self.half_len:self.npts]
        egf_positive = -np.diff(cf_positve)/self.delta
        data_revese = st.data[::-1]
        cf_negitve = data_revese[self.half_len:self.npts]
        egf_negitve = -np.diff(cf_negitve)/self.delta
        egf_data = np.hstack((egf_negitve[::-1],[0],egf_positive))

        ### Bandpass the seismogram and write into sac format ###
        st.data = egf_data
        file_filted = st.filter('bandpass',freqmin=1/float(T_band[1]),freqmax=1/float(T_band[0]),corners=4, zerophase=True)
        t = np.linspace(st.stats.sac.b,-st.stats.sac.b-st.stats.sac.delta,st.stats.sac.npts)

        #normalized
        #self.data = st.data/max(abs(st.data))
        self.data = st.data
        self.Fs = 1/st.stats.delta
        return self.data

    def fft(self,color,imageflag,label):
        if label == 'Orig':
            Data = self.data
        elif label == 'Posi':
            Data = self.DataP
        elif label == 'Nega':
            Data = self.DataN
        elif label == 'Posi_win':
            Data = self.winp
        elif label == 'Nega_win':
            Data = self.winn
        n = len(Data)
        k = np.arange(n)
        T = n / self.Fs
        frq = k / T
        frq1 = frq[range(int(n / 2))]
        YY = np.fft.fft(Data)
        Y = np.fft.fft(Data) / n
        Y1 = Y[range(int(n / 2))]
        Y1[0] = Y1[1]
        
        if imageflag:
            plt.plot(frq1,abs(Y1),color)
            plt.xlim(0,1)
            #plt.ylim(0,0.04)
        return frq1,abs(Y1)

def get_fft_std(dayrange,flag):
    fft = []
    Daynum = 0
    for iday in range(dayrange[0],dayrange[1]):
        file = glob.glob('*day'+str(iday)+'*SAC')
        Daynum += 1
        print(file,Daynum)
        sac = sacs()
        sac.readsac(file[0])
        sac.sides()
        [frq,Y] = sac.fft('k',False,flag)
        if Daynum == 1:
            fft = Y
        else:
            fft = np.vstack((fft,Y))
    std = get_std(fft)
    std_smooth = smooth(std, 5)
    return std_smooth,frq,fft


if __name__ == "__main__":
    #Days()
    Single()
