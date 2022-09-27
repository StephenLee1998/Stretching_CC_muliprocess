# This Code ultilize the Stretching Method to measure the
# velocity perturbation from moving-window stack of ncfs in multiprocess
# Sep 27 , 2022   By li_chao , NJU

import obspy
import os
import glob
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from multiprocessing import Pool
from matplotlib import pyplot as plt

def cal_correlation(data1,data2):
        D1 = pd.Series(data1)
        D2 = pd.Series(data2)
        cc = D1.corr(D2)
        return(cc)

class Stretch_CC:
    # Class for bandpass , window and stretch
    def strech(self,data,data_standard,flag):
        # Define the stretching percent range and step
        ratio = np.arange(-0.1,0.1,0.0005)
        # Select the side
        if flag == 'pos':
            T=self.T1
        elif flag == 'neg':
            T=self.T2
        t=self.t
        cc = []
        # Stretch
        for R in ratio:
            Ts = t*(1+R)
            f = interp1d(Ts, data,kind='linear')
            data_strech = f(T)
            CC = cal_correlation(data_strech,data_standard)
            cc.append(CC)
        cc_max = max(cc)
        idcc = np.where(cc == cc_max)[0][0]
        return ratio[idcc],cc_max

    def window(self,data_standard):
        n_zero  = self.half_len
        dt = self.delta
        Trange = self.Trange

        i1 = n_zero + int(self.Trange[1]/dt)
        i2 = n_zero + int(self.Trange[0]/dt)
        i3 = n_zero - int(self.Trange[1]/dt)
        i4 = n_zero - int(self.Trange[0]/dt)
    
        self.T1 = np.arange((i1-self.half_len)*self.delta, (i2-self.half_len)*self.delta ,self.delta)
        self.T2 = np.arange((i4-self.half_len)*self.delta, (i3-self.half_len)*self.delta ,self.delta)
    
        win_positive = data_standard[i1:i2]
        win_negative = data_standard[i4:i3]
    
        #strech(data,standard_positive,T,t)
        return win_positive,win_negative
    
    def bpfilter(self,data,day,Trange,out_path,v_win):
        ### Read the data needed
        tr = obspy.read(data)
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
        file_filted = st.filter('bandpass',freqmin=float(1/float(Trange[1])),freqmax=float(1/float(Trange[0])),corners=4, zerophase=True)
        out = out_path + '/BP_' + str(day) + '.sac'
        file_filted.write(out)
        
        self.t = np.linspace(st.stats.sac.b,-st.stats.sac.b-st.stats.sac.delta,st.stats.sac.npts)
        self.Trange = v_win

        return st.data,self.t


def _streching(daily,station_pair):
    # get the date information
    day = daily.split('_')[2].split('d.')[0].split('day')[1]
    year = daily.split('.dat')[0].split('d_')[1]

    # define the window and period 
    t1 = 35
    t2 = 15
    #Periods : support of multi period bands
    Periods = [[2,8]]

    # define the path and Standard file 
    out_path = '/work/li_chao/work/Axial_Seamount/Out/BP_'+ station_pair
    out_png = '/work/li_chao/work/Axial_Seamount/process_code/Distant_Station/DirectWave/PB_2_8/'+ station_pair+'/'
    
    standard_file = '/work/li_chao/work/Axial_Seamount/NCFs/Standard/DS/Z-Z/ascii/'+station_pair.split('_')[0]+'-'+station_pair.split('_')[1]+'.SAC'
    
    # Stretching in Periods' loop
    for Trange in Periods:
        T_label = str(Trange[0]) + '-' + str(Trange[1])
        out_dir = out_path+'/BP_'+T_label+'s'
        if os.path.isdir(out_dir)== False:
            os.mkdir(out_dir)
    
        filename = T_label+'_ccmax_'
    
        sacs = Strech_CC()
        # Bandpass the moving-window ncf and Standard ncf
        [data,T] = sacs.bpfilter(daily,day,Trange,out_dir,[t1,t2])
        [data_standard,T] = sacs.bpfilter(standard_file,day,Trange,out_dir,[t1,t2])
        # get the window for measurment 
        [d_positive,d_negative] = sacs.window(data_standard)
        # Stretch the positive and negative sides of ncf
        [ratio_positive,cc_p] = sacs.strech(data,d_positive,'pos')
        [ratio_negative,cc_n] = sacs.strech(data,d_negative,'neg')
        print(daily)
    
        # wirte the streching ratio 
        with open(out_png+'/'+filename+'positive.txt','a') as f1:
            f1.write(str(day)+' '+str(100*ratio_positive)+' '+str(cc_p)+' '+str(year)+'\n')
            
        with open(out_png+'/'+filename+'negative.txt','a') as f2:
            f2.write(str(day)+' '+str(100*ratio_negative)+' '+str(cc_n)+' '+str(year)+'\n')


def Thread():
    # Define the station pairs to be calculated 
    station_pair = ['AXCC1_AXBA1','AXEC2_AXBA1','AXEC1_AXBA1','AXAS1_AXBA1','AXID1_AXBA1','AXAS2_AXBA1','AXEC3_AXBA1']

    for pair in station_pair:
        # The input path
        orig_path = '/work/li_chao/work/Axial_Seamount/NCFs/Distant_Station/'+pair+'/Z-Z/ascii'
        # The output path 
        out_files = '/work/li_chao/work/Axial_Seamount/process_code/Distant_Station/DirectWave/PB_2_8/'+ pair.split('_')[0]+'_'+pair.split('_')[1]+'/'
        os.system('mkdir %s' % (out_files))

        os.chdir(orig_path)
        os.system('rm *txt')
        daily_data = glob.glob("*.dat.SAC")
        
    
        p = Pool(14)   # Create process pool
        
        for daily in daily_data:
            p.apply_async(_streching , args = (daily,pair)) # Stretching in multiprocess
        
        p.close()
        p.join()
            
            

if __name__ == "__main__":
    Thread()
    
