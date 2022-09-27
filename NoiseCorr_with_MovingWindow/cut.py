#! /usr/env/python

import glob
import os
import numpy as np
import subprocess
day = np.arange(331,366)
path1 = '/project/li_chao/software/Yao/NoiseCorr-2016Jul-v4.2/data_obs/AXEC2/2016'
path2 = '/project/li_chao/software/Yao/NoiseCorr-2016Jul-v4.2/data_obs/AXCC1/2016'

for i in range(len(day)):
    dayid = str(day[i])
    os.chdir(path1)
    a = glob.glob("*%s*" % (dayid))
    npts1 = os.system("saclst npts f %s" % (a[0]))
    os.chdir(path2)
    b = glob.glob("*%s*" % (dayid))
    npts2 = os.system("saclst npts f %s" % (a[0]))
    if npts1 > npts2 :
        npts1 = npts2
    print(a,b)
    f1 = path1+'/'+a[0]
    f2 = path2+'/'+b[0]
    s = ''
    s += 'cut b n %s' % (npts1)
    s += 'r %s %s\n' %(f1,f2)
    s += 'w over\n'
    s += 'q\n'
    subprocess.Popen(['sac'],stdin=subprocess.PIPE).communicate(s.encode())
     