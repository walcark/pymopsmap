import numpy as np
import sys
import random

types=("basalt","basalt_andesite","andesite","dacite","rhyolite")
mr_data=np.loadtxt('jgrd54020-sup-0015-2016JD026328-ds15.txt',skiprows=27)
mi_data=np.loadtxt('jgrd54020-sup-0016-2016JD026328-ds16.txt',skiprows=27)

for i_type in range(len(types)):
  f=open('refrac_%s'%types[i_type],'w')
  for i_wvl in range(mr_data.shape[0]):
    f.write('%f %f %f\n'%(mr_data[i_wvl,0]*0.001,mr_data[i_wvl,1+i_type],mi_data[i_wvl,1+i_type]))
  f.close()
