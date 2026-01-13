# Plot the optical properties of single particles in MOPSMAP data set as function of size parameter
# The netcdf files need to be given as commandline arguments

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys
import string
import calc_scat_matr

fig=plt.figure(figsize=(8,12))
rcParams.update({'font.size': 10})

ax1=plt.subplot(411)
ax2=plt.subplot(412,sharex=ax1)
ax3=plt.subplot(413,sharex=ax1)
ax4=plt.subplot(414,sharex=ax1)

for i_file in range(1,len(sys.argv)):

  input_file=Dataset(sys.argv[i_file], mode='r')

  sizepara = input_file.variables['sizepara'][:]
  qext = input_file.variables['qext'][:]
  qsca = input_file.variables['qsca'][:]

  lmax = input_file.variables['lmax'][:]
  a1 = input_file.variables['a1'][:]
  a2 = input_file.variables['a2'][:]
  a3 = input_file.variables['a3'][:]
  a4 = input_file.variables['a4'][:]
  b1 = input_file.variables['b1'][:]
  b2 = input_file.variables['b2'][:]

  input_file.close()

  label=sys.argv[i_file].split('/')[-1]
  label=label[0:(string.find(label,".",-5,-1))]
  #label=label[15:20]

  ax1.semilogx(sizepara,qext,label=label)
  ax2.semilogx(sizepara,qsca/qext)

  a1_index=[]
  for i in range(len(sizepara)):
    a1_index.append(1+i+sum(lmax[0:i]))

  ax3.semilogx(sizepara,a1[a1_index]/3.)

  depol=[]
  for i in range(len(sizepara)):
    l_start=i+sum(lmax[0:i])
    l_end=l_start+lmax[i]
    a1_temp=a1[l_start:l_end]
    a2_temp=a2[l_start:l_end]
    a3_temp=a3[l_start:l_end]
    a4_temp=a4[l_start:l_end]
    b1_temp=b1[l_start:l_end]
    b2_temp=b2[l_start:l_end]
    angle=np.array([0,90,180])
    matr=calc_scat_matr.calc_scat_matr(a1_temp,a2_temp,a3_temp,a4_temp,b1_temp,b2_temp,angle)
    depol.append(1.0-matr[1,-1]/matr[0,-1])

  ax4.semilogx(sizepara,depol)

ax1.set_xlim([1,1000]) # size parameter range of the plot
ax4.set_xlabel('size parameter $x_c$')
ax1.legend(fontsize=8)

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

ax1.set_ylim([0,5.0])
ax1.set_ylabel('$q_{ext}$')

ax2.set_ylim([0,1.05])
ax2.set_ylabel('$\omega_0$')

ax3.set_ylim([0,1.001])
ax3.set_ylabel('$g$')

ax4.set_ylim([0,1.001])
ax4.set_ylabel('$\delta_l$')

plt.tight_layout()
plt.savefig('plot.png',dpi=300)
