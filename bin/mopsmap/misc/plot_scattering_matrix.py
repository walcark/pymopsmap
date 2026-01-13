from matplotlib import rcParams,use
# required for working on server
use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import sys

data=np.loadtxt(sys.argv[1])

if data.shape[1]==3:
  only_phase_function=True
else:
  only_phase_function=False

# prepare plot
fig=plt.figure(figsize=(14,8))
rcParams.update({'font.size': 12})

if only_phase_function==False:
  ax1=plt.subplot(321)
  ax2=plt.subplot(322)
  ax3=plt.subplot(323,sharex=ax1)
  ax4=plt.subplot(324,sharex=ax2)
  ax5=plt.subplot(325,sharex=ax1)
  ax6=plt.subplot(326,sharex=ax2)

  ax1.grid()
  ax2.grid()
  ax3.grid()
  ax4.grid()
  ax5.grid()
  ax6.grid()
  ax1.set_xlim([0,180])
  ax2.set_xlim([0,180])
  ax2.set_ylim([0,1])
  ax3.set_ylim([-1,1])
  ax4.set_ylim([-1,1])
  ax5.set_ylim([-1,1])
  ax6.set_ylim([-1,1])
  ax5.set_xlabel('scattering angle')
  ax6.set_xlabel('scattering angle')
  ax1.set_title ('a1')
  ax1.set_ylabel('a1')
  ax2.set_title ('a2/a1')
  ax2.set_ylabel('a2/a1')
  ax3.set_title ('a3/a1')
  ax3.set_ylabel('a3/a1')
  ax4.set_title ('a4/a1')
  ax4.set_ylabel('a4/a1')
  ax5.set_title ('b1/a1')
  ax5.set_ylabel('b1/a1')
  ax6.set_title ('b2/a1')
  ax6.set_ylabel('b2/a1')
else:
  ax1=plt.subplot(111)
  ax1.grid()
  ax1.set_xlim([0,180])
  ax1.set_title ('a1')
  ax1.set_ylabel('a1')
  ax1.set_xlabel('scattering angle')

wavelengths=np.unique(data[:,0])
n_wvl=wavelengths.shape[0]
n_angles=data.shape[0]/n_wvl
if only_phase_function==False:
  data=data.reshape((n_wvl,n_angles,8))
else:
  data=data.reshape((n_wvl,n_angles,3))

# plot only four wavelengths if there are many wavelength
if data.shape[0]>6:
  i_wvl_max=data.shape[0]-1
  i_wvl_list=(0,i_wvl_max/3,i_wvl_max*2/3,i_wvl_max)
else:
  i_wvl_list=range(data.shape[0])

for i_wvl in i_wvl_list:

  ax1.semilogy(data[i_wvl,:,1],data[i_wvl,:,2],label="%8.3f $\mu$m"%wavelengths[i_wvl])
  if only_phase_function==False:
    ax2.plot(data[i_wvl,:,1],data[i_wvl,:,3]/data[i_wvl,:,2])
    ax3.plot(data[i_wvl,:,1],data[i_wvl,:,4]/data[i_wvl,:,2])
    ax4.plot(data[i_wvl,:,1],data[i_wvl,:,5]/data[i_wvl,:,2])
    ax5.plot(data[i_wvl,:,1],data[i_wvl,:,6]/data[i_wvl,:,2])
    ax6.plot(data[i_wvl,:,1],data[i_wvl,:,7]/data[i_wvl,:,2])

ax1.legend(fontsize=12)
fig.tight_layout()
plt.savefig(sys.argv[1]+'.png',dpi=100)
