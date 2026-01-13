import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

import mopsmap_python_interface

fig=plt.figure(figsize=(8,10))
rcParams.update({'font.size': 18})

ax1=plt.subplot(311)
ax2=plt.subplot(312,sharex=ax1)
ax3=plt.subplot(313,sharex=ax1)

wvl=(0.4,0.6,1.0)

size_equ='cs'

n=(6.67,0.898,0.1016,0.000472)
r_mod=(0.0212,0.07,0.39,1.9)
sigma=(2.24,1.95,2.00,2.15)
r_min=0.01
r_max=np.logspace(np.log10(0.5), np.log10(60), num=200)
m=('file refr_water_soluble','file refr_mineral','file refr_mineral','file refr_mineral')
nonabs_fraction=0.0
shape=('spheres','spheroid distr_file ar_kandler','spheroid distr_file ar_kandler','spheroid distr_file ar_kandler')

num_theta=2

# do the calculations
results_ext=np.zeros((len(wvl),len(r_max)))
results_ssa=np.zeros((len(wvl),len(r_max)))
results_g=np.zeros((len(wvl),len(r_max)))
results_a=np.zeros((len(wvl),len(r_max)))
results_v=np.zeros((len(wvl),len(r_max)))

for i_wvl in range(len(wvl)):
  for i_rmax in range(len(r_max)):

    output = mopsmap_python_interface.call_mopsmap(wvl[i_wvl],size_equ,n,r_mod,sigma,r_min,r_max[i_rmax],m,nonabs_fraction,shape,num_theta)

    results_ext[i_wvl,i_rmax]=output['ext_coeff'][0]
    results_ssa[i_wvl,i_rmax]=output['ssa'][0]
    results_g[i_wvl,i_rmax]=output['g'][0]
    results_a[i_wvl,i_rmax]=output['a'][0]
    results_v[i_wvl,i_rmax]=output['v'][0]

# make the plots
ax1.semilogx(r_max,results_ext[0,:]/results_ext[0,-1],label="$\\alpha_{ext}$ at 400nm", color='b')
ax1.semilogx(r_max,results_ext[1,:]/results_ext[1,-1],label="$\\alpha_{ext}$ at 600nm", color='g')
ax1.semilogx(r_max,results_ext[2,:]/results_ext[2,-1],label="$\\alpha_{ext}$ at 1000nm", color='r')
ax1.semilogx(r_max,results_a[0,:]/results_a[0,-1],label="a - 'area'", color='k')
ax1.semilogx(r_max,results_v[0,:]/results_v[0,-1],label="M - 'mass'", color='grey')
ax1.set_ylim((0,1.0))
ax1.axvline(x=1.25,color='y')
ax1.axvline(x=5.0,color='y')
ax1.axvline(x=10.0,color='y')
ax1.legend()
ticks = np.arange(0.0, 1.1, 0.1)
ax1.set_yticks(ticks)
ax1.set_ylabel('normalized value')
print results_a[0,:]/results_a[0,-1], results_v[0,:]/results_v[0,-1]

ax2.semilogx(r_max,results_ssa[0,:],color='b')
ax2.semilogx(r_max,results_ssa[1,:],color='g')
ax2.semilogx(r_max,results_ssa[2,:],color='r')
ticks = np.arange(0.76, 1.04, 0.04)
ax2.set_yticks(ticks)
ax2.set_ylim((0.76,1.0))
ax2.axvline(x=1.25,color='y')
ax2.axvline(x=5.0,color='y')
ax2.axvline(x=10.0,color='y')
ax2.set_ylabel('single scat. alb. $\omega_0$')
print results_ssa[0,:],results_ssa[1,:],results_ssa[2,:]

ax3.semilogx(r_max,results_g[0,:],color='b')
ax3.semilogx(r_max,results_g[1,:],color='g')
ax3.semilogx(r_max,results_g[2,:],color='r')
ticks = np.arange(0.6, 0.84, 0.04)
ax3.set_yticks(ticks)
ax3.set_ylim((0.6,0.8))
ax3.axvline(x=1.25,color='y')
ax3.axvline(x=5.0,color='y')
ax3.text(1.25*1.005,0.602,'PM2.5',color='y',fontsize=13)
ax3.text(5.0*1.005,0.602,'PM10',color='y',fontsize=13)
ax3.axvline(x=10.0,color='y')
ax3.set_ylabel('asymmetry para. $g$')
print results_g[0,:],results_g[1,:],results_g[2,:]

ax3.set_xlabel('cut-off radius $r_{max}$ [$\mu$m]')
ax3.set_xlim((r_max[0],r_max[-1]))

ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
ax2.grid(which='minor', alpha=0.2)
ax2.grid(which='major', alpha=0.5)
ax3.grid(which='minor', alpha=0.2)
ax3.grid(which='major', alpha=0.5)

fig.tight_layout()
plt.savefig('plot.pdf',dpi=300)
