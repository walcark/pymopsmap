import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

import mopsmap_python_interface

fig=plt.figure(figsize=(8,6))
rcParams.update({'font.size': 13})

ax1=plt.subplot(211)
ax2=plt.subplot(212,sharex=ax1)

# parameters for the calculations
wvl=(0.45,0.525,0.635)

size_equ='cs'

r_mod=np.logspace(np.log10(0.001), np.log10(2), num=150)
n=1
sigma=1.6
r_min=0.001
r_max=5.0

m='file refr_mineral'
nonabs_fraction=0.5

num_theta=1801

# do the mopsmap calculations
results_sphere=np.zeros((len(wvl),len(r_mod),num_theta))
results_spheroid=np.zeros((len(wvl),len(r_mod),num_theta))

for i_rmod in range(len(r_mod)):

  shape='sphere'
  output = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod[i_rmod],sigma,r_min,r_max,m,nonabs_fraction,shape,num_theta)
  results_sphere[:,i_rmod,:]=output['a1'][:,:]

  shape='spheroid distr_file ar_kandler'
  output = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod[i_rmod],sigma,r_min,r_max,m,nonabs_fraction,shape,num_theta)
  results_spheroid[:,i_rmod,:]=output['a1'][:,:]

# do the angular integration
angle_degr=output['angle']
angle_rad=angle_degr*np.pi/180.0

weight_ts_ideal=np.zeros(angle_rad.shape[0])
weight_ts_para=np.zeros(angle_rad.shape[0])
weight_bs_ideal=np.zeros(angle_rad.shape[0])
weight_bs_para=np.zeros(angle_rad.shape[0])

for i in range(len(angle_degr)):
  weight_ts_ideal[i]=np.sin(angle_rad[i])
  if angle_degr[i]>=90:
    weight_bs_ideal[i]=np.sin(angle_rad[i])

  if angle_degr[i]>10 and angle_degr[i]<171:
    weight_ts_para[i]=1.01*np.sin(angle_rad[i])**1.19
    weight_bs_para[i]=max(0,1.01*np.sin(angle_rad[i])**1.19*min(1,(angle_degr[i]-70.25)/39.99))

for i in range(len(angle_degr)):
  print angle_degr[i],weight_ts_ideal[i],weight_ts_para[i],weight_bs_ideal[i],weight_bs_para[i]

results_integrated_sphere=np.zeros((len(wvl),len(r_mod),4))
results_integrated_spheroid=np.zeros((len(wvl),len(r_mod),4))
for i_wvl in range(len(wvl)):
  for i_rmod in range(len(r_mod)):

    results_integrated_sphere[i_wvl,i_rmod,0]=np.trapz(weight_ts_ideal*results_sphere[i_wvl,i_rmod,:],x=angle_rad)
    results_integrated_sphere[i_wvl,i_rmod,1]=np.trapz(weight_ts_para *results_sphere[i_wvl,i_rmod,:],x=angle_rad)
    results_integrated_sphere[i_wvl,i_rmod,2]=np.trapz(weight_bs_ideal*results_sphere[i_wvl,i_rmod,:],x=angle_rad)
    results_integrated_sphere[i_wvl,i_rmod,3]=np.trapz(weight_bs_para *results_sphere[i_wvl,i_rmod,:],x=angle_rad)

    results_integrated_spheroid[i_wvl,i_rmod,0]=np.trapz(weight_ts_ideal*results_spheroid[i_wvl,i_rmod,:],x=angle_rad)
    results_integrated_spheroid[i_wvl,i_rmod,1]=np.trapz(weight_ts_para *results_spheroid[i_wvl,i_rmod,:],x=angle_rad)
    results_integrated_spheroid[i_wvl,i_rmod,2]=np.trapz(weight_bs_ideal*results_spheroid[i_wvl,i_rmod,:],x=angle_rad)
    results_integrated_spheroid[i_wvl,i_rmod,3]=np.trapz(weight_bs_para *results_spheroid[i_wvl,i_rmod,:],x=angle_rad)

rayleigh_ts_factor=np.zeros(len(wvl))
rayleigh_bs_factor=np.zeros(len(wvl))
for i_wvl in range(len(wvl)):
  rayleigh_ts_factor[i_wvl]=results_integrated_sphere[i_wvl,0,0]/results_integrated_sphere[i_wvl,0,1]
  rayleigh_bs_factor[i_wvl]=results_integrated_sphere[i_wvl,0,2]/results_integrated_sphere[i_wvl,0,3]

# make the plots
colors=['b','g','r']
for i_wvl in range(len(wvl)):
  ax1.semilogx(r_mod,results_integrated_sphere[i_wvl,:,0]/results_integrated_sphere[i_wvl,:,1]/rayleigh_ts_factor[i_wvl],label="$\lambda$ = %3.0f nm, spheres"%(wvl[i_wvl]*1000), color=colors[i_wvl], linestyle="-")
  ax1.semilogx(r_mod,results_integrated_spheroid[i_wvl,:,0]/results_integrated_spheroid[i_wvl,:,1]/rayleigh_ts_factor[i_wvl],label="$\lambda$ = %3.0f nm, spheroids"%(wvl[i_wvl]*1000), color=colors[i_wvl], linestyle=":")
  ax2.semilogx(r_mod,results_integrated_sphere[i_wvl,:,2]/results_integrated_sphere[i_wvl,:,3]/rayleigh_bs_factor[i_wvl], color=colors[i_wvl],linestyle="-")
  ax2.semilogx(r_mod,results_integrated_spheroid[i_wvl,:,2]/results_integrated_spheroid[i_wvl,:,3]/rayleigh_bs_factor[i_wvl], color=colors[i_wvl],linestyle=":")
ax1.legend()
ticks = np.arange(0.0, 3.1, 0.2)
ax1.set_yticks(ticks)
ax1.set_ylim((1.0,2.0))
ax2.set_ylim((0.8,1.2))
ax2.set_xlim((0.01,1.0))
ax2.set_xlabel('median diameter [$\mu$m]')
ax2.set_xlabel('mode radius $r_{mod}$ [$\mu$m]')
ax1.set_ylabel('$C_{ts}$')
ax2.set_ylabel('$C_{bs}$')

ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
ax2.grid(which='minor', alpha=0.2)
ax2.grid(which='major', alpha=0.5)

plt.savefig('plot.pdf',dpi=300)
