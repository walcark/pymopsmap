#!/usr/bin/python
# -*- coding: utf8 -*-

import mopsmap_python_interface
import numpy as np

# different types of wavelength input
wvl=(0.355,0.532,1.064)
wvl=0.344
wvl=np.linspace(0.4,1.0,num=7)

# size equivalence
size_equ='vol'
size_equ='cs'

# two log-normal modes
n=(1.*10**8,0.5*10**8) # concentration in m-3
r_mod=(0.1,0.4)
sigma=(2.6,2.2)
r_min=0.01
r_max=10

# refractive index
m="1.54 0.002"
m=("1.54 0.000","1.56 0.01")

nonabs_fraction=0.5
nonabs_fraction=0.

shape="spheroid oblate 1.7"
shape="spheroid distr_file '../data/ar_kandler'"

num_theta=2

results=mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod,sigma,r_min,r_max,m,nonabs_fraction,shape,num_theta)

print '*** wavelength-independent parameters ***'
print '  effective radius:         %8.4f µm'%(np.mean(results['r_eff']))
print '  number concentration:       %8.4e m^-3'%(np.mean(results['n']))
print '  cross section density:      %8.4e m^-1'%(np.mean(results['a']))
print '  volume density:             %8.4e '%(np.mean(results['v']))
print '  mass concentration:         %8.4e g m^-3'%(np.mean(results['m']))

for i_wvl in range(len(results['wvl'])):
  if i_wvl>0:
    print
    print '*** wavelength range from %8.4f µm to %8.4f µm ***'%(results['wvl'][i_wvl-1],results['wvl'][i_wvl])
    print '  extinction angstrom:      %8.4f'%results['ext_angstrom'][i_wvl]
    print '  scattering angstrom:      %8.4f'%results['sca_angstrom'][i_wvl]
    print '  absorption angstrom:      %8.4f'%results['abs_angstrom'][i_wvl]
    print '  backscatter angstrom:     %8.4f'%results['back_angstrom'][i_wvl]
  print
  print '*** wavelength %8.4f µm ***'%results['wvl'][i_wvl]
  print '  extinction coeff.:          %8.4e m^-1'%(results['ext_coeff'][i_wvl])
  print '  single scattering albedo: %8.4f'%(results['ssa'][i_wvl])
  print '  asymmetry parameter:      %8.4f'%(results['g'][i_wvl])
  print '  phasefunc. at 0°:       %10.4f'%(results['a1'][i_wvl,0])
  print '  phasefunc. at 180°:       %8.4f'%(results['a1'][i_wvl,-1])
  print '  2nd matrix elem. at 180°: %8.4f'%(results['a2'][i_wvl,-1])
  print '  backscatter coef.:          %8.4e m^-1 sr^-1'%results['back_coeff'][i_wvl]
  print '  lidar ratio:              %8.4f sr'%results['S'][i_wvl]
  print '  lin. depol. ratio         %8.4f'%results['delta_l'][i_wvl]
