import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

import mopsmap_python_interface

delta=0.01  # relative perturbation for each parameter

# parameters of the reference ensemble of prolate spheroids with a given aspect ratio
wvl=0.355
wvl=0.532
num_theta=2
size_equ='cs'
n=1

input_ref={}
input_ref['r_mod']=0.1
input_ref['sigma']=2.6
input_ref['r_min']=0.001
input_ref['r_max']=20
input_ref['mr']=1.53
input_ref['mi']=0.0063
input_ref['nonabs_fraction']=0
input_ref['ar']=2.0

# do reference calculation
input=input_ref.copy()
results_ref = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,input['r_mod'],input['sigma'],input['r_min'],input['r_max'],'%f %f'%(input['mr'],input['mi']),input['nonabs_fraction'],'spheroid prolate %f'%input['ar'],num_theta)

print '   d ext    |   d SSA    |   d g      | d back_coef|   d S      |  d delta_l |   d n      |    d a     |   d v      |   value     | lower value| upper_value|  parameter'
print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e|     ---          ---          ---     | ref_value '%(results_ref['ext_coeff'][0],results_ref['ssa'][0],results_ref['g'][0],results_ref['back_coeff'][0],results_ref['S'][0],results_ref['delta_l'][0],results_ref['n'][0],results_ref['a'][0],results_ref['v'][0])


# perturb each parameter and calculate Jacobians
for key in sorted(input_ref):

  # but do not consider parameter which were set to zero
  if input[key]>0:

    # first slightly decrease the parameter value
    input=input_ref.copy()
    input[key]=input[key]*(1.0-delta)
    results_lower = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,input['r_mod'],input['sigma'],input['r_min'],input['r_max'],'%f %f'%(input['mr'],input['mi']),input['nonabs_fraction'],'spheroid prolate %f'%input['ar'],num_theta)
    value_lower=input[key]

    # then slightly increase the parameter value
    input=input_ref.copy()
    input[key]=input[key]*(1.0+delta)
    results_upper = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,input['r_mod'],input['sigma'],input['r_min'],input['r_max'],'%f %f'%(input['mr'],input['mi']),input['nonabs_fraction'],'spheroid prolate %f'%input['ar'],num_theta)
    value_upper=input[key]


    dext=(results_upper['ext_coeff'][0]-results_lower['ext_coeff'][0])/(value_upper-value_lower)
    dssa=(results_upper['ssa'][0]-results_lower['ssa'][0])/(value_upper-value_lower)
    dg=(results_upper['g'][0]-results_lower['g'][0])/(value_upper-value_lower)
    dback_coeff=(results_upper['back_coeff'][0]-results_lower['back_coeff'][0])/(value_upper-value_lower)
    dS=(results_upper['S'][0]-results_lower['S'][0])/(value_upper-value_lower)
    ddelta_l=(results_upper['delta_l'][0]-results_lower['delta_l'][0])/(value_upper-value_lower)
    dn=(results_upper['n'][0]-results_lower['n'][0])/(value_upper-value_lower)
    da=(results_upper['a'][0]-results_lower['a'][0])/(value_upper-value_lower)
    dv=(results_upper['v'][0]-results_lower['v'][0])/(value_upper-value_lower)

    print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e|%12.5e %12.5e %12.5e | d %s'%(dext,dssa,dg,dback_coeff,dS,ddelta_l,dn,da,dv,input_ref[key],value_lower,value_upper,key)
