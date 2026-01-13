# This python script calls MOPSMAP with random microphysical parameters and determines the mean and the range of the modeled parameters
# The microphysical parameters are choosen with constant probability within the given range, e.g. input_ref['sigma']=2.6 and input_dev['sigma']=0.1 means that sigma is a random number between 2.5 and 2.7
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import random

import mopsmap_python_interface

# number of sampled ensembles (high n_sample reduces the statistical uncertainty of the result)
n_sample=1000

wvl=0.532
num_theta=2
size_equ='cs'
n=1

input_ref={}
input_dev={}

input_ref['sigma']=2.6
input_dev['sigma']=0.1

input_ref['r_mod']=0.1
input_dev['r_mod']=0.01

input_ref['r_min']=0.001
input_dev['r_min']=0.0

input_ref['r_max']=20
input_dev['r_max']=0.0

input_ref['mr']=1.53
input_dev['mr']=0.03

input_ref['mi']=0.0063
input_dev['mi']=0.002

input_ref['nonabs_fraction']=0
input_dev['nonabs_fraction']=0

input_ref['ar']=2.0
input_dev['ar']=0.5

# do reference calculation
input=input_ref.copy()
results_ref = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,input['r_mod'],input['sigma'],input['r_min'],input['r_max'],'%f %f'%(input['mr'],input['mi']),input['nonabs_fraction'],'spheroid prolate %f'%input['ar'],num_theta)

print '    ext     |    SSA     |    g       | back_coef  |    S       |   delta_l  |    n       |     a      |     v       | '
print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e | ref_value '%(results_ref['ext_coeff'][0],results_ref['ssa'][0],results_ref['g'][0],results_ref['back_coeff'][0],results_ref['S'][0],results_ref['delta_l'][0],results_ref['n'][0],results_ref['a'][0],results_ref['v'][0])

# do monte carlo sampling
ext=np.zeros(n_sample)
ssa=np.zeros(n_sample)
g=np.zeros(n_sample)
back_coeff=np.zeros(n_sample)
S=np.zeros(n_sample)
delta_l=np.zeros(n_sample)
n1=np.zeros(n_sample)
a=np.zeros(n_sample)
v=np.zeros(n_sample)

for i_sample in range(n_sample):

  input=input_ref.copy()
  for key in sorted(input_ref):
    random_number=2.0*(random.random()-0.5) # range from -1 to +1
    input[key]=input[key]+random_number*input_dev[key]

  results_tmp = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,input['r_mod'],input['sigma'],input['r_min'],input['r_max'],'%f %f'%(input['mr'],input['mi']),input['nonabs_fraction'],'spheroid prolate %f'%input['ar'],num_theta)

  ext[i_sample]=results_tmp['ext_coeff'][0]
  ssa[i_sample]=results_tmp['ssa'][0]
  g[i_sample]=results_tmp['g'][0]
  back_coeff[i_sample]=results_tmp['back_coeff'][0]
  S[i_sample]=results_tmp['S'][0]
  delta_l[i_sample]=results_tmp['delta_l'][0]
  n1[i_sample]=results_tmp['n'][0]
  a[i_sample]=results_tmp['a'][0]
  v[i_sample]=results_tmp['v'][0]

print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e | mean '%(np.mean(ext),np.mean(ssa),np.mean(g),np.mean(back_coeff),np.mean(S),np.mean(delta_l),np.mean(n1),np.mean(a),np.mean(v))
print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e | min '%(np.min(ext),np.min(ssa),np.min(g),np.min(back_coeff),np.min(S),np.min(delta_l),np.min(n1),np.min(a),np.min(v))
print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e | max '%(np.max(ext),np.max(ssa),np.max(g),np.max(back_coeff),np.max(S),np.max(delta_l),np.max(n1),np.max(a),np.max(v))
print '%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e | std '%(np.std(ext),np.std(ssa),np.std(g),np.std(back_coeff),np.std(S),np.std(delta_l),np.std(n1),np.std(a),np.std(v))
