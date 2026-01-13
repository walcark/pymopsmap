import numpy as np
import subprocess
import collections

# needs to be adjusted to your installation
path_optical_dataset='../optical_dataset/'
path_mopsmap_executable='../mopsmap'

# Interface for (multi-modal) log-normal distributions

# Input parameters:
# wvl: wavelength (a single number or a list of numbers)
# size_equ: size equivalence
# n,r_mod,sigma: parameters of log-normal size modes (single parameter each or lists with same lengths describing the parameters of each mode)
# r_min,r_max: minimum and maximum radii (if single parameter it is applied to all modes, otherwise if list it is applied to each mode)
# m: refractive index, given as a string which is added in input file after 'refrac' (if single string it is applied to all modes, otherwise if list it is applied to each mode)
# nonabs_fraction: ratio of non-absorbing particles (if single parameter it is applied to all modes, otherwise if list it is applied to each mode)
# shape: particle shape, given as a string which is added in input file after 'shape' (if single string it is applied to all modes, otherwise if list it is applied to each mode)
# num_theta: number of scattering angles in output

# The results are returned as a dictionary with the keywords as listed near the end of this file

def call_mopsmap(wvl,size_equ,n,r_mod,sigma,r_min,r_max,m,nonabs_fraction,shape,num_theta):

  # create a input file for the Fortran code and a wavelength file
  mopsmap_input_file = open('tmp_mopsmap.inp', 'w')
  mopsmap_wvl_file = open('tmp_mopsmap.wvl', 'w')

  # write wavelength file
  wvl=np.array(wvl,ndmin=1)
  mopsmap_input_file.write("wavelength file tmp_mopsmap.wvl \n")
  for i_wvl in range(wvl.shape[0]):
    mopsmap_wvl_file.write('%10.8f \n'%wvl[i_wvl])
  mopsmap_wvl_file.close()

  # write size_equ
  mopsmap_input_file.write('size_equ %s\n'%size_equ)

  # write modes
  n=np.array(n,ndmin=1)
  r_mod=np.array(r_mod,ndmin=1)
  sigma=np.array(sigma,ndmin=1)
  r_min=np.array(r_min,ndmin=1)
  r_max=np.array(r_max,ndmin=1)
  if n.shape!=r_mod.shape or n.shape!=sigma.shape:
    print "shapes of n, r_mod, and sigma do not agree"
    raise SystemExit()
  if n.shape[0]>1 and r_min.shape[0]==1:
    r_min=np.resize(r_min,n.shape[0])
    r_min[:]=r_min[0]
  if n.shape[0]>1 and r_max.shape[0]==1:
    r_max=np.resize(r_max,n.shape[0])
    r_max[:]=r_max[0]

  if isinstance(m, basestring):
    m=[m,]
  if n.shape[0]>1 and len(m)==1:
    for i in range(1,n.shape[0]):
      m.append(m[0])

  if isinstance(nonabs_fraction, (int, long, float, complex)):
    nonabs_fraction=[nonabs_fraction,]
  if n.shape[0]>1 and len(nonabs_fraction)==1:
    for i in range(1,n.shape[0]):
      nonabs_fraction.append(nonabs_fraction[0])

  if isinstance(shape, basestring):
    shape=[shape,]
  if n.shape[0]>1 and len(shape)==1:
    for i in range(1,n.shape[0]):
      shape.append(shape[0])

  for i_mode in range(n.shape[0]):
    mopsmap_input_file.write('mode %d size log_normal %f %f %f %f %f\n'%(i_mode+1,r_mod[i_mode],sigma[i_mode],n[i_mode],r_min[i_mode],r_max[i_mode]))
    mopsmap_input_file.write('mode %d refrac %s\n'%(i_mode+1,m[i_mode]))
    mopsmap_input_file.write('mode %d refrac nonabs_fraction %f\n'%(i_mode+1,nonabs_fraction[i_mode]))
    mopsmap_input_file.write('mode %d shape %s\n'%(i_mode+1,shape[i_mode]))

  mopsmap_input_file.write('scatlib \'%s\'\n'%path_optical_dataset)
  mopsmap_input_file.write('output integrated\n')
  mopsmap_input_file.write('output scattering_matrix\n')
  mopsmap_input_file.write('output volume_scattering_function\n')
  mopsmap_input_file.write('output num_theta %i\n'%num_theta)
  mopsmap_input_file.write('output lidar\n')
  mopsmap_input_file.write('output digits 15\n')
  mopsmap_input_file.write('output ascii_file tmp_mopsmap\n')

  mopsmap_input_file.close()

  # after writing the input file now start mopsmap
  p = subprocess.Popen([path_mopsmap_executable, 'tmp_mopsmap.inp'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
  stdout1,stderr1 = p.communicate()

  if stdout1 or stderr1:
    if stdout1:
      print stdout1
    if stderr1:
      print stderr1
    raise SystemExit()

  # read the mopsmap output files into numpy arrays
  output_integrated=np.loadtxt('tmp_mopsmap.integrated',ndmin=1,dtype=[('wvl', 'f8'),('ext_coeff', 'f8'), ('ssa','f8'),('g','f8'),('r_eff','f8'),('n','f8'),('a','f8'),('v','f8'),('m','f8'),('ext_angstrom','f8'),('sca_angstrom','f8'),('abs_angstrom','f8')])
  output_matrix=np.loadtxt('tmp_mopsmap.scattering_matrix',ndmin=1,dtype=[('wvl', 'f8'),('angle', 'f8'), ('a1','f8'),('a2','f8'),('a3','f8'),('a4','f8'),('b1','f8'),('b2','f8')])
  output_vol_scat=np.loadtxt('tmp_mopsmap.volume_scattering_function',ndmin=1,dtype=[('wvl', 'f8'),('angle', 'f8'), ('a1_vol','f8')])
  output_lidar=np.loadtxt('tmp_mopsmap.lidar',ndmin=1,dtype=[('wvl', 'f8'),('ext_coeff', 'f8'), ('back_coeff','f8'), ('S','f8'), ('delta_l','f8'),('ext_angstrom','f8'),('back_angstrom','f8')])

  # store the results in an easier-to-use way
  num_wvl=output_integrated['wvl'].shape[0]
  num_angles=output_matrix['angle'].shape[0]/num_wvl

  results={}
  results['wvl']=output_integrated['wvl']
  results['ext_coeff']=output_integrated['ext_coeff']
  results['ssa']=output_integrated['ssa']
  results['g']=output_integrated['g']
  results['r_eff']=output_integrated['r_eff']
  results['n']=output_integrated['n']
  results['a']=output_integrated['a']
  results['v']=output_integrated['v']
  results['m']=output_integrated['m']
  results['ext_angstrom']=output_integrated['ext_angstrom']
  results['sca_angstrom']=output_integrated['sca_angstrom']
  results['abs_angstrom']=output_integrated['abs_angstrom']
  results['angle']=output_matrix['angle'][0:num_angles]
  results['a1']=output_matrix['a1'].reshape((num_wvl,num_angles))
  results['a2']=output_matrix['a2'].reshape((num_wvl,num_angles))
  results['a3']=output_matrix['a3'].reshape((num_wvl,num_angles))
  results['a4']=output_matrix['a4'].reshape((num_wvl,num_angles))
  results['b1']=output_matrix['b1'].reshape((num_wvl,num_angles))
  results['b2']=output_matrix['b2'].reshape((num_wvl,num_angles))
  results['a1_vol']=output_vol_scat['a1_vol'].reshape((num_wvl,num_angles))
  results['back_coeff']=output_lidar['back_coeff']
  results['S']=output_lidar['S']
  results['delta_l']=output_lidar['delta_l']
  results['back_angstrom']=output_lidar['back_angstrom']

  return results
