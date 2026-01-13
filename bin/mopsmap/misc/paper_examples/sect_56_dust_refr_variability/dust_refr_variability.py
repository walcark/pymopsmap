import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

import mopsmap_python_interface

rcParams.update({'font.size': 13})

# general settings
wvl=0.355

size_equ='cs'

n=np.array((6.667,0.898,0.1016,0.000472))
n=n*10**12
r_mod=(0.0212,0.07,0.39,1.9)
sigma=(2.24,1.95,2.00,2.15)
shape=("sphere","spheroid distr_file 'ar_kandler'","spheroid distr_file 'ar_kandler'","spheroid distr_file 'ar_kandler'")

opac_mi=0.0166 # at 355nm

num_theta=1801

# read the mi distribution
fobj = open("kandler_TEB_550_sm_Appendix_S1_dustonly.txt", "r")

fobj.readline() # first line is a comment
mi_bins=[]
header=fobj.readline().split()
for i in range(len(header)):
  if header[i]=="k":
    mi_bins.append((float(header[i-2]),float(header[i+2])))

d_bins=[]
mi_dist=[]
for line in fobj:
  line=line.split()
  if len(line)>0:

    wvl_tmp=float(line[1])/1000.0

    if abs(wvl_tmp-wvl)<0.00001:

      d_min=float(line[2])
      d_max=float(line[6])
      d_bins.append((d_min,d_max))

      mi_dist_tmp=[int(x) for x in line[7:24]]
      mi_dist_tmp=[float(x)/sum(mi_dist_tmp) for x in mi_dist_tmp]
      mi_dist.append(mi_dist_tmp)

fobj.close()

# calculate average mi
mi_avg=0
for i_d in range(len(d_bins)):
  mi_avg1=0
  for i_mi in range(len(mi_dist[i_d])):
    mi_avg1+=mi_dist[i_d][i_mi]*0.5*(mi_bins[i_mi][0]+mi_bins[i_mi][1])
  mi_avg+=mi_avg1
mi_avg=mi_avg/len(d_bins)

# do the accurate calculation
nonabs_fraction=0
n_tmp=[]
r_mod_tmp=[]
sigma_tmp=[]
r_min_tmp=[]
r_max_tmp=[]
m_tmp=[]
shape_tmp=[]

# waso
n_tmp.append(n[0])
r_mod_tmp.append(r_mod[0])
sigma_tmp.append(sigma[0])
r_min_tmp.append(d_bins[0][0]*0.5)
r_max_tmp.append(d_bins[-1][1]*0.5)
m_tmp.append('1.53 0.005')
shape_tmp.append(shape[0])

for i_d in range(len(d_bins)):
  for i_mi in range(len(mi_dist[i_d])):
    for i_mode in range(1,len(n)):
      n_tmp.append(n[i_mode]*mi_dist[i_d][i_mi])
      r_mod_tmp.append(r_mod[i_mode])
      sigma_tmp.append(sigma[i_mode])
      r_min_tmp.append(d_bins[i_d][0]*0.5)
      r_max_tmp.append(d_bins[i_d][1]*0.5)
      m_tmp.append('1.53 %f'%(0.5*(mi_bins[i_mi][0]+mi_bins[i_mi][1])))
      shape_tmp.append(shape[1])

result_accurate = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n_tmp,r_mod_tmp,sigma_tmp,r_min_tmp,r_max_tmp,m_tmp,nonabs_fraction,shape_tmp,num_theta)

# do calculation with average mi
nonabs_fraction=0
m_tmp=[]
m_tmp.append('1.53 0.005')
for i in range(3):
  m_tmp.append('1.53 %f'%mi_avg)
r_min=d_bins[0][0]*0.5
r_max=d_bins[-1][1]*0.5
result_average = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod,sigma,r_min,r_max,m_tmp,nonabs_fraction,shape,num_theta)

# do calculation with parameterization 50% nonabsorbing and 50% absorbing particles
nonabs_fraction=0.5
result_para50 = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod,sigma,r_min,r_max,m_tmp,nonabs_fraction,shape,num_theta)

# do calculation with parameterization 25% nonabsorbing and 75% absorbing particles
nonabs_fraction=0.25
result_para25 = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod,sigma,r_min,r_max,m_tmp,nonabs_fraction,shape,num_theta)

# do calculation with parameterization 75% nonabsorbing and 25% absorbing particles
nonabs_fraction=0.75
result_para75 = mopsmap_python_interface.call_mopsmap(wvl,size_equ,n,r_mod,sigma,r_min,r_max,m_tmp,nonabs_fraction,shape,num_theta)

# print out some results
print "************************************************"
print wvl,mi_avg,opac_mi
print 'accurate ', result_accurate['r_eff'][0], result_accurate['g'][0], result_accurate['ext_coeff'][0], result_accurate['ssa'][0], result_accurate['S'][0], result_accurate['delta_l'][0]
print 'average  ', result_average['r_eff'][0], result_average['g'][0], result_average['ext_coeff'][0], result_average['ssa'][0], result_average['S'][0], result_average['delta_l'][0]
print 'para25   ', result_para25['r_eff'][0], result_para25['g'][0], result_para25['ext_coeff'][0], result_para25['ssa'][0], result_para25['S'][0], result_para25['delta_l'][0]
print 'para50   ', result_para50['r_eff'][0], result_para50['g'][0], result_para50['ext_coeff'][0], result_para50['ssa'][0], result_para50['S'][0], result_para50['delta_l'][0]
print 'para75   ', result_para75['r_eff'][0], result_para75['g'][0], result_para75['ext_coeff'][0], result_para75['ssa'][0], result_para75['S'][0], result_para75['delta_l'][0]

print np.mean(((result_average['a1_vol'][0]-result_accurate['a1_vol'][0])/result_accurate['a1_vol'][0])**2)**0.5
print np.mean(((result_para25['a1_vol'][0]-result_accurate['a1_vol'][0])/result_accurate['a1_vol'][0])**2)**0.5
print np.mean(((result_para50['a1_vol'][0]-result_accurate['a1_vol'][0])/result_accurate['a1_vol'][0])**2)**0.5
print np.mean(((result_para75['a1_vol'][0]-result_accurate['a1_vol'][0])/result_accurate['a1_vol'][0])**2)**0.5

# make the plot
plt.semilogy(result_average['angle'],result_average['a1_vol'][0],label='average $m_i$',color="k")
plt.semilogy(result_para25['angle'],result_para25['a1_vol'][0],color="b",linestyle="--",lw=0.5,label='average $m_i$ with nonabs. fraction $\mathcal{X}$=0.25')
plt.semilogy(result_para50['angle'],result_para50['a1_vol'][0],label='average $m_i$ with nonabs. fraction $\mathcal{X}$=0.50',color="b")
plt.semilogy(result_para75['angle'],result_para75['a1_vol'][0],color="b",linestyle="-",lw=0.5,label='average $m_i$ with nonabs. fraction $\mathcal{X}$=0.75')
plt.semilogy(result_accurate['angle'],result_accurate['a1_vol'][0],label='measured $m_i$ distribution',color="r",alpha=0.8)

handles,labels = plt.gca().get_legend_handles_labels()
handles = [handles[4],handles[0],handles[1],handles[2],handles[3]]
labels = [labels[4],labels[0],labels[1],labels[2],labels[3]]
plt.legend(handles,labels,loc=1)

plt.ylabel('volume scattering function $\widetilde{a}_1$ (arb. scale)')
plt.grid()
plt.xlim((0.,180.))
plt.xlabel('scattering angle $\\theta$')
plt.savefig('plot.pdf',dpi=300)
