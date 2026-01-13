from matplotlib import rcParams,use
# required for working on server
use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import sys

data=np.loadtxt(sys.argv[1],ndmin=2)

#if data.ndim<2:
#  print "too few wavelengths for a plot"
#  raise SystemExit()

# prepare plot
fig=plt.figure(figsize=(14,8))
rcParams.update({'font.size': 12})

ax1=plt.subplot(311)
ax2=plt.subplot(312,sharex=ax1)
ax3=plt.subplot(313,sharex=ax1)

ax1.grid()
ax2.grid()
ax3.grid()
ax1.set_xlim([data[0,0]*0.99,data[-1,0]*1.01])
ax3.set_xlabel('wavelength [$\mu$m]')
ax1.set_ylabel ('extinction coefficent [$m^{-1}$]')
ax2.set_ylabel ('single scattering albedo')
ax3.set_ylabel ('asymmetry parameter')

ax1.plot(data[:,0],data[:,1],'r',marker='o')
ax2.plot(data[:,0],data[:,2],'r',marker='o')
ax3.plot(data[:,0],data[:,3],'r',marker='o')

fig.tight_layout()
plt.savefig(sys.argv[1]+'.png',dpi=100)
