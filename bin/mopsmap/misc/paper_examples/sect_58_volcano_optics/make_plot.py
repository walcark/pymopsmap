# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

volcanoes=('GRI','KEL','SAK','EYJ','SPU','RED','SOU','MSH','CHA')
volcano_names=(u'Grimsvötn','Mount Kelud','Mount Sakuraj.',u'Eyjafjallajökull','Mount Spurr','Mount Redoubt','Soufriere Hills','Mount St. Helens',u'Chaitén')
colors=('#002fd5','#00a7ff','#00ce3d','#00ce3d','#00ce3d','#ffad00','#ffad00','#ffad00','#ff0000')
lt=('-','-','-','--','-.','-','--','-.','-')

fig=plt.figure(figsize=(8,8))
rcParams.update({'font.size': 12})

ax1=plt.subplot(211)
ax2=plt.subplot(212,sharex=ax1)

ticks = np.arange(0.72, 1.00, 0.02)
ax1.set_yticks(ticks)
ax1.set_ylabel('single scattering albedo $\omega_0$')
ax1.set_ylim((0.8,1.027))

ticks = np.arange(0.6, 1.00, 0.01)
ax2.set_yticks(ticks)
ax2.set_ylabel('asymmetry parameter $g$')

ticks = np.arange(300, 1600, 100)
ax2.set_xticks(ticks)
ax2.set_xlim((299,1501))
ax2.set_xlabel('wavelength $\lambda$ [nm]')

ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
ax2.grid(which='minor', alpha=0.2)
ax2.grid(which='major', alpha=0.5)

for i_vol in range(len(volcanoes)):
  data=np.loadtxt("results/"+volcanoes[i_vol]+".integrated")
  reff=data[0,4]
  ax1.plot(data[:,0]*1000.,data[:,2],color=colors[i_vol],linestyle=lt[i_vol],label="%s"%(volcano_names[i_vol]))
  ax2.plot(data[:,0]*1000.,data[:,3],color=colors[i_vol],linestyle=lt[i_vol])

ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0215), ncol=3 )
plt.savefig('plot.pdf',dpi=300)
