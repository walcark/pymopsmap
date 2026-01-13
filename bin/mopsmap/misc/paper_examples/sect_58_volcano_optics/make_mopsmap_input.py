import numpy as np
import sys

ar_grid=(1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0)
ar_dist=np.zeros(len(ar_grid))

print "output integrated"
print "output scattering_matrix"
print "output lidar"
print "output ascii_file 'results/%s'"%sys.argv[2]
print "wavelength range 0.3 1.5 0.025"
print "refrac file %s"%sys.argv[3]
print "refrac nonabs_fraction 0.5"
print "scatlib '/home/josef/uni/optical_dataset'"

# read the size-shape data
data=np.loadtxt(sys.argv[1],skiprows=27)

i_mode=1
for i in range(data.shape[0]):

  r=data[i,0]*0.5
  if r>47.5:       # larger size are not covered by the data set
    r=47.5

  ar=1./data[i,14]
  if ar>5:         # larger aspect ratios are not covered by the data set
    ar=5

  if i==0:         # at the beginning
    r_last=r

  if r==r_last:    # if same radius just add particle
    for i_ar in range(1,len(ar_grid)):
      if ar<=ar_grid[i_ar]:
        weight=(ar_grid[i_ar]-ar)/(ar_grid[i_ar]-ar_grid[i_ar-1])
        ar_dist[i_ar]+=1.0-weight
        ar_dist[i_ar-1]+=weight
        break

  if r!=r_last or i==(data.shape[0]-1):
    for i_ar in range(len(ar_grid)):    # write the radius with the ar-distr
      if ar_dist[i_ar]!=0:
        print "mode %d size %f %f"%(i_mode,r_last,ar_dist[i_ar])
        print "mode %d shape spheroid prolate %f"%(i_mode,ar_grid[i_ar])
        i_mode+=1

    ar_dist=np.zeros(len(ar_grid))      # set ar_dist to zero and add new particle
    for i_ar in range(1,len(ar_grid)):
      if ar<=ar_grid[i_ar]:
        weight=(ar_grid[i_ar]-ar)/(ar_grid[i_ar]-ar_grid[i_ar-1])
        ar_dist[i_ar]+=1.0-weight
        ar_dist[i_ar-1]+=weight
        break

  r_last=r
