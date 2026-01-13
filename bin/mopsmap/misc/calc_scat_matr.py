# more or less direct translation of matr.f90 (adapted from T-matrix code of Michael Mishchenko) to Python

import numpy as np

def calc_scat_matr(a1,a2,a3,a4,b1,b2,angles):

  matr=np.zeros([6,angles.shape[0]])

  for i_angle in range(angles.shape[0]):

    u=np.cos(angles[i_angle]/180.0*np.pi)

    f11=0.0; f2=0.0; f3=0.0; f44=0.0; f12=0.0; f34=0.0
    p1=0.0; p2=0.0; p3=0.0; p4=0.0
    pp1=1.0
    pp2=0.25*(1.0+u)*(1.0+u)
    pp3=0.25*(1.0-u)*(1.0-u)
    pp4=6.0**0.5*0.25*(u*u-1.0)

    for l in range(a1.shape[0]):

      f11=f11+a1[l]*pp1
      f44=f44+a4[l]*pp1

      pl1=2.0*l+1.0
      p=(pl1*u*pp1-l*p1)/(l+1.0)
      p1=pp1
      pp1=p

      if l<2:
        continue

      f2=f2+(a2[l]+a3[l])*pp2
      f3=f3+(a2[l]-a3[l])*pp3
      f12=f12+b1[l]*pp4
      f34=f34+b2[l]*pp4

      pl2=l*(l+1.0)*u
      pl3=(l+1.0)*(l*l-4.0)
      pl4=1.0/(l*((l+1.0)*(l+1.0)-4.0))

      p=(pl1*(pl2-4.0)*pp2-pl3*p2)*pl4
      p2=pp2
      pp2=p
      p=(pl1*(pl2+4.0)*pp3-pl3*p3)*pl4
      p3=pp3
      pp3=p
      p=(pl1*u*pp4-(l*l-4.0)**0.5*p4)/((l+1.0)*(l+1.0)-4.0)**0.5
      p4=pp4
      pp4=p

    matr[:,i_angle]=np.array([f11,(f2+f3)*0.5,(f2-f3)*0.5,f44,f12,f34])

  return matr
