! Copyright (C) 1998-2018 Michael Mishchenko (NASA), Josef Gasteiger (University Vienna, LMU Munich)
!
! This file is part of MOPSMAP, see <https://mopsmap.net>.
!
! MOPSMAP is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! MOPSMAP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with MOPSMAP. If not, see <http://www.gnu.org/licenses/>.

! *** This subroutine MATR was originally programmed by Michael Mishchenko (NASA).
! *** It is part of Mishchenko's T-matrix program, which is available (last checked: June 2018) from
! ***     https://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html
! *** It was converted from Fortran77 to a newer Fortran version by Josef Gasteiger.

subroutine matr(coeff,lmax,n_theta,angle,mat)

use module_input, only: dp,pi

implicit none

integer :: lmax,n_theta

real(dp) :: coeff(6,0:lmax)
real(dp) :: angle(n_theta)
real(dp) :: mat(6,n_theta)

integer :: i_theta,l
real(dp) :: u,l_dble
real(dp) :: f11,f2,f3,f44,f12,f34
real(dp) :: p1,p2,p3,p4
real(dp) :: pp1,pp2,pp3,pp4
real(dp) :: pl1,pl2,pl3,pl4,p


do i_theta=1,n_theta

  u=cos(angle(i_theta)/180d0*pi)

  f11=0d0; f2=0d0; f3=0d0; f44=0d0; f12=0d0; f34=0d0
  p1=0d0; p2=0d0; p3=0d0; p4=0d0
  pp1=1d0
  pp2=0.25d0*(1d0+u)*(1d0+u)
  pp3=0.25d0*(1d0-u)*(1d0-u)
  pp4=dsqrt(6d0)*0.25d0*(u*u-1d0)

  do l=0,lmax

    l_dble=l

    f11=f11+coeff(1,l)*pp1
    f44=f44+coeff(4,l)*pp1

    pl1=2d0*l_dble+1d0
    p=(pl1*u*pp1-l_dble*p1)/(l_dble+1d0)
    p1=pp1
    pp1=p

    if(l.lt.2) cycle

    f2=f2+(coeff(2,l)+coeff(3,l))*pp2
    f3=f3+(coeff(2,l)-coeff(3,l))*pp3
    f12=f12+coeff(5,l)*pp4
    f34=f34+coeff(6,l)*pp4

    pl2=l_dble*(l_dble+1d0)*u
    pl3=(l_dble+1d0)*(l_dble*l_dble-4d0)
    pl4=1d0/(l_dble*((l_dble+1d0)*(l_dble+1d0)-4d0))

    p=(pl1*(pl2-4d0)*pp2-pl3*p2)*pl4
    p2=pp2
    pp2=p
    p=(pl1*(pl2+4d0)*pp3-pl3*p3)*pl4
    p3=pp3
    pp3=p
    p=(pl1*u*pp4-dsqrt(l_dble*l_dble-4d0)*p4)/dsqrt((l_dble+1d0)*(l_dble+1d0)-4d0)
    p4=pp4
    pp4=p

  end do

  mat(:,i_theta)=(/f11,(f2+f3)*0.5d0,(f2-f3)*0.5d0,f44,f12,f34/)

end do

end subroutine
