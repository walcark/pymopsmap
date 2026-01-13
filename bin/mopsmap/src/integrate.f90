! Copyright (C) 2007-2018 Josef Gasteiger (University Vienna, LMU Munich)
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

subroutine integrate(n,x,y,x_min,x_max,results)

! ***************
! *** integration with trapezoidal rule
! *** x_min needs to be between x(1) and x(2)
! *** x_max needs to be between x(n-1) and x(n)
! ***************

use module_input, only: dp

implicit none

integer :: n,i
real(dp) :: x(n),y(n)
real(dp) :: x_min,x_max
real(dp) :: results

! *** check if xrange fits the x-array
if(n.lt.2.or.x_min.lt.x(1).or.x_min.gt.x(2).or.x_max.lt.x(n-1).or.x_max.gt.x(n)) then
  print*,"integration error"
  print*,"x_min=",x_min,"x_max=",x_max
  print*,"n=",n,"x=",x
  stop
end if

if(n.eq.2) then
  results=(x_max-x_min)*((y(2)+y(1))/2d0+(y(2)-y(1))/(x(2)-x(1))*((x_max+x_min)-(x(2)+x(1)))/2d0)
else
  ! *** first interval
  results=(x(2)-x_min)*((y(2)+y(1))/2d0+(y(2)-y(1))/(x(2)-x(1))*((x(2)+x_min)-(x(2)+x(1)))/2d0)
  ! *** last interval
  results=results+(x_max-x(n-1))*((y(n)+y(n-1))/2d0+(y(n)-y(n-1))/(x(n)-x(n-1))*((x_max+x(n-1))-(x(n)+x(n-1)))/2d0)
  ! *** the other intervals
  do i=3,n-1
    results=results+(y(i)+y(i-1))/2d0*(x(i)-x(i-1))
  end do
end if

end subroutine
