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

subroutine interpolate_linear(num_array,array,value,i_upper,i_lower,weight_upper,weight_lower)

  use module_input, only: dp

  implicit none

  ! *** input
  integer :: num_array
  real(dp) :: array(num_array)
  real(dp) :: value

  ! *** output
  integer :: i_upper,i_lower
  real(dp) :: weight_upper,weight_lower

  integer :: i

  if(num_array.lt.1) then

    print*,"error 'num_array<1' in subroutine interpolate"
    stop

  elseif(num_array.eq.1) then

    if(abs(1d0-value/array(1)).gt.1d-6) then
      print*,"error 'value!=array(1)' in subroutine interpolate"
      stop
    end if
    i_upper=1
    i_lower=1
    weight_upper=1d0
    weight_lower=1d0

  else

    ! *** check whether array is monotonically increasing
    do i=2,num_array
      if(array(i).lt.array(i-1)) then
        print*,"error 'array not monotonically increasing' in subroutine interpolate"
        stop
      end if
    end do

    ! find the nearest grid points
    do i=1,num_array
      if(array(i).gt.value) exit
    end do

    if(array(num_array).lt.value) then
      print*,"error 'value>array(num_array)' in subroutine interpolate"
      stop
    end if
    if(array(1).gt.value) then
      print*,"error 'value<array(1)' in subroutine interpolate"
      stop
    end if

    i=min(i,num_array)
    i_upper=i
    i_lower=i-1

    weight_upper=(value-array(i-1))/(array(i)-array(i-1))
    weight_lower=(array(i)-value)/(array(i)-array(i-1))

    if(abs(weight_upper).lt.1d-9) weight_upper=0
    if(abs(weight_lower).lt.1e-9) weight_lower=0

  end if

end subroutine
