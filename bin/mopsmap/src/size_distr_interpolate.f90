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

real(dp) function size_distr_interpolate(radius,len,distr,type)

! *******************************************************
! *** distr(:,1) contains radii, distr(:,2) the SD values
! *******************************************************

  use module_input, only: dp,pi
  use module_aux_var

  implicit none

  real(dp), parameter :: lnlog=0.434294481903252d0 ! *** 1/ln(10)

  real(dp) :: radius
  integer :: len
  real(dp) :: distr(len,2)
  character(len=30) :: type

  integer :: i_upper,i_lower
  real(dp) :: weight_upper,weight_lower


  ! *** size distribution is assumed constant outside range (improves integration over range)
  if(radius.le.distr(1,1)) then
    size_distr_interpolate=distr(1,2)
  elseif(radius.ge.distr(len,1)) then
    size_distr_interpolate=distr(len,2)
  else
    call interpolate_linear(len,distr(:,1),radius,i_upper,i_lower,weight_upper,weight_lower)
    size_distr_interpolate=distr(i_upper,2)*weight_upper+distr(i_lower,2)*weight_lower
  end if

  ! *** convert data to dndr
  if(type.eq."dndr") then
    size_distr_interpolate=size_distr_interpolate*1d0
  elseif(type.eq."dndlnr") then
    size_distr_interpolate=size_distr_interpolate/radius
  elseif(type.eq."dndlogr") then
    size_distr_interpolate=size_distr_interpolate/radius*lnlog
  elseif(type.eq."dadr") then
    size_distr_interpolate=size_distr_interpolate/(radius**2*pi)
    size_distr_interpolate=size_distr_interpolate*1d12         ! consider that radius is in micrometer but cross section area in square meters
  elseif(type.eq."dadlnr") then
    size_distr_interpolate=size_distr_interpolate/(radius**3*pi)
    size_distr_interpolate=size_distr_interpolate*1d12
  elseif(type.eq."dadlogr") then
    size_distr_interpolate=size_distr_interpolate/(radius**3*pi)*lnlog
    size_distr_interpolate=size_distr_interpolate*1d12
  elseif(type.eq."dvdr") then
    size_distr_interpolate=size_distr_interpolate/(radius**3*pi*4d0/3d0)
    size_distr_interpolate=size_distr_interpolate*1d18         ! consider that radius is in micrometer but volume in cubic meters
  elseif(type.eq."dvdlnr") then
    size_distr_interpolate=size_distr_interpolate/(radius**4*pi*4d0/3d0)
    size_distr_interpolate=size_distr_interpolate*1d18
  elseif(type.eq."dvdlogr") then
    size_distr_interpolate=size_distr_interpolate/(radius**4*pi*4d0/3d0)*lnlog
    size_distr_interpolate=size_distr_interpolate*1d18
  else
    print*,"Error: Unknown size_distr_table_type"
    stop
  end if

end function
