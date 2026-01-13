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

subroutine add_contribution(i_wavelength,i_contribution)

use module_dbindex, only: avail_np,avail_eps,avail_mreal,avail_mimag
use module_input, only: dp,pi,mode,size_equ
use module_results
use module_contributions
use netcdf

implicit none

! *** arrays for the integration
real(dp),allocatable :: r(:),n(:),r2_n(:),r3_n(:),csca_n(:),cext_n(:)
real(dp),allocatable :: csca_n_a1(:),csca_n_a2(:),csca_n_a3(:),csca_n_a4(:),csca_n_b1(:),csca_n_b2(:)
integer,allocatable :: l0(:)

! *** data from the netcdf file
real(dp),allocatable :: a1(:),a2(:),a3(:),a4(:),b1(:),b2(:)
real(dp),allocatable :: qext(:),qsca(:),sizepara(:)
integer,allocatable :: lmax(:)
real :: rat

! *** variables for handling of netcdf files
integer :: status,ncid
integer :: sizepara_dimid,index_coeff_dimid
integer :: sizepara_varid
integer :: a1_id,a2_id,a3_id,a4_id,b1_id,b2_id,lmax_id,qext_id
integer :: qsca_id,rat_id
integer :: size_len,index_coeff_len

! *** misc.
integer :: i_size_lower,i_size_upper
integer :: index_coeff_lower,index_coeff_upper
integer :: lmax_max

integer :: i_wavelength,i_contribution
real(dp) :: result_tmp
real(dp) :: log_distr
real(dp) :: mod_gamma
real(dp) :: size_distr_interpolate
integer :: i_size,l

integer :: i_size_min
real(dp) :: rmin_tmp

character (len=1000) :: filename


! *******************************************
! *** open netcdf file required for contribution
! *******************************************
call get_nc_filename(avail_np(contributions(i_contribution)%i_np),avail_eps(contributions(i_contribution)%i_eps),avail_mreal(contributions(i_contribution)%i_mreal),avail_mimag(contributions(i_contribution)%i_mimag),2,filename)
status=nf90_open(trim(filename), 0, ncid)
if(status.ne.nf90_noerr) then
  print *,"Error opening "//trim(filename)//": "//trim(nf90_strerror(status))
  stop
end if

! *************************************************
! *** find out the ids of the variables in the file
! *************************************************
status=nf90_inq_dimid(ncid,'sizepara',sizepara_dimid)
status=nf90_inq_dimid(ncid,'index_coeff',index_coeff_dimid)
status=nf90_inq_varid(ncid,'sizepara',sizepara_varid)
status=nf90_inq_varid(ncid,'a1',a1_id)
status=nf90_inq_varid(ncid,'a2',a2_id)
status=nf90_inq_varid(ncid,'a3',a3_id)
status=nf90_inq_varid(ncid,'a4',a4_id)
status=nf90_inq_varid(ncid,'b1',b1_id)
status=nf90_inq_varid(ncid,'b2',b2_id)
status=nf90_inq_varid(ncid,'lmax',lmax_id)
status=nf90_inq_varid(ncid,'qext',qext_id)
status=nf90_inq_varid(ncid,'qsca',qsca_id)
status=nf90_inq_varid(ncid,'rat',rat_id)

status=nf90_inquire_dimension(ncid,sizepara_dimid,len=size_len)
status=nf90_inquire_dimension(ncid,index_coeff_dimid,len=index_coeff_len)

! ************************************
! *** find out relevant sizeparameters
! ************************************
allocate(sizepara(size_len))
status=nf90_get_var(ncid,sizepara_varid,sizepara)
status=nf90_get_var(ncid,rat_id,rat)

! *** convert sizeparameters to other size definitions if required
if(size_equ.eq.1) sizepara=sizepara*rat
if(size_equ.eq.2) sizepara=sizepara*rat**3

! *** determine which sizeparameters are needed for the integration from rmin to rmax
i_size_lower=1
do i_size=1,size_len
  if(contributions(i_contribution)%rmin.gt.sizepara(i_size)*wavelength_array(i_wavelength)/(2d0*pi)) then
    i_size_lower=i_size
  else
    exit
  end if
end do

i_size_upper=size_len
do i_size=size_len,1,-1
  if(contributions(i_contribution)%rmax.lt.sizepara(i_size)*wavelength_array(i_wavelength)/(2d0*pi)) then
    i_size_upper=i_size
  else
    exit
  end if
end do

i_size_upper=max(i_size_upper,i_size_lower+1)  ! *** at least two size parameters are necessary for integration

if(contributions(i_contribution)%rmax.gt.sizepara(i_size_upper)*wavelength_array(i_wavelength)/(2d0*pi)) then
  print'(A,ES10.3,A,ES10.3,A)',"Error: Maximum specified particle size ",contributions(i_contribution)%rmax," at wavelength ", wavelength_array(i_wavelength)," is not covered by "//trim(filename)
  stop
end if

if(contributions(i_contribution)%rmin.lt.sizepara(i_size_lower)*wavelength_array(i_wavelength)/(2d0*pi)) then
  print'(A,ES10.3,A,ES10.3,A)',"Error: Minimum specified particle size ",contributions(i_contribution)%rmin," at wavelength ", wavelength_array(i_wavelength)," is not covered by "//trim(filename)
  stop
end if

deallocate(sizepara)

! ********************************************
! *** find out required expansion coefficients
!*********************************************
allocate(lmax(size_len))
status=nf90_get_var(ncid,lmax_id,lmax)

index_coeff_lower=i_size_lower+sum(lmax(1:i_size_lower-1))
index_coeff_upper=i_size_upper+sum(lmax(1:i_size_upper))

deallocate(lmax)

! *****************************************************
! *** allocate the variables for the data from the file
! *****************************************************
size_len=i_size_upper-i_size_lower+1
index_coeff_len=index_coeff_upper-index_coeff_lower+1

allocate(lmax(size_len))
allocate(a1(index_coeff_len))
allocate(a2(index_coeff_len))
allocate(a3(index_coeff_len))
allocate(a4(index_coeff_len))
allocate(b1(index_coeff_len))
allocate(b2(index_coeff_len))

allocate(sizepara(size_len))
allocate(qext(size_len))
allocate(qsca(size_len))

! *******************************
! *** read the data from the file
! *******************************
status=nf90_get_var(ncid,lmax_id,lmax,(/i_size_lower/),(/size_len/))
status=nf90_get_var(ncid,a1_id,a1,(/index_coeff_lower/),(/index_coeff_len/))
status=nf90_get_var(ncid,a2_id,a2,(/index_coeff_lower/),(/index_coeff_len/))
status=nf90_get_var(ncid,a3_id,a3,(/index_coeff_lower/),(/index_coeff_len/))
status=nf90_get_var(ncid,a4_id,a4,(/index_coeff_lower/),(/index_coeff_len/))
status=nf90_get_var(ncid,b1_id,b1,(/index_coeff_lower/),(/index_coeff_len/))
status=nf90_get_var(ncid,b2_id,b2,(/index_coeff_lower/),(/index_coeff_len/))

status=nf90_get_var(ncid,sizepara_varid,sizepara,(/i_size_lower/),(/size_len/))
status=nf90_get_var(ncid,qext_id,qext,(/i_size_lower/),(/size_len/))
status=nf90_get_var(ncid,qsca_id,qsca,(/i_size_lower/),(/size_len/))

status=nf90_get_var(ncid,rat_id,rat)

status=nf90_close(ncid)

lmax_max=maxval(lmax)
if(lmax_max.gt.max_num_coeff) then
  print *,"Error: max_num_coeff too small for "//trim(filename)
  stop
end if

if(lmax_max.gt.num_coeff) num_coeff=lmax_max

! ***********************************
! *** allocate arrays for integration
! ***********************************
allocate(n(size_len))
allocate(r(size_len))
allocate(r2_n(size_len))
allocate(r3_n(size_len))
allocate(csca_n(size_len))
allocate(cext_n(size_len))
n=0;r=0;r2_n=0;r3_n=0;csca_n=0;cext_n=0

allocate(l0(size_len))
allocate(csca_n_a1(size_len))
allocate(csca_n_a2(size_len))
allocate(csca_n_a3(size_len))
allocate(csca_n_a4(size_len))
allocate(csca_n_b1(size_len))
allocate(csca_n_b2(size_len))

! ******************************************************************************************
! *** convert cross-section-equivalent size parameters to other size definitions if required
! ******************************************************************************************
if(size_equ.eq.1) sizepara=sizepara*rat
if(size_equ.eq.2) sizepara=sizepara*rat**3

! ****************************************************************************************
! *** fill the arrays required for the integration (expansion coefficients are done later)
! ****************************************************************************************
do i_size=1,size_len

  r(i_size)=sizepara(i_size)*wavelength_array(i_wavelength)/(2d0*pi)

  select case (contributions(i_contribution)%size_type)

    case ('mono', 'bin')
      n(i_size)=contributions(i_contribution)%n/(contributions(i_contribution)%rmax-contributions(i_contribution)%rmin)

    case ('log_normal')
      n(i_size)=contributions(i_contribution)%n*log_distr(r(i_size),contributions(i_contribution)%size_distr_parameter(1),contributions(i_contribution)%size_distr_parameter(2))

    case ('mod_gamma')
      n(i_size)=contributions(i_contribution)%n*mod_gamma(r(i_size),contributions(i_contribution)%size_distr_parameter(1),contributions(i_contribution)%size_distr_parameter(2),contributions(i_contribution)%size_distr_parameter(3))

    case ('distr_file', 'distr_list')
      n(i_size)=contributions(i_contribution)%n*size_distr_interpolate(r(i_size),mode_aux(contributions(i_contribution)%i_mode)%size_distr_table_len,mode_aux(contributions(i_contribution)%i_mode)%size_distr_table,mode(contributions(i_contribution)%i_mode)%size_distr_table_type)

    case default
      print'(A,I6)', "Error: Unknown size_type "//trim(contributions(i_contribution)%size_type)//" for contribution ",i_contribution
      stop

  end select

  r2_n(i_size)=r(i_size)**2*n(i_size)
  r3_n(i_size)=r(i_size)**3*n(i_size)

  cext_n(i_size)=qext(i_size)*pi*r2_n(i_size)
  csca_n(i_size)=qsca(i_size)*pi*r2_n(i_size)

  ! *** correction because of qext and qsca being normalized by cross section
  if(size_equ.eq.1) then
    cext_n(i_size)=cext_n(i_size)/rat**2
    csca_n(i_size)=csca_n(i_size)/rat**2
  elseif(size_equ.eq.2) then
    cext_n(i_size)=cext_n(i_size)/rat**6
    csca_n(i_size)=csca_n(i_size)/rat**6
  end if

  l0(i_size)=i_size+sum(lmax(1:i_size-1)) ! *** l0(i_size) is the index in the a1-array corresponding to l=0 for i_size

end do

! *************************************************************
! *** now integrate over size distribution of current contribution
! *** and add the results to sums over all contributions
! *************************************************************
call integrate(size_len,r,n,contributions(i_contribution)%rmin,contributions(i_contribution)%rmax,result_tmp)
n_sum=n_sum+result_tmp

call integrate(size_len,r,r2_n,contributions(i_contribution)%rmin,contributions(i_contribution)%rmax,result_tmp)
r2_n_sum=r2_n_sum+result_tmp
if(size_equ.eq.0) then
  cs_sum=cs_sum+pi*result_tmp
elseif(size_equ.eq.1) then
  cs_sum=cs_sum+pi*result_tmp*rat**(-2)
elseif(size_equ.eq.2) then
  cs_sum=cs_sum+pi*result_tmp*rat**(-6)
end if

call integrate(size_len,r,r3_n,contributions(i_contribution)%rmin,contributions(i_contribution)%rmax,result_tmp)
r3_n_sum=r3_n_sum+result_tmp
if(size_equ.eq.0) then
  vol_sum=vol_sum+4d0/3d0*pi*result_tmp*rat**3
  mass_sum=mass_sum+4d0/3d0*pi*contributions(i_contribution)%density*result_tmp*rat**3
elseif(size_equ.eq.1) then
  vol_sum=vol_sum+4d0/3d0*pi*result_tmp
  mass_sum=mass_sum+4d0/3d0*pi*contributions(i_contribution)%density*result_tmp
elseif(size_equ.eq.2) then
  vol_sum=vol_sum+4d0/3d0*pi*result_tmp*rat**(-6)
  mass_sum=mass_sum+4d0/3d0*pi*contributions(i_contribution)%density*result_tmp*rat**(-6)
end if

call integrate(size_len,r,cext_n,contributions(i_contribution)%rmin,contributions(i_contribution)%rmax,result_tmp)
cext_sum=cext_sum+result_tmp

call integrate(size_len,r,csca_n,contributions(i_contribution)%rmin,contributions(i_contribution)%rmax,result_tmp)
csca_sum=csca_sum+result_tmp

! *******************************************************
! *** finally the expansion coefficients are integrated
! *******************************************************
do l=0,lmax_max

  csca_n_a1=0;csca_n_a2=0;csca_n_a3=0;csca_n_a4=0;csca_n_b1=0;csca_n_b2=0
  i_size_min=size_len
  do i_size=1,size_len
    if (lmax(i_size).ge.l) then
      csca_n_a1(i_size)=csca_n(i_size)*a1(l0(i_size)+l)
      csca_n_a2(i_size)=csca_n(i_size)*a2(l0(i_size)+l)
      csca_n_a3(i_size)=csca_n(i_size)*a3(l0(i_size)+l)
      csca_n_a4(i_size)=csca_n(i_size)*a4(l0(i_size)+l)
      csca_n_b1(i_size)=csca_n(i_size)*b1(l0(i_size)+l)
      csca_n_b2(i_size)=csca_n(i_size)*b2(l0(i_size)+l)
      i_size_min=min(i_size_min,max(1,i_size-1))
    end if
  end do

  rmin_tmp=max(r(i_size_min),contributions(i_contribution)%rmin) ! *** rmin adjusted to where expansion coefficients start to deviate from zero

  call integrate(size_len-i_size_min+1,r(i_size_min:size_len),csca_n_a1(i_size_min:size_len),rmin_tmp,contributions(i_contribution)%rmax,result_tmp)
  a1_sum(l)=a1_sum(l)+result_tmp
  call integrate(size_len-i_size_min+1,r(i_size_min:size_len),csca_n_a2(i_size_min:size_len),rmin_tmp,contributions(i_contribution)%rmax,result_tmp)
  a2_sum(l)=a2_sum(l)+result_tmp
  call integrate(size_len-i_size_min+1,r(i_size_min:size_len),csca_n_a3(i_size_min:size_len),rmin_tmp,contributions(i_contribution)%rmax,result_tmp)
  a3_sum(l)=a3_sum(l)+result_tmp
  call integrate(size_len-i_size_min+1,r(i_size_min:size_len),csca_n_a4(i_size_min:size_len),rmin_tmp,contributions(i_contribution)%rmax,result_tmp)
  a4_sum(l)=a4_sum(l)+result_tmp
  call integrate(size_len-i_size_min+1,r(i_size_min:size_len),csca_n_b1(i_size_min:size_len),rmin_tmp,contributions(i_contribution)%rmax,result_tmp)
  b1_sum(l)=b1_sum(l)+result_tmp
  call integrate(size_len-i_size_min+1,r(i_size_min:size_len),csca_n_b2(i_size_min:size_len),rmin_tmp,contributions(i_contribution)%rmax,result_tmp)
  b2_sum(l)=b2_sum(l)+result_tmp

end do

deallocate(a1,a2,a3,a4,b1,b2,lmax)
deallocate(sizepara,qext,qsca,csca_n,cext_n)
deallocate(csca_n_a1,csca_n_a2,csca_n_a3,csca_n_a4,csca_n_b1,csca_n_b2,r,n,r2_n,r3_n)

end subroutine add_contribution
