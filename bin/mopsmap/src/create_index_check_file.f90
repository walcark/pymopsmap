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

subroutine check_scatdata_file(filename)

use netcdf

implicit none

character(*) :: filename

integer :: status,ncid
integer :: np_id,eps_id,rat_id
integer :: np
real :: eps,rat,rat2

integer :: i,dim_id,var_id,dim_len
integer,allocatable :: array(:)
real,allocatable :: array1(:),array2(:)

! *** file exists?
status=nf90_open(trim(filename), 0, ncid)
if(status.ne.nf90_noerr) print*,"could not find "//trim(filename)

! *** check for right rat
status=nf90_inq_varid(ncid,'np',np_id)
status=nf90_get_var(ncid,np_id,np)
status=nf90_inq_varid(ncid,'eps',eps_id)
status=nf90_get_var(ncid,eps_id,eps)

status=nf90_inq_varid(ncid,'rat',rat_id)
if(status.ne.nf90_noerr) print*,"variable rat not found in "//trim(filename)

status=nf90_get_var(ncid,rat_id,rat)

if(np.eq.-1.or.np.eq.-3.or.np.eq.-4.or.np.eq.-5) then
  if(eps.gt.1) then
    rat2=sqrt(1.-1./eps**2)
    rat2=0.25*(2.*eps**(2./3.)+eps**(-4./3.)*LOG((1.+rat2)/(1.-rat2))/rat2)
    rat2=1./sqrt(rat2)
  else
    rat2=sqrt(1.-eps**2)
    rat2=0.5*(eps**(2./3.)+eps**(-1./3.)*asin(rat2)/rat2)
    rat2=1./sqrt(rat2)
  end if
  if(abs(rat-rat2).gt.0.00001) print*,"rat not correct in "//trim(filename)
end if


! *** check whether lmax<0 exists
status=nf90_inq_dimid(ncid,'sizepara',dim_id)
status=nf90_inquire_dimension(ncid,dim_id,len=dim_len)
allocate(array(dim_len))
status=nf90_inq_varid(ncid,'lmax',var_id)
if(status.eq.0) then
  status=nf90_get_var(ncid,var_id,array)
  do i=1,dim_len
    if(array(i).lt.0) print*,"lmax<0 in "//trim(filename)
  end do
end if
deallocate(array)

! *** some tests for T-matrix
if(np.ne.-1.or.eps.eq.1) goto 123

status=nf90_inq_dimid(ncid,'sizepara',dim_id)
status=nf90_inquire_dimension(ncid,dim_id,len=dim_len)
allocate(array1(dim_len))
status=nf90_inq_varid(ncid,'qext',var_id)
status=nf90_get_var(ncid,var_id,array1)
allocate(array2(dim_len))
status=nf90_inq_varid(ncid,'qsca',var_id)
status=nf90_get_var(ncid,var_id,array2)

! *** check for omega > 1.001
do i=1,dim_len
  if(array2(i)/array1(i).gt.1.001) print'(A,F8.4,A)',"omega=",array2(i)/array1(i), " in "//trim(filename)
end do

deallocate(array1,array2)

123 status=nf90_close(ncid)

end subroutine
