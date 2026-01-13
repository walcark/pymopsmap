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

subroutine read_dbindex

use netcdf
use module_dbindex

implicit none

integer :: status,ncid
integer :: mreal_dimid,mimag_dimid,np_dimid,eps_dimid
integer :: mreal_varid,mimag_varid,np_varid,eps_varid,max_sizepara_varid,rat_varid
character(len=1000) :: index_file

call get_nc_filename(-1,0d0,0d0,0d0,3,index_file)

status=nf90_open(index_file, 0, ncid)
if(status.ne.nf90_noerr) then
  print'(A)',"Error '"//trim(nf90_strerror(status))//"' during opening of index file '"//trim(index_file)//"'"
  stop
end if

status=nf90_inq_dimid(ncid,'np', np_dimid)
status=nf90_inq_dimid(ncid,'eps', eps_dimid)
status=nf90_inq_dimid(ncid,'mreal', mreal_dimid)
status=nf90_inq_dimid(ncid,'mimag', mimag_dimid)

status=nf90_inq_varid(ncid,'np', np_varid)
status=nf90_inq_varid(ncid,'eps', eps_varid)
status=nf90_inq_varid(ncid,'mreal', mreal_varid)
status=nf90_inq_varid(ncid,'mimag', mimag_varid)
status=nf90_inq_varid(ncid,'max_sizepara', max_sizepara_varid)
status=nf90_inq_varid(ncid,'rat', rat_varid)

status=nf90_inquire_dimension(ncid,np_dimid,len=num_avail_np)
status=nf90_inquire_dimension(ncid,eps_dimid,len=num_avail_eps)
status=nf90_inquire_dimension(ncid,mreal_dimid,len=num_avail_mreal)
status=nf90_inquire_dimension(ncid,mimag_dimid,len=num_avail_mimag)

allocate(avail_np(num_avail_np))
allocate(avail_eps(num_avail_eps))
allocate(avail_mreal(num_avail_mreal))
allocate(avail_mimag(num_avail_mimag))
allocate(max_sizepara(num_avail_np,num_avail_eps,num_avail_mreal,num_avail_mimag))
allocate(rat(num_avail_np,num_avail_eps))

status=nf90_get_var(ncid,np_varid,avail_np)
status=nf90_get_var(ncid,eps_varid,avail_eps)
status=nf90_get_var(ncid,mreal_varid,avail_mreal)
status=nf90_get_var(ncid,mimag_varid,avail_mimag)
status=nf90_get_var(ncid,max_sizepara_varid,max_sizepara)
status=nf90_get_var(ncid,rat_varid,rat)

status=nf90_close(ncid)

! check whether prolate and oblate eps are ordered as expected
if(avail_eps(1).ne.0.2d0.or.avail_eps(16).ne.1d0.or.avail_eps(31).ne.5d0.or.num_avail_eps.ne.31) then
  print*,"Error: Content of variable eps of data set index file differs from what is expected by the code. Please contact the code authors for clarification."
  stop
end if

allocate(avail_eps_width(num_avail_eps))
avail_eps_width(1:15)=(/0.2d0,0.4d0,0.4d0,0.4d0,0.4d0,0.3d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0,0.2d0/) ! width with respect to aspect ratio not eps
avail_eps_width(16)=0.2d0
avail_eps_width(31:17:-1)=avail_eps_width(1:15)

! check whether first two elements of avail_np are -1 and -6; this is assumed in subroutine make_contributions
if(avail_np(1).ne.-1.or.avail_np(2).ne.-7) then
  print*,"Error: Content of variable np of data set index file differs from what is expected by the code. Please contact the code authors for clarification."
  stop
end if

end subroutine
