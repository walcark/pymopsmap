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

program create_index

use netcdf
use module_input, only: scatlib_path,dp

implicit none

integer,parameter :: num_np=8
integer,parameter :: num_eps=31
integer,parameter :: num_mreal=36
integer,parameter :: num_mimag=21

integer :: np(num_np)
real(dp) :: eps(num_eps)
real(dp) :: mreal(num_mreal)
real(dp) :: mimag(num_mimag)
real(dp) :: max_sizepara(num_np,num_eps,num_mreal,num_mimag)
real(dp) :: min_sizepara(num_np,num_eps,num_mreal,num_mimag)
real(dp) :: rat(num_np,num_eps)
real(dp) :: tmp

integer :: num_particles(num_np,num_eps,num_mreal,num_mimag)

logical :: file_found(num_np,num_eps,num_mreal,num_mimag)

logical :: is_avail_np(num_np)
logical :: is_avail_eps(num_eps)
logical :: is_avail_mreal(num_mreal)
logical :: is_avail_mimag(num_mimag)

integer :: num_np_avail
integer :: num_eps_avail
integer :: num_mreal_avail
integer :: num_mimag_avail

integer :: np_avail(num_np)
real(dp) :: eps_avail(num_eps)
real(dp) :: mreal_avail(num_mreal)
real(dp) :: mimag_avail(num_mimag)
real(dp) :: max_sizepara_avail(num_np,num_eps,num_mreal,num_mimag)
real(dp) :: min_sizepara_avail(num_np,num_eps,num_mreal,num_mimag)
real(dp) :: rat_avail(num_np,num_eps)

integer :: i_np,i_eps,i_mreal,i_mimag

character(len=1000) :: filename

integer :: ncid,status

integer :: sizepara_dimid,sizepara_varid
integer :: sizepara_len
real(dp),allocatable :: sizepara(:)

integer :: np_dimid,eps_dimid,mreal_dimid,mimag_dimid
integer :: np_varid,eps_varid,mreal_varid,mimag_varid
integer :: var_dim(4)
integer :: max_sizepara_id
integer :: rat_id,tmp_id

scatlib_path="./"

np=(/-1,-7,-109,-110,-111,-118,-126,-127/)
eps(17:31)=(/1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0,3.4d0,3.8d0,4.2d0,4.6d0,5.0d0/)
eps(16)=1.
eps(15:1:-1)=1d0/eps(17:31)
mreal=(/0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.00d0,1.04d0,1.08d0,1.12d0,1.16d0,1.20d0,1.24d0,1.28d0,1.32d0,1.36d0,1.40d0,1.44d0,1.48d0,1.52d0,1.56d0,1.60d0,1.64d0,1.68d0,1.76d0,1.84d0,1.92d0,2.00d0,2.20d0,2.40d0,2.60d0,2.80d0,3.00d0/)
mimag=(/0.0d0,0.00053750d0,0.00107500d0,0.00152028d0,0.00215000d0,0.00304056d0,0.00430000d0,0.00608112d0,0.00860000d0,0.01216224d0,0.01720000d0,0.02432447d0,0.03440000d0,0.04864895d0,0.06880000d0,0.09729789d0,0.13760000d0,0.27520000d0,0.55040000d0,1.10080000d0,2.20160000d0/)

max_sizepara=0d0
max_sizepara_avail=0d0
min_sizepara=1000d0
min_sizepara_avail=1000d0
rat=1d0
rat_avail=1d0
num_particles=0

file_found=.false.

is_avail_np=.false.
is_avail_eps=.false.
is_avail_mreal=.false.
is_avail_mimag=.false.


do i_np=1,num_np
  do i_eps=1,num_eps
    do i_mreal=1,num_mreal
      do i_mimag=1,num_mimag

        call get_nc_filename(np(i_np),eps(i_eps),mreal(i_mreal),mimag(i_mimag),1,filename)

        status=nf90_open(trim(filename), 0, ncid)

        if(status.eq.0) then

          ! *** extract maximum sizeparameters
          status=nf90_inq_dimid(ncid,'sizepara',sizepara_dimid)
          status=nf90_inq_varid(ncid,'sizepara',sizepara_varid)

          status=nf90_inquire_dimension(ncid,sizepara_dimid,len=sizepara_len)
          allocate(sizepara(sizepara_len))
          status=nf90_get_var(ncid,sizepara_varid,sizepara)

          max_sizepara(i_np,i_eps,i_mreal,i_mimag)=sizepara(sizepara_len)
          min_sizepara(i_np,i_eps,i_mreal,i_mimag)=sizepara(1)

          deallocate(sizepara)

          ! *** extract rat
          status=nf90_inq_varid(ncid,'rat',rat_id)
          status=nf90_get_var(ncid,rat_id,tmp)

          rat(i_np,i_eps)=tmp

          ! *** check if eps, mreal, mimag agree with filename
          status=nf90_inq_varid(ncid,'eps',tmp_id)
          status=nf90_get_var(ncid,tmp_id,tmp)
          if(np(i_np).gt.-10.and.abs(tmp-eps(i_eps)).gt.0.0001) print*,"eps in file and filename do not agree in ",trim(filename)

          status=nf90_inq_varid(ncid,'mreal',tmp_id)
          status=nf90_get_var(ncid,tmp_id,tmp)
          if(abs(tmp-mreal(i_mreal)).gt.0.0001) print*,"mreal in file and filename do not agree in ",trim(filename)

          status=nf90_inq_varid(ncid,'mimag',tmp_id)
          status=nf90_get_var(ncid,tmp_id,tmp)
          if(abs(tmp-mimag(i_mimag)).gt.0.0001) print*,"mimag in file and filename do not agree in ",trim(filename)

          status=nf90_close(ncid)

          call check_scatdata_file(filename)

          num_particles(i_np,i_eps,i_mreal,i_mimag)=num_particles(i_np,i_eps,i_mreal,i_mimag)+sizepara_len

          file_found(i_np,i_eps,i_mreal,i_mimag)=.true.

          is_avail_eps(i_eps)=.true.
          is_avail_np(i_np)=.true.
          is_avail_mreal(i_mreal)=.true.
          is_avail_mimag(i_mimag)=.true.

        end if

      end do
    end do
  end do
end do

! *** Create arrays with the available parameters
num_np_avail=0
do i_np=1,num_np

  if(is_avail_np(i_np)) then

    num_np_avail=num_np_avail+1
    np_avail(num_np_avail)=np(i_np)

    num_eps_avail=0
    do i_eps=1,num_eps

      if(is_avail_eps(i_eps)) then

        num_eps_avail=num_eps_avail+1
        eps_avail(num_eps_avail)=eps(i_eps)

        rat_avail(num_np_avail,num_eps_avail)=rat(i_np,i_eps)

        num_mreal_avail=0
        do i_mreal=1,num_mreal

          if(is_avail_mreal(i_mreal)) then

            num_mreal_avail=num_mreal_avail+1
            mreal_avail(num_mreal_avail)=mreal(i_mreal)

            num_mimag_avail=0
            do i_mimag=1,num_mimag

              if(is_avail_mimag(i_mimag)) then

                num_mimag_avail=num_mimag_avail+1
                mimag_avail(num_mimag_avail)=mimag(i_mimag)

                max_sizepara_avail(num_np_avail,num_eps_avail,num_mreal_avail,num_mimag_avail)=max_sizepara(i_np,i_eps,i_mreal,i_mimag)
                min_sizepara_avail(num_np_avail,num_eps_avail,num_mreal_avail,num_mimag_avail)=min_sizepara(i_np,i_eps,i_mreal,i_mimag)

              end if
            end do
          end if
        end do
      end if
    end do
  end if
end do

! *** check whether maximum sizeparameter of spheroids from TMM increases with decreasing deviation from sphere
if(num_eps_avail.eq.31) then
  do i_mreal=1,num_mreal_avail
    do i_mimag=1,num_mimag_avail
      do i_eps=2,15
        if(max_sizepara_avail(1,i_eps,i_mreal,i_mimag).lt.max_sizepara_avail(1,i_eps-1,i_mreal,i_mimag)) then
          print*,"max_sizepara problem at ",eps_avail(i_eps),mreal_avail(i_mreal),mimag_avail(i_mimag)
          max_sizepara_avail(1,i_eps-1,i_mreal,i_mimag)=max_sizepara_avail(1,i_eps,i_mreal,i_mimag)
        end if
      end do
      do i_eps=30,17,-1
        if(max_sizepara_avail(1,i_eps,i_mreal,i_mimag).lt.max_sizepara_avail(1,i_eps+1,i_mreal,i_mimag)) then
          print*,"max_sizepara problem at ",eps_avail(i_eps),mreal_avail(i_mreal),mimag_avail(i_mimag)
          max_sizepara_avail(1,i_eps+1,i_mreal,i_mimag)=max_sizepara_avail(1,i_eps,i_mreal,i_mimag)
        end if
      end do
    end do
  end do
end if

! *** create the index file
status=nf90_create('index.nc', nf90_classic_model, ncid)

status=nf90_def_dim(ncid,'np', num_np_avail, np_dimid)
status=nf90_def_dim(ncid,'eps', num_eps_avail, eps_dimid)
status=nf90_def_dim(ncid,'mreal', num_mreal_avail, mreal_dimid)
status=nf90_def_dim(ncid,'mimag', num_mimag_avail, mimag_dimid)

status=nf90_def_var(ncid,'np', nf90_int, (/np_dimid/), np_varid)
status=nf90_def_var(ncid,'eps', nf90_double, (/eps_dimid/), eps_varid)
status=nf90_def_var(ncid,'mreal', nf90_double, (/mreal_dimid/), mreal_varid)
status=nf90_def_var(ncid,'mimag', nf90_double, (/mimag_dimid/), mimag_varid)

var_dim(1)=np_dimid
var_dim(2)=eps_dimid
var_dim(3)=mreal_dimid
var_dim(4)=mimag_dimid

status=nf90_def_var(ncid,'max_sizepara', nf90_double, var_dim, max_sizepara_id)

status=nf90_def_var(ncid,'rat', nf90_double, var_dim(1:2), rat_id)

status=nf90_enddef(ncid)

status=nf90_put_var(ncid,np_varid,np_avail(1:num_np_avail))
status=nf90_put_var(ncid,eps_varid,eps_avail(1:num_eps_avail))
status=nf90_put_var(ncid,mreal_varid,mreal_avail(1:num_mreal_avail))
status=nf90_put_var(ncid,mimag_varid,mimag_avail(1:num_mimag_avail))
status=nf90_put_var(ncid,max_sizepara_id,max_sizepara_avail(1:num_np_avail,1:num_eps_avail,1:num_mreal_avail,1:num_mimag_avail))
status=nf90_put_var(ncid,rat_id,rat_avail(1:num_np_avail,1:num_eps_avail))

status=nf90_close(ncid)

print*,"*******************"
print*,"Number of particles"
print*,"*******************"
print*
print'(A,I8)',"    total: ",sum(num_particles(:,:,:,:))
print*
print*,"grouped by NP:"
do i_np=1,num_np
  if(is_avail_np(i_np)) print'(I9,A,I7,A,I6,A)',np(i_np),": ",sum(num_particles(i_np,:,:,:))," (",count(file_found(i_np,:,:,:))," files)"
end do
print*
print*,"grouped by EPS:"
do i_eps=1,num_eps
  if(is_avail_eps(i_eps)) print'(F9.5,A,I7,A,I6,A)',eps(i_eps),": ",sum(num_particles(:,i_eps,:,:))," (",count(file_found(:,i_eps,:,:))," files)"
end do
print*
print*,"grouped by MREAL:"
do i_mreal=1,num_mreal
  if(is_avail_mreal(i_mreal)) print'(F9.5,A,I7,A,I6,A)',mreal(i_mreal),": ",sum(num_particles(:,:,i_mreal,:))," (",count(file_found(:,:,i_mreal,:))," files)"
end do
print*
print*,"grouped by MIMAG:"
do i_mimag=1,num_mimag
  if(is_avail_mimag(i_mimag)) print'(F9.5,A,I7,A,I6,A)',mimag(i_mimag),": ",sum(num_particles(:,:,:,i_mimag))," (",count(file_found(:,:,:,i_mimag))," files)"
end do
print*

end program
