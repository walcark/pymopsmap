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

subroutine init_shape

use module_input
use module_dbindex
use module_aux_var

implicit none

integer :: i_mode,i_eps

integer :: i,i1,i2
real(dp) :: w1,w2

real(dp) :: log_distr_ar

real(dp) :: ratio_prolate
real(dp) :: eps,weight
integer :: np

character(len=1000) :: filename, arg_str

logical :: file_exists

do i_mode=1,num_modes

  allocate(mode_aux(i_mode)%eps_distr(num_avail_eps))
  mode_aux(i_mode)%eps_distr=0d0
  allocate(mode_aux(i_mode)%irregular_np_distr(num_avail_np))
  mode_aux(i_mode)%irregular_np_distr=0d0

  if(mode(i_mode)%shape_type.eq."sphere") then

    call interpolate_linear(num_avail_eps,avail_eps,1d0,i1,i2,w1,w2)
    mode_aux(i_mode)%eps_distr(i1)=w1
    mode_aux(i_mode)%eps_distr(i2)=w2

  elseif(mode(i_mode)%shape_type.eq."spheroid") then

    if((mode(i_mode)%shape_parameter(1)-0d0).lt.1d-4) then ! oblate
      call interpolate_linear(num_avail_eps,avail_eps,mode(i_mode)%shape_parameter(2),i1,i2,w1,w2)
      mode_aux(i_mode)%eps_distr(i1)=w1
      mode_aux(i_mode)%eps_distr(i2)=w2
    elseif((mode(i_mode)%shape_parameter(1)-1d0).lt.1d-4) then ! prolate
      call interpolate_linear(num_avail_eps,-1./avail_eps,-mode(i_mode)%shape_parameter(2),i1,i2,w1,w2)
      mode_aux(i_mode)%eps_distr(i1)=w1
      mode_aux(i_mode)%eps_distr(i2)=w2
    else
      print*,"error 2 in init_shape"
      stop
    end if

  elseif(mode(i_mode)%shape_type.eq."spheroid_log_normal") then

    if(mode(i_mode)%shape_parameter(2).lt.1) then ! oblate

      do i_eps=17,31
        mode_aux(i_mode)%eps_distr(i_eps)=avail_eps_width(i_eps)*log_distr_ar(avail_eps(i_eps),mode(i_mode)%shape_parameter(3),mode(i_mode)%shape_parameter(4))
      end do
      mode_aux(i_mode)%eps_distr(17:31)=mode(i_mode)%shape_parameter(1)*(1.-mode(i_mode)%shape_parameter(2))*mode_aux(i_mode)%eps_distr(17:31)/sum(mode_aux(i_mode)%eps_distr(17:31))

    end if

    if(mode(i_mode)%shape_parameter(2).gt.0) then ! prolate

      do i_eps=1,15
        mode_aux(i_mode)%eps_distr(i_eps)=avail_eps_width(i_eps)*log_distr_ar(1./avail_eps(i_eps),mode(i_mode)%shape_parameter(3),mode(i_mode)%shape_parameter(4))
      end do
      mode_aux(i_mode)%eps_distr(1:15)=mode(i_mode)%shape_parameter(1)*mode(i_mode)%shape_parameter(2)*mode_aux(i_mode)%eps_distr(1:15)/sum(mode_aux(i_mode)%eps_distr(1:15))

    end if

    mode_aux(i_mode)%eps_distr(16)=1.-mode(i_mode)%shape_parameter(1)

  elseif(mode(i_mode)%shape_type.eq."spheroid_distr_file") then

    inquire(file=mode(i_mode)%shape_distr_file,exist=file_exists)
    if(.not.file_exists) then
      print'(A)',"Error: Could not find distr_file "//trim(mode(i_mode)%shape_distr_file)
      stop
    end if

    open(1,file=mode(i_mode)%shape_distr_file)

    read(1,*,err=103,end=103) ratio_prolate ! read ratio of prolate spheroids

    mode_aux(i_mode)%eps_distr(:)=0

    do while(.true.)

      read(1,*,err=103,end=103) eps,weight

      if(eps.lt.1.or.eps.gt.5) then
        print'(A)',"Error while parsing distr_file "//trim(mode(i_mode)%shape_distr_file)//": Aspect ratio out of range. Only aspect ratios from 1 to 5 are allowed."
        stop
      end if

      if(ratio_prolate.lt.1) then ! oblate
        call interpolate_linear(num_avail_eps,avail_eps,eps,i1,i2,w1,w2)
        mode_aux(i_mode)%eps_distr(i1)=mode_aux(i_mode)%eps_distr(i1)+w1*(1-ratio_prolate)*weight
        mode_aux(i_mode)%eps_distr(i2)=mode_aux(i_mode)%eps_distr(i2)+w2*(1-ratio_prolate)*weight
      end if

      if(ratio_prolate.gt.0) then ! prolate
        call interpolate_linear(num_avail_eps,-1./avail_eps,-eps,i1,i2,w1,w2)
        mode_aux(i_mode)%eps_distr(i1)=mode_aux(i_mode)%eps_distr(i1)+w1*ratio_prolate*weight
        mode_aux(i_mode)%eps_distr(i2)=mode_aux(i_mode)%eps_distr(i2)+w2*ratio_prolate*weight
      end if

    end do

    103 close(1)

    if(sum(mode_aux(i_mode)%eps_distr(:)).le.0) then

      print'(A)',"Error while parsing distr_file "//trim(mode(i_mode)%shape_distr_file)
      stop

    end if

    mode_aux(i_mode)%eps_distr(:)=mode_aux(i_mode)%eps_distr(:)/sum(mode_aux(i_mode)%eps_distr(:))

  ! *** if only a single irregular shape is used
  elseif(mode(i_mode)%shape_type.eq."irregular") then

    ! *** search for index in data set
    do i=1,num_avail_np
      if((abs(avail_np(i)-mode(i_mode)%shape_parameter(1))).lt.1d-4) i1=i
    end do

    if(i1.gt.0) then
      mode_aux(i_mode)%irregular_np_distr(i1)=1
    else
      print'(A,F10.5,A)',"Error NP ",mode(i_mode)%shape_parameter(1)," not found in data set."
      stop
    end if

  elseif(mode(i_mode)%shape_type.eq."irregular_distr_file") then

    continue  ! *** work is done below

  else

    print*,"error 1 in init_shape"
    stop

  end if

  ! *** read distribution of irregular shapes from file if required
  if(mode(i_mode)%shape_irregular_overlay.or.mode(i_mode)%shape_type.eq."irregular_distr_file") then

    if(mode(i_mode)%shape_irregular_overlay) then
      filename=mode(i_mode)%shape_irregular_overlay_file
    else
      filename=mode(i_mode)%shape_distr_file
    end if

    inquire(file=filename,exist=file_exists)
    if(.not.file_exists) then
      print'(A)',"Error: Could not find irregular shape file "//trim(mode(i_mode)%shape_irregular_overlay_file)
      stop
    end if

    open(1,file=filename)
    do while(.true.)
      np=0
      read(1,*,err=104,end=104) arg_str,weight
      if(arg_str.eq."A".or.arg_str.eq."-109") np=-109
      if(arg_str.eq."B".or.arg_str.eq."-110") np=-110
      if(arg_str.eq."C".or.arg_str.eq."-111") np=-111
      if(arg_str.eq."D".or.arg_str.eq."-118") np=-118
      if(arg_str.eq."E".or.arg_str.eq."-126") np=-126
      if(arg_str.eq."F".or.arg_str.eq."-127") np=-127
      i1=0
      do i=1,num_avail_np
        if(avail_np(i).eq.np) i1=i
      end do
      if(i1.gt.0) then
        mode_aux(i_mode)%irregular_np_distr(i1)=weight
      else
        print'(A)',"NP not found while parsing irregular shape file "//trim(mode(i_mode)%shape_irregular_overlay_file)
        stop
      end if
    end do

    104 close(1)

    if(sum(mode_aux(i_mode)%irregular_np_distr(:)).le.0) then

      print'(A)',"Error while parsing irregular shape file "//trim(mode(i_mode)%shape_irregular_overlay_file)
      stop

    end if

    ! *** normalize sum to one

    mode_aux(i_mode)%irregular_np_distr(:)=mode_aux(i_mode)%irregular_np_distr(:)/sum(mode_aux(i_mode)%irregular_np_distr)

  end if


end do


end subroutine
