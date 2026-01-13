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

subroutine init_size

use module_input
use module_aux_var

implicit none

integer :: i,j,i_mode

real(dp) :: r,n

logical :: file_exists

character(len=:), allocatable :: size_distr_list_tmp

do i_mode=1,num_modes

  if(mode(i_mode)%size_type.eq."distr_file".or.mode(i_mode)%size_type.eq."bin_file") then

    inquire(file=mode(i_mode)%size_distr_file,exist=file_exists)
    if(.not.file_exists) then
      print'(A)',"Error: Could not find size distr_file "//trim(mode(i_mode)%size_distr_file)
      stop
    end if

    ! *** first count number of lines
    mode_aux(i_mode)%size_distr_table_len=0
    open(1,file=mode(i_mode)%size_distr_file)

    do while(.true.)

      read(1,*,err=101,end=101) r,n

      mode_aux(i_mode)%size_distr_table_len=mode_aux(i_mode)%size_distr_table_len+1

    end do

    101 rewind(1)

    allocate(mode_aux(i_mode)%size_distr_table(mode_aux(i_mode)%size_distr_table_len,2))

    do i=1,mode_aux(i_mode)%size_distr_table_len

      read(1,*,err=102,end=102) mode_aux(i_mode)%size_distr_table(i,:)

      if(i.gt.1.and.mode_aux(i_mode)%size_distr_table(i,1).le.mode_aux(i_mode)%size_distr_table(i-1,1)) then
        print'(A)',"Error: Sizes not in ascending order in size distr_file "//trim(mode(i_mode)%size_distr_file)
        stop
      end if

    end do

    102 close(1)

  end if

  if(mode(i_mode)%size_type.eq."distr_list") then

    size_distr_list_tmp=mode(i_mode)%size_distr_list

    mode_aux(i_mode)%size_distr_table_len=0
    do while(.true.)

      read(size_distr_list_tmp,*,end=103,err=103)
      i=index(size_distr_list_tmp," ")+1
      size_distr_list_tmp=adjustl(trim(size_distr_list_tmp(i:)))

      read(size_distr_list_tmp,*,end=103,err=103)
      i=index(size_distr_list_tmp," ")+1
      size_distr_list_tmp=adjustl(trim(size_distr_list_tmp(i:)))

      mode_aux(i_mode)%size_distr_table_len=mode_aux(i_mode)%size_distr_table_len+1

      if(i.eq.1.or.len(size_distr_list_tmp).lt.2) goto 104

    end do

    103 print'(A)','Error while reading size distr_list.'
    stop

    105 print'(A)','Error while reading size distr_list. Please check if for each radius a size distribution value is given. '
    stop

    104 size_distr_list_tmp=mode(i_mode)%size_distr_list

    allocate(mode_aux(i_mode)%size_distr_table(mode_aux(i_mode)%size_distr_table_len,2))

    do i=1,mode_aux(i_mode)%size_distr_table_len

      read(size_distr_list_tmp,*,end=105,err=105) mode_aux(i_mode)%size_distr_table(i,1),mode_aux(i_mode)%size_distr_table(i,2)
      j=index(size_distr_list_tmp," ")+1
      size_distr_list_tmp=adjustl(trim(size_distr_list_tmp(j:)))
      j=index(size_distr_list_tmp," ")+1
      size_distr_list_tmp=adjustl(trim(size_distr_list_tmp(j:)))

      if(i.gt.1.and.mode_aux(i_mode)%size_distr_table(i,1).le.mode_aux(i_mode)%size_distr_table(i-1,1)) then
        print'(A)',"Error: Sizes not in ascending order in size distr_list"
        stop
      end if

    end do

  end if

  if(mode(i_mode)%size_type.eq."distr_file".or.mode(i_mode)%size_type.eq."bin_file".or.mode(i_mode)%size_type.eq."distr_list") then

    ! *** if option diameter is used
    if(diameter.eqv..True.) then

      ! *** convert values from diameter to radius (which is used internally)
      mode_aux(i_mode)%size_distr_table(:,1)=mode_aux(i_mode)%size_distr_table(:,1)*0.5d0

      ! *** size_distr_table_type dndd etc is changed to dndr
      i=len_trim(mode(i_mode)%size_distr_table_type)
      if(mode(i_mode)%size_distr_table_type(i:i).eq."d") mode(i_mode)%size_distr_table_type(i:i)="r"

      ! *** convert values dXdd tp dXdr (in case of dlnd or dlogd no conversion is needed)
      if(mode(i_mode)%size_distr_table_type.eq."dndr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*2d0
      if(mode(i_mode)%size_distr_table_type.eq."dadr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*2d0
      if(mode(i_mode)%size_distr_table_type.eq."dvdr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*2d0

    end if

    ! *** mode(i_mode)%n will be multiplied with the distribution from the file
    mode(i_mode)%n=1d0
    mode(i_mode)%rmin=mode_aux(i_mode)%size_distr_table(1,1)
    mode(i_mode)%rmax=mode_aux(i_mode)%size_distr_table(mode_aux(i_mode)%size_distr_table_len,1)

  end if

end do


end subroutine
