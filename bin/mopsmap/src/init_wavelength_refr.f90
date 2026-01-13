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

subroutine init_wavelength_refr

use module_input
use module_dbindex
use module_aux_var

implicit none

real(dp) :: line(3)
integer :: i_wavelength,i_mode,i_line,n_lines

integer :: i1,i2
real(dp) :: w1,w2

real(dp),allocatable :: tmp_table(:,:)

logical :: file_exists


! **************************
! *** initialize wavelengths
! **************************
num_wavelengths=0
if(use_wavelength_from_refr_file.or.use_wavelength_file) then

  ! check if wavelengths from refractive index file should be used
  if(use_wavelength_from_refr_file) then

    wavelength_file=mode(1)%refr_file

    if(wavelength_file.eq."") then

      print'(A)',"Error: Option 'wavelength from_refrac_file' needs to be combined with option 'refrac file'"
      stop

    end if

  end if

  ! *** read wavelengths from file
  inquire(file=wavelength_file,exist=file_exists)
  if(.not.file_exists) then
    print'(A)',"Error: Could not find wavelength file "//trim(wavelength_file)
    stop
  end if

  open(1,file=wavelength_file)

  do while(.true.)

    read(1,*,end=100) line(1)
    num_wavelengths=num_wavelengths+1

  end do

  100 rewind(1)
  allocate(wavelength_array(num_wavelengths))

  do i_wavelength=1,num_wavelengths

    read(1,*) line(1)
    wavelength_array(i_wavelength)=line(1)

  end do

  close(1)

  if(num_wavelengths.lt.1) then

    print'(A)',"Error: No wavelengths found in wavelength file "//trim(wavelength_file)
    stop

  end if

elseif(use_wavelength_range) then

  num_wavelengths=int(1+(1d0+1d-9)*(wavelength_max-wavelength_min)/wavelength_step)
  allocate(wavelength_array(num_wavelengths))

  do i_wavelength=1,num_wavelengths

    wavelength_array(i_wavelength)=wavelength_min+wavelength_step*(i_wavelength-1)

  end do

  ! *** necessary to exclude error due to wavelength outside refractive index range
  if(abs(wavelength_array(num_wavelengths)-wavelength_max)/wavelength_max.lt.1d-9) wavelength_array(num_wavelengths)=wavelength_max

elseif(use_wavelength_list) then

  ! *** can be improved
  allocate(wavelength_array(1000))
  do while(.true.)

    num_wavelengths=num_wavelengths+1

    read(wavelength_list,*,end=103,err=103) line(1)

    wavelength_array(num_wavelengths)=line(1)

    i1=index(wavelength_list," ")+1
    wavelength_list=adjustl(trim(wavelength_list(i1:)))

    if(i1.eq.1.or.len(wavelength_list).lt.2) goto 101

  end do

  103 print'(A)','Error while reading wavelength list.'
  stop

  101 continue

else

  num_wavelengths=1
  allocate(wavelength_array(num_wavelengths))
  wavelength_array(1)=wavelength

end if

! *** sort wavelengths in ascending order ('bubble sort')
do i1=num_wavelengths-1,1,-1
  do i2=1,i1
    if(wavelength_array(i2).gt.wavelength_array(i2+1)) then
      w1=wavelength_array(i2)
      wavelength_array(i2)=wavelength_array(i2+1)
      wavelength_array(i2+1)=w1
    end if
  end do
end do

! *** some checks
if(num_wavelengths.eq.0) then

  print'(A)',"Error: No wavelengths given in input file!"
  stop

end if

! ******************************************************
! *** initialize refractive indices at these wavelengths
! ******************************************************
allocate(mode_aux(num_modes))
do i_mode=1,num_modes

  allocate(mode_aux(i_mode)%mreal(num_wavelengths))
  allocate(mode_aux(i_mode)%mimag(num_wavelengths))

  if(mode(i_mode)%refr_file.ne."") then

    if(i_mode.gt.1.and.mode(i_mode)%refr_file.eq.mode(i_mode-1)%refr_file) then

      ! *** just copy the refractive index from the previous mode if it uses the same file
      mode_aux(i_mode)%mreal=mode_aux(i_mode-1)%mreal
      mode_aux(i_mode)%mimag=mode_aux(i_mode-1)%mimag

    else

      ! *** read refractive index from file
      inquire(file=mode(i_mode)%refr_file,exist=file_exists)
      if(.not.file_exists) then
        print'(A)',"Error: Could not find refractive index file "//trim(mode(i_mode)%refr_file)
        stop
      end if

      n_lines=0
      open(1,file=mode(i_mode)%refr_file)
      do while(.true.)
        read(1,*,end=102) line ! line is an array of 3 floats
        n_lines=n_lines+1
      end do
      102 rewind(1)
      allocate(tmp_table(n_lines,3))
      do i_line=1,n_lines
        read(1,*) tmp_table(i_line,:) ! line is an array of 3 floats
        if(i_line.gt.1.and.tmp_table(i_line,1).le.tmp_table(i_line-1,1)) then
          print'(A)',"Error: Wavelengths in refractive index file "//trim(mode(i_mode)%refr_file)//" are not in ascending order!"
          stop
        end if
      end do
      close(1)

      if(tmp_table(1,1).gt.wavelength_array(1).or.tmp_table(n_lines,1).lt.wavelength_array(num_wavelengths)) then
        print'(A,F10.5,A,F10.5,A)',"Error: Refractive index file "//trim(mode(i_mode)%refr_file)//" does not cover required wavelength range from ",wavelength_array(1)," to ",wavelength_array(num_wavelengths)
        stop
      end if

      do i_wavelength=1,num_wavelengths
        call interpolate_linear(n_lines,tmp_table(:,1),wavelength_array(i_wavelength),i1,i2,w1,w2)
        mode_aux(i_mode)%mreal(i_wavelength)=w1*tmp_table(i1,2)+w2*tmp_table(i2,2)
        mode_aux(i_mode)%mimag(i_wavelength)=w1*tmp_table(i1,3)+w2*tmp_table(i2,3)
      end do

      deallocate(tmp_table)

    end if

  else

    ! *** use given wavelength-independent refractive index
    mode_aux(i_mode)%mreal(:)=mode(i_mode)%mreal
    mode_aux(i_mode)%mimag(:)=mode(i_mode)%mimag

  end if

end do

! *** if nonabs_fraction is larger than zero increase the imaginary part of the absorbing fraction so that the average imaginary part remains unchanged
do i_mode=1,num_modes

  mode_aux(i_mode)%mimag(:)=mode_aux(i_mode)%mimag(:)/(1d0-mode(i_mode)%nonabs_fraction)

end do

end subroutine
