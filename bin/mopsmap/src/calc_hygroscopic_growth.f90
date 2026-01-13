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

subroutine calc_hygroscopic_growth

use module_input
use module_aux_var

implicit none

real(dp) :: growth_factor
integer :: i_mode,i_wavelength,i_line,n_lines
logical :: file_exists

real(dp) :: line(3)
real(dp),allocatable :: tmp_table(:,:)
integer :: i1,i2
real(dp) :: w1,w2

real(dp) :: water_mreal(num_wavelengths)
real(dp) :: water_mimag(num_wavelengths)

if(relative_humidity.gt.0) then

  if(relative_humidity.gt.99) then
    print'(A)',"Error: Relative humidity > 99% not supported."
    stop
  end if

  ! *** read the water refractive index
  inquire(file=water_refr_file,exist=file_exists)
  if(.not.file_exists) then
    print'(A)',"Error: Could not find water refractive index file "//trim(water_refr_file)
    stop
  end if

  n_lines=0
  open(1,file=water_refr_file)
  do while(.true.) ! first count the number of lines
    read(1,*,end=101) line
    n_lines=n_lines+1
  end do
  101 rewind(1)
  allocate(tmp_table(n_lines,3))
  do i_line=1,n_lines
    read(1,*) tmp_table(i_line,:) ! then save the data in an array
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
    water_mreal(i_wavelength)=w1*tmp_table(i1,2)+w2*tmp_table(i2,2)
    water_mimag(i_wavelength)=w1*tmp_table(i1,3)+w2*tmp_table(i2,3)
  end do

  ! *** go through all modes
  do i_mode=1,num_modes

    ! *** check if kappa was provided in the input file
    if(mode(i_mode)%kappa.lt.0) then
      print'(A,I4)',"Error: No kappa value provided for mode ",i_mode
      stop
    end if

    ! *** calculate growth factor
    growth_factor=(1d0+mode(i_mode)%kappa*relative_humidity/(100d0-relative_humidity))**(1d0/3d0)

    ! *** modify size distribution
    mode(i_mode)%rmin=mode(i_mode)%rmin*growth_factor
    mode(i_mode)%rmax=mode(i_mode)%rmax*growth_factor

    if(mode(i_mode)%size_type.eq.'log_normal') then
      mode(i_mode)%size_distr_parameter(1)=mode(i_mode)%size_distr_parameter(1)*growth_factor
    end if

    if(mode(i_mode)%size_type.eq.'mod_gamma') then
      print'(A,I4)',"Error: The combination of hygroscopic growth with the mod_gamma size distribution is not implemented."
      stop
    end if

    if(allocated(mode_aux(i_mode)%size_distr_table)) then

      ! *** multiply the radii with the growth_factor
      mode_aux(i_mode)%size_distr_table(:,1)=mode_aux(i_mode)%size_distr_table(:,1)*growth_factor

      ! *** adjust the actual size_distr values, first for the change in dr
      if(mode(i_mode)%size_distr_table_type.eq."dndr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)/growth_factor
      if(mode(i_mode)%size_distr_table_type.eq."dadr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)/growth_factor
      if(mode(i_mode)%size_distr_table_type.eq."dvdr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)/growth_factor
      ! *** then for the change in da and dv
      if(mode(i_mode)%size_distr_table_type.eq."dadr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*growth_factor**2
      if(mode(i_mode)%size_distr_table_type.eq."dadlnr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*growth_factor**2
      if(mode(i_mode)%size_distr_table_type.eq."dadlogr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*growth_factor**2
      if(mode(i_mode)%size_distr_table_type.eq."dvdr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*growth_factor**3
      if(mode(i_mode)%size_distr_table_type.eq."dvdlnr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*growth_factor**3
      if(mode(i_mode)%size_distr_table_type.eq."dvdlogr") mode_aux(i_mode)%size_distr_table(:,2)=mode_aux(i_mode)%size_distr_table(:,2)*growth_factor**3

    end if

    ! *** modify refractive index
    do i_wavelength=1,num_wavelengths

      mode_aux(i_mode)%mreal(i_wavelength)=(mode_aux(i_mode)%mreal(i_wavelength)+water_mreal(i_wavelength)*(growth_factor**3-1d0))/growth_factor**3
      mode_aux(i_mode)%mimag(i_wavelength)=(mode_aux(i_mode)%mimag(i_wavelength)+water_mimag(i_wavelength)*(growth_factor**3-1d0))/growth_factor**3

    end do

    ! *** modify density
    mode(i_mode)%density=(mode(i_mode)%density+1.00d0*(growth_factor**3-1d0))/growth_factor**3

  end do

end if

end subroutine
