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

subroutine calc_optics

use module_input
use module_aux_var
use module_results

implicit none

integer :: i_theta,i_wavelength

real(dp) :: angle
logical :: file_exists

! ********************************
! *** initialize scattering angles
! ********************************
if(output_theta_filename=="") then

  allocate(theta_array(output_n_theta))
  do i_theta=1,output_n_theta
    theta_array(i_theta)=180.*(i_theta-1d0)/(output_n_theta-1d0)
  end do

else

  inquire(file=output_theta_filename,exist=file_exists)
  if(.not.file_exists) then
    print'(A)',"Error: Could not find theta file "//trim(output_theta_filename)
    stop
  end if

  output_n_theta=0
  open(1,file=output_theta_filename)
  do while(.true.)
    read(1,*,end=102) angle
    output_n_theta=output_n_theta+1
  end do
  102 rewind(1)
  allocate(theta_array(output_n_theta))
  do i_theta=1,output_n_theta
    read(1,*) angle
    theta_array(i_theta)=angle
  end do
  close(1)

  ! *** for lidar output set last angle to 180°
  if((output_lidar).and.abs(theta_array(output_n_theta)-180d0).gt.1d-6) then
    print'(A)',"Error: For lidar output, the last angle in theta file "//trim(output_theta_filename)//" has to be set to 180°"
    stop
  end if

end if

! ************************************
! *** calculations for each wavelength
! ************************************
allocate(optics(num_wavelengths))
do i_wavelength=1,num_wavelengths

  call make_contributions(i_wavelength) ! *** first determine the required contributions
  call mix_contributions(i_wavelength)  ! *** then integrate the optical properties over each contribution and add up

end do


end subroutine
