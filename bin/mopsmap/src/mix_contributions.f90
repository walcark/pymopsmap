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

subroutine mix_contributions(i_wavelength)

use module_input, only: output_n_theta,write_status
use module_results
use module_contributions

implicit none

integer :: i_wavelength,i_contribution
real :: time

! *** initialize summation variables
csca_sum=0;cext_sum=0;n_sum=0
a1_sum=0;a2_sum=0;a3_sum=0;a4_sum=0;b1_sum=0;b2_sum=0
r2_n_sum=0;r3_n_sum=0;vol_sum=0;cs_sum=0;mass_sum=0
num_coeff=0

! *** summation of the contributions
do i_contribution=1,current_n_contributions

  if(contributions(i_contribution)%n.gt.0) call add_contribution(i_wavelength,i_contribution)

  if(write_status) then
    call cpu_time(time)
    open(120,file='tmp_status.txt',status='replace')
    if(i_wavelength.lt.num_wavelengths.or.i_contribution.lt.current_n_contributions) then
      write(120,'(A,F10.4)')  'Wavelength',wavelength_array(i_wavelength)
      write(120,'(A,I7,A,I7)')'Wavelength', i_wavelength, ' out of ',num_wavelengths
      write(120,'(A,I7,A,I7)')'Component ',i_contribution, ' out of ',current_n_contributions
    end if
    write(120,'(A,F10.2,A)') 'CPU time',time,' seconds'
    close(120)
  end if

end do

! *** store optical properties at this wavelength
optics(i_wavelength)%cext=cext_sum
optics(i_wavelength)%csca=csca_sum

optics(i_wavelength)%num_coeff=num_coeff

allocate(optics(i_wavelength)%coeff(6,0:num_coeff))
optics(i_wavelength)%coeff(1,:)=a1_sum(0:num_coeff)/csca_sum
optics(i_wavelength)%coeff(2,:)=a2_sum(0:num_coeff)/csca_sum
optics(i_wavelength)%coeff(3,:)=a3_sum(0:num_coeff)/csca_sum
optics(i_wavelength)%coeff(4,:)=a4_sum(0:num_coeff)/csca_sum
optics(i_wavelength)%coeff(5,:)=b1_sum(0:num_coeff)/csca_sum
optics(i_wavelength)%coeff(6,:)=b2_sum(0:num_coeff)/csca_sum

optics(i_wavelength)%n=n_sum
optics(i_wavelength)%reff=r3_n_sum/r2_n_sum

optics(i_wavelength)%cs=cs_sum
optics(i_wavelength)%vol=vol_sum
optics(i_wavelength)%mass=mass_sum

! *** calculate scattering matrix from expansion coefficients
allocate(optics(i_wavelength)%mat(6,output_n_theta))
call matr(optics(i_wavelength)%coeff,num_coeff,output_n_theta,theta_array,optics(i_wavelength)%mat)

end subroutine
