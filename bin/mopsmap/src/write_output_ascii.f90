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

subroutine write_output_ascii

use module_input
use module_aux_var, only: num_wavelengths, wavelength_array
use module_results

implicit none

character(len=1000) :: tmp

character(len=7) :: number_format

character(len=1000), allocatable :: tmp_array(:)

integer :: n_col,i_col

integer :: i_wavelength,i_theta,i_coeff,n_coeff_tmp

real(dp) :: angst_ext,angst_sca,angst_abs,angst_back

if(output_ascii_to_file) then

  if(output_integrated) open(31,file=trim(output_ascii_filename)//".integrated")
  if(output_phase_function) open(32,file=trim(output_ascii_filename)//".phase_function")
  if(output_scattering_matrix) open(33,file=trim(output_ascii_filename)//".scattering_matrix")
  if(output_volume_scattering_function) open(34,file=trim(output_ascii_filename)//".volume_scattering_function")
  if(output_coeff) open(35,file=trim(output_ascii_filename)//".coeff")
  if(output_lidar) open(36,file=trim(output_ascii_filename)//".lidar")

end if

write(number_format,'(A2,I2.2,A1,I2.2)') 'ES',output_digits+7,'.',output_digits-1

! write the header information about the output

if(output_header) then

  if(output_integrated) then
    n_col=12
    allocate(tmp_array(n_col))
    tmp_array(1)="wavelength [µm]"
    tmp_array(2)="extinction coefficient [m^-1]"
    tmp_array(3)="single scattering albedo [1]"
    tmp_array(4)="asymmetry parameter [1]"
    if(size_equ.eq.0) tmp_array(5)="cross section equivalent effective radius [µm]"
    if(size_equ.eq.1) tmp_array(5)="volume equivalent effective radius [µm]"
    if(size_equ.eq.2) tmp_array(5)="volume cross section ratio equivalent effective radius [µm]"
    tmp_array(6)="number of particles per atmospheric volume [m^-3]"
    tmp_array(7)="particle cross section per atmospheric volume [m^-1]"
    tmp_array(8)="particle volume per atmospheric volume [1]"
    tmp_array(9)="particle mass per atmospheric volume [g m^-3]"
    tmp_array(10)="angstrom exponent of the extinction coefficient for interval from previous wavelength to current wavelength [1]"
    tmp_array(11)="angstrom exponent of the scattering coefficient for interval from previous wavelength to current wavelength [1]"
    tmp_array(12)="angstrom exponent of the absorption coefficient for interval from previous wavelength to current wavelength [1]"
    do i_col=1,n_col
      if(output_ascii_to_file) then
        write(31,'(A2,I2,A,A)') "# ",i_col,". column: ",trim(tmp_array(i_col))
      else
        write(*,'(A2,I2,A,A)') "# ",i_col+1,". column of integrated: ",trim(tmp_array(i_col))
      end if
    end do
    deallocate(tmp_array)
  end if

  if(output_phase_function) then
    n_col=3
    allocate(tmp_array(n_col))
    tmp_array(1)="wavelength [µm]"
    tmp_array(2)="scattering angle [°]"
    tmp_array(3)="normalized phase function [1]"
    do i_col=1,n_col
      if(output_ascii_to_file) then
        write(32,'(A2,I2,A,A)') "# ",i_col,". column: ",trim(tmp_array(i_col))
      else
        write(*,'(A2,I2,A,A)') "# ",i_col+1,". column of phase_function: ",trim(tmp_array(i_col))
      end if
    end do
    deallocate(tmp_array)
  end if

  if(output_scattering_matrix) then
    n_col=8
    allocate(tmp_array(n_col))
    tmp_array(1)="wavelength [µm]"
    tmp_array(2)="scattering angle [°]"
    tmp_array(3)="scattering matrix element a1 [1]"
    tmp_array(4)="scattering matrix element a2 [1]"
    tmp_array(5)="scattering matrix element a3 [1]"
    tmp_array(6)="scattering matrix element a4 [1]"
    tmp_array(7)="scattering matrix element b1 [1]"
    tmp_array(8)="scattering matrix element b2 [1]"
    do i_col=1,n_col
      if(output_ascii_to_file) then
        write(33,'(A2,I2,A,A)') "# ",i_col,". column: ",trim(tmp_array(i_col))
      else
        write(*,'(A2,I2,A,A)') "# ",i_col+1,". column of scattering_matrix: ",trim(tmp_array(i_col))
      end if
    end do
    deallocate(tmp_array)
  end if

  if(output_volume_scattering_function) then
    n_col=3
    allocate(tmp_array(n_col))
    tmp_array(1)="wavelength [µm]"
    tmp_array(2)="scattering angle [°]"
    tmp_array(3)="volume scattering function [m^-1 sr^-1]"
    do i_col=1,n_col
      if(output_ascii_to_file) then
        write(34,'(A2,I2,A,A)') "# ",i_col,". column: ",trim(tmp_array(i_col))
      else
        write(*,'(A2,I2,A,A)') "# ",i_col+1,". column of volume_scattering_function: ",trim(tmp_array(i_col))
      end if
    end do
    deallocate(tmp_array)
  end if

  if(output_coeff) then
    n_col=8
    allocate(tmp_array(n_col))
    tmp_array(1)="wavelength [µm]"
    tmp_array(2)="index l [1]"
    tmp_array(3)="expansion coefficient alpha_1^l [1]"
    tmp_array(4)="expansion coefficient alpha_2^l [1]"
    tmp_array(5)="expansion coefficient alpha_3^l [1]"
    tmp_array(6)="expansion coefficient alpha_4^l [1]"
    tmp_array(7)="expansion coefficient beta_1^l [1]"
    tmp_array(8)="expansion coefficient beta_2^l [1]"
    do i_col=1,n_col
      if(output_ascii_to_file) then
        write(35,'(A2,I2,A,A)') "# ",i_col,". column: ",trim(tmp_array(i_col))
      else
        write(*,'(A2,I2,A,A)') "# ",i_col+1,". column of coeff: ",trim(tmp_array(i_col))
      end if
    end do
    deallocate(tmp_array)
  end if

  if(output_lidar) then
    n_col=9
    allocate(tmp_array(n_col))
    tmp_array(1)="wavelength [µm]"
    tmp_array(2)="extinction coefficient [m^-1]"
    tmp_array(3)="backscatter coefficient [m^-1 sr^-1]"
    tmp_array(4)="lidar ratio [sr]"
    tmp_array(5)="linear depolarization ratio [1]"
    tmp_array(6)="angstrom exponent of the extinction coefficient for interval from previous wavelength to current wavelength [1]"
    tmp_array(7)="angstrom exponent of the backscatter coefficient for interval from previous wavelength to current wavelength [1]"
    tmp_array(8)="extinction to mass conversion factor [g m^-2]"
    tmp_array(9)="mass to backscatter conversion factor [m^2 g^-1 sr^-1]"
    do i_col=1,n_col
      if(output_ascii_to_file) then
        write(36,'(A2,I2,A,A)') "# ",i_col,". column: ",trim(tmp_array(i_col))
      else
        write(*,'(A2,I2,A,A)') "# ",i_col+1,". column of lidar: ",trim(tmp_array(i_col))
      end if
    end do
    deallocate(tmp_array)
  end if

end if

! loop for writing the results
do i_wavelength=1,num_wavelengths

  ! first calculate Angstrom exponents
  if(i_wavelength.eq.1) then

    angst_ext=sqrt(-1d0*num_wavelengths)
    angst_sca=sqrt(-1d0*num_wavelengths)
    angst_abs=sqrt(-1d0*num_wavelengths)
    angst_back=sqrt(-1d0*num_wavelengths)

  else

    angst_ext=-log(optics(i_wavelength)%cext/optics(i_wavelength-1)%cext)/log(wavelength_array(i_wavelength)/wavelength_array(i_wavelength-1))
    angst_sca=-log(optics(i_wavelength)%csca/optics(i_wavelength-1)%csca)/log(wavelength_array(i_wavelength)/wavelength_array(i_wavelength-1))
    angst_abs=-log((optics(i_wavelength)%cext-optics(i_wavelength)%csca)/(optics(i_wavelength-1)%cext-optics(i_wavelength-1)%csca))/log(wavelength_array(i_wavelength)/wavelength_array(i_wavelength-1))
    angst_back=-log((optics(i_wavelength)%csca*optics(i_wavelength)%mat(1,output_n_theta))/(optics(i_wavelength-1)%csca*optics(i_wavelength-1)%mat(1,output_n_theta)))/log(wavelength_array(i_wavelength)/wavelength_array(i_wavelength-1))

  end if

  if(output_integrated) then

    write(tmp,'(12'//number_format//')') wavelength_array(i_wavelength),&
      optics(i_wavelength)%cext*1d-12,&
      optics(i_wavelength)%csca/optics(i_wavelength)%cext,&
      optics(i_wavelength)%coeff(1,1)/3d0,&
      optics(i_wavelength)%reff,&
      optics(i_wavelength)%n,&
      optics(i_wavelength)%cs*1d-12,&
      optics(i_wavelength)%vol*1d-18,&
      optics(i_wavelength)%mass*1d-12,&
      angst_ext,&
      angst_sca,&
      angst_abs

    if(output_ascii_to_file) then
      write(31,*) trim(tmp)
    else
      write(*,'(A)') 'integrated '//trim(tmp)
    end if

  end if

  if(output_phase_function) then

    do i_theta=1,output_n_theta

      write(tmp,'(8'//number_format//')') wavelength_array(i_wavelength),&
        theta_array(i_theta),&
        optics(i_wavelength)%mat(1,i_theta)

      if(output_ascii_to_file) then
        write(32,*) trim(tmp)
      else
        write(*,'(A)') 'phase_function '//trim(tmp)
      end if

    end do

  end if

  if(output_scattering_matrix) then

    do i_theta=1,output_n_theta

      write(tmp,'(8'//number_format//')') wavelength_array(i_wavelength),&
        theta_array(i_theta),&
        optics(i_wavelength)%mat(:,i_theta)

      if(output_ascii_to_file) then
        write(33,*) trim(tmp)
      else
        write(*,'(A)') 'scattering_matrix '//trim(tmp)
      end if

    end do

  end if

  if(output_volume_scattering_function) then

    do i_theta=1,output_n_theta

      write(tmp,'(3'//number_format//')') wavelength_array(i_wavelength),&
        theta_array(i_theta),&
        optics(i_wavelength)%mat(1,i_theta)*optics(i_wavelength)%csca/(4d0*pi)*1d-12

      if(output_ascii_to_file) then
        write(34,*) trim(tmp)
      else
        write(*,'(A)') 'volume_scattering_function '//trim(tmp)
      end if

    end do

  end if

  if(output_coeff) then

    if(output_n_coeff.ge.0) then

      n_coeff_tmp=min(output_n_coeff,optics(i_wavelength)%num_coeff)

    else

      n_coeff_tmp=optics(i_wavelength)%num_coeff

    end if

    do i_coeff=0,n_coeff_tmp-1

      write(tmp,'('//number_format//',I5,6'//number_format//')') wavelength_array(i_wavelength),&
        i_coeff,&
        optics(i_wavelength)%coeff(:,i_coeff)

      if(output_ascii_to_file) then
        write(35,*) trim(tmp)
      else
        write(*,'(A)') 'coeff '//trim(tmp)
      end if

    end do

  end if

  if(output_lidar) then

    write(tmp,'(9'//number_format//')') wavelength_array(i_wavelength),&
      optics(i_wavelength)%cext*1d-12,&
      optics(i_wavelength)%csca*optics(i_wavelength)%mat(1,output_n_theta)/(4d0*pi)*1d-12,&
      4d0*pi*optics(i_wavelength)%cext/(optics(i_wavelength)%csca*optics(i_wavelength)%mat(1,output_n_theta)),&
      (1d0-optics(i_wavelength)%mat(2,output_n_theta)/optics(i_wavelength)%mat(1,output_n_theta))/(1d0+optics(i_wavelength)%mat(2,output_n_theta)/optics(i_wavelength)%mat(1,output_n_theta)),&
      angst_ext,&
      angst_back,&
      optics(i_wavelength)%mass/optics(i_wavelength)%cext,&
      optics(i_wavelength)%csca*optics(i_wavelength)%mat(1,output_n_theta)/(4d0*pi)/optics(i_wavelength)%mass

    if(output_ascii_to_file) then
      write(36,*) trim(tmp)
    else
      write(*,'(A)') 'lidar '//trim(tmp)
    end if

  end if

end do

end subroutine
