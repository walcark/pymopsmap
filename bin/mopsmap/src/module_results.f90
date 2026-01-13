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

module module_results

  use module_input, only: dp

  implicit none

  integer,parameter :: max_num_coeff=6000

  real(dp),allocatable :: theta_array(:)

  ! *** intermediate variables
  real(dp) :: csca_sum
  real(dp) :: cext_sum
  real(dp) :: n_sum
  real(dp) :: a1_sum(0:max_num_coeff)
  real(dp) :: a2_sum(0:max_num_coeff)
  real(dp) :: a3_sum(0:max_num_coeff)
  real(dp) :: a4_sum(0:max_num_coeff)
  real(dp) :: b1_sum(0:max_num_coeff)
  real(dp) :: b2_sum(0:max_num_coeff)
  real(dp) :: r2_n_sum
  real(dp) :: r3_n_sum
  real(dp) :: vol_sum
  real(dp) :: cs_sum
  real(dp) :: mass_sum
  integer :: num_coeff

  ! *** variables with the final results
  type optics_single_wavelength

    real(dp) :: cext
    real(dp) :: csca

    integer :: num_coeff
    real(dp),allocatable :: coeff(:,:)

    real(dp),allocatable :: mat(:,:)

    real(dp) :: n
    real(dp) :: reff
    real(dp) :: cs
    real(dp) :: vol
    real(dp) :: mass

  end type
  type(optics_single_wavelength), allocatable :: optics(:)

end module
