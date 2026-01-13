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

module module_aux_var

  use module_input, only: dp

  implicit none

  integer :: num_wavelengths

  real(dp), allocatable :: wavelength_array(:)

  type single_mode_aux

    real(dp), allocatable :: mreal(:)
    real(dp), allocatable :: mimag(:)

    real(dp), allocatable :: eps_distr(:)
    real(dp), allocatable :: irregular_np_distr(:)

    real(dp), allocatable :: size_distr_table(:,:)
    integer :: size_distr_table_len

  end type
  type(single_mode_aux), allocatable :: mode_aux(:)

end module
