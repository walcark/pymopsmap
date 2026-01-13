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

module module_input

  implicit none

  integer,parameter :: dp=selected_real_kind(15) ! double precision float kind
  real(dp),parameter :: pi=3.14159265358979d0

  integer,parameter :: max_size_distr_parameters=3
  integer,parameter :: max_shape_parameters=4

  integer :: num_modes = -1
  integer :: size_equ = 0
  logical :: diameter = .false.

  real(dp) :: relative_humidity = -1

  character(len=1000) :: water_refr_file = "../data/refr_water_segelstein"

  character(len=1000) :: scatlib_path = "optical_dataset"

! microphysical parameters
  type single_mode

    character(len=:), allocatable :: size_type
    real(dp) :: n = -1
    real(dp) :: rmin = -1
    real(dp) :: rmax = -1
    real(dp) :: size_distr_parameter(max_size_distr_parameters) = -1
    character(len=:), allocatable :: size_distr_file
    character(len=:), allocatable :: size_distr_table_type
    character(len=:), allocatable :: size_distr_list

    character(len=:), allocatable :: shape_type
    real(dp) :: shape_parameter(max_shape_parameters) = -1
    character(len=:), allocatable :: shape_distr_file

    logical :: shape_irregular_overlay = .false.
    real(dp) :: shape_irregular_overlay_xmin = -1
    real(dp) :: shape_irregular_overlay_xmax = -1
    character(len=:), allocatable :: shape_irregular_overlay_file

    real(dp) :: mreal = -1
    real(dp) :: mimag = -1
    character(len=:), allocatable :: refr_file
    real(dp) :: nonabs_fraction = 0

    real(dp) :: kappa = -1

    real(dp) :: density = 2.6

  end type
  type(single_mode), allocatable :: mode(:)

! wavelengths
  real(dp) :: wavelength = -1
  real(dp) :: wavelength_min = -1
  real(dp) :: wavelength_max = -1
  real(dp) :: wavelength_step = -1
  character(len=:), allocatable :: wavelength_file
  character(len=:), allocatable :: wavelength_list
  logical :: use_wavelength_range = .false.
  logical :: use_wavelength_file = .false.
  logical :: use_wavelength_from_refr_file = .false.
  logical :: use_wavelength_list = .false.

! output options
  logical :: output_netcdf = .false.
  logical :: output_netcdf_reff = .false.
  character(len=:), allocatable :: output_netcdf_filename
  logical :: output_ascii_to_file = .false.
  character(len=:), allocatable :: output_ascii_filename

  logical :: output_integrated = .false.
  logical :: output_phase_function = .false.
  logical :: output_scattering_matrix = .false.
  logical :: output_volume_scattering_function = .false.
  logical :: output_coeff = .false.
  logical :: output_lidar = .false.
  integer :: output_digits = 6
  integer :: output_n_theta = 181
  character(len=:), allocatable :: output_theta_filename
  integer :: output_n_coeff = -1
  logical :: output_header = .false.

  logical :: write_debug = .false.
  logical :: write_status = .false.

end module
