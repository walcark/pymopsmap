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

module module_contributions

  use module_input, only: dp,max_size_distr_parameters
  use module_aux_var

  implicit none

  integer,parameter :: n_contributions_alloc_step=1000 ! *** each time not enough contributions are allocated, this number of contributions is allocated in addition

  integer :: current_n_contributions = 0 ! *** counter for number of used contributions

  type single_contribution
    integer :: i_mode = -1
    integer :: i_np = -1
    integer :: i_eps = -1
    integer :: i_mreal = -1
    integer :: i_mimag = -1
    real(dp) :: density = -1
    character(len=30) :: size_type = ""
    real(dp) :: n = -1
    real(dp) :: rmin = -1
    real(dp) :: rmax = -1
    real(dp) :: size_distr_parameter(max_size_distr_parameters) = -1
  end type
  type(single_contribution), allocatable :: contributions(:)
  type(single_contribution), allocatable :: contributions_tmp(:)

contains

! copy contribution i_contribution as a new contribution
subroutine copy_contribution(i_contribution)

  implicit none

  integer :: i_contribution

  current_n_contributions=current_n_contributions+1

  if(current_n_contributions.gt.size(contributions)) then
    allocate(contributions_tmp(size(contributions)+n_contributions_alloc_step))
    contributions_tmp(1:size(contributions))=contributions
    deallocate(contributions)
    call move_alloc(contributions_tmp,contributions)
  end if

  contributions(current_n_contributions)=contributions(i_contribution)

end subroutine

! make a new contribution based on parameters of mode i_mode
subroutine new_contribution(i_mode)

  use module_input, only: mode

  implicit none

  integer :: i_mode

  current_n_contributions=current_n_contributions+1

  if(current_n_contributions.gt.size(contributions)) then
    allocate(contributions_tmp(size(contributions)+n_contributions_alloc_step))
    contributions_tmp(1:size(contributions))=contributions
    deallocate(contributions)
    call move_alloc(contributions_tmp,contributions)
  end if

  contributions(current_n_contributions)%i_mode=i_mode
  contributions(current_n_contributions)%i_np=1
  contributions(current_n_contributions)%i_eps=0
  contributions(current_n_contributions)%i_mreal=0
  contributions(current_n_contributions)%i_mimag=0
  contributions(current_n_contributions)%density=mode(i_mode)%density
  contributions(current_n_contributions)%size_type=trim(mode(i_mode)%size_type)
  contributions(current_n_contributions)%n=mode(i_mode)%n
  contributions(current_n_contributions)%rmin=mode(i_mode)%rmin
  contributions(current_n_contributions)%rmax=mode(i_mode)%rmax
  contributions(current_n_contributions)%size_distr_parameter=mode(i_mode)%size_distr_parameter

end subroutine

end module
