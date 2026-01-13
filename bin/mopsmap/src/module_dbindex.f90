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

module module_dbindex

  use module_input, only: dp

  implicit none

  integer :: num_avail_mreal
  integer :: num_avail_mimag
  integer :: num_avail_np
  integer :: num_avail_eps

  real(dp),allocatable :: avail_mreal(:)
  real(dp),allocatable :: avail_mimag(:)
  integer,allocatable :: avail_np(:)
  real(dp),allocatable :: avail_eps(:)
  real(dp),allocatable :: rat(:,:)
  real(dp),allocatable :: max_sizepara(:,:,:,:)

  real(dp),allocatable :: avail_eps_width(:)

end module
