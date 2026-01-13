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

real(dp) function log_distr_ar(ar,ar0,sigma)

  use module_input, only: dp

  real(dp) :: ar,ar0,sigma

  log_distr_ar=0.39894228d0/((ar-1d0)*sigma)*exp(-((log(ar-1d0)-log(ar0-1d0))**2)/(2d0*sigma**2))

  ! *** 0.39894228 is 1/sqrt(2*pi)

end function
