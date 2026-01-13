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

real(dp) function log_distr(radius,r_mod,sigma)

! ****************************************************
! *** return dn(r)/dr of a log-normal distribution ***
! ****************************************************

  use module_input, only: dp

  real(dp) :: radius,r_mod,sigma

  log_distr=1d0/(2.506628275d0*(radius*log(sigma)))*exp(-((log(radius)-log(r_mod))**2)/(2d0*(log(sigma))**2))

end function
