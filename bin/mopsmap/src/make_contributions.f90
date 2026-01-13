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

subroutine make_contributions(i_wavelength)

use module_input
use module_dbindex
use module_aux_var
use module_contributions

implicit none

integer :: i,j,i_bin
integer :: i_wavelength,i_mode,i_np,i_eps,i_contribution,n_contributions_tmp
real(dp) :: rmin_mode,rmax_mode,rmin_dda,rmax_dda
real(dp) :: rmin_bin,rmax_bin,rmin_comp,rmax_comp,bin_fraction
logical :: two_ranges

integer :: i1,i2
real(dp) :: w1,w2
integer :: i3,i4
real(dp) :: w3,w4


if (allocated(contributions)) deallocate(contributions)
allocate(contributions(n_contributions_alloc_step))
current_n_contributions=0  ! *** counter for number of used contributions (increased by subroutines new_contribution and add_contribution)

! ****************************************************************************************************************************************
! *** loop over all user-defined modes to create shape-specific contributions and a separation between spheres/spheroids and irregular shapes
! ****************************************************************************************************************************************
do i_mode=1,num_modes

  if(mode(i_mode)%n.gt.0) then  ! *** only consider modes with number density >0

    rmin_mode=mode(i_mode)%rmin ! *** minimum and maximum size of mode
    rmax_mode=mode(i_mode)%rmax

    ! *** first the case when only irregular shapes are requested
    if(mode(i_mode)%shape_type(1:9).eq.'irregular') then

      do i_np=1,num_avail_np

        if(mode_aux(i_mode)%irregular_np_distr(i_np).gt.0) then

          call new_contribution(i_mode)

          contributions(current_n_contributions)%i_np=i_np
          contributions(current_n_contributions)%i_eps=0
          contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*mode_aux(i_mode)%irregular_np_distr(i_np) ! *** this is the number concentration

        end if

      end do

    end if

    ! *** in case the irregular shape overlay is requested
    if(mode(i_mode)%shape_irregular_overlay) then

      rmin_dda=wavelength_array(i_wavelength)/2d0/pi*mode(i_mode)%shape_irregular_overlay_xmin
      rmax_dda=wavelength_array(i_wavelength)/2d0/pi*mode(i_mode)%shape_irregular_overlay_xmax

      if(rmin_dda.lt.rmax_mode.and.rmax_dda.gt.rmin_mode) then ! *** add irregular contribution if there is an overlap between dda range and mode range

        rmin_dda=max(rmin_dda,rmin_mode)
        rmax_dda=min(rmax_dda,rmax_mode)

        do i_np=1,num_avail_np

          if(mode_aux(i_mode)%irregular_np_distr(i_np).gt.0) then

            call new_contribution(i_mode)

            contributions(current_n_contributions)%i_np=i_np
            contributions(current_n_contributions)%i_eps=0
            contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*mode_aux(i_mode)%irregular_np_distr(i_np) ! *** this is the number concentration
            contributions(current_n_contributions)%rmin=rmin_dda
            contributions(current_n_contributions)%rmax=rmax_dda

          end if

        end do

      else  ! *** no overlap

        rmin_dda=-1
        rmax_dda=-1

      end if

    else ! *** no irregular shape overlay requested

      rmin_dda=-1
      rmax_dda=-1

    end if ! *** end of irregular shape overlay case

    ! *** loop over aspect ratios
    do i_eps=1,num_avail_eps

      if(mode_aux(i_mode)%eps_distr(i_eps).gt.0) then  ! *** only make a contribution if aspect ratio occurs

        two_ranges=.false.
        if(rmin_dda.gt.0) then
          if(rmin_dda.le.rmin_mode.and.rmax_dda.ge.rmax_mode) rmin_mode=rmax_mode   ! *** nothing required since mode range is already covered by irregular shapes
          if(rmin_dda.le.rmin_mode.and.rmax_dda.le.rmax_mode.and.rmax_dda.ge.rmin_mode) rmin_mode=rmax_dda  ! *** lower part of range already covered by irregular shapes
          if(rmin_dda.ge.rmin_mode.and.rmin_dda.le.rmax_mode.and.rmax_dda.ge.rmax_mode) rmax_mode=rmin_dda  ! *** upper part of range already covered by irregular shapes
          if(rmin_dda.ge.rmin_mode.and.rmin_dda.le.rmax_mode.and.rmax_dda.ge.rmin_mode.and.rmax_dda.le.rmax_mode) two_ranges=.true. ! *** range covered by irregular shapes divides range into two separate parts
        end if

        ! *** create contribution if required
        if(rmax_mode.gt.rmin_mode) then

          call new_contribution(i_mode)

          if(abs(avail_eps(i_eps)-1d0).lt.0.001) then
            contributions(current_n_contributions)%i_np=1  ! spheres
          else
            contributions(current_n_contributions)%i_np=2  ! spheroids
          end if
          contributions(current_n_contributions)%i_eps=i_eps
          contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*mode_aux(i_mode)%eps_distr(i_eps) ! *** this is the number concentration
          contributions(current_n_contributions)%rmin=rmin_mode
          contributions(current_n_contributions)%rmax=rmax_mode
          if(two_ranges.eqv..true.) then
            contributions(current_n_contributions)%rmax=rmin_dda   ! then contributions(current_n_contributions) represents the first range
          end if

          if(two_ranges.eqv..true.) then
            call new_contribution(i_mode)                  ! this is the second range
            if(abs(avail_eps(i_eps)-1d0).lt.0.001) then
              contributions(current_n_contributions)%i_np=1  ! spheres
            else
              contributions(current_n_contributions)%i_np=2  ! spheroids
            end if
            contributions(current_n_contributions)%i_eps=i_eps
            contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*mode_aux(i_mode)%eps_distr(i_eps) ! *** this is the number concentration
            contributions(current_n_contributions)%rmin=rmax_dda
            contributions(current_n_contributions)%rmax=rmax_mode
          end if

        end if

      end if

    end do ! *** end of loop over aspect ratios

  end if

end do  ! *** end of loop over user-defined modes

! *****************************************************************************************
! *** now separate the contributions into refractive indices available in the optical data set
! *****************************************************************************************
n_contributions_tmp=current_n_contributions
do i_contribution=1,n_contributions_tmp

  if(mode_aux(contributions(i_contribution)%i_mode)%mreal(i_wavelength).gt.avail_mreal(num_avail_mreal).or.mode_aux(contributions(i_contribution)%i_mode)%mreal(i_wavelength).lt.avail_mreal(1)) then
    print'(A,F10.5,A,F10.5,A,F10.5,A,F10.5,A)',"Error: The real part of the refractive index, ",mode_aux(contributions(i_contribution)%i_mode)%mreal(i_wavelength)," at wavelength ",wavelength_array(i_wavelength),", is not within the range covered by the scattering data set (from ",avail_mreal(1)," to ",avail_mreal(num_avail_mreal),")!"
    stop
  end if

  if(mode_aux(contributions(i_contribution)%i_mode)%mimag(i_wavelength).gt.avail_mimag(num_avail_mimag).or.mode_aux(contributions(i_contribution)%i_mode)%mimag(i_wavelength).lt.avail_mimag(1)) then
    print'(A,F10.5,A,F10.5,A,F10.5,A,F10.5,A)',"Error: The imaginary part of the refractive index, ",mode_aux(contributions(i_contribution)%i_mode)%mimag(i_wavelength)," at wavelength ",wavelength_array(i_wavelength),", is not within the range covered by the scattering data set (from ",avail_mimag(1)," to ",avail_mimag(num_avail_mimag),")!"
    stop
  end if

  call interpolate_linear(num_avail_mreal,avail_mreal,mode_aux(contributions(i_contribution)%i_mode)%mreal(i_wavelength),i1,i2,w1,w2)
  call interpolate_linear(num_avail_mimag,avail_mimag,mode_aux(contributions(i_contribution)%i_mode)%mimag(i_wavelength),i3,i4,w3,w4)

  if(w1*w3.gt.0) then
    call copy_contribution(i_contribution)
    contributions(current_n_contributions)%i_mreal=i1
    contributions(current_n_contributions)%i_mimag=i3
    contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*w1*w3*(1d0-mode(contributions(i_contribution)%i_mode)%nonabs_fraction)
  end if

  if(w2*w3.gt.0) then
    call copy_contribution(i_contribution)
    contributions(current_n_contributions)%i_mreal=i2
    contributions(current_n_contributions)%i_mimag=i3
    contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*w2*w3*(1d0-mode(contributions(i_contribution)%i_mode)%nonabs_fraction)
  end if

  if(w1*w4.gt.0) then
    call copy_contribution(i_contribution)
    contributions(current_n_contributions)%i_mreal=i1
    contributions(current_n_contributions)%i_mimag=i4
    contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*w1*w4*(1d0-mode(contributions(i_contribution)%i_mode)%nonabs_fraction)
  end if

  if(w2*w4.gt.0) then
    call copy_contribution(i_contribution)
    contributions(current_n_contributions)%i_mreal=i2
    contributions(current_n_contributions)%i_mimag=i4
    contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*w2*w4*(1d0-mode(contributions(i_contribution)%i_mode)%nonabs_fraction)
  end if

  ! *** now the non-absorbing fraction
  if(mode(contributions(i_contribution)%i_mode)%nonabs_fraction.gt.0) then

    if(w1.gt.0) then
      call copy_contribution(i_contribution)
      contributions(current_n_contributions)%i_mreal=i1
      contributions(current_n_contributions)%i_mimag=1
      contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*w1*mode(contributions(i_contribution)%i_mode)%nonabs_fraction
    end if

    if(w2.gt.0) then
      call copy_contribution(i_contribution)
      contributions(current_n_contributions)%i_mreal=i2
      contributions(current_n_contributions)%i_mimag=1
      contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*w2*mode(contributions(i_contribution)%i_mode)%nonabs_fraction
    end if

  end if

  contributions(i_contribution)%n=-1 ! the original contributions gets deactivated

end do

! ***********************************************************************
! *** in case of 'bin_file' each bin needs to become a separate contribution
! ***********************************************************************
n_contributions_tmp=current_n_contributions
do i_contribution=1,n_contributions_tmp

  if(contributions(i_contribution)%size_type.eq."bin_file".and.contributions(i_contribution)%n.gt.0) then

    do i_bin=1,mode_aux(contributions(i_contribution)%i_mode)%size_distr_table_len-1

      rmin_bin=mode_aux(contributions(i_contribution)%i_mode)%size_distr_table(i_bin,1)
      rmax_bin=mode_aux(contributions(i_contribution)%i_mode)%size_distr_table(i_bin+1,1)
      rmin_comp=contributions(i_contribution)%rmin
      rmax_comp=contributions(i_contribution)%rmax

      if(rmin_bin.lt.rmax_comp.and.rmax_bin.gt.rmin_comp) then
        call copy_contribution(i_contribution)
        contributions(current_n_contributions)%rmin=max(rmin_bin,rmin_comp)
        contributions(current_n_contributions)%rmax=min(rmax_bin,rmax_comp)
        bin_fraction=(contributions(current_n_contributions)%rmax-contributions(current_n_contributions)%rmin)/(rmax_bin-rmin_bin) ! *** this is neccessary because sometimes not the complete bin is covered by the current contribution (e.g. if TMM is supplemented by IGOM)
        contributions(current_n_contributions)%n=contributions(current_n_contributions)%n*mode_aux(contributions(i_contribution)%i_mode)%size_distr_table(i_bin,2)*bin_fraction
        contributions(current_n_contributions)%size_type='bin'
      end if

    end do

    contributions(i_contribution)%n=-1 ! the original contributions gets deactivated

  end if

end do

! *** check if a contribution was created
if(current_n_contributions.lt.1) then
  print*,"Error: No contributions created."
  stop
end if

! *** write out the created contributions
if(write_debug) then
  open(120,file='tmp_contributions.txt')
  write(120,*) wavelength_array(i_wavelength)
  do i=1,current_n_contributions
    if(contributions(i)%n.gt.-100d0) then
      write(120,'(I7,I5,4F11.5)',advance='no') i,avail_np(contributions(i)%i_np),avail_eps(contributions(i)%i_eps),avail_mreal(contributions(i)%i_mreal),avail_mimag(contributions(i)%i_mimag),contributions(i)%density
      write(120,'(A1,A,F14.5,2F11.5)',advance='no') " ",trim(contributions(i)%size_type),contributions(i)%n,contributions(i)%rmin,contributions(i)%rmax
      do j=1,max_size_distr_parameters
        write(120,'(F11.5)',advance='no') contributions(i)%size_distr_parameter(j)
      end do
      if(contributions(i)%size_type.eq."distr_file") then
        write(120,'(A,A)',advance='no') "; ",trim(mode(contributions(i)%i_mode)%size_distr_table_type)
        do j=1,mode_aux(contributions(i)%i_mode)%size_distr_table_len
          write(120,'(A,2F11.5,A,A)',advance='no') " ; ",mode_aux(contributions(i)%i_mode)%size_distr_table(j,1),mode_aux(contributions(i)%i_mode)%size_distr_table(j,2)
        end do
      end if
      write(120,*)
    end if
  end do
  close(120)
end if

end subroutine
