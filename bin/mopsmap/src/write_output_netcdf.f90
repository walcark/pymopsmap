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

subroutine write_output_netcdf

use netcdf
use module_input
use module_aux_var
use module_results

implicit none

character(len=1000) :: text

real,allocatable :: theta(:,:,:,:)
real,allocatable :: phase(:,:,:,:)
integer,allocatable :: ntheta(:,:,:)
real :: ext(1,num_wavelengths)
real :: ssa(1,num_wavelengths)

real,allocatable :: pmom(:,:,:,:)

integer :: i
integer :: status,ncid
integer :: wavelength_dimid,phamat_dimid,hum_dimid,theta_dimid,mom_dimid
integer :: wavelength_varid,hum_varid,theta_id,phase_id,ntheta_id,ext_id,ssa_id,rho_id,pmom_id,nmom_id
integer :: refre_id,refim_id

integer :: i_wavelength,i_matrelem
integer :: num_matrelem


if(.not.output_netcdf) return

num_matrelem=6

if(output_n_coeff.le.0) then

  output_n_coeff=129

end if

allocate(theta(output_n_theta,num_matrelem,1,num_wavelengths))
allocate(phase(output_n_theta,num_matrelem,1,num_wavelengths))
allocate(ntheta(num_matrelem,1,num_wavelengths))
allocate(pmom(output_n_coeff,num_matrelem,1,num_wavelengths))

! ****************************
! *** prepare output variables
! ****************************
do i_wavelength=1,num_wavelengths

  do i_matrelem=1,num_matrelem

    theta(:,i_matrelem,1,i_wavelength)=real(theta_array(output_n_theta:1:-1))
    ntheta(i_matrelem,1,i_wavelength)=output_n_theta
    select case(i_matrelem)
      case(1)
        phase(:,i_matrelem,1,i_wavelength)=real(optics(i_wavelength)%mat(1,output_n_theta:1:-1))
      case(2)
        phase(:,i_matrelem,1,i_wavelength)=real(optics(i_wavelength)%mat(5,output_n_theta:1:-1))
      case(3)
        phase(:,i_matrelem,1,i_wavelength)=real(optics(i_wavelength)%mat(3,output_n_theta:1:-1))
      case(4)
        phase(:,i_matrelem,1,i_wavelength)=-real(optics(i_wavelength)%mat(6,output_n_theta:1:-1))
      case(5)
        phase(:,i_matrelem,1,i_wavelength)=real(optics(i_wavelength)%mat(2,output_n_theta:1:-1))
      case(6)
        phase(:,i_matrelem,1,i_wavelength)=real(optics(i_wavelength)%mat(4,output_n_theta:1:-1))
    end select

    ! not sure if this is the correct order for libRadtran
    pmom(:,i_matrelem,1,i_wavelength)=real(optics(i_wavelength)%coeff(i_matrelem,0:(output_n_coeff-1)))

  end do

  if (pmom(1,1,1,i_wavelength).ne.1) then
    if (abs(pmom(1,1,1,i_wavelength)-1).gt.0.001) then
      print*, 'Warning: First legendre moment of the phase function not equal to 1 for i_wavelength ',i_wavelength
    end if
    pmom(1,1,1,i_wavelength)=1
  end if

  ext(1,i_wavelength)=real(1000d0*optics(i_wavelength)%cext/(optics(i_wavelength)%mass))
  ssa(1,i_wavelength)=real(optics(i_wavelength)%csca/optics(i_wavelength)%cext)

end do

! ******************
! *** write the file
! ******************
status=nf90_create(output_netcdf_filename, nf90_clobber, ncid)

status=nf90_def_dim(ncid,'nlam',num_wavelengths,wavelength_dimid)
status=nf90_def_dim(ncid,'nmommax',output_n_coeff,mom_dimid)
status=nf90_def_dim(ncid,'nphamat',num_matrelem,phamat_dimid)
if(output_netcdf_reff) then
  status=nf90_def_dim(ncid,'nreff',1,hum_dimid)
else
  status=nf90_def_dim(ncid,'nhum',1,hum_dimid)
end if
status=nf90_def_dim(ncid,'nthetamax',output_n_theta,theta_dimid)

status=nf90_def_var(ncid,'wavelen',nf90_double,(/wavelength_dimid/),wavelength_varid)
if(output_netcdf_reff) then
  status=nf90_def_var(ncid,'reff',nf90_double,(/hum_dimid/),hum_varid)
else
  status=nf90_def_var(ncid,'hum',nf90_double,(/hum_dimid/),hum_varid)
end if
status=nf90_def_var(ncid,'theta',nf90_float,(/theta_dimid,phamat_dimid,hum_dimid,wavelength_dimid/),theta_id)
status=nf90_def_var(ncid,'phase',nf90_float,(/theta_dimid,phamat_dimid,hum_dimid,wavelength_dimid/),phase_id)
status=nf90_def_var(ncid,'pmom',nf90_float,(/mom_dimid,phamat_dimid,hum_dimid,wavelength_dimid/),pmom_id)
status=nf90_def_var(ncid,'nmom',nf90_int,(/phamat_dimid,hum_dimid,wavelength_dimid/),nmom_id)
status=nf90_def_var(ncid,'ntheta',nf90_int,(/phamat_dimid,hum_dimid,wavelength_dimid/),ntheta_id)
status=nf90_def_var(ncid,'ext',nf90_double,(/hum_dimid,wavelength_dimid/),ext_id)
status=nf90_def_var(ncid,'ssa',nf90_double,(/hum_dimid,wavelength_dimid/),ssa_id)
status=nf90_def_var(ncid,'rho',nf90_double,(/hum_dimid/),rho_id)

if(num_modes.eq.1) then
  status=nf90_def_var(ncid,'refre',nf90_double,(/wavelength_dimid/),refre_id)
  status=nf90_def_var(ncid,'refim',nf90_double,(/wavelength_dimid/),refim_id)
  if(mode(1)%size_type.eq."mono") then
    write(text,'(A,F8.3)') "radius ",(mode(1)%rmin+mode(1)%rmax)*0.5
    status=nf90_put_att(ncid,nf90_global,"size",trim(text))
  elseif(mode(1)%size_type.eq."log_normal") then
    write(text,'(A,F7.3,A,F7.3,A,F7.3,A,F7.3)') "log-normal: r_mode ",mode(1)%size_distr_parameter(1), "; sigma ",mode(1)%size_distr_parameter(2), "; r_min ",mode(1)%rmin, "; r_max ",mode(1)%rmax
    status=nf90_put_att(ncid,nf90_global,"size_distr",trim(text))
  elseif(mode(1)%size_type.eq."mod_gamma") then
    write(text,'(A,F7.3,A,F7.3,A,F7.3,A,F7.3,A,F7.3)') "mod_gamma: alpha ",mode(1)%size_distr_parameter(1), "; B ",mode(1)%size_distr_parameter(2), "; gamma ",mode(1)%size_distr_parameter(3), "; r_min ",mode(1)%rmin, "; r_max ",mode(1)%rmax
    status=nf90_put_att(ncid,nf90_global,"size_distr",trim(text))
  end if
  if(mode(1)%shape_type.eq."sphere") then
    write(text,'(A)') "sphere"
    status=nf90_put_att(ncid,nf90_global,"shape",trim(text))
  elseif(mode(1)%shape_type.eq."spheroid") then
    if(mode(1)%shape_parameter(1).eq.0) then
      write(text,'(A,F7.3)') "oblate spheroid with aspect ratio", mode(1)%shape_parameter(2)
    else
      write(text,'(A,F7.3)') "prolate spheroid with aspect ratio", mode(1)%shape_parameter(2)
    end if
    status=nf90_put_att(ncid,nf90_global,"shape",trim(text))
  elseif(mode(1)%shape_type.eq."spheroid_log_normal") then
    write(text,'(A,F7.3,A,F7.3,A,F7.3,A,F7.3)') "spheroids with log-normal aspect ratio distribution: ratio of spheroids ", mode(1)%shape_parameter(1), "; ratio of prolate spheroids (of all spheroids) ", mode(1)%shape_parameter(2), "; epsilon'_0 ", mode(1)%shape_parameter(3), "; sigma ", mode(1)%shape_parameter(4)
    status=nf90_put_att(ncid,nf90_global,"shape_distr",trim(text))
  end if
end if

text="Netcdf file created using mopsmap"
status=nf90_put_att(ncid,nf90_global,"file_info",trim(text))
status=nf90_put_att(ncid,nf90_global,"version",20090626)
text="Mie, T-matrix method, geometric optics"
status=nf90_put_att(ncid,nf90_global,"parameterization",trim(text))

status=nf90_enddef(ncid)

status=nf90_put_var(ncid,wavelength_varid,wavelength_array)
if(output_netcdf_reff) then
  status=nf90_put_var(ncid,hum_varid,(/dble(optics(1)%reff)/))
else
  status=nf90_put_var(ncid,hum_varid,(/0d0/))
end if
status=nf90_put_var(ncid,theta_id,theta)
status=nf90_put_var(ncid,phase_id,phase)
status=nf90_put_var(ncid,pmom_id,real(pmom))
status=nf90_put_var(ncid,nmom_id,reshape((/(output_n_coeff,i=1,num_matrelem*num_wavelengths)/), shape=(/num_matrelem,1,num_wavelengths/)))
status=nf90_put_var(ncid,ntheta_id,reshape((/(output_n_theta,i=1,num_matrelem*num_wavelengths)/), shape=(/num_matrelem,1,num_wavelengths/)))
status=nf90_put_var(ncid,ext_id,dble(ext))
status=nf90_put_var(ncid,ssa_id,dble(ssa))
status=nf90_put_var(ncid,rho_id,(/dble(optics(1)%mass/optics(1)%vol)/))
if(num_modes.eq.1) then
  status=nf90_put_var(ncid,refre_id,mode_aux(1)%mreal(1:num_wavelengths))
  status=nf90_put_var(ncid,refim_id,mode_aux(1)%mimag(1:num_wavelengths))
end if
status=nf90_close(ncid)

end subroutine
