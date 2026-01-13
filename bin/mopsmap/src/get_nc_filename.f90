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

subroutine get_nc_filename(np,eps,mreal,mimag,with_dir,filename)

  ! effect of variable with_dir:
  ! with_dir=0: only filename
  ! with_dir=1: filename with directory within the scattering library
  ! with_dir=2: filename with complete path
  ! with_dir=3: complete path to index file

  use module_input, only: scatlib_path,dp

  integer :: np
  real(dp) :: eps,mreal,mimag
  integer :: with_dir
  character (*) :: filename

  character (len=6) :: cmreal
  character (len=8) :: cmimag
  character (len=5) :: caspect
  character (len=3) :: cnp
  character (len=15) :: shapename


  select case (np)
  case(-1)
    write(caspect,'(f5.3)') eps
    shapename='spheroid'
  case(-2)
    write(caspect,'(f5.3)') eps
    shapename='cylinder'
  case(-6)
    write(caspect,'(f5.3)') eps
    shapename='spheroid_igom'
  case(-7)
    write(caspect,'(f5.3)') eps
    shapename='spheroid_merged'
  case(2:100)
    write(caspect,'(f5.3)') eps
    write(cnp,'(I2)') np
    shapename='chebys'//cnp
    if(shapename(7:7).eq." ") shapename(7:7)="0"
  case(-109)
    shapename='shapeA'
  case(-110)
    shapename='shapeB'
  case(-111)
    shapename='shapeC'
  case(-118)
    shapename='shapeD'
  case(-126)
    shapename='shapeE'
  case(-127)
    shapename='shapeF'
  case default
    print *,'Falscher Parameter NP bei Teilchenform in subroutine get_nc_filename'
    stop
  end select

  write(cmreal,'(f6.4)') mreal
  write(cmimag,'(f8.6)') mimag
  if(caspect(1:1).eq." ") caspect(1:1)="0"
  filename=trim(shapename)//'_'//caspect//'_'//cmreal//'_'//cmimag//'.nc'
  if(np.le.-100) filename=trim(shapename)//'_'//cmreal//'_'//cmimag//'.nc'
  if(np.eq.-1.and.abs(eps-1.).lt.0.001) filename='sphere_'//cmreal//'_'//cmimag//'.nc'

! add directory
  if(with_dir.eq.1.or.with_dir.eq.2) then
    if((np.eq.-1).and.(abs(eps-1.).lt.0.001)) then
      filename=trim("spheres/"//filename)
    elseif(np.eq.-1) then
      filename=trim("spheroids/"//filename)
    elseif(np.eq.-6) then
      filename=trim("spheroids_igom/"//filename)
    elseif(np.eq.-7) then
      filename=trim("spheroids_merged/"//filename)
    elseif(np.le.-100) then
      filename=trim("irregular/"//filename)
    end if
  end if

! add complete path
  if(with_dir.eq.2) filename=trim(scatlib_path)//"/"//filename
  if(with_dir.eq.3) filename=trim(scatlib_path)//"/index.nc"

end subroutine
