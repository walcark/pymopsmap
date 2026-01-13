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

subroutine read_input

use module_input

implicit none

character(len=1000) :: filename,line,option,arg_str,arg_str2,arg_str3,err_str

real(dp) :: arg_real(6)

integer :: arg_int

integer :: i_mode,i_mode_min,i_mode_max

integer :: i_line

logical :: file_exists

integer :: i


! Read name of input file which is the first command line argument
i=command_argument_count()
if(i.gt.0) then
  call get_command_argument(1,filename)
else
  filename="input_mopsmap.txt"
end if

! check if input file exists
inquire(file=filename,exist=file_exists)
if(.not.file_exists) then
  print'(A)',"Error: Could not find input file "//trim(filename)
  stop
end if

! open and parse the input file
open(11,file=filename,status='old',err=1112)

! first determine how many modes are in the input file
num_modes=1
do while(.true.)
  read(11,'(A)',end=1113) line
  line=adjustl(line)
  if(line(1:4).eq."mode") then
    read(line,*,err=1110,end=1110) option,i_mode
    if(i_mode.gt.num_modes) num_modes=i_mode
  end if
end do

1113 rewind(11)

! allocate the required number of modes
allocate(mode(num_modes))

i_line=0

do while(.true.)

  ! count the lines
  i_line=i_line+1

  ! each line in input file is read into variable line which is processed below
  read(11,'(A)',end=1111) line

  ! prepare part of the error message (in case it is needed later)
  write(err_str,'(A,I6,A)') "Error parsing "//trim(filename)//" on line ",i_line," '"//trim(line)//"':"

  ! test if input line is longer than what can be stored in variable line (assuming that this is the case if only few spaces are at the end of the fixed-length variable 'line')
  if(len(line)-len(trim(line)).lt.10) then
    print'(A,I6,A)',trim(err_str)//" Line too long! It should be less than ",len(line)-10," characters long."
    stop
  end if

  ! remove spaces at the beginning of the line
  line=adjustl(line)

  ! if line is empty or starts with '#' continue with next line
  if(line.eq."") cycle
  if(line(1:1).eq."#") cycle

  ! check if line is for a specific mode
  ! if yes, save the mode number in i_mode_min and i_mode_max and remove first part from variable line
  if(line(1:4).eq."mode") then

    read(line,*,err=1110,end=1110) option,i_mode_min
    i_mode_max=i_mode_min

    line=line(5:) ! remove "mode"
    line=adjustl(line) ! remove leading spaces
    i=index(line,' ') ! search for first space (after the mode number)
    line=line(i:)  ! remove the mode number
    line=adjustl(line) ! remove leading spaces --> after that line should be without mode and mode number

  else

    i_mode_min=1
    i_mode_max=num_modes

  end if

  ! read the first 'word' from variable line
  ! this 'word' specifies option described in the user manual
  read(line,*,err=1110,end=1110) option

  select case(option)

    case('size_equ')

      read(line,*,err=1110,end=1110) option,arg_str

      select case(arg_str)

        case('cs')

          size_equ=0

        case('vol')

          size_equ=1

        case('vol_cs_ratio')

          size_equ=2

        case default

          print'(A)',trim(err_str)//" Unknown argument "//trim(arg_str)//" for size_equ"
          stop

      end select

    case('diameter')

      diameter=.True.

    case('size')

      read(line,*,err=1110,end=1110) option,arg_str

      if(arg_str.eq."log_normal") then

        read(line,*,err=1110,end=1110) option,arg_str,arg_real(1:5)
        if(arg_real(1).le.0) then
          print'(A)',trim(err_str)//" r_mod <= 0, i.e. r_mod is not in valid range."
          stop
        end if
        if(arg_real(2).le.1) then
          print'(A)',trim(err_str)//" sigma <= 1, i.e. sigma not in valid range."
          stop
        end if
        if(arg_real(3).lt.0) then
          print'(A)',trim(err_str)//" n < 0, i.e. n not in valid range."
          stop
        end if
        if(arg_real(4).gt.arg_real(5)) then
          print'(A)',trim(err_str)//" r_min > r_max"
          stop
        end if
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%size_type="log_normal"
          mode(i_mode)%n=arg_real(3)
          mode(i_mode)%rmin=arg_real(4)
          mode(i_mode)%rmax=arg_real(5)
          mode(i_mode)%size_distr_parameter(1:2)=arg_real(1:2)
        end do

      elseif(arg_str.eq."mod_gamma") then

        read(line,*,err=1110,end=1110) option,arg_str,arg_real(1:6)
        if(arg_real(1).lt.0) then
          print'(A)',trim(err_str)//" n < 0, i.e. n not in valid range."
          stop
        end if
        if(arg_real(2).gt.arg_real(3)) then
          print'(A)',trim(err_str)//" r_min > r_max"
          stop
        end if
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%size_type="mod_gamma"
          mode(i_mode)%n=arg_real(1)
          mode(i_mode)%rmin=arg_real(2)
          mode(i_mode)%rmax=arg_real(3)
          mode(i_mode)%size_distr_parameter(1:3)=arg_real(4:6)
        end do

      elseif(arg_str.eq."distr_file") then

        read(line,*,err=1110,end=1110) option,arg_str,arg_str2,arg_str3
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%size_type="distr_file"
          mode(i_mode)%size_distr_table_type=arg_str2
          mode(i_mode)%size_distr_file=arg_str3
          mode(i_mode)%n=1d0
        end do

      elseif(arg_str.eq."bin_file") then

        read(line,*,err=1110,end=1110) option,arg_str,arg_str2
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%size_type="bin_file"
          mode(i_mode)%size_distr_file=arg_str2
          mode(i_mode)%n=1d0
        end do

      elseif(arg_str.eq."distr_list") then

        read(line,*,err=1110,end=1110) option,arg_str,arg_str2
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%size_type="distr_list"
          mode(i_mode)%size_distr_table_type=arg_str2
          mode(i_mode)%n=1d0
          i=index(line,trim(arg_str2))+len(trim(arg_str2))
          mode(i_mode)%size_distr_list=adjustl(trim(line(i:)))
        end do

      else

        read(line,*,err=1110,end=1110) option,arg_real(1:2)
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%size_type="mono"
          mode(i_mode)%n=arg_real(2)
          mode(i_mode)%rmin=arg_real(1)*(1d0-1d-12) ! a very narrow size interval is used
          mode(i_mode)%rmax=arg_real(1)*(1d0+1d-12)
        end do

     end if

    case('shape')

      read(line,*,err=1110,end=1110) option,arg_str

      select case(arg_str)

        case('sphere','spheres')

        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%shape_type="sphere"
        end do

        case('spheroid','spheroids')

          read(line,*,err=1110,end=1110) option,arg_str,arg_str2

          if(arg_str2.eq."log_normal") then

            read(line,*,err=1110,end=1110) option,arg_str,arg_str2,arg_real(1:4)
            if(arg_real(1).lt.0.or.arg_real(1).gt.1) then
              print'(A)',trim(err_str)//" zeta_1 < 0 or zeta_1 > 1, i.e. zeta_1 is not in valid range."
              stop
            end if
            if(arg_real(2).lt.0.or.arg_real(2).gt.1) then
              print'(A)',trim(err_str)//" zeta_2 < 0 or zeta_2 > 1, i.e. zeta_2 is not in valid range."
              stop
            end if
            if(arg_real(3).le.1) then
              print'(A)',trim(err_str)//" epsilon_0 <= 1, i.e. epsilon_0 is not in valid range."
              stop
            end if
            if(arg_real(4).le.0) then
              print'(A)',trim(err_str)//" sigma_ar <= 0, i.e. sigma_ar is not in valid range."
              stop
            end if
            do i_mode=i_mode_min,i_mode_max
              mode(i_mode)%shape_type="spheroid_log_normal"
              mode(i_mode)%shape_parameter(1:4)=arg_real(1:4)
            end do

          elseif(arg_str2.eq."distr_file") then

            read(line,*,err=1110,end=1110) option,arg_str,arg_str2,arg_str3
            do i_mode=i_mode_min,i_mode_max
              mode(i_mode)%shape_type="spheroid_distr_file"
              mode(i_mode)%shape_distr_file=arg_str3
            end do

          else

            read(line,*,err=1110,end=1110) option,arg_str,arg_str2,arg_real(1)
            if(arg_str2.ne."oblate".and.arg_str2.ne."prolate") then
              print'(A)',trim(err_str)//" Unexpected shape argument "//trim(arg_str2)
              stop
            end if
            do i_mode=i_mode_min,i_mode_max
              mode(i_mode)%shape_type="spheroid"
              if(arg_str2.eq."oblate") mode(i_mode)%shape_parameter(1)=0.
              if(arg_str2.eq."prolate") mode(i_mode)%shape_parameter(1)=1.
              mode(i_mode)%shape_parameter(2)=arg_real(1)
              if(mode(i_mode)%shape_parameter(2).lt.1d0) then
                print'(A)',trim(err_str)//" An aspect ratio smaller than 1 was given but the aspect ratios needs to be at least 1."
                stop
              end if
              if(mode(i_mode)%shape_parameter(2).gt.5d0) then
                print'(A)',trim(err_str)//" Aspect ratio greater than 5 but only aspect ratios <=5 are covered by the data set."
                stop
              end if
            end do

          end if

        case('irregular')

          read(line,*,err=1110,end=1110) option,arg_str,arg_str2

          if(arg_str2.eq."distr_file") then

            read(line,*,err=1110,end=1110) option,arg_str,arg_str2,arg_str3
            do i_mode=i_mode_min,i_mode_max
              mode(i_mode)%shape_type="irregular_distr_file"
              mode(i_mode)%shape_distr_file=arg_str3
            end do

          else

            read(line,*,err=1110,end=1110) option,arg_str,arg_str2
            do i_mode=i_mode_min,i_mode_max
              mode(i_mode)%shape_type="irregular"
              if(arg_str2.eq."A".or.arg_str2.eq."-109") mode(i_mode)%shape_parameter(1)=-109
              if(arg_str2.eq."B".or.arg_str2.eq."-110") mode(i_mode)%shape_parameter(1)=-110
              if(arg_str2.eq."C".or.arg_str2.eq."-111") mode(i_mode)%shape_parameter(1)=-111
              if(arg_str2.eq."D".or.arg_str2.eq."-118") mode(i_mode)%shape_parameter(1)=-118
              if(arg_str2.eq."E".or.arg_str2.eq."-126") mode(i_mode)%shape_parameter(1)=-126
              if(arg_str2.eq."F".or.arg_str2.eq."-127") mode(i_mode)%shape_parameter(1)=-127

              if(mode(i_mode)%shape_parameter(1).eq.-1) then
                print'(A)',trim(err_str)//" Unknown irregular shape "//trim(arg_str2)
                stop
              end if

            end do

          end if

        case('irregular_overlay')

          read(line,*,err=1110,end=1110) option,arg_str,arg_str2,arg_real(1:2)
          do i_mode=i_mode_min,i_mode_max
            mode(i_mode)%shape_irregular_overlay=.true.
            mode(i_mode)%shape_irregular_overlay_file=arg_str2
            mode(i_mode)%shape_irregular_overlay_xmin=arg_real(1)
            mode(i_mode)%shape_irregular_overlay_xmax=arg_real(2)
          end do

        case default

          print'(A)',trim(err_str)//" Unknown shape "//trim(arg_str)
          stop

      end select

    case('refr','refrac')

      read(line,*,err=1110,end=1110) option,arg_str
      if(arg_str.eq."file") then
        read(line,*,err=1110,end=1110) option,arg_str,arg_str2
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%refr_file=arg_str2
        end do
      elseif(arg_str.eq."nonabs_fraction") then
        read(line,*,err=1110,end=1110) option,arg_str,arg_real(1)
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%nonabs_fraction=arg_real(1)
        end do
      else
        read(line,*,err=1110,end=1110) option,arg_real(1:2)
        do i_mode=i_mode_min,i_mode_max
          mode(i_mode)%mreal=arg_real(1)
          mode(i_mode)%mimag=arg_real(2)
        end do
      end if

    case('density')

      read(line,*,err=1110,end=1110) option,arg_real(1)
      do i_mode=i_mode_min,i_mode_max
        mode(i_mode)%density=arg_real(1)
      end do

    case('kappa')

      read(line,*,err=1110,end=1110) option,arg_real(1)
      do i_mode=i_mode_min,i_mode_max
        mode(i_mode)%kappa=arg_real(1)
      end do

    case('rH','rh')

      read(line,*,err=1110,end=1110) option,arg_real(1)
      relative_humidity=arg_real(1)

    case('water_refrac_file')

      read(line,*,err=1110) option,arg_str
      water_refr_file=arg_str

    case('wavelength')

      read(line,*,err=1110,end=1110) option,arg_str
      if(arg_str.eq."range") then
        read(line,*,err=1110,end=1110) option,arg_str,arg_real(1:3)
        wavelength_min=arg_real(1)
        wavelength_max=arg_real(2)
        wavelength_step=arg_real(3)
        use_wavelength_range=.true.
        use_wavelength_file=.false.
        use_wavelength_from_refr_file=.false.
        use_wavelength_list=.false.
      elseif(arg_str.eq."file") then
        read(line,*,err=1110,end=1110) option,arg_str,arg_str2
        wavelength_file=arg_str2
        use_wavelength_range=.false.
        use_wavelength_file=.true.
        use_wavelength_from_refr_file=.false.
        use_wavelength_list=.false.
      elseif(arg_str.eq."from_refrac_file") then
        use_wavelength_range=.false.
        use_wavelength_file=.false.
        use_wavelength_from_refr_file=.true.
        use_wavelength_list=.false.
      elseif(arg_str.eq."list") then
        i=index(line,"list")+5
        wavelength_list=adjustl(trim(line(i:)))
        use_wavelength_range=.false.
        use_wavelength_file=.false.
        use_wavelength_from_refr_file=.false.
        use_wavelength_list=.true.
      else
        read(line,*,err=1110,end=1110) option,arg_real(1)
        wavelength=arg_real(1)
        use_wavelength_range=.false.
        use_wavelength_file=.false.
        use_wavelength_from_refr_file=.false.
        use_wavelength_list=.false.
      end if

    case('scatlib')

      read(line,*,err=1110,end=1110) option,arg_str
      scatlib_path=arg_str

    case('output')

      read(line,*,err=1110,end=1110) option,arg_str

      select case(arg_str)

        case('netcdf')

          read(line,*,err=1110,end=1110) option,arg_str,arg_str2
          output_netcdf=.true.
          output_netcdf_filename=arg_str2
          read(line,*,iostat=i) option,arg_str,arg_str2,arg_str3
          if(i.eq.0) then
            if(trim(arg_str3).eq.'hum') then
              output_netcdf_reff=.false.
            elseif(trim(arg_str3).eq.'reff') then
              output_netcdf_reff=.true.
            else
              output_netcdf_reff=.true.
            end if
          else
            output_netcdf_reff=.true.
          end if

        case('ascii_file')

          output_ascii_to_file=.true.
          read(line,*,iostat=i) option,arg_str,arg_str2
          if(i.eq.0) then
            output_ascii_filename=arg_str2
          else
            output_ascii_filename=filename ! name of input file
          end if

        case('num_coeff')

          read(line,*,err=1110,end=1110) option,arg_str,arg_int
          output_n_coeff=arg_int

        case('num_theta')

          read(line,*,err=1110,end=1110) option,arg_str,arg_int
          output_n_theta=arg_int
          if(output_n_theta.lt.2) then
            print'(A)',trim(err_str)//" num_theta<2 not allowed"
            stop
          end if

        case('theta_file')

          read(line,*,err=1110,end=1110) option,arg_str,arg_str2
          output_theta_filename=arg_str2

        case('integrated')

          output_integrated=.true.

        case('phase_function')

          output_phase_function=.true.

        case('scattering_matrix')

          output_scattering_matrix=.true.

        case('volume_scattering_function')

          output_volume_scattering_function=.true.

        case('coeff')

          output_coeff=.true.

        case('lidar')

          output_lidar=.true.

        case('digits')

          read(line,*,err=1110,end=1110) option,arg_str,arg_int
          output_digits=arg_int
          if(output_digits.lt.3.or.output_digits.gt.15) then
            print'(A)',trim(err_str)//" digits needs to be in the range from 3 to 15."
            stop
          end if

        case('header')

          output_header=.true.

        case default

          print'(A)',trim(err_str)//" Unknown option output "//trim(arg_str)
          stop

      end select

    case('debug')

      write_debug=.true.

    case('status')

      write_status=.true.

    case default

      print'(A)',trim(err_str)//" Unknown option "//trim(option)
      stop

  end select

end do

1110 print'(A)',trim(err_str)//" Option not recognized"
stop

1112 print'(A)',"Error opening "//trim(filename)
stop

1111 continue

! if the mode input option was not given, only one mode is assumed
if(num_modes.lt.1) num_modes=1

! check whether size, shape and refractive index have been given for all modes
do i_mode=1,num_modes

  if(mode(i_mode)%refr_file.eq."".and.(mode(i_mode)%mreal.eq.-1.or.mode(i_mode)%mimag.eq.-1)) then
    print'(A,I3)',"Error parsing "//trim(filename)//": No refractive index specified for mode ",i_mode
    stop
  end if

  if(mode(i_mode)%shape_type.eq."") then
    print'(A,I3)',"Error parsing "//trim(filename)//": No particle shape specified for mode ",i_mode
    stop
  end if

  if(mode(i_mode)%n.lt.0) then
    print'(A,I3)',"Error parsing "//trim(filename)//": No particle size specified for mode ",i_mode
    stop
  end if

  if(mode(i_mode)%shape_irregular_overlay.and.mode(i_mode)%shape_type(1:9).eq."irregular" ) then
    print'(A,I3)',"Error parsing "//trim(filename)//": Combination of irregular_overlay with irregular not allowed for mode ",i_mode
    stop
  end if

end do

! check whether an output option has been given
if(.not.(output_integrated.or.output_phase_function.or.output_scattering_matrix.or.output_volume_scattering_function.or.output_coeff.or.output_lidar.or.output_netcdf)) then

  write(*,'(A)') "Error parsing "//trim(filename)//": No output option specified in input file!"
  stop

end if

! check whether a wavelength has been defined
if(.not.(wavelength.gt.0.or.use_wavelength_range.or.use_wavelength_file.or.use_wavelength_from_refr_file.or.use_wavelength_list)) then

  write(*,'(A)') "Error parsing "//trim(filename)//": No wavelength specified in input file!"
  stop

end if

! convert sizes from diameter to radius if diameter-keyword was given in the input file (in the subsequent code always radius is assumed)
if(diameter.eqv..True.) then

  do i_mode=1,num_modes

    mode(i_mode)%rmin=0.5*mode(i_mode)%rmin
    mode(i_mode)%rmax=0.5*mode(i_mode)%rmax

    if(mode(i_mode)%size_type=="log_normal") then
      mode(i_mode)%size_distr_parameter(1)=0.5*mode(i_mode)%size_distr_parameter(1)
    elseif(mode(i_mode)%size_type=="mod_gamma") then
      continue
    elseif(mode(i_mode)%size_type=="mono") then
      continue
    elseif(mode(i_mode)%size_type=="distr_file") then
      continue
    elseif(mode(i_mode)%size_type=="bin_file") then
      continue
    else
      print'(A,I3)',"Error unknown size_type converting from diameter to radius for mode ", i_mode
      stop
    end if

  end do

end if

end subroutine
