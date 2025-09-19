! This file is part of Bottom RedOx Model (BROM, v.1.1).
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Shamil Yakubov,
!                     Elizaveta Protsenko, Phil Wallhead, Anfisa Berezina, 
!                     Beatriz Arellano-Nava
!-----------------------------------------------------------------------


    module io_netcdf

      use input_mod
      use types_mod
  
      use netcdf
      use fabm_omp, only: type_fabm_model => type_fabm_omp_model
      use fabm_types, only: attribute_length, rk
  
  
      implicit none
      
      private                             !everything is private unless made public
      !public functions
      public :: scan_forcing_dimensions, load_variable_year, load_variable_year_1d
      public :: init_netcdf, save_netcdf, close_netcdf

      !-------------For saving data in a netcdf file
      ! NetCDF file handle
      integer, save :: ncid = -1

      ! Dimension IDs
      integer, save :: dim_time = -1, dim_depth = -1, dim_depth_interface = -1

      ! Coordinate variable IDs
      integer, save :: var_time = -1, var_z = -1, var_z1 = -1

      ! Physics variable IDs
      integer, save :: var_T = -1, var_S = -1, var_Kz = -1
      integer, save :: var_swrad = -1, var_hice = -1
      integer, save :: var_airsea_co2 = -1   ! time-only CO2 air-sea flux

      ! FABM state variables
      integer, allocatable, save :: parameter_id(:)         ! tracer concentrations
      integer, allocatable, save :: parameter_fick_id(:)    ! tracer-specific flux across interfaces
      integer, allocatable, save :: parameter_sink_id(:)    ! tracer-specific sinking fluxes

      ! FABM diagnostics
      integer, allocatable, save :: parameter_id_diag(:)

      ! Internal bookkeeping
      logical, save :: nc_ready = .false.
      integer, save :: time_index = 0
  
  
  
      contains
  !=======================================================================================================================
        !------------------------------------------------------------
        ! Scan metadata (depth, year/day_of_year coverage)
        !------------------------------------------------------------
        subroutine scan_forcing_dimensions(filename, years, year_start_idx, year_last_idx, days_in_year, nrecs_in_year, &
                                           start_year, first_day, last_day, repeat_forcing_year, use_hice, use_swradWm2, z_w)

          !------------------- args -------------------
          character(*), intent(in)  :: filename
          integer,      allocatable, intent(out) :: years(:), year_start_idx(:), year_last_idx(:)
          integer,      allocatable, intent(out) :: days_in_year(:), nrecs_in_year(:)
          real(rk),     allocatable, intent(out) :: z_w(:)

          integer,      intent(in) :: start_year, first_day, last_day
          integer,      intent(in) :: repeat_forcing_year, use_hice, use_swradWm2
          !------------------- locals -----------------
          class(type_input), allocatable :: nc
          real(rk), allocatable :: tmp(:), full_depth(:)
          integer,  allocatable :: year_in(:), doy_in(:)
          integer :: nt, nz, i, k, kstart, kk, rem                                          
  
          
          ! Open file
          if (allocated(nc)) deallocate(nc)
          allocate(nc); nc = type_input(filename)

          ! depth axis
          if (.not. nc%var_exists("depth")) stop 'FATAL: Missing depth in forcing file'
          nz = nc%get_1st_dim_length("depth")
          allocate(z_w(nz))
          z_w = abs(nc%get_column("depth"))

          ! detect time vars
          if (.not. nc%var_exists("year") .or. (.not.(nc%var_exists("doy") .or. nc%var_exists("day_of_year")))) then
            stop 'FATAL: Require numeric year and doy/day_of_year in forcing file'
          end if

          ! Read year(:) and convert to integer
          tmp = nc%get_column("year")
          nt  = size(tmp)
          if (nt <= 0) stop 'FATAL (io_netcdf): number of years <= 0.'

          allocate(year_in(nt))
          year_in = int(tmp)

          ! Read doy/day_of_year and convert to integer
          if (nc%var_exists("doy")) then
            tmp = nc%get_column('doy')
          else
            tmp = nc%get_column('day_of_year')
          end if
          allocate(doy_in(nt))
          doy_in = int(tmp)

          ! Build unique-year bins
          k = 1
          do i = 2, nt
            if (year_in(i) /= year_in(i-1)) k = k + 1
          end do

          if (allocated(years)) deallocate(years, year_start_idx, year_last_idx)
          allocate(years(k), year_start_idx(k), year_last_idx(k))

          k = 1
          years(k) = year_in(1)
          year_start_idx(k) = 1
          do i = 2, nt
            if (year_in(i) /= year_in(i-1)) then
              year_last_idx(k) = i - 1
              k = k + 1
              years(k) = year_in(i)
              year_start_idx(k) = i
            end if
          end do
          year_last_idx(k) = nt

          ! Per-year summaries
          if (allocated(days_in_year)) deallocate(days_in_year, nrecs_in_year)
          allocate(days_in_year(size(years)), nrecs_in_year(size(years)))
          do k = 1, size(years)
            days_in_year(k)  = max(1, maxval(doy_in(year_start_idx(k):year_last_idx(k))))
            nrecs_in_year(k) = year_last_idx(k) - year_start_idx(k) + 1
          end do

          ! Print summary
          write(*,'(a)') 'Forcing years detected:  year   days   records'
          do k = 1, size(years)
            write(*,'(i8,2x,i5,2x,i8)') years(k), days_in_year(k), nrecs_in_year(k)
          end do

          ! Coverage checks
          kstart = findloc(years, start_year, dim=1)
          if (kstart == 0) then
            write(*,*) 'FATAL (io_netcdf): start_year ', start_year, ' not found in forcing.'
            stop
          end if

          if (first_day < 1 .or. first_day > days_in_year(kstart)) then
            write(*,*) 'FATAL (io_netcdf): first_day=', first_day, ' outside year ', years(kstart), &
                      ' [1..', days_in_year(kstart), '].'
            stop
          end if

          if (repeat_forcing_year == 1) then
            if (last_day < first_day) then
              write(*,*) 'FATAL (io_netcdf): last_day < first_day (', last_day, ' < ', first_day, ').'
              stop
            end if
          else
            rem = last_day - first_day + 1
            if (rem <= 0) then
              write(*,*) 'FATAL (io_netcdf): last_day < first_day (', last_day, ' < ', first_day, ').'
              stop
            end if

            rem = rem - (days_in_year(kstart) - (first_day - 1))
            kk  = kstart + 1
            do while (rem > 0)
              if (kk > size(years)) then
                write(*,*) 'FATAL (io_netcdf): forcing ends before last_day is reached. Need ', rem, &
                          ' more day(s) after year ', years(kk-1), '.'
                stop
              end if
              rem = rem - days_in_year(kk)
              kk  = kk + 1
            end do
          end if

          ! Variable presence (required & optional)
          if (.not. nc%var_exists('temperature')) stop 'FATAL (io_netcdf): Missing variable "temperature"'
          if (.not. nc%var_exists('salinity'))    stop 'FATAL (io_netcdf): Missing variable "salinity"'
          if (.not. nc%var_exists('Kz'))          stop 'FATAL (io_netcdf): Missing variable "Kz"'
          if (use_hice == 1) then
            if (.not. nc%var_exists('hice')) stop 'FATAL (io_netcdf): Missing variable "hice"'
          end if
          if (use_swradWm2 == 1) then
            if (.not. nc%var_exists('swradWm2')) stop 'FATAL (io_netcdf): Missing variable "swradWm2"'
          end if

          deallocate(nc)

      end subroutine scan_forcing_dimensions


  !=======================================================================================================================
      !------------------------------------------------------------
      ! Load 2D var (depth × time) for a given year
      ! Assumes var(:,:) with depth dimension present.
      !   Handles both (depth,time) and (time,depth) layouts.
      ! Inputs:
      !   varname  - name of variable (e.g. "temperature")
      !   kyear    - index into years(:), i.e. which year bin
      !
      ! Output:
      !   out(nz, nrec) - variable slice for that year
      !                   (depth × records_in_year)
      !------------------------------------------------------------
      subroutine load_variable_year(filename, varname, year, years, year_start_idx, year_last_idx, out)
          character(*), intent(in)  :: filename, varname
          integer, intent(in)          :: year
          integer, intent(in)          :: years(:), year_start_idx(:), year_last_idx(:)
          real(rk), allocatable, intent(out) :: out(:,:)
          type(type_input) :: nc

          real(rk), allocatable :: full(:,:)
          integer :: k, nrec, nz

          nc = type_input(filename)

          ! find which bin this year belongs to
          k = findloc(years, year, dim=1)
          if (k <= 0) stop 'FATAL: year not in forcing file'

          nrec = year_last_idx(k) - year_start_idx(k) + 1
          nz = nc%get_1st_dim_length("depth")
          full = nc%get_array(trim(varname))
          
          if (allocated(out)) deallocate(out)
          allocate(out(nz, nrec))

          if (size(full,1) == nz) then
            ! Layout (depth, time)
            if (size(full,2) < year_last_idx(k)) then
              stop 'FATAL (io_netcdf): not enough time records in variable '//trim(varname)
            end if 
            if (size(full,2) == nz)  write(*,*) 'WARNING: forcing variable "', trim(varname), '" has shape (nz,nz). Assuming layout (depth,time).'
            out(:,:) = full(:, year_start_idx(k):year_last_idx(k))

          else if (size(full,2) == nz) then
            ! Layout (time, depth)   
            if (size(full,1) < year_last_idx(k)) then
              stop 'FATAL (io_netcdf): not enough time records in variable '//trim(varname)
            end if         
            out(:,:) = transpose(full(year_start_idx(k):year_last_idx(k), :))        
          else
            stop 'FATAL (io_netcdf): variable "'//trim(varname)//'" has unexpected shape'
          end if
      end subroutine load_variable_year

      !------------------------------------------------------------
      ! Load 1D var (time only) for a given year
      !   Assumes var(:) with only a time dimension.
      ! Inputs:
      !   varname  - name of variable (e.g. "swradWm2")
      !   year     - actual year number
      !
      ! Output:
      !   out(nrec) - time series for that year
      !------------------------------------------------------------
      subroutine load_variable_year_1d(filename, varname, year, years, year_start_idx, year_last_idx, out)
          character(*), intent(in)  :: filename, varname
          integer, intent(in)          :: year
          integer, intent(in)          :: years(:), year_start_idx(:), year_last_idx(:)
          real(rk), allocatable, intent(out) :: out(:)
    
          type(type_input) :: nc
          real(rk), allocatable :: full(:)
          integer :: k, nrec

          nc = type_input(filename)

          k = findloc(years, year, dim=1)
          if (k == 0) stop 'FATAL (io_netcdf): requested year not found in forcing file'

          nrec = year_last_idx(k) - year_start_idx(k) + 1
          full = nc%get_column(trim(varname))
          ! Defensive check: does the file have enough records?
          if (size(full) < year_last_idx(k)) then
            write(*,*) 'FATAL (io_netcdf): variable "', trim(varname), &
                      '" has only', size(full), 'records, need at least', year_last_idx(k)
            stop
          end if

          ! Allocate and slice to this year
          if (allocated(out)) deallocate(out)
          allocate(out(nrec))
          out(:) = full(year_start_idx(k):year_last_idx(k))
      end subroutine load_variable_year_1d
  
  !=======================================================================================================================


      subroutine init_netcdf(filename, k_max, z, z1, model, use_hice, start_year)
  
        !Input variables
        character(*), intent(in) :: filename
        integer,      intent(in) :: k_max, use_hice, start_year
        real(rk),     intent(in) :: z(:), z1(:)
        class(type_fabm_model), pointer, intent(in) :: model

        integer :: ntr, ndiag, i, ilast
        character(len=128) :: vname
        character(len=16)  :: yearstr

        print *, 'Initialising NetCDF version: ', trim(nf90_inq_libvers())
        write(yearstr,'(I4)') start_year

        ! Create file
        call check_err(nf90_create(trim(filename), nf90_clobber, ncid), "creating file")

        ! Dimensions
        call check_err(nf90_def_dim(ncid, "time", nf90_unlimited, dim_time), "defining time dimension")
        call check_err(nf90_def_dim(ncid, "depth", k_max, dim_depth), "defining depth dimension")
        call check_err(nf90_def_dim(ncid, "depth_interface", k_max+1, dim_depth_interface), "defining depth_interface dimension")

        ! Coordinates
        call check_err(nf90_def_var(ncid, "time", nf90_double, (/dim_time/), var_time), "defining time variable")
        call check_err(nf90_put_att(ncid, var_time, "units", "days since "//trim(yearstr)//"-01-01 00:00:00"))
        call check_err(nf90_put_att(ncid, var_time, "standard_name", "time"))
        call check_err(nf90_put_att(ncid, var_time, "long_name", "time"))
        call check_err(nf90_put_att(ncid, var_time, "axis", "T"))

        call check_err(nf90_def_var(ncid, "depth", nf90_double, (/dim_depth/), var_z), "defining depth variable")
        call check_err(nf90_put_att(ncid, var_z, "units", "m"))
        call check_err(nf90_put_att(ncid, var_z, "positive", "down"))
        call check_err(nf90_put_att(ncid, var_z, "long_name", "depth"))
        call check_err(nf90_put_att(ncid, var_z, "axis", "Z"))

        call check_err(nf90_def_var(ncid, "depth_interface", nf90_double, (/dim_depth_interface/), var_z1), "defining depth_interface variable")
        call check_err(nf90_put_att(ncid, var_z1, "units", "m"))
        call check_err(nf90_put_att(ncid, var_z1, "positive", "down"))
        call check_err(nf90_put_att(ncid, var_z1, "long_name", "depth at interfaces"))

        ! Physics variables
        call check_err(nf90_def_var(ncid, "temperature", nf90_double, (/dim_depth, dim_time/), var_T), "defining temperature")
        call check_err(nf90_put_att(ncid, var_T, "units", "degree_Celsius"))
        call check_err(nf90_put_att(ncid, var_T, "long_name", "sea water temperature"))

        call check_err(nf90_def_var(ncid, "salinity", nf90_double, (/dim_depth, dim_time/), var_S), "defining salinity")
        !call check_err(nf90_put_att(ncid, var_S, "units", "1e-3"))
        call check_err(nf90_put_att(ncid, var_S, "long_name", "sea water salinity"))

        call check_err(nf90_def_var(ncid, "Kz", nf90_double, (/dim_depth_interface, dim_time/), var_Kz), "defining Kz")
        call check_err(nf90_put_att(ncid, var_Kz, "units", "m2 s-1"))
        call check_err(nf90_put_att(ncid, var_Kz, "long_name", "vertical eddy diffusivity at interfaces"))


        call check_err(nf90_def_var(ncid, "swrad", nf90_double, (/dim_time/), var_swrad), "defining swrad")
        call check_err(nf90_put_att(ncid, var_swrad, "units", "W m-2"))
        call check_err(nf90_put_att(ncid, var_swrad, "long_name", "surface downward shortwave radiation"))

        if (use_hice /= 0) then
          call check_err(nf90_def_var(ncid, "hice", nf90_double, (/dim_time/), var_hice), "defining hice")
          call check_err(nf90_put_att(ncid, var_hice, "units", "m"))
          call check_err(nf90_put_att(ncid, var_hice, "long_name", "sea ice thickness"))
        end if

        ! Scalars
        call check_err(nf90_def_var(ncid, "air_sea_flux_CO2", nf90_double, (/dim_time/), var_airsea_co2), "defining air_sea_flux")
        call check_err(nf90_put_att(ncid, var_airsea_co2, "units", "mol m-2 d-1"))
        call check_err(nf90_put_att(ncid, var_airsea_co2, "long_name", "air-sea flux of CO2"))

        ! FABM tracers
        ntr = size(model%interior_state_variables)
        allocate(parameter_id(ntr), parameter_fick_id(ntr), parameter_sink_id(ntr))
        do i = 1, ntr
          vname = model%interior_state_variables(i)%path
          ilast = index(vname, '/', .true.)
          if (ilast > 0) vname = vname(ilast+1:)

          ! main tracer concentration
          call check_err(nf90_def_var(ncid, trim(vname), nf90_double, (/dim_depth, dim_time/), parameter_id(i)), "defining tracer "//trim(vname))
          if (len_trim(model%interior_state_variables(i)%units) > 0) &
              call check_err(nf90_put_att(ncid, parameter_id(i), "units", trim(model%interior_state_variables(i)%units)))
          call check_err(nf90_put_att(ncid, parameter_id(i), "long_name", trim(vname)))
          
          ! interface fluxes
          call check_err(nf90_def_var(ncid, "fick_"//trim(vname), nf90_double, (/dim_depth_interface, dim_time/), parameter_fick_id(i)), "defining fick_"//trim(vname))
          call check_err(nf90_put_att(ncid, parameter_fick_id(i), "units", "mol m-2 d-1"))

          call check_err(nf90_def_var(ncid, "sink_"//trim(vname), nf90_double, (/dim_depth_interface, dim_time/), parameter_sink_id(i)), "defining sink_"//trim(vname))
          call check_err(nf90_put_att(ncid, parameter_sink_id(i), "units", "mol m-2 d-1"))
        end do

        ! FABM diagnostics
        ndiag = size(model%interior_diagnostic_variables)
        allocate(parameter_id_diag(ndiag))
        do i = 1, ndiag
          if (model%interior_diagnostic_variables(i)%save) then
              vname = model%interior_diagnostic_variables(i)%path
              ilast = index(vname, '/', .true.)
              if (ilast > 0) vname = vname(ilast+1:)

              call check_err(nf90_def_var(ncid, trim(vname), nf90_double, (/dim_depth, dim_time/), parameter_id_diag(i)), "defining diagnostic "//trim(vname))
              if (len_trim(model%interior_diagnostic_variables(i)%units) > 0) &
                call check_err(nf90_put_att(ncid, parameter_id_diag(i), "units", trim(model%interior_diagnostic_variables(i)%units)))
              call check_err(nf90_put_att(ncid, parameter_id_diag(i), "long_name", trim(vname)))
          else
              parameter_id_diag(i) = -1
          end if
        end do

        call check_err(nf90_enddef(ncid), "ending NetCDF definition")

        ! Write depth coordinates once
        call check_err(nf90_put_var(ncid, var_z,  z), "writing depth")
        call check_err(nf90_put_var(ncid, var_z1, z1), "writing depth_interface")

        ! Init counters
        time_index = 0
        nc_ready   = .true.

  
      end subroutine init_netcdf
  !=======================================================================================================================
  
    
  
  
  
  !=======================================================================================================================
      subroutine save_netcdf(k_max, cc, t, s, kz, &
                             model, swradWm2,use_hice, hice, &
                             fick_per_day, sink_per_day, air_sea_flux_co2, &
                             time_output)
  
        ! Arguments
        integer, intent(in)                  :: k_max, use_hice
        real(rk), intent(in)                 :: time_output        ! time in days since start
        real(rk), dimension(:,:), intent(in) :: cc              ! tracer concentrations (depth x tracer)        
        real(rk), dimension(:,:), intent(in) :: fick_per_day    ! fluxes across interfaces (depth+1 x tracer)
        real(rk), dimension(:,:), intent(in) :: sink_per_day    ! sinking fluxes across interfaces (depth+1 x tracer)
        real(rk), dimension(:),   intent(in) :: t, s, kz        ! physical profiles
        class(type_fabm_model), pointer      :: model
        real(rk), intent(in)                 :: swradWm2, hice  ! surface forcing scalars
        real(rk), intent(in) :: air_sea_flux_co2                ! air-sea flux of CO2
                            
        ! Locals
        integer :: ip
        real(rk), allocatable :: diag(:)
        real(rk) :: tmp(1)    ! buffer for scalar writes (NetCDF requires rank-1 array, not bare scalar)
                            
        if (.not. nc_ready) stop "FATAL: save_netcdf called before init_netcdf"
      
        ! Increment time record counter
        time_index = time_index + 1
      
        !------------------------------------------------------------
        ! Write time coordinate
        ! Note: tmp(1) is used because nf90_put_var expects a vector,
        ! even for length-1 writes.
        !------------------------------------------------------------
        tmp(1) = time_output
        call check_err(nf90_put_var(ncid, var_time, tmp, start=(/time_index/), count=(/1/)), "writing time")
      
        !------------------------------------------------------------
        ! Write physics profiles (depth-resolved)
        !------------------------------------------------------------
        call put_column(var_T,   t, "temperature")
        call put_column(var_S,   s, "salinity")
        call put_column(var_Kz,  kz, "Kz")
      
        !------------------------------------------------------------
        ! Write scalar forcings
        !------------------------------------------------------------      
        if (use_hice == 1) then
           tmp(1) = hice
           call check_err(nf90_put_var(ncid, var_hice, tmp, start=(/time_index/), count=(/1/)), "writing hice")
        end if
      
        tmp(1) = swradWm2
        call check_err(nf90_put_var(ncid, var_swrad, tmp, start=(/time_index/), count=(/1/)), "writing swrad")

        !------------------------------------------------------------
        ! Write air-sea flux of CO2
        !------------------------------------------------------------
        tmp(1) = air_sea_flux_co2
        call check_err(nf90_put_var(ncid, var_airsea_co2, tmp, start=(/time_index/), count=(/1/)), "writing air_sea_flux_CO2")


      
        !------------------------------------------------------------
        ! Write tracer concentrations and fluxes
        !------------------------------------------------------------
        do ip = 1, size(model%interior_state_variables)
           call check_err(nf90_put_var(ncid, parameter_id(ip), cc(:,ip), start=(/1,time_index/), count=(/k_max,1/)), "writing tracer "//trim(model%interior_state_variables(ip)%name))
           call check_err(nf90_put_var(ncid, parameter_fick_id(ip), fick_per_day(:,ip), start=(/1,time_index/), count=(/k_max+1,1/)), "writing fick "//trim(model%interior_state_variables(ip)%name))
           call check_err(nf90_put_var(ncid, parameter_sink_id(ip), sink_per_day(:,ip), start=(/1,time_index/), count=(/k_max+1,1/)), "writing sink "//trim(model%interior_state_variables(ip)%name))
        end do
      
        !------------------------------------------------------------
        ! Write diagnostics (only those flagged with %save)
        !------------------------------------------------------------
        do ip = 1, size(model%interior_diagnostic_variables)
           if (model%interior_diagnostic_variables(ip)%save) then
              diag = model%get_interior_diagnostic_data(ip)
              ! Guard against extreme values that NetCDF cannot handle
              if (maxval(abs(diag)) < 1.0e37_rk) then
                 call check_err(nf90_put_var(ncid, parameter_id_diag(ip), diag, start=(/1,time_index/), count=(/size(diag),1/)), "writing diagnostic "//trim(model%interior_diagnostic_variables(ip)%name))
              end if
           end if
        end do
      
        !------------------------------------------------------------
        ! Flush buffers to disk (safer for long runs, avoids data loss if crash)
        !------------------------------------------------------------
        call check_err(nf90_sync(ncid), "syncing NetCDF")
        
        contains
          !------------------------------------------------------------
          ! Function to write a depth profile (vector) for this timestep
          !------------------------------------------------------------
          subroutine put_column(varid, arr, label)
            integer, intent(in) :: varid
            real(rk), intent(in) :: arr(:)
            character(*), intent(in) :: label
            integer :: st_local
            st_local = nf90_put_var(ncid, varid, arr, start=(/1,time_index/), count=(/size(arr),1/))
            call check_err(st_local, "writing "//trim(label))
          end subroutine put_column
  
      end subroutine save_netcdf
  !======================================================================================================================= 
   
  
  
  
  !=======================================================================================================================
      subroutine close_netcdf()
        use netcdf
        implicit none
        if (nc_ready) then
           call check_err(nf90_close(ncid), "closing NetCDF file")
      
           if (allocated(parameter_id))       deallocate(parameter_id)
           if (allocated(parameter_fick_id))  deallocate(parameter_fick_id)
           if (allocated(parameter_sink_id))  deallocate(parameter_sink_id)
           if (allocated(parameter_id_diag))  deallocate(parameter_id_diag)
      
           nc_ready = .false.
           ncid     = -1
           write(*,'(a)') "NetCDF file closed successfully"
        end if
      end subroutine close_netcdf
  !=======================================================================================================================
  
  
  
  
  
  !===============Helpers=================================================================================================
  
    subroutine lower_bound_int(a, x, idx)
        integer, intent(in) :: a(:), x
        integer, intent(out):: idx
        integer :: lo, hi, mid, n
        n = size(a); lo = 1; hi = n+1
        do while (lo < hi)
          mid = (lo + hi)/2
          if (a(mid) < x) then
            lo = mid + 1
          else
            hi = mid
          end if
        end do
        idx = min(lo, n)
    end subroutine

    subroutine upper_bound_int(a, x, idx)
        integer, intent(in) :: a(:), x
        integer, intent(out):: idx
        integer :: lo, hi, mid, n
        n = size(a); lo = 1; hi = n+1
        do while (lo < hi)
          mid = (lo + hi)/2
          if (a(mid) <= x) then
            lo = mid + 1
          else
            hi = mid
          end if
        end do
        idx = min(lo, n)
    end subroutine

    pure logical function is_leap_gregorian(y) result(isleap)
      integer, intent(in) :: y
      isleap = (mod(y,4)==0 .and. (mod(y,100)/=0 .or. mod(y,400)==0))
    end function

  !=======================================================================================================================
  
  
  
    
  !=======================================================================================================================
      subroutine check_err(status, msg)
        integer, intent(in) :: status
        character(*), intent(in), optional :: msg
      
        if (status /= NF90_NOERR) then
           print *, "NetCDF error: ", trim(nf90_strerror(status))
           if (present(msg)) print *, "  while: ", trim(msg)
           stop "Stopped due to NetCDF error"
        end if
      end subroutine check_err
  !=======================================================================================================================
  
  
  
  
  !=======================================================================================================================
      subroutine svan(s, t, po, sigma)
  !/*
  !c------specific volume anomaly based on 1980 equation of state for
  !c      seawater and 1978 practical salinity scale
  !c      pressure          PO     decibars
  !c      temperature        T     degree celsius (IPTS-68)
  !c      salinitty          S     (PSS-78)
  !c      spec.vol.anom.  SVAN     1.0E-8 m**3/Kg
  !c      density anom.   SIGMA    Kg/m**3
  !*/
  
      real(rk)  p,sig,sr,r1,r2,r3,s,t,po,sigma
      real(rk)  a,b,c,d,e,a1,b1,aw,bw,k,ko,kw,k35,v350p,sva,dk,gam,pk,dr35p,dvan
  
      real(rk) r3500, r4 ,dr350
      data  r3500 /1028.1063/, r4/4.8314E-4/,dr350/28.106331/
  
      p=po/10.
      sr=sqrt(abs(s))
  
      r1= ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*t &
             -9.095290E-3)*t+6.793952E-2)*t-28.263737
      r2= (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t &
            +8.24493E-1
      r3= (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3
  
      sig=(r4*s + r3*sr + r2)*s +r1
  
      v350p=1.0/r3500
      sva=-sig*v350p/(r3500+sig)
      sigma= sig + dr350
  
      if (p.eq.0.0) return
  
      e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
         bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
      b = bw + e*s
  
      d= 1.91075E-4
      c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
      aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
      a = (d*sr + c)*s + aw
  
      b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
      a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
      kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t &
             +148.4206)*t-1930.06
      ko = (b1*sr + a1)*s + kw
  
      dk = (b*p+a)*p+ko
      k35 = (5.03217E-5*p+3.359406)*p+21582.27
      gam=p/k35
      pk=1.0-gam
      sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))
  
      v350p= v350p*pk
      dr35p=gam/v350p
      dvan= sva/(v350p*(v350p+sva))
      sigma = dr350 + dr35p -dvan
      return
      end subroutine svan
  !=======================================================================================================================
  
  
      end module io_netcdf
  