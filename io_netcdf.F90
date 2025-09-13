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
      public :: open_forcing_file         ! open NetCDF, detect layout, build year bins
      public :: load_variable_year        ! fetch vertical profiles per year
      public :: load_variable_year_1d        ! fetch vertical profiles per year
      public :: find_year_index
      public :: days_in_year, years, year_start_idx, year_last_idx, nrecs_in_year
      public init_netcdf, save_netcdf, close_netcdf

      logical, save :: forcing_init = .false.
      class(type_input), allocatable, save :: nc

      integer, save :: nt = 0, nz = 0
      character(len=16), save :: depth_name = 'depth'
      logical, save :: dims_are_depth_time = .true.   ! flip if the file uses (time,depth)

      ! time axes (we keep only year+doy; 'time' strings are converted into these)
      integer,  allocatable, save :: year_in(:), doy_in(:)

      ! per-year bins: year k covers indices [year_start_idx(k) .. year_last_idx(k)-1]
      integer, allocatable, save :: years(:)           ! list of years present
      integer, allocatable, save :: year_start_idx(:)  ! first forcing index of each year
      integer, allocatable, save :: year_last_idx(:)    ! last  forcing index of each year
      integer, allocatable, save :: days_in_year(:)    ! number of days in each year
      integer, allocatable, save :: nrecs_in_year(:)   ! number of records in each year



      !netCDF file id
      integer               :: nc_id
      integer, allocatable  :: parameter_id(:)
      integer, allocatable  :: parameter_fick_id(:)
      integer, allocatable  :: parameter_sink_id(:)
      integer, allocatable  :: parameter_id_diag(:)
  
      integer               :: i_id, z_id, z2_id, time_id, swradWm2_id, hice_id
      integer               :: pH_id, T_id, S_id, Kz_id, Kz_sol_id, Kz_par_id, w_sol_id, w_par_id, gas_air_sea_id
      integer               :: pCO2_id, Om_Ca_id, Om_Ar_id
  
      logical               :: first
  
  
  
      contains
  !=======================================================================================================================
      subroutine open_forcing_file(z_w)

        use io_ascii, only: get_brom_name, get_brom_par

        character(len=:), allocatable :: ncfile
        integer :: start_year, first_day, last_day, repeat_forcing_year
        integer :: use_hice, use_swradWm2    ! Optional variables to load

        logical :: has_year, has_doy, has_day_of_year
        integer :: i, k, kstart, kk, rem
        real(rk), allocatable :: tmp(:)
        real(rk), allocatable, intent(out) :: z_w(:)
         
        ! Read settings
        ncfile = get_brom_name("ncinfile_name")
        start_year = get_brom_par("start_year")
        first_day = get_brom_par("first_day")
        last_day = get_brom_par("last_day")
        repeat_forcing_year = get_brom_par("repeat_forcing_year")
        use_hice    = get_brom_par("use_hice")
        use_swradWm2= get_brom_par("use_swradWm2")

        ! Open file
        if (allocated(nc)) deallocate(nc)
        allocate(nc); nc = type_input(ncfile)

        if (.not. nc%var_exists(depth_name)) stop 'FATAL (io_netcdf): Missing depth dimension in netcdf file'
        nz = nc%get_1st_dim_length(depth_name); if (nz<=0) stop 'FATAL (io_netcdf): Depth dimension <=0'
        if (allocated(z_w)) deallocate(z_w)
        allocate(z_w(nz))
        z_w = nc%get_column(depth_name)

        ! --- detect layout ---
        has_year        = nc%var_exists('year')
        has_doy         = nc%var_exists('doy')
        has_day_of_year = nc%var_exists('day_of_year')

        if (.not. has_year .or. (.not. has_doy .and. .not. has_day_of_year)) then
          stop 'FATAL (io_netcdf): Require numeric "year" and "doy" (or "day_of_year")'
        end if

        ! Read year(:) as real -> int
        tmp = nc%get_column('year')
        nt  = size(tmp); if (nt<=0) stop 'FATAL (io_netcdf): nt <= 0 (year)'
        if (allocated(year_in)) deallocate(year_in, doy_in)
        allocate(year_in(nt)); year_in = int(tmp)

        ! Read doy/day_of_year as real -> int
        if (has_doy) then
          tmp = nc%get_column('doy')
        else
          tmp = nc%get_column('day_of_year')
        end if
        allocate(doy_in(nt));  doy_in  = int(tmp)

        ! --- Build per-year bins: years(k) covers indices [year_start_idx(k) .. year_last_idx(k)] ---
        k = 1
        do i = 2, nt
          if (year_in(i) /= year_in(i-1)) k = k + 1
        end do
        if (allocated(years)) deallocate(years, year_start_idx, year_last_idx)
        allocate(years(k), year_start_idx(k), year_last_idx(k))

        k = 1; years(k) = year_in(1); year_start_idx(k) = 1
        do i = 2, nt
          if (year_in(i) /= year_in(i-1)) then
            year_last_idx(k) = i-1
            k = k + 1
            years(k) = year_in(i); year_start_idx(k) = i
          end if
        end do
        year_last_idx(k) = nt

        ! Per-year summaries (unique days present; total records)
        if (allocated(days_in_year)) deallocate(days_in_year, nrecs_in_year)
        allocate(days_in_year(size(years)), nrecs_in_year(size(years)))
        do k=1, size(years)
          days_in_year(k)  = max(1, maxval(doy_in(year_start_idx(k):year_last_idx(k))))
          nrecs_in_year(k) = year_last_idx(k) - year_start_idx(k) + 1
        end do

        ! Print Summary
        write(*,'(a)') 'Forcing years detected:  year   days   records'
        do k=1, size(years)
          write(*,'(i8,2x,i5,2x,i8)') years(k), days_in_year(k), nrecs_in_year(k)
        end do

        ! --- Sanity check: coverage for simulation period ---
        kstart = find_year_index(start_year)
        if (kstart == 0) then
          write(*,*) 'FATAL (io_netcdf): start_year ', start_year, ' not found in forcing.'
          stop
        end if

        ! day bounds in the start year
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
          ! No repeat: we must have enough consecutive years in the file to cover last_day
          rem = last_day - first_day + 1
          if (rem <= 0) then
            write(*,*) 'FATAL (io_netcdf): last_day < first_day (', last_day, ' < ', first_day, ').'
            stop
          end if

          ! consume remaining days across bins until satisfied
          rem = rem - (days_in_year(kstart) - (first_day - 1))
          kk = kstart + 1
          do while (rem > 0)
            if (kk > size(years)) then
              write(*,*) 'FATAL (io_netcdf): forcing ends before last_day is reached. Need ', rem, &
                        ' more day(s) after year ', years(kk-1), '.'
              stop
            end if
            rem = rem - days_in_year(kk)
            kk = kk + 1
          end do
        end if        

        ! -------- sanity check: variables present (required & optional-enabled) --------
        if (.not. nc%var_exists('temperature')) stop 'FATAL (io_netcdf): Missing variable "temperature" in netcdf file'
        if (.not. nc%var_exists('salinity')) stop 'FATAL (io_netcdf): Missing variable "salinity" in netcdf file'
        if (.not. nc%var_exists('Kz')) stop 'FATAL (io_netcdf): Missing variable "Kz" in netcdf file'
        if (use_hice.eq.1) then
          if (.not. nc%var_exists('hice')) stop 'FATAL (io_netcdf): Missing variable "hice" in netcdf file'
        end if
        if (use_swradWm2.eq.1) then
          if (.not. nc%var_exists('swradWm2')) stop 'FATAL (io_netcdf): Missing variable "swradWm2" in netcdf file'
        end if

        forcing_init = .true.

      end subroutine open_forcing_file
    



      !------------------------------------------------------------
      ! Load a single year's forcing data for one variable.
      !   Assumes var(:,:) with depth dimension present.
      !   Handles both (depth,time) and (time,depth) layouts.
      ! Inputs:
      !   varname  - name of variable (e.g. "temperature")
      !   kyear    - index into years(:), i.e. which year bin
      !
      ! Output:
      !   out(nz, nrec) - variable slice for that year
      !                   (depth × records_in_year)
      !------------------------------------------------------------
      subroutine load_variable_year(varname, year, out)
        use types_mod, only: rk
        implicit none
        character(*), intent(in)  :: varname
        integer,      intent(in)  :: year      ! actual year number
        real(rk),     allocatable, intent(out) :: out(:,:)
      
        real(rk), allocatable :: full(:,:)
        integer :: k, nrec
      
        if (.not. forcing_init) stop 'FATAL (io_netcdf): call open_forcing_file first'
      
        ! find which bin this year belongs to
        k = find_year_index(year)
        if (k == 0) stop 'FATAL (io_netcdf): requested year not found in forcing file'
      
        nrec = year_last_idx(k) - year_start_idx(k) + 1
      
        full = nc%get_array(trim(varname))
      
        if (size(full,1) == nz) then
          ! Layout (depth, time)
          allocate(out(nz, nrec))
          out(:,:) = full(:, year_start_idx(k):year_last_idx(k))
      
        else if (size(full,2) == nz) then
          ! Layout (time, depth)
          allocate(out(nz, nrec))
          out(:,:) = transpose(full(year_start_idx(k):year_last_idx(k), :))
      
        else
          stop 'FATAL (io_netcdf): variable "'//trim(varname)//'" has unexpected shape'
        end if
      end subroutine load_variable_year 


      !------------------------------------------------------------
      ! Load a single year's forcing data for a time-only variable.
      !   Assumes var(:) with only a time dimension.
      ! Inputs:
      !   varname  - name of variable (e.g. "swradWm2")
      !   year     - actual year number
      !
      ! Output:
      !   out(nrec) - time series for that year
      !------------------------------------------------------------
      subroutine load_variable_year_1d(varname, year, out)
        use types_mod, only: rk
        implicit none
        character(*), intent(in)  :: varname
        integer,      intent(in)  :: year      ! actual year number
        real(rk),     allocatable, intent(out):: out(:)

        real(rk), allocatable :: full(:)
        integer :: k, nrec

        if (.not. forcing_init) stop 'FATAL (io_netcdf): call open_forcing_file first'

        ! find which bin this year belongs to
        k = find_year_index(year)
        if (k == 0) stop 'FATAL (io_netcdf): requested year not found in forcing file'

        nrec = year_last_idx(k) - year_start_idx(k) + 1

        full = nc%get_column(trim(varname))

        ! Slice to year
        allocate(out(nrec))
        out(:) = full(year_start_idx(k):year_last_idx(k))
      end subroutine load_variable_year_1d
  
  
  
  
  !=======================================================================================================================
      subroutine init_netcdf(fn, k_max, model, use_swradWm2, use_hice, year)
  
      !Input variables
      character(len=*), intent(in)     :: fn
      integer, intent(in)              :: k_max, use_swradWm2, use_hice, year
      class (type_fabm_model), pointer :: model
  
      !Local variables
      integer                          :: z_dim_id, z2_dim_id, time_dim_id
      integer                          :: ip, iret, ilast
      integer, parameter               :: time_len = NF90_UNLIMITED
      character(len=4)                 :: yearstr
      integer                          :: dim1d
      integer                          :: dim_ids(2), dim_ids2(2), dim_ids0(1)
  
    write(*,*) "k_max = ", k_max
  
      first = .true.
      print *, 'NetCDF version: ', trim(nf90_inq_libvers())
      nc_id = -1
      call check_err(nf90_create(fn, NF90_CLOBBER, nc_id))
  
      !Define the dimensions
      call check_err(nf90_def_dim(nc_id, "z", k_max, z_dim_id))
      call check_err(nf90_def_dim(nc_id, "z2", k_max+1, z2_dim_id))
      call check_err(nf90_def_dim(nc_id, "time", time_len, time_dim_id))
  
      !Define coordinates
      dim1d = z_dim_id
      call check_err(nf90_def_var(nc_id, "z", NF90_REAL, dim1d, z_id))
      call check_err(nf90_put_att(nc_id, z_id, "positive", "down"))
      call check_err(nf90_put_att(nc_id, z_id, "long_name", "depth at layer midpoints"))
      call check_err(nf90_put_att(nc_id, z_id, "units", "metres"))
      call check_err(nf90_put_att(nc_id, z_id, "axis", "Z"))
      dim1d = z2_dim_id
      call check_err(nf90_def_var(nc_id, "z2", NF90_REAL, dim1d, z2_id))
      call check_err(nf90_put_att(nc_id, z2_id, "positive", "down"))
      call check_err(nf90_put_att(nc_id, z2_id, "long_name", "depth at layer interfaces"))
      call check_err(nf90_put_att(nc_id, z2_id, "units", "metres"))
      call check_err(nf90_put_att(nc_id, z2_id, "axis", "Z"))
      dim1d = time_dim_id
      call check_err(nf90_def_var(nc_id, "time", NF90_REAL, dim1d, time_id))
      write(yearstr,'(i4)') year
      call check_err(nf90_put_att(nc_id, time_id, "long_name", "time"))
      call check_err(nf90_put_att(nc_id, time_id, "units", "days since "//trim(yearstr)//"-01-01 00:00:00"))
      call check_err(nf90_put_att(nc_id, time_id, "axis", "T"))
  
      write(*,*) "Init params"
      write(*,*) "z_dim_id: ", z_id
      write(*,*) "z2_dim_id: ", z2_dim_id
  
      !Define state variables
      dim_ids = (/z_dim_id, time_dim_id/)
      dim_ids2 = (/z2_dim_id, time_dim_id/)
      dim_ids0 = (/time_dim_id/)
      allocate(parameter_id(size(model%interior_state_variables)))
      allocate(parameter_fick_id(size(model%interior_state_variables)))
      allocate(parameter_sink_id(size(model%interior_state_variables)))
      do ip=1,size(model%interior_state_variables)
          ilast = index(model%interior_state_variables(ip)%path,'/',.true.)
          call check_err(nf90_def_var(nc_id, model%interior_state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id(ip)))  ! was ilast+1:
          call check_err(nf90_def_var(nc_id, 'fick:'//model%interior_state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids2, parameter_fick_id(ip)))
          call check_err(nf90_def_var(nc_id, 'sink:'//model%interior_state_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids2, parameter_sink_id(ip)))
          call check_err(set_attributes(ncid=nc_id, id=parameter_id(ip), units=model%interior_state_variables(ip)%units, &
              long_name=model%interior_state_variables(ip)%long_name, missing_value=model%interior_state_variables(ip)%missing_value))
          call check_err(set_attributes(ncid=nc_id, id=parameter_fick_id(ip), units='mmol/m^2/day', &
              long_name='fick:'//model%interior_state_variables(ip)%long_name,missing_value=model%interior_state_variables(ip)%missing_value))
          call check_err(set_attributes(ncid=nc_id, id=parameter_sink_id(ip), units='mmol/m^2/day', &
              long_name='sink:'//model%interior_state_variables(ip)%long_name,missing_value=model%interior_state_variables(ip)%missing_value))
          call check_err(nf90_put_att(nc_id, parameter_fick_id(ip), "positive", "down"))
          call check_err(nf90_put_att(nc_id, parameter_sink_id(ip), "positive", "down"))
      end do
  
      !Define diagnostic variables
      allocate(parameter_id_diag(size(model%interior_diagnostic_variables)))
      !do ip=1,size(model%interior_diagnostic_variables)
      !    write(*,*) model%interior_diagnostic_variables(ip)%name
      !enddo
      do ip=1,size(model%interior_diagnostic_variables)
          if (model%interior_diagnostic_variables(ip)%save) then
              ilast = index(model%interior_diagnostic_variables(ip)%path,'/',.true.)
              call check_err(nf90_def_var(nc_id, model%interior_diagnostic_variables(ip)%path(ilast+1:), NF90_REAL, dim_ids, parameter_id_diag(ip)))
              call check_err(set_attributes(ncid=nc_id, id=parameter_id_diag(ip), units=model%interior_diagnostic_variables(ip)%units, &
                  long_name=model%interior_diagnostic_variables(ip)%long_name,missing_value=model%interior_diagnostic_variables(ip)%missing_value))
          end if
      end do
  
      !Define forcing variables used in the run
      call check_err(nf90_def_var(nc_id, "T", NF90_REAL, dim_ids, T_id))
      call check_err(nf90_put_att(nc_id, T_id, "long_name", "temperature"))
      call check_err(nf90_put_att(nc_id, T_id, "units", "degC"))
      call check_err(nf90_def_var(nc_id, "S", NF90_REAL, dim_ids, S_id))
      call check_err(nf90_put_att(nc_id, S_id, "long_name", "salinity"))
      call check_err(nf90_def_var(nc_id, "Kz", NF90_REAL, dim_ids2, Kz_id))
      call check_err(nf90_put_att(nc_id, Kz_id, "long_name", "vertical eddy diffusivity"))
      call check_err(nf90_put_att(nc_id, Kz_id, "units", "m2/s"))
      call check_err(nf90_def_var(nc_id, "Kz_sol", NF90_REAL, dim_ids2, Kz_sol_id))
      call check_err(nf90_put_att(nc_id, Kz_sol_id, "long_name", "total vertical diffusivity of a solute"))
      call check_err(nf90_put_att(nc_id, Kz_sol_id, "units", "m2/s"))
      call check_err(nf90_def_var(nc_id, "Kz_par", NF90_REAL, dim_ids2, Kz_par_id))
      call check_err(nf90_put_att(nc_id, Kz_par_id, "long_name", "total vertical diffusivity of a particulate"))
      call check_err(nf90_put_att(nc_id, Kz_par_id, "units", "m2/s"))
      call check_err(nf90_def_var(nc_id, "w_sol", NF90_REAL, dim_ids2, w_sol_id))
      call check_err(nf90_put_att(nc_id, w_sol_id, "long_name", "total advective velocity of a solute"))
      call check_err(nf90_put_att(nc_id, w_sol_id, "units", "m/s"))
      call check_err(nf90_def_var(nc_id, "w_par", NF90_REAL, dim_ids2, w_par_id))
      call check_err(nf90_put_att(nc_id, w_par_id, "long_name", "total advective velocity of a particulate"))
      call check_err(nf90_put_att(nc_id, w_par_id, "units", "m/s"))
      call check_err(nf90_def_var(nc_id, "gas_air_sea", NF90_REAL, dim_ids, gas_air_sea_id))  ! DIC air-sea flux
      call check_err(nf90_put_att(nc_id, gas_air_sea_id, "long_name", "gas_air_sea"))
      call check_err(nf90_put_att(nc_id, gas_air_sea_id, "units", "mmol/m2/d"))
      if (use_swradWm2.eq.1) then
          call check_err(nf90_def_var(nc_id, "swradWm2", NF90_REAL, time_dim_id, swradWm2_id))
          call check_err(nf90_put_att(nc_id, swradWm2_id, "units", "W/m2"))
      end if
      if (use_hice.eq.1) then
          call check_err(nf90_def_var(nc_id, "hice", NF90_REAL, time_dim_id, hice_id))
          call check_err(nf90_put_att(nc_id, hice_id, "units", "m"))
      end if
  
      call check_err(nf90_enddef(nc_id))
  
      end subroutine init_netcdf
  !=======================================================================================================================
  
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine save_netcdf(k_max, julianday, cc, t, s, kz, kzti, wti, &
          model, z, hz, swradWm2, use_swradWm2, hice, use_hice, fick_per_day, sink_per_day, &
          ip_sol, ip_par, i_day, id, idt, output_step, input_step, gas_air_sea, multiyears_physics ) !i_sec_pr) ! i_day here = i_day + 1
  
      !Input variables
      integer, intent(in)                    :: k_max, julianday, use_swradWm2, input_step
      integer, intent(in)                    :: use_hice, ip_sol, ip_par, i_day, id, idt, output_step, multiyears_physics
      real(rk), dimension(:,:), intent(in)   :: cc, t, s, kz, kzti, wti, fick_per_day, sink_per_day, gas_air_sea
      class (type_fabm_model), pointer :: model
      real(rk), dimension(:), intent(in)     :: z, hz, swradWm2, hice
  
      !Local variables
      integer, dimension(1)                  :: start_z, count_z, start_z2, count_z2, start_time, count_time, start_x, count_x
      integer, dimension(2)                  :: start_cc, count_cc, start_flux, count_flux
      real(rk)                               :: temp_matrix(k_max), dum(1), z2(k_max+1), day_part
      !Note: The input arguments to nf90_put_var MUST be vectors, even if the length is 1
      !      Removing the dimension(1) or (1) from dum above triggers a spurious error "not finding nf90_put_var"
      integer                                :: ip, i, i_sec, istep_out, i_day_share
        !integer,parameter         :: timestepkind = selected_int_kind(12)   !how to make int(8)
        !integer(kind=timestepkind):: i_sec  !declared as int(8) it will not work with netcdf functions 
  !    day_part=real(id)/real(idt)
  !    i_sec=((i_day-1)*86400 + int((86400*(id/100))/(idt/100)))/output_step ! as in Horten
  !    istep_out = int((julianday)*86400/input_step)
  
      day_part=real(id)/real(idt)
      i_day_share = 86400/output_step  !a multipier allowing to decrease max int number in i_sec
!      i_sec=int(((i_day-1)*86400 + int(86400*id/idt))/output_step) ! time count for saving  (in array numbers)
      i_sec=int(((i_day-1)*i_day_share + int(i_day_share*id/idt))) ! time count for saving  (in array numbers)      
      istep_out = max(1, int(((julianday-1)*86400 + int(86400*id/idt))/input_step)) ! time count to select data from arrays, i.e. temp, salt
  
   !Define nf90_put_var arguments "start" and "count" for z, z2, time, (cc,t,s) and (fick,kz)
      start_z = 1
      count_z = k_max
      start_z2 = 1
      count_z2 = k_max+1
      start_time = i_sec
      count_time = 1
      start_cc = (/1, i_sec/) ! i_sec
      count_cc = (/k_max, 1/)
      start_flux = (/1, i_sec/) ! i_sec
      count_flux = (/k_max+1, 1/)
      !At first call only, output depth variable mz = -1*z
      if (first) then
          call check_err(nf90_put_var(nc_id, z_id, z, start_z, count_z))
          z2(1:k_max) = z(1:k_max) - 0.5_rk*hz(1:k_max)
          z2(k_max+1) = z(k_max) + 0.5_rk*hz(k_max)
          call check_err(nf90_put_var(nc_id, z2_id, z2, start_z2, count_z2))
          !Note: nc_id, z_id and z2_id are available to the entire module and defined in init_netcdf
          first = .false.
      end if
      dum(1) = (real(i_day)+day_part)
      !For all calls output cc, fick, diagnostics and forcings (t,s,kz)
      if (nc_id.ne.-1) then
          call check_err(nf90_put_var(nc_id, time_id, dum, start_time, count_time))
          do ip=1,size(model%interior_state_variables)
              call check_err(nf90_put_var(nc_id, parameter_id(ip), cc(:,ip), start_cc, count_cc))
              call check_err(nf90_put_var(nc_id, gas_air_sea_id, gas_air_sea(:,ip), start_cc, count_cc))
              call check_err(nf90_put_var(nc_id, parameter_fick_id(ip), fick_per_day(:,ip), start_flux, count_flux))
              call check_err(nf90_put_var(nc_id, parameter_sink_id(ip), sink_per_day(:,ip), start_flux, count_flux))
          end do
          do ip=1,size(model%interior_diagnostic_variables)
              if (model%interior_diagnostic_variables(ip)%save) then
                  temp_matrix = model%get_interior_diagnostic_data(ip)
                  if (maxval(abs(temp_matrix)).lt.1.0E37) then
                      !This is to avoid "NetCDF numeric conversion error" for some NetCDF configurations
                      call check_err(nf90_put_var(nc_id, parameter_id_diag(ip), temp_matrix, start_cc, count_cc))
                  end if
              end if
          end do
          if(multiyears_physics.gt.0) then
            call check_err(nf90_put_var(nc_id, T_id, t(:,i_sec), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, S_id, s(:,i_sec), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, Kz_id, kz(:,i_sec), start_flux, count_flux))
          else
            call check_err(nf90_put_var(nc_id, T_id, t(:,istep_out), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, S_id, s(:,istep_out), start_cc, count_cc))
            call check_err(nf90_put_var(nc_id, Kz_id, kz(:,istep_out), start_flux, count_flux))
          endif

  
          call check_err(nf90_put_var(nc_id, Kz_sol_id, kzti(:,ip_sol), start_flux, count_flux))
          call check_err(nf90_put_var(nc_id, Kz_par_id, kzti(:,ip_par), start_flux, count_flux))
          call check_err(nf90_put_var(nc_id, w_sol_id, wti(:,ip_sol), start_flux, count_flux))
          call check_err(nf90_put_var(nc_id, w_par_id, wti(:,ip_par), start_flux, count_flux))
  
          if (use_swradWm2.eq.1) then
              dum(1) = swradWm2(i_day) !julianday
              call check_err(nf90_put_var(nc_id, swradWm2_id, dum, start_time, count_time))
          end if
  
          if (use_hice.eq.1) then
              dum(1) = hice(i_day) !julianday
              call check_err(nf90_put_var(nc_id, hice_id, dum, start_time, count_time))
          end if
          call check_err(nf90_sync(nc_id))
      end if
  
      end subroutine save_netcdf
  !=======================================================================================================================
  
  
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine close_netcdf()
  
      implicit none
      if (nc_id.ne.-1) then
          call check_err(nf90_close(nc_id))
          deallocate(parameter_id)
          deallocate(parameter_fick_id)
          deallocate(parameter_sink_id)
          deallocate(parameter_id_diag)
          write (*,'(a)') "finished"
      end if
      nc_id = -1
  
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

    integer function find_year_index(yeartofind) result(k)
      integer, intent(in) :: yeartofind
      integer :: j
      k = 0
      do j=1, size(years)
        if (years(j) == yeartofind) then
          k = j; return
        end if
      end do
    end function

    pure logical function is_leap_gregorian(y) result(isleap)
      integer, intent(in) :: y
      isleap = (mod(y,4)==0 .and. (mod(y,100)/=0 .or. mod(y,400)==0))
    end function











      integer function set_attributes(ncid,id,                         &
                                      units,long_name,                 &
                                      valid_min,valid_max,valid_range, &
                                      scale_factor,add_offset,         &
                                      FillValue,missing_value,         &
                                      C_format,FORTRAN_format)
      !
      ! !DESCRIPTION:
      !  This routine is used to set a number of attributes for
      !  variables. The routine makes heavy use of the {\tt optional} keyword.
      !  The list of recognized keywords is very easy to extend. We have
      !  included a sub-set of the COARDS conventions.
      !
      ! !USES:
      !  IMPLICIT NONE
      !
      ! !INPUT PARAMETERS:
      integer, intent(in)                     :: ncid,id
      character(len=*), optional              :: units,long_name
      real, optional                          :: valid_min,valid_max
      real, optional                          :: valid_range(2)
      real, optional                          :: scale_factor,add_offset
      double precision, optional              :: FillValue,missing_value
      character(len=*), optional              :: C_format,FORTRAN_format
      !
      ! !REVISION HISTORY:
      !  Original author(s): Karsten Bolding & Hans Burchard
      !
      ! !LOCAL VARIABLES:
      integer                                 :: iret
      real                                    :: vals(2)
      !
      !
      !-----------------------------------------------------------------------
      !
      if (present(units)) then
          iret = nf90_put_att(ncid,id,'units',trim(units))
      end if
  
      if (present(long_name)) then
          iret = nf90_put_att(ncid,id,'long_name',trim(long_name))
      end if
  
      if (present(C_format)) then
          iret = nf90_put_att(ncid,id,'C_format',trim(C_format))
      end if
  
      if (present(FORTRAN_format)) then
          iret = nf90_put_att(ncid,id,'FORTRAN_format',trim(FORTRAN_format))
      end if
  
      if (present(valid_min)) then
          vals(1) = valid_min
          iret = nf90_put_att(ncid,id,'valid_min',vals(1:1))
      end if
  
      if (present(valid_max)) then
          vals(1) = valid_max
          iret = nf90_put_att(ncid,id,'valid_max',vals(1:1))
      end if
  
      if (present(valid_range)) then
          vals(1) = valid_range(1)
          vals(2) = valid_range(2)
          iret = nf90_put_att(ncid,id,'valid_range',vals(1:2))
      end if
  
      if (present(scale_factor)) then
          vals(1) = scale_factor
          iret = nf90_put_att(ncid,id,'scale_factor',vals(1:1))
      end if
  
      if (present(add_offset)) then
          vals(1) = add_offset
          iret = nf90_put_att(ncid,id,'add_offset',vals(1:1))
      end if
  
      if (present(FillValue)) then
          vals(1) = FillValue
          iret = nf90_put_att(ncid,id,'_FillValue',vals(1:1))
      end if
  
      if (present(missing_value)) then
          vals(1) = missing_value
          iret = nf90_put_att(ncid,id,'missing_value',vals(1:1))
      end if
  
      set_attributes = 0
  
      return
  
      end function set_attributes
  !=======================================================================================================================
  
  
  
  
  
  
  !=======================================================================================================================
      subroutine check_err(status)
  
      integer, intent (in) :: status
  
      if (status .ne. NF90_NOERR) then
          print *, trim(nf90_strerror(status))
          stop
      endif
  
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
  