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
!                     Elizaveta Protsenko, Phil Wallhead,
!                     Anfisa Berezina, Matvey Novikov, Beatriz Arellano-Nava
!-----------------------------------------------------------------------
!
! Code was adapted for FABM-1.0.3
! More info https://github.com/fabm-model/fabm/wiki/FABM-1.0
!-----------------------------------------------------------------------
    module brom_transport
                                                                                                                                                                                                                                                                        

        ! AB check fabm subroutines
        use fabm_omp
        use fabm_types  !, only: attribute_length, rk ! check attribute_length
        use io_netcdf
        use io_ascii
        use mtridiagonal, only: init_tridiagonal,clean_tridiagonal
        use ids         !Provides access to variable indices id_O2 etc.
    
    
        implicit none
        private
        public init_brom_transport, do_brom_transport, clear_brom_transport
    
        real(rk), parameter :: pi=3.141592653589793_rk
    
        !FABM model with all data and procedures related to biogeochemistry
        class (type_fabm_omp_model), pointer :: model


        type (type_horizontal_standard_variable), parameter :: id_hice = type_horizontal_standard_variable(name='hice',units='m') ! horizontal - 2D
        type (type_horizontal_standard_variable), parameter :: id_aice = type_horizontal_standard_variable(name='aice',units='-')
        type (type_interior_standard_variable) :: & ! check if used!!
            volume_of_cell = type_interior_standard_variable(name='volume_of_cell')
    
        !Solution parameters to be read from brom.yaml (see brom.yaml for details)
        integer   :: k_min, k_wat_bbl, k_bbl_sed !z-axis related
        integer   :: k_points_below_water, k_max, k_storm !z-axis related
        integer   :: par_max                     !no. BROM variables
        integer   :: diff_method, bioturb_across_SWI  !vertical diffusivity related
        integer   :: h_relax             !horizontal transport  (relaxation) switches
        integer   :: use_swradWm2, use_hice ! use input for light, ice, calculate Kz
        integer   :: input_type, port_initial_state !I/O related
        integer   :: bio_model ! basic ecosystem model: 0- for BROM_bio (default) 1- for OxyDep
        real(rk)  :: water_layer_thickness
        real(rk)  :: K_O2s, gargett_a0, gargett_q, mult_Kz, Kz_storm
    
        ! Time input and output
        real(rk)  :: dt
        integer   :: start_year, first_day, last_day, repeat_forcing_year  !time related ! 
        integer   :: year_index, calendar_year, days_in_yr
        integer, allocatable :: years(:), year_start_idx(:), year_last_idx(:)
        integer, allocatable :: days_in_year(:), nrecs_in_year(:)
        integer   :: freq_turb, freq_sed  !time related ! ?? freq_sed, freq_turb
        integer   :: i_day, sim_day, output_step ! 
    
        character(len=64) :: forcing_filename, icfile_name, outfile_name, output_filename
        character :: hmix_file
    
        !Forcings to be provided to FABM: These must have the TARGET attribute
        real(rk), allocatable, target, dimension(:)   :: hice, aice, swradWm2
        real(rk), allocatable, target, dimension(:)   :: surf_flux, bott_flux, bott_source, Izt, pressure, cell_thickness
        real(rk), allocatable, target, dimension(:,:) :: t, s
        real(rk), allocatable, target, dimension(:,:) :: vv, dVV, cc, cc_out, dcc, dcc_R, wbio ! add the description cc - all params, dcc - volumes of solids
        real(rk), allocatable, target, dimension(:,:) :: wbio_2d
        ! vv -> 2d
    
        !Surface and bottom forcings, used within brom-transport only
        real(rk), allocatable, dimension(:,:)    :: cc_top, cc_bottom
    
        !Horizontal mixing forcings, used within brom-transport only
        real(rk), allocatable, dimension(:,:,:)  :: cc_hmix ! relaxation
    
        !Grid parameters and forcings for water column only
        real(rk), allocatable, dimension(:)        :: z_w, dz_w, hz_w
        real(rk), allocatable, dimension(:,:)      :: t_w, s_w, kz_w
    
        !Grid parameters and forcings for full column including water and sediments
        real(rk), allocatable, dimension(:)        :: z, dz, hz, z1, z_s1
        real(rk), allocatable, dimension(:)        :: bc_top, bc_bottom, kz_bio, alpha, phi, phi1, phi_inv, tortuosity, w_b, u_b, wat_content
        real(rk), allocatable, dimension(:,:)      :: kztCFL, wCFL
        real(rk), allocatable, dimension(:,:)      :: kz, kzti, fick, fick_per_day, sink, sink_per_day, bcpar_top, bcpar_bottom
        real(rk), allocatable, dimension(:,:)      :: kz_mol, pF1, pF2, wti, pWC
        integer, allocatable, dimension(:)         :: is_solid, is_gas, k_wat, k_sed, k_sed1, k_bbl1
        integer, allocatable, dimension(:)         :: bctype_top, bctype_bottom, hmixtype
        integer, allocatable, dimension(:,:)       :: kz_molCFL
        character(len=attribute_length), allocatable, dimension(:)    :: par_name
    
                ! variables for building the grid
        integer                                   :: sel, iday
    !    integer                                   :: k_sed(k_max-k_bbl_sed), k_sed1(k_max+1-k_bbl_sed), k_bbl1(k_bbl_sed-k_wat_bbl)
        real(rk)                                  :: z_wat_bbl, z_bbl_sed, kz_gr !, z1(k_max+1), z_s1(k_max+1), phi1(k_max+1)
        real(rk)                                  :: hz_sed_min, dbl_thickness, kz_mol0
        real(rk)                                  :: a1_bioirr, a2_bioirr
        real(rk)                                  :: kz_bioturb_max, z_const_bioturb, z_decay_bioturb
        real(rk)                                  :: phi_0, phi_inf, z_decay_phi, w_binf, rho_def, wat_con_0, wat_con_inf
    
        integer                                   :: inj_changing !for changing with time injection
        real(rk)                                  :: inj_square   ! square of the layer with injection
        !Constant forcings that can be read as parameters from brom.yaml
        real(rk) :: wind_speed, pco2_atm, mu0_musw, dphidz_SWI, area_col 
     
        ! Injection of something as a function or years
        real(rk)     :: inj_smth(400)          

        real(rk)     :: latitude, Io   !Variables used to calculate surface irradiance from latitude
           ! Environment
        real(rk),target :: doy_frac    !Day of year plus the fraction of the day
        real(rk),target :: day_frac    !Fraction of time during the day

        real(rk), allocatable, dimension(:)        :: rho
    
        !Counters
        integer   :: i, k, ip, ip_sol, ip_par, kj, ij !, i_dummy AB
        real(rk)  :: i_dummy ! AB

        contains
    
    !=======================================================================================================================
    !
    !  name: brom_transport.init_brom_transport
    !  @param
    !  @return
    !
        subroutine init_brom_transport()
    
            !Initialises the offline vertical transport model BROM-transport
        
            use ids         !Provides access to variable indices id_O2 etc
        
            implicit none
            integer :: bc_units_convert
        
            !Reading brom.yaml
            call init_common()
        
            !Get grid and numerical solution parameters from from brom.yaml
            dt = get_brom_par("dt")
            start_year = get_brom_par("start_year")
            first_day = get_brom_par("first_day")   ! First day of start_year to start the simulation (If 1st of January then first_day=1)   
            last_day = get_brom_par("last_day")

            forcing_filename = get_brom_name("forcing_filename")
            repeat_forcing_year = get_brom_par("repeat_forcing_year")

            freq_turb = get_brom_par("freq_turb")
            freq_sed  = get_brom_par("freq_sed ")   

            water_layer_thickness = get_brom_par("water_layer_thickness")
            k_storm = get_brom_par("k_storm")
            hz_sed_min = get_brom_par("hz_sed_min")
            k_points_below_water = get_brom_par("k_points_below_water")
            area_col = get_brom_par("area_col")   
        
            ! for free length output (assumed to be a day fraction)
            output_step = get_brom_par("output_step")
        
            bio_model = get_brom_par("bio_model")
            diff_method = get_brom_par("diff_method")
            bioturb_across_SWI = get_brom_par("bioturb_across_SWI")
            input_type = get_brom_par("input_type")
            use_swradWm2 = get_brom_par("use_swradWm2")
            use_hice = get_brom_par("use_hice")
            port_initial_state = get_brom_par("port_initial_state")
            icfile_name = get_brom_name("icfile_name")
            outfile_name = get_brom_name("outfile_name")
            output_filename = get_brom_name("output_filename")
            K_O2s = get_brom_par("K_O2s")
            h_relax =  get_brom_par("h_relax")
            bc_units_convert = get_brom_par("bc_units_convert")
            ! light connected parameters (if not available in the forcing file)
            latitude = get_brom_par("latitude")
            Io = get_brom_par("Io")                    !W m-2 maximum surface downwelling irradiance at latitudes <= 23.5N,S
        
            ! vertical grid params    
            dbl_thickness = get_brom_par("dbl_thickness")

            !Set constant forcings
            wind_speed = get_brom_par("wind_speed")    ! 10m wind speed [m s-1]
            pco2_atm   = get_brom_par("pco2_atm")      ! CO2 partical pressure [ppm]
        
            !Molecular diffusivity of solutes (single constant value, infinite dilution)
            kz_mol0 = get_brom_par("kz_mol0")
            mu0_musw = get_brom_par("mu0_musw")
        
            !Bioturbation
            kz_bioturb_max = get_brom_par("kz_bioturb_max")
            z_const_bioturb = get_brom_par("z_const_bioturb")
            z_decay_bioturb = get_brom_par("z_decay_bioturb")
        
            !Bioirrigation
            a1_bioirr = get_brom_par("a1_bioirr")
            a2_bioirr = get_brom_par("a2_bioirr")
        
            !Porosity
            phi_0 = get_brom_par("phi_0")
            phi_inf = get_brom_par("phi_inf")
            z_decay_phi = get_brom_par("z_decay_phi")
            wat_con_0 = get_brom_par("wat_con_0")
            wat_con_inf = get_brom_par("wat_con_inf")
        
            !Vertical advection in the sediments
            w_binf = get_brom_par("w_binf")
        
            mult_Kz = get_brom_par("mult_Kz")
            Kz_storm = get_brom_par("Kz_storm")

            k_min = 1
            !Initialize FABM model from fabm.yaml
            model => fabm_create_omp_model()
            par_max = size(model%interior_state_variables)

            !----------------Open forcing data----------------------------------------------------------------------------------------------------------------
            if (input_type.eq.0) then !Input sinusoidal seasonal changes (hypothetical)
                stop 'FATAL (brom-transport): input_type=0 (sinusoidal forcing) not supported at the moment.'
            end if
            if (input_type.eq.1) then !Input physics from ascii
                stop 'FATAL (brom-transport): input_type=1 (ASCII forcing) not supported at the moment.'
            end if
            if (input_type.eq.2) then 
                !-----------------------------------------------------------------
                ! Scanning depth and time dimensions from netCDF forcing file
                !   - opens forcing file and figures out years and days per year
                !   - loads depth dimension (z_w)
                !-----------------------------------------------------------------
                call scan_forcing_dimensions(forcing_filename, years, year_start_idx, year_last_idx, days_in_year, nrecs_in_year, &
                                             start_year, first_day, last_day, repeat_forcing_year, use_hice, use_swradWm2, z_w)
                write(*,*) "NetCDF forcing file successfully opened (depth axis and time dimensions)"
                !Note: This uses the netCDF file to set z_w = depth at layer midpoints and checks that needed forcing variables are present. 
            end if

            !Determine total number of vertical grid points (layers) now that k_wat_bbl is determined
            k_wat_bbl = size(z_w)
            k_max = k_wat_bbl + k_points_below_water

            !Determine number of days in the first year
            year_index = findloc(years, start_year, dim=1)
            days_in_yr = days_in_year(year_index)

            !------------------------------------------------------------
            ! Time check: ensure first_day is valid for start_year
            !------------------------------------------------------------
            if (first_day < 1 .or. first_day > days_in_yr) then
                write(*,*) "WARNING: first_day =", first_day, &
                        " is outside valid range [1,", days_in_yr, "] for start_year =", start_year
                write(*,*) "Resetting first_day = 1 (January 1st)."
                first_day = 1
            end if
        
            !Allocate full grid variables now that k_max is knownk
            allocate(z(k_max))
            allocate(dz(k_max))
            allocate(hz(k_max))
            allocate(cc_hmix(par_max,k_max,days_in_yr))
            allocate(kz_mol(k_max+1,par_max))
            allocate(kz_bio(k_max+1))
            allocate(pF1(k_max,par_max))
            allocate(pF2(k_max+1,par_max))
            allocate(pWC(k_max+1,par_max)) ! water content??
            allocate(alpha(k_max))
            allocate(phi(k_max))
            allocate(wat_content(k_max))
            allocate(phi1(k_max+1))
            allocate(phi_inv(k_max))
            allocate(tortuosity(k_max+1))
            allocate(w_b(k_max+1))
            allocate(u_b(k_max+1))
            allocate(wti(k_max+1,par_max)) ! vertical velocity
            allocate(cc(k_max,par_max))
            allocate(cc_out(k_max,par_max))
            allocate(dcc(k_max,par_max))
            allocate(dcc_R(k_max,par_max))
            allocate(fick(k_max+1,par_max))
            allocate(fick_per_day(k_max+1,par_max))
            allocate(wbio(k_max,par_max))    !sinking vertical velocity (m/s, negative for sinking)
            allocate(wbio_2d(k_max,par_max))    !sinking vertical velocity (m/s, negative for sinking)
            allocate(sink(k_max+1,par_max))  !sinking flux (mmol/m2/s, positive downward)
            allocate(sink_per_day(k_max+1,par_max))
            allocate(vv(k_max,1))
            allocate(dVV(k_max,1))
            allocate(Izt(k_max))
            allocate(pressure(k_max))
            allocate(cell_thickness(k_max))
            allocate(kzti(k_max+1,par_max))
            allocate(kztCFL(k_max-1,par_max))
            allocate(wCFL(k_max-1,par_max))        
            allocate(z1(k_max+1))
            allocate(z_s1(k_max+1))
            allocate(kz_molCFL(k_max-1,par_max))

            call compute_thicknesses(z_w, dz_w, hz_w)

            if (k_points_below_water==0) then
                ! Classical only water-column case (no BBL, no sediments).
                k_max=k_wat_bbl
                z=z_w
                dz=dz_w
                hz=hz_w
                k_bbl_sed=k_wat_bbl !needed for Irradiance calculations
                write(*,*) "Constructed only water column grid (", k_max, " layers)."
            else
                ! Adds BBL and sediment layers
                ! The deepest layer in the water column is adjusted to insert BBL layers,
                ! Sediment layers are added below until reaching the bottom.
                call build_vert_grid(z, dz, hz, z_w, dz_w, hz_w, k_wat_bbl, k_max, k_bbl_sed)
                write(*,*) "Constructed water+BBL+sediment grid (", k_max, " layers; BBL ends at ", k_bbl_sed, ")."
                allocate(k_bbl1(k_bbl_sed-k_wat_bbl))
                allocate(k_wat(k_bbl_sed))
                allocate(k_sed(k_max-k_bbl_sed))
                allocate(k_sed1(k_max+1-k_bbl_sed))
                k_wat = (/(k,k=1,k_bbl_sed)/)       !Index vector for all points in the water column
                k_sed = (/(k,k=k_bbl_sed+1,k_max)/) !Index vector for all points in the sediments
                k_sed1 = (/(k,k=k_bbl_sed+1,k_max+1)/) !Indices of layer interfaces in the sediments (including the SWI)
            endif

            !------------------------------------------------------------
            ! Compute key interface depths and index vectors for water, BBL, and sediment
            if (k_points_below_water.gt.0) then  
                z_wat_bbl = z(k_wat_bbl+1) - 0.5_rk*hz(k_wat_bbl+1)  ! Depth of water–BBL interface (top of BBL)
                z_bbl_sed = z(k_bbl_sed+1) - 0.5_rk*hz(k_bbl_sed+1)  ! Depth of BBL–sediment interface (sediment-water interface
                ! Define index vectors for different domains        
                k_sed  = (/(k,k=k_bbl_sed+1,k_max)/)     ! Indices of layer midpoints in the sediments        
                k_sed1 = (/(k,k=k_bbl_sed+1,k_max+1)/)   ! Indices of layer interfaces in the sediments (including SWI)        
                k_bbl1 = (/(k,k=k_wat_bbl+1,k_bbl_sed)/) ! Indices of layer interfaces in the BBL (including its top)
            else
                z_wat_bbl = z(k_wat_bbl)    ! If no BBL/sediment points exist take interface depths simply as the last water column midpoints
                z_bbl_sed = z(k_bbl_sed)
            endif

            ! Depth of layer interfaces (z1):
            ! for each layer, compute the TOP interface as midpoint – half thickness
            z1(1:k_max) = z(:) - 0.5_rk*hz(:)
            z1(k_max+1) = z(k_max) + 0.5_rk*hz(k_max)
            

            ! Depth of interfaces relative to SWI (zero at sediment-water interface)
            z_s1 = z1 - z_bbl_sed
            !------------------------------------------------------------

            !------------------------------------------------------------
            ! Loading initial variables and building the forcing arrays for the full vertical grid            
            call load_variable_year(forcing_filename, "temperature", start_year, years, year_start_idx, year_last_idx, t_w)
            call load_variable_year(forcing_filename, 'salinity', start_year, years, year_start_idx, year_last_idx, s_w)
            call load_variable_year(forcing_filename, "Kz", start_year, years, year_start_idx, year_last_idx, kz_w)
            if (use_swradWm2 == 1) then
                call load_variable_year_1d(forcing_filename, 'swradWm2', start_year, years, year_start_idx, year_last_idx, swradWm2)
            else
                if (allocated(swradWm2)) deallocate(swradWm2)
                allocate(swradWm2(days_in_yr))
                call build_swrad_year(Io, latitude, days_in_yr, swradWm2)
            end if
            if (use_hice.eq.1) then
                call load_variable_year_1d(forcing_filename, 'hice', start_year, years, year_start_idx, year_last_idx, hice)
                call load_variable_year_1d(forcing_filename, 'aice', start_year, years, year_start_idx, year_last_idx, aice)
            else 
                if (allocated(hice)) deallocate(hice)
                allocate(hice(days_in_yr))
                hice = 0.0_rk
            end if
            write(*,'(A,I6,A)') "Forcing data for year ", start_year, " loaded successfully."
            if (k_points_below_water>0) then
                ! Construct full-depth annual forcing arrays (T, S, Kz) for the model
                ! Below the water column, repeats the bottom value (constant T, S).
                call build_year_forcing(k_max, k_wat_bbl, k_bbl_sed, &
                                        z, z1, z_bbl_sed, dbl_thickness, &
                                        t_w, s_w, kz_w, t, s, kz)                               
            end if   

            write(*,*) "Initialized depth-dependent environment data (T, S and Kz) in BBL and sediments."

            !------------------------------------------------------------
        

            !Allocate biological variables now that par_max is known
            allocate(surf_flux(par_max))     !surface flux (tracer unit * m/s, positive for tracer entering column)
            allocate(bott_flux(par_max))     !bottom flux (tracer unit * m/s, positive for tracer entering column)
            allocate(bott_source(par_max))   !surface flux (tracer unit * m/s, positive for tracer entering column)
            allocate(bc_top(par_max))
            allocate(bc_bottom(par_max))
            allocate(bctype_top(par_max))
            allocate(bctype_bottom(par_max))
            allocate(bcpar_top(par_max,3))
            allocate(bcpar_bottom(par_max,3))
            allocate(par_name(par_max))
            allocate(cc_top(par_max,days_in_yr))
            allocate(cc_bottom(par_max,days_in_yr))
            allocate(is_solid(par_max))
            allocate(is_gas(par_max))
            allocate(hmixtype(par_max))
            allocate(rho(par_max))
        
        
            !Retrieve the parameter names from the model structure
            do ip=1,par_max
                par_name(ip) = model%interior_state_variables(ip)%name
            end do
            !Make the named parameter indices id_O2 etc.
            call get_ids(par_name)
        
        if (id_O2.lt.1) id_O2=id_oxy       ! in BROM we have O2 and in OxyDep we have OXY
        if (id_POML.lt.1) id_POML=id_POM   ! in BROM we have POML and in OxyDep we have POM
        
            !Get boudary condition parameters from brom.yaml:
            !bctype = 0, 1, 2, 3 for no flux (default), Dirichlet constant, Dirichlet sinusoid, and Dirichlet netcdf input respectively
            do ip=1,par_max
                bctype_top(ip) = get_brom_par('bctype_top_' // trim(par_name(ip)),0.0_rk)
                if (bctype_top(ip).eq.1) then
                    bc_top(ip) = get_brom_par('bc_top_' // trim(par_name(ip)))
                    write(*,*) "Constant Dirichlet upper boundary condition for " // trim(par_name(ip))
                    write(*,'(a, es10.3)') " = ", bc_top(ip)
                else if (bctype_top(ip).eq.2) then     !Model: bc_top = a1top + a2top*sin(omega*(day_of_year-a3top))
                    bcpar_top(ip,1) = get_brom_par('a1top_' // trim(par_name(ip)))
                    bcpar_top(ip,2) = get_brom_par('a2top_' // trim(par_name(ip)))
                    bcpar_top(ip,3) = get_brom_par('a3top_' // trim(par_name(ip)))
                    write(*,*) "Sinusoidal Dirichlet upper boundary condition for " // trim(par_name(ip))
                    write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_top(ip,1), " + ", &
                        bcpar_top(ip,2), "*sin(omega*(day_of_year -", bcpar_top(ip,3), "))"
                else if (bctype_top(ip).eq.3) then     !Read from netcdf
                    write(*,*) "NetCDF specified Dirichlet upper boundary condition for " // trim(par_name(ip))
                else if (bctype_top(ip).eq.4) then     !Read from ascii or calc from calintiiy
                    write(*,*) "Upper boundary condition from NODC " // trim(par_name(ip))
                end if
                bctype_bottom(ip) = get_brom_par('bctype_bottom_' // trim(par_name(ip)),0.0_rk)
                if (bctype_bottom(ip).eq.1) then
                    bc_bottom(ip) = get_brom_par('bc_bottom_' // trim(par_name(ip)))
                    write(*,*) "Constant Dirichlet lower boundary condition for " // trim(par_name(ip))
                    write(*,'(a, es10.3)') " = ", bc_bottom(ip)
                else if (bctype_bottom(ip).eq.2) then  !Model: bc_bottom = a1bottom + a2bottom*sin(omega*day_of_year-a3bottom))
                    bcpar_bottom(ip,1) = get_brom_par('a1bottom_' // trim(par_name(ip)))
                    bcpar_bottom(ip,2) = get_brom_par('a2bottom_' // trim(par_name(ip)))
                    bcpar_bottom(ip,3) = get_brom_par('a3bottom_' // trim(par_name(ip)))
                    write(*,*) "Sinusoidal Dirichlet lower boundary condition for " // trim(par_name(ip))
                    write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_bottom(ip,1), " + ", &
                    bcpar_bottom(ip,2), "*sin(omega*(day_of_year -", bcpar_bottom(ip,3), "))"
                else if (bctype_bottom(ip).eq.3) then  !Read from netcdf
                    write(*,*) "NetCDF specified Dirichlet lower boundary condition for " // trim(par_name(ip))
                end if
            end do
        
            write(*,*) "All other boundary conditions use surface and bottom fluxes from FABM"
        
            
            !Initialize tridiagonal matrix if necessary
            if (diff_method.gt.0) then
                call init_tridiagonal(k_max)
                write(*,*) "Initialized tridiagonal matrix"
            end if
        
            !Set model domain
            call model%set_domain(k_max)
        
        
            !Initial volumes of layers:
            vv(1:k_max,1) = 1.0_rk
        
            !Make they (full) pressure variable to pass to FABM
            !This is used by brom_eqconst.F90 to compute equilibrium constants for pH calculations
            !and by brom_carb.F90 to compute equilibrium constants for saturation states (subroutine CARFIN)
            pressure(:) = z(:) + 10.0_rk
            cell_thickness(:) = hz(:)
            doy_frac = 0.0_rk
        
            !Point FABM to array slices with biogeochemical state.
            do ip=1,par_max
                call model%link_interior_state_data(ip, cc(:,ip))
            end do
        
            !Link temperature and salinity data to FABM (needs to be redone every time a new day starts)
            call model%link_interior_data(fabm_standard_variables%temperature, t(:,1))
            call model%link_interior_data(fabm_standard_variables%practical_salinity, s(:,1))
        
            !Link other data needed by FABM
            call model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux, Izt)  !W m-2
            call model%link_interior_data(fabm_standard_variables%pressure, pressure)                              !dbar
            call model%link_interior_data(fabm_standard_variables%depth, z)                                    
            call model%link_interior_data(fabm_standard_variables%cell_thickness, cell_thickness)
            call model%link_horizontal_data(fabm_standard_variables%wind_speed, wind_speed)
            call model%link_horizontal_data(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air, pco2_atm)
            call model%link_horizontal_data(fabm_standard_variables%latitude, latitude)
            call model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, swradWm2(1))
            call model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, doy_frac)
            if (use_hice.eq.1) then
                call model%link_horizontal_data(type_horizontal_standard_variable(name='hice'), hice(1))
                call model%link_horizontal_data(type_horizontal_standard_variable(name='aice'), aice(1))
            endif    
            call model%link_interior_data(volume_of_cell, vv(:,1))
        
        
            !Check FABM is ready
            call model%start()    
        
            !Allow FABM models to use their default initialization (this sets cc)
            call model%initialize_interior_state(1, k_max)
        
            !Read initial values from ascii file if req'd
            if (port_initial_state.eq.1) call porting_initial_state_variables(trim(icfile_name), start_year, &
                                                                            first_day, k_max, par_max, par_name, cc, vv)
        
            if (port_initial_state.eq.2) then
            call porting_initial_state_variables(trim(icfile_name), start_year, &
                                                    first_day, k_max, par_max, par_name, cc, vv)
            endif
        
            !Check biological parameters for non-zero values
            if (bio_model.lt.1) then  !case BROM
                do ip=1,par_max
                    if (ip.eq.id_Phy.or.ip.eq.id_Het.or.ip.eq.id_Baae.or.ip.eq.id_Baan.or.ip.eq.id_Bhae.or.ip.eq.id_Bhan) then                    
                        do k=1,k_max
                            if(cc(k,ip).le.0.0_rk) cc(k,ip)= 1.0E-7 !-11
                        enddo
                    endif
                enddo
            endif
        
            if (bio_model.ge.1) then    !case OxyDep
                do ip=1,par_max
                    if (ip.eq.id_Phy.or.ip.eq.id_Het) then
                        
                        do k=1,k_max
                            if(cc(k,ip).le.0.0_rk) cc(k,ip)= 1.0E-7 !-11
                        enddo
                    
                    endif
                enddo
            endif
        
        

        
            !Establish which variables will be treated as solid phase in the sediments, based on the biological velocity (sinking/floating) from FABM.
        
            wbio = 0.0_rk

            call model%get_vertical_movement(1,k_max, wbio)
            !do k=1,k_wat_bbl
            !    call model%get_vertical_movement(1,k, wbio(k,:))
            !enddo
            wbio = -1.0_rk * wbio !FABM returns NEGATIVE wbio for sinking; sign change here means that wbio is POSITIVE for sinking
            is_solid = 0
            is_gas = 0
            ip_sol = 0
            ip_par = 0
            write(*,*) "The following variables are assumed to join the solid phase in the sediments"
            do ip=1,par_max
                if (wbio(k_wat_bbl,ip).gt.0.0_rk) then
                    is_solid(ip) = 1 !Any variables that SINKS (wbio>0) in the bottom cell of the water column will become "solid" in the sediments
        !            write(*,*) trim(par_name(ip))
                    if (ip_par.eq.0) ip_par = ip !Set ip_par = index of first particulate variable
                elseif (wbio(k_wat_bbl,ip).lt.0.0_rk) then
                    is_gas(ip) = 1 !Any variables that floats (wbio>1)
        !            write(*,*) trim(par_name(ip))    
                else
                    if (ip_sol.eq.0) ip_sol = ip !Set ip_par = index of first solute variable
                end if
            end do
            
            is_solid(26:28) = 1

            do ip=1,par_max
                if (is_gas(ip).eq.1) then
                    write(*,*) "Gaseous variable: ", trim(par_name(ip))   
                endif    
            enddo
        
            write(*,*) "The following variables are particulate matter:"
            !Density of particles
            rho_def = get_brom_par("rho_def")
            rho = 0.0_rk
            do ip=1,par_max
                if (is_solid(ip).eq.1) then
                    rho(ip) = get_brom_par('rho_' // trim(par_name(ip)),rho_def)
                    write(*,*) "Assumed density of ", trim(par_name(ip)), " = ", rho(ip)
                end if
            end do


            !------------------------------------------------------------
            ! Sediment physics: porosity, tortuosity, diffusivities, compaction
            !------------------------------------------------------------

            ! Initialize arrays with safe defaults
            phi         = 1.0_rk     ! Porosity at layer midpoints (fraction of pore space)
            phi_inv     = 1.0_rk/phi ! Inverse porosity (for quick conversion factors)
            wat_content = 1.0_rk     ! Water content of sediment (fraction)
            phi1        = 1.0_rk     ! Porosity at layer interfaces
            pF1         = 1.0_rk     ! Factor to convert concentration units (total volume -> phase volume)
            pF2         = 1.0_rk     ! Porosity-related restriction factor for fluxes across interfaces
            pWC         = 1.0_rk     ! Factor to convert concentration for porewater/solid phase (analytical data)
            tortuosity  = 1.0_rk     ! Tortuosity at interfaces (effective diffusion path length factor)
            kz_mol      = 0.0_rk     ! Effective molecular diffusivity in sediments [m2/s]
            kz_bio      = 0.0_rk     ! Bioturbation diffusivity [m2/s]
            alpha       = 0.0_rk     ! Bioirrigation rate [s^-1]
            w_b         = 0.0_rk     ! Background advective velocity of solids (burial rate) [m/s]
            u_b         = 0.0_rk     ! Background advective velocity of porewater [m/s]


            !---------------- Porosity profiles -------------------------
            phi(k_sed) = phi_inf + (phi_0 - phi_inf) * &
                        exp(-1.0_rk * (z(k_sed) - z_bbl_sed) / z_decay_phi)

            dphidz_SWI = -1.0_rk * (phi_0 - phi_inf) / z_decay_phi
            !water content (wat_content) (assumed constant in time)
            wat_content(k_sed) = wat_con_inf + (wat_con_0 - wat_con_inf) * &
                                exp(-1.0_rk * (z(k_sed) - z_bbl_sed) / z_decay_phi)
            !Porosity on layer interfaces (phi1)
            phi1(k_sed) = phi_inf + (phi_0 - phi_inf) * &
                        exp(-1.0_rk * (z(k_sed) - 0.5_rk*hz(k_sed) - z_bbl_sed) / z_decay_phi)
            
            phi1(k_max+1) = phi_inf + (phi_0 - phi_inf) * &
                            exp(-1.0_rk * (z(k_max) + 0.5_rk*hz(k_max) - z_bbl_sed) / z_decay_phi)

            !---------------- Conversion factors ------------------------
            !Porosity factors used in diffusivity calculations (pF1, pF2)
            !(assumed constant in time but will vary between solutes vs. solids)
            !These allow us to use a single equation to model diffusivity updates in the water column and sediments, for both solutes and solids:
            ! dC/dt = d/dz(pF2*kzti*d/dz(pF1*C)) where C has units [mass per unit total volume (water+sediments)]
            do ip=1,par_max
                if (is_solid(ip) == 0) then
                    ! Solutes
                    pF1(k_sed,ip) = 1.0_rk / phi(k_sed)
                    pF2(k_sed1,ip) = phi1(k_sed1)
                    pWC(k_sed,ip)  = 1.0_rk / wat_content(k_sed)

                else if (is_solid(ip) == 1) then
                    ! Solids
                    pF1(k_sed,ip) = 1.0_rk / (1.0_rk - phi(k_sed))
                    pF2(k_sed1,ip) = 1.0_rk - phi1(k_sed1)
                    pWC(k_sed,ip)  = 1.0_rk / (1.0_rk - wat_content(k_sed))
                end if
            end do

            !---------------- Tortuosity (Boudreau 1996, Eq. 4.120) ------
            tortuosity(:) = sqrt(1.0_rk - 2.0_rk*log(phi1(:)))

            !---------------- Molecular diffusivity ---------------------
            do ip=1,par_max
                if (is_solid(ip) == 0) then
                    kz_mol(1:k_max+1,ip) = kz_mol0
                    kz_mol(k_sed1,ip)    = mu0_musw * kz_mol0 / (tortuosity(k_sed1)**2)

                    ! CFL check for stability
                    kz_molCFL(:,ip) = (kz_mol(2:k_max,ip)*dt/freq_turb) / (dz(1:k_max-1)**2)
                    if (diff_method == 0 .and. maxval(kz_molCFL(:,ip)) > 0.5_rk) then
                        write(*,*) "WARNING (BROM): CFL > 0.5 for solute ", trim(par_name(ip))
                    end if
                end if
            end do

            !---------------- Bioturbation diffusivity ------------------
            do k=k_bbl_sed+(2-bioturb_across_SWI),k_max+1
                if (z_s1(k) < z_const_bioturb) then
                    kz_bio(k) = kz_bioturb_max
                else
                    kz_bio(k) = kz_bioturb_max * exp(-1.0_rk * (z_s1(k)-z_const_bioturb) / z_decay_bioturb)
                end if
            end do

            !---------------- Bioirrigation -----------------------------
            alpha(k_sed) = a1_bioirr * exp(-1.0_rk*a2_bioirr*(z(k_sed)-z_bbl_sed))

            !---------------- Compaction velocities ---------------------
            w_b(k_sed1) = ((1.0_rk - phi_inf) / (1.0_rk - phi1(k_sed1))) * w_binf  !Boudreau (1997), Eqn 3.67; Holzbecher (2002) Eqn 3
            u_b(k_sed1) = (phi_inf / phi1(k_sed1)) * w_binf                        !Boudreau (1997), Eqn 3.68; Holzbecher (2002) Eqn 12

            

            !------------------------------------------------------------
            ! Write depth-dependent physical properties to ASCII file
            !------------------------------------------------------------
            open(unit=12, file='DepthProperties.dat', status='replace')
            ! Header line
            write(12,'(a)') ' k   z[m]        hz[m]       phi        phi1       tortuosity   kz_mol[m2/s]   kz_bio[m2/s]   alpha[1/s]     w_b[m/s]      u_b[m/s]'
            
            do k=1,k_max
            write(12,'(i4,1x,f10.4,1x,f10.4,1x,f10.5,1x,f10.5,1x,f10.5,1x,es12.4,1x,es12.4,1x,es12.4,1x,es12.4,1x,es12.4)') &
                    k, z(k), hz(k), phi(k), phi1(k), tortuosity(k), &
                    kz_mol(k,1), kz_bio(k), alpha(k), w_b(k), u_b(k)
            end do

            close(12)
            
            !------------------------------------------------------------
            !                  Horizontal relaxation
            !------------------------------------------------------------
            
            write(*,*) "Constructed porosity, diffusivity, and bioturbation profiles (BBL & sediments)"


        
            !Get horizontal relaxation parameters from brom.yaml:
            !Complete hydrophysical forcings
            cc_hmix=0.0_rk
            !hmixtype = 0, 1 or 2  for no horizontal relaxation (default), box model mixing respectively
                do ip=1,par_max
                    hmixtype(ip) = get_brom_par('hmix_' // trim(par_name(ip)),0.0_rk)
                    hmixtype(ip)=hmixtype(ip)
        
                    if (hmixtype(ip).eq.1) then
                        write(*,*) "Horizontal relaxation assumed for " // trim(par_name(ip))
                    end if
                    if (hmixtype(ip).eq.2) then
            !            hmix_file = get_brom_name("hmix_filename_" // trim(par_name(ip)(13:))) !niva_oxydep_NUT
                        write(*,*) "Horizontal relaxation (ASCII) assumed for " // trim(par_name(ip))
                    !else if (hmixtype(ip).eq.2) then  ! read relaxation files in two cases, also for top bondary condition  #bctype_top(ip).eq.4.or.
                        open(20, file= get_brom_name("hmix_filename_" // trim(par_name(ip))))!'' // hmix_file
                        write(*,*) ("hmix_filename_" // trim(par_name(ip)))
                        do k=1,k_wat_bbl
                            do i_day=1,days_in_yr
                                read(20, *) i_dummy,i_dummy,cc_hmix(ip,k,i_day) ! data for relaxation(par_max,k_max,days_in_yr)
                            end do
                        end do
                        close(20)
        
                    else if (bctype_top(ip).eq.4.) then
                        open(40, file = get_brom_name("bc_top_filename_" // trim(par_name(ip)))) !to boundary condition file
                        do i_day=1,days_in_yr
                            read(40, *) cc_top(ip,i_day)
                        end do
                        close(40)
                    end if
                end do

            !------------------------------------------------------------


            !!_____Patch for STORM_______________________________!
            !if(k_storm.gt.0) then
            !    kz(1:(k_storm),:)=Kz_storm
            !endif


            !convert bottom boundary values from 'mass/pore water ml' for dissolved and 'mass/mass' for solids into 'mass/total volume'
            if (bc_units_convert.eq.1) then
                do ip=1,par_max
                    if (bctype_bottom(ip).eq.1) bc_bottom(ip) = bc_bottom(ip)/pF1(k_max,ip)
                enddo
            end if  
        

            open(8,FILE = 'burying_rate.dat')
            !Initialize output
            call init_netcdf(trim(output_filename), k_max, z, z1, model, use_hice, start_year)
    
        end subroutine init_brom_transport
    !=======================================================================================================================
    
    
     
    
    !=======================================================================================================================
        subroutine do_brom_transport()
    
        !Executes the offline vertical transport model BROM-transport
    
        use calculate, only:  calculate_phys, calculate_sed, calculate_sed_eya

        implicit none
    
        integer      :: day_of_year, model_year
        integer      :: substep, steps_per_day !time related
        real(rk)     :: time_output, sim_sec, next_output_sec

        integer      :: surf_flux_with_diff              !1 to include surface fluxes in diffusion update, 0 to include in bgc update
        integer      :: bott_flux_with_diff              !1 to include bottom fluxes in diffusion update, 0 to include in bgc update
        integer      :: constant_w_sed                   !1 to assume constant burial (advection) velocities in the sediments
        integer      :: dynamic_w_sed                    !1 to assume dynamic burial (advection) velocities in the sediments depending on dVV(k_bbl_sed)
        integer      :: show_maxmin, show_kztCFL, show_wCFL, show_nan, show_nan_kztCFL, show_nan_wCFL     !options for runtime output to screen
        integer      :: sediments_units_convert !options for conversion of concentrations units in the sediment

        integer      :: trawling_switch                  ! Switch to run a trawling experiment
        integer      :: k_trawling,k_erosion,start_trawling,k_suspension, closest_k_depth  ! Auxiliary integer variables for bottom trawling experiments
        integer      :: k_inj,inj_switch,inj_num,start_inj,stop_inj    !#number of layer and column to inject into, start day, stop day number
        real(rk)     :: cnpar                            !"Implicitness" parameter for GOTM vertical diffusion (set in brom.yaml)
        real(rk)     :: cc0                              !Resilient concentration (same for all variables)
        real(rk)     :: omega                            !angular frequency of sinusoidal forcing = 2*pi/365 rads/day
        real(rk)     :: O2stat                           !oxygen status of sediments (factor modulating the bioirrigation rate)
        real(rk)     :: a1_bioirr                        !to detect whether or not bioirrigation is activated
        real(rk)     :: tau_relax                        ! relaxation time (AB)
        real(rk)     :: fresh_PM_poros                   ! porosity of fresh precipitated PM (i.e. dVV)
        real(rk)     :: w_binf                   ! baseline rate of burying
        real(rk)     :: bu_co                   ! "Burial coeficient" for setting velosity exactly to the SWI proportional to the
                                               !   settling velocity in the water column (0<bu_co<1), 0 - for no setting velosity, (nd)
        
        real(rk)     :: air_sea_flux_CO2        ! Air-sea flux of CO2 (per day?)

        real(rk)     :: injection_rate                   ! injection rate
        real(rk)     :: injection_rate_ini               ! injection rate initial constant or multiplier if injection rate is a function
        real(rk)     :: start_inj_part_day     ! share of day to start injection first day (0=midnight, 0.75= 18h00m etc.)
        real(rk)     :: depth_trawling         ! Penetration depth of the trawling device in the sediments [m].
        real(rk)     :: depth_erosion          ! Depth of the eroded layer during trawling [m].
        real(rk)     :: thickness_suspension   ! Thickness of the water layer where resuspension occurs [m].
        real(rk)     :: vol_trawling           ! Volume of a set of layers where substances are mixed after trawling
        real(rk)     :: mass_trawling          ! mass of state variable (cc) in a set of layers during trawling.
        real(rk)     :: sumdz                  ! Depth within the sediments to interpolate with layers beneath
        real(rk)     :: interp_slope           ! auxiliar variable to interpolate bottom layers of sediments onto upper layers after a trawling event
        
        character(len=attribute_length), allocatable, dimension(:)    :: inj_var_name
    
        omega = 2.0_rk*pi/365.0_rk
    
        !Get parameters for the time-stepping and vertical diffusion / sedimentation
        cnpar = get_brom_par("cnpar")
        dynamic_w_sed = get_brom_par("dynamic_w_sed")
        constant_w_sed = get_brom_par("constant_w_sed")
        fresh_PM_poros = get_brom_par("fresh_PM_poros")
        w_binf = get_brom_par("w_binf")
        bu_co = get_brom_par("bu_co")
        cc0 = get_brom_par("cc0")
        a1_bioirr = get_brom_par("a1_bioirr")
        surf_flux_with_diff = get_brom_par("surf_flux_with_diff")
        bott_flux_with_diff = get_brom_par("bott_flux_with_diff")
        show_maxmin = get_brom_par("show_maxmin")
        show_kztCFL = get_brom_par("show_kztCFL")
        show_wCFL = get_brom_par("show_wCFL")
        show_nan = get_brom_par("show_nan")
        show_nan_kztCFL = get_brom_par("show_nan_kztCFL")
        show_nan_wCFL = get_brom_par("show_nan_wCFL")        
        sediments_units_convert = get_brom_par("sediments_units_convert")
        tau_relax = get_brom_par("tau_relax")
        injection_rate_ini = get_brom_par("injection_rate_ini")
        k_inj = get_brom_par("k_injection")
        inj_switch = get_brom_par("injection_switch")
        start_inj = get_brom_par("start_inj")
        start_inj_part_day = get_brom_par("start_inj_part_day")
        stop_inj = get_brom_par("stop_inj")
        trawling_switch = get_brom_par("trawling_switch")
        depth_erosion = get_brom_par("depth_erosion")
        depth_trawling = get_brom_par("depth_trawling")
        thickness_suspension = get_brom_par("thickness_suspension")
        start_trawling = get_brom_par("start_trawling")
        
        kzti = 0.0_rk
        sink=0.0_rk
        wti=0.0_rk
        dVV = 0.0_rk

        write(*,*) "Starting the simulation"

        !--------------------------------------------------------------------------------------------------
        ! Main Loop over simulation days
        !   - Advances the simulation one day at a time
        !   - Handles calendar rollover and forcing reload per year
        !   - Updates FABM forcing data (T, S, Kz, swrad, ice, etc.) for the current day
        !   - Sets Dirichlet boundary conditions for this day
        !   - Calls the subloop over time steps within the day
        !   - Handles daily diagnostics and optional ASCII output
        !--------------------------------------------------------------------------------------------------

        ! Initialize time variables
        steps_per_day = int(86400._rk/dt)   !number of time-steps per day  
        day_of_year = first_day - 1
        model_year = 1
        year_index   = findloc(years, start_year, dim=1)
        calendar_year = start_year
        sim_sec        = 0.0_rk
        next_output_sec = 0.0_rk !real(output_step, rk)

        do sim_day = first_day,last_day
            ! Advance to next day
            day_of_year = day_of_year + 1

            ! Handle year rollover
            if (day_of_year > days_in_year(year_index)) then
                year_index  = year_index + 1
                if (year_index > size(years)) then
                    if (repeat_forcing_year == 1) then
                        year_index = 1
                    else
                        stop "FATAL: ran out of forcing data"
                    end if
                end if
                day_of_year = 1
                model_year = model_year + 1
                calendar_year = calendar_year + 1
                days_in_yr = days_in_year(year_index)
            end if

            !-------- Reload data for the new year if needed ------------------------------------------------------------------------
            if (model_year > 1 .and. day_of_year == 1 .and. repeat_forcing_year == 0) then
                ! Loads forcing data every new year
                call load_variable_year(forcing_filename, "temperature", calendar_year, years, year_start_idx, year_last_idx, t_w)
                call load_variable_year(forcing_filename, 'salinity', calendar_year, years, year_start_idx, year_last_idx, s_w)
                call load_variable_year(forcing_filename, "Kz", calendar_year, years, year_start_idx, year_last_idx, kz_w)
                if (use_swradWm2 == 1) then
                    call load_variable_year_1d(forcing_filename, 'swradWm2', calendar_year, years, year_start_idx, year_last_idx, swradWm2)
                else
                    if (allocated(swradWm2)) deallocate(swradWm2)
                    allocate(swradWm2(days_in_yr))
                    call build_swrad_year(Io, latitude, days_in_yr, swradWm2)
                end if
                if (use_hice.eq.1) then
                    call load_variable_year_1d(forcing_filename, 'hice', calendar_year, years, year_start_idx, year_last_idx, hice)
                    call load_variable_year_1d(forcing_filename, 'aice', calendar_year, years, year_start_idx, year_last_idx, aice)
                else 
                    if (allocated(hice)) deallocate(hice)
                    allocate(hice(days_in_yr))
                    hice = 0.0_rk
                end if
                write(*,'(A,I6,A)') "Forcing data for year ", calendar_year, " loaded successfully."

                if (k_points_below_water>0) then
                    ! Construct full-depth annual forcing arrays (T, S, Kz) for the model
                    ! Below the water column, repeats the bottom value (constant T, S).
                    call build_year_forcing(k_max, k_wat_bbl, k_bbl_sed, &
                                            z, z1, z_bbl_sed, dbl_thickness, &
                                            t_w, s_w, kz_w, t, s, kz)                          
                end if                
            end if
                
            ! ------------- Update FABM with corresponding day t,s,kz values ---------------------------------------------------
            call model%link_interior_data(fabm_standard_variables%temperature, t(:,day_of_year))
            call model%link_interior_data(fabm_standard_variables%practical_salinity, s(:,day_of_year))
            call model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, swradWm2(day_of_year))
            if (use_hice.eq.1) then
                call model%link_horizontal_data(type_horizontal_standard_variable(name='hice'), hice(day_of_year))
                call model%link_horizontal_data(type_horizontal_standard_variable(name='aice'), aice(day_of_year))
            end if
            !--------------------------------------------------------------------------------------------------------------------    

        

            !-------------------Set time-varying Dirichlet boundary conditions for current day---------------------------------
            do ip=1,par_max                
                !Sinusoidal variations
                if (bctype_top(ip).eq.2) bc_top(ip) = bcpar_top(ip,1) + &
                    bcpar_top(ip,2)*sin(omega*(day_of_year-bcpar_top(ip,3)))
                if (bctype_bottom(ip).eq.2) bc_bottom(ip) = bcpar_bottom(ip,1) + &
                    bcpar_bottom(ip,2)*sin(omega*(day_of_year-bcpar_bottom(ip,3)))
    
                !Variations read from netcdf
                if (bctype_top(ip).eq.3) bc_top(ip) = cc_top(ip,day_of_year)
                if (bctype_bottom(ip).eq.3) bc_bottom(ip) = cc_bottom(ip,day_of_year)
    
                !Variations read from ascii file and/or calculated as a function of something
                if (bctype_top(ip).eq.4) then
                    bc_top(ip) = cc_hmix(ip,1,day_of_year)
                end if
    
                !SO4 in mmol/m3, SO4/Salt from Morris, A.W. and Riley, J.P.(1966) quoted in Dickson et al.(2007)
                if (bctype_top(ip).eq.5) bc_top(ip)=(0.1400_rk/96.062_rk)*(s(1,day_of_year)/1.80655_rk)*1.e6_rk !.and.ip.eq.id_SO4
                !if (bctype_top(ip).eq.4.and.ip.eq.id_Alk) bc_top(ip)=0.068*s(1,1,day_of_year)
                !         !Alk in mmol/m3, Alk/Salt from Murray, 2014
            enddo
            !---------------------------------------------------------------------------------------------------------

            !---------------Console output per day----------------------------------------------------------------------------
            if (k_points_below_water.gt.0) then
                write (*,'(a, i4, a, i4, 3(a, f10.4))') " model year:", model_year, "; dayofyear:", day_of_year, &
                    "; w_sed 0 (cm/yr):", wti(k_bbl_sed,1)*365.*8640000., &
                    "; w_sed 1 (cm/yr):", wti(k_bbl_sed+1,1)*365.*8640000.        
            else
                write (*,'(a, i4, a, i4, a, f9.4)') " model year:", model_year, "; dayofyear:", day_of_year,"; w_sed (cm/yr):", wti(k_bbl_sed,1)*days_in_year(year_index)*8640000.0_rk
            endif
            !----------------------------------------------------------------------------------------------------------------

                    
            !------------------------------------------------------------------
            ! Subloop within the current simulation day (time steps of length dt)
            !   - Updates day fraction and passes it to FABM
            !   - Performs operator splitting:
            !       * Vertical diffusion
            !       * Bioirrigation
            !       * Biogeochemical reactions (FABM sources)
            !       * Particle sinking and sediment transport
            !       * Horizontal relaxation (if enabled)
            !       * Injection and trawling experiments (if enabled)
            !   - Advances tracer concentrations and volume changes
            !   - Checks CFL conditions and NaNs
            !   - Triggers NetCDF output when output interval is reached
            !------------------------------------------------------------------
            do substep=1,steps_per_day 
                !Note: The numerical approach here is Operator Splitting with tracer transport processes assumed to be
                !numerically more demanding than the biogeochemistry (hence freq_turb, freq_sed >= 1) (Butenschon et al., 2012)
                day_frac = real(substep, rk)/real(steps_per_day, rk)
                doy_frac = real(day_of_year-1, rk) + day_frac
                ! doy_frac=real(day_of_year)      
                call model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, doy_frac)           

                !---------------------------------------------------------------------------------
                !--------- PHYSICAL PROCESSES ----------------------------------------------------

                !---------vertical diffusion--------------------------------
                call calculate_phys(k_max, par_max, model, cc, kzti, fick, &
                    dcc, bctype_top, bctype_bottom, bc_top, bc_bottom, &
                    surf_flux, bott_flux, bott_source, k_bbl_sed, dz, hz, kz(:,day_of_year), &
                    kz_mol, kz_bio, id_O2, K_O2s, dt, freq_turb, &
                    diff_method, cnpar, surf_flux_with_diff,bott_flux_with_diff, &
                    bioturb_across_SWI, pF1, pF2, phi_inv, is_solid, cc0)

    
                !_______bioirrigation_____________!
                if (a1_bioirr.gt.0.0_rk) then
                    dcc = 0.0_rk                    
                    !Oxygen status of sediments set by O2 level just above sediment surface
                    O2stat = cc(k_bbl_sed,id_O2) / (cc(k_bbl_sed,id_O2) + K_O2s)
                    do ip=1,par_max
                        if (is_solid(ip).eq.0) then
                            !Calculate tendencies dcc    
                            !Schluter et al. (2000), Meile et al. (2001)
                            !Note use of factor pF1 = 1/phi to convert cc from [mass per unit total volume]
                            !to [mass per unit volume pore water]
                            dcc(k_sed,ip) = O2stat*alpha(k_sed)*phi(k_sed) * &
                                (cc(k_bbl_sed,ip) - pF1(k_sed,ip)*cc(k_sed,ip))
    
                            !Bottom cell of water column receives -1 * sum of all exchange fluxes (conservation of mass)
                            dcc(k_bbl_sed,ip) = -1.0_rk * sum(dcc(k_sed,ip)*hz(k_sed)) / hz(k_bbl_sed)
                            !Bottom cell of water column receives -1 * sum of all exchange fluxes (conservation of mass)
    
                            !Update concentrations
                            !Simple Euler time step for all k in the sediments (index vector k_sed)
                            cc(k_sed,ip) = cc(k_sed,ip) + dt*dcc(k_sed,ip)
                            cc(k_bbl_sed,ip) = cc(k_bbl_sed,ip) + dt*dcc(k_bbl_sed,ip)
                            if (bctype_bottom(ip).gt.0) cc(k_max,ip) = bc_bottom(ip) !Reassert Dirichlet BC if required
                            cc(k_sed,ip) = max(cc0, cc(k_sed,ip)) !Impose resilient concentration
                        end if
                    end do
                end if
    
                !_____water_biogeochemistry_______!
                call model%prepare_inputs()  ! This ensures that all variables that the source routines depend on are computed.
                dcc = 0.0_rk
                call model%get_interior_sources(1, k_max, dcc)
    
                !Add surface and bottom fluxes if treated here
                if (surf_flux_with_diff.eq.0) then
                    surf_flux = 0.0_rk
                    call model%get_surface_sources(surf_flux(:))
                    fick(k_min,:) = surf_flux(:)
                    do ip=1,par_max
                        dcc(k_min,ip) = dcc(k_min,ip) + surf_flux(ip) / hz(k_min)
                    end do
                end if
    
                if (bott_flux_with_diff.eq.0) then
                    bott_flux = 0.0_rk
                    bott_source = 0.0_rk
                    call model%get_bottom_sources(bott_flux(:),bott_source(:))
                    sink(k_max+1,:) = bott_flux(:)
                    do ip=1,par_max
                        dcc(k_max,ip) = dcc(k_max,ip) + bott_flux(ip) / hz(k_max)
                    end do
                end if
    
                call model%finalize_outputs()  ! This ensures that diagnostics unrelated to source routines are computed.
    
                if (dynamic_w_sed.eq.1) dcc_R(:,:) = dcc(:,:) !Record biological reaction terms for use in calculate_sed

                ! Updating concentrations with reaction changes provided by FABM (dcc)
                ! Using explicit Forward Euler method
                cc(:,:) = cc(:,:) + dt*dcc(:,:)
    
                !Reassert Dirichlet BCs
                do ip=1,par_max
                    if (bctype_top(ip).gt.0) then
                        cc(k_min,ip) = bc_top(ip)
                    end if
                    if (bctype_bottom(ip).gt.0) then
                        cc(k_max,ip) = bc_bottom(ip)
                    end if
                enddo
    
                cc(:,:) = max(cc0, cc(:,:)) !Impose resilient concentration
    
           !Compute vertical velocity in water column (sinking/floating) using FABM.
           wbio = 0.0_rk
           call model%get_vertical_movement(1,k_max, wbio)
           wbio = -wbio !FABM returns NEGATIVE wbio for sinking; sign change here means that wbio is POSITIVE for sinking
    
            !_______Particles sinking_________!
            dcc = 0.0_rk
            call calculate_sed_eya(k_max, par_max, model, cc, wti, &
                    sink, dcc, dVV, bctype_top, bctype_bottom, bc_top, &
                    bc_bottom, hz, dz, k_bbl_sed, wbio, w_b, u_b, &
                    dt, freq_sed, dynamic_w_sed, constant_w_sed, is_solid, &
                    rho, phi1, fick, k_sed1, K_O2s, kz_bio, &
                    id_O2, dphidz_SWI, cc0, bott_flux, bott_source, w_binf, bu_co, is_gas)             

    
           !_______Calculate changes of volumes of cells_________!
            dVV(:,1)=0.0           
            k=k_bbl_sed
            do ip=1,par_max !Sum over contributions from each particulate variable
                if (is_solid(ip).eq.1) then
                    !change of Volume of a cell as a function of biology, salt precipitation and sinking
                    dVV(k,1)= dVV(k,1)+(sink(k-1,ip))/rho(ip) ! m3/m2/s = m/s
                end if
            end do
            if (fresh_PM_poros.gt.0.0_rk) then
                dVV(k,1)= dVV(k,1)+(sink(k-1,id_POMR))/rho(id_POMR)*(1.0_rk/(1.0_rk-fresh_PM_poros)-1.0_rk) 
            endif
    
            if (k_points_below_water.gt.0) then
                if (substep.eq.1) write (8,'(a, i8,a, i4, a, i4,  a, f6.3,7(a, e9.3))') &
                    "sim_day:", sim_day, " year:", model_year, "; jday:", day_of_year, &
                    " ; w_sed_bl(cm/yr):", wti(k_bbl_sed+2,1)*365.*8640000., &
                  " ; dVV(k_bbl_sed): ", dVV(k_bbl_sed,1), &
                   " ;dVV_POML: " , sink(k_bbl_sed-1,id_POML)/rho(id_POML), &
                   " ;dVV_POMR: " , sink(k_bbl_sed-1,id_POMR)/rho(id_POML), &
                   " ;dVV_Mn4: " , sink(k_bbl_sed-1,id_Mn4)/rho(id_Mn4), &
                   " ;dVV_Fe3: " , sink(k_bbl_sed-1,id_Fe3)/rho(id_Fe3), & 
                   " ;sink_Mn4: " , sink(k_bbl_sed-1,id_Mn4), &
                   " ;rho_Mn4: " , rho(id_Mn4)
        endif
    
           !________Horizontal relaxation_________!
            if (h_relax.eq.1) then
                dcc = 0.0_rk                
                    do ip=1,par_max
                        if  (hmixtype(ip).ge.1) then
                            do k=1,k_wat_bbl
                            !Calculate tendency dcc (water column only)
                                dcc(k,ip) = (cc_hmix(ip,k,day_of_year)-cc(k,ip))/tau_relax
                            !Update concentration (water column only)
                                cc(k,ip) = cc(k,ip) + dt*dcc(k,ip)
                            end do
                            cc(:,ip) = max(cc0, cc(:,ip)) !Impose resilient concentration
                        end if
                    end do
            endif
    
    
    
            !________Injection____________________!
            ! Source of "substance" in XXXX mmol/sec, should be devided to the volume of the grid cell, i.e. dz(k)*area_col
            if (inj_switch.ne.0)  then
                if (inj_changing.ne.0) then
                    do ip = 1, par_max
                        if (par_name(ip).eq.get_brom_name("inj_var_name")) exit
                        inj_num = ip+1
                    end do
                    cc(k_inj,inj_num)=cc(k_inj,inj_num) &
                        +  dt*inj_smth(model_year+1)/(area_col*dz(k_inj))
                else     
                if (sim_day.ge.start_inj.and.sim_day.lt.stop_inj) then
                    if (inj_switch.eq.1)  then
                        do ip = 1, par_max
                            if (par_name(ip).eq.get_brom_name("inj_var_name")) exit
                            inj_num = ip+1
                        end do
                        cc(k_inj,inj_num)=cc(k_inj,inj_num) &
                            +  dt*injection_rate_ini/(area_col*dz(k_inj))
                    else
                    if (sim_day.ne.start_inj.or.(real(substep)/real(steps_per_day)).gt.start_inj_part_day) then
                        inj_switch=0 !!!!!! we do it only once
                        !print *, "injection num", inj_num
                        !"cc(i_inj,k_inj,inj_num)=cc(i_inj,k_inj,inj_num) & !+86400.0_rk*dt
                        !"         +86400.0_rk*dt/freq_float &
                        !"         *injection_rate/(dx(i_inj)*dy*dz(k_inj))
                        !cc(:,k_inj,inj_num)=cc(:,k_inj,inj_num)+86400.0_rk*dt*injection_rate/(dx(i_inj)*dx(i_inj)*dz(k_inj))
                        do k=2,k_inj
                        injection_rate = dz(k_inj-k-1)*injection_rate_ini/(z(k_inj)-z(1))    !"triangle weight"  distr. in the layer 0-k_inj    
                        cc(k,inj_num)=cc(k,inj_num) &
                        !       dt/freq_float &  ! w/o  for MgOH2
                            +  dt*injection_rate/(area_col*dz(k))
                        enddo
                    end if
                    endif
                endif
                end if
            endif


            !________Bottom Trawling____________________!
            if (trawling_switch.ne.0)  then
                if (sim_day.eq.start_trawling.and.substep.eq.1) then
                    write (*,*)  "Trawling bottom "
                    ! Calculating indices for the layers where erosion and resuspension occurs
                    k_suspension = find_closest_index(z, z(k_bbl_sed)-thickness_suspension)
                    k_erosion = find_closest_index(z, z(k_bbl_sed+1)+depth_erosion)
                    k_trawling = find_closest_index(z, z(k_bbl_sed+1)+depth_trawling-depth_erosion)                    

                    do ip = 1, par_max
                        ! Calculating the mass of each substance eroded from the sediments
                        mass_trawling=0.0_rk
                        do k=k_bbl_sed+1,k_erosion
                            mass_trawling = mass_trawling + (cc(k,ip)*(area_col*hz(k)))
                        enddo
                        
                        ! Calculating volume of layer where the eroded mass (mass_trawling) will be resuspended  
                        vol_trawling = 0.0_rk    
                        do k=k_suspension, k_bbl_sed
                            vol_trawling = vol_trawling + (area_col * hz(k))
                        enddo                    
                        
                        
                        ! Calculating changes of concentratations in the water column due to resuspension          
                        do k=k_suspension, k_bbl_sed
                            cc(k,ip) = cc(k,ip) + (mass_trawling / vol_trawling)    
                        enddo
                

                        ! Set concentration in the upper sediment layer equal to the top of the remaining sediment layer
                        cc(k_bbl_sed + 1, ip) = cc(k_erosion + 1, ip)                 
                        sumdz = 0.0_rk
                        ! Moving the remaining sediment layers upwards
                        do k=k_bbl_sed + 2, k_max - 1                        
                            sumdz = sumdz + dz(k-1)
                            if ((z(k_erosion + 1)+sumdz).lt.z(k_max)) then
                                closest_k_depth = find_closest_index(z, z(k_erosion + 1)+sumdz)
                                ! If the depth is close enough to the depth of a layer beneath, the concentration is just moved upwards
                                if (abs(z(closest_k_depth) - (z(k_erosion + 1) + sumdz)) < hz_sed_min) then                            
                                    cc(k, ip) = cc(closest_k_depth, ip)    
                                else 
                                !g Otherwise interpolate the value
                                    if (z(closest_k_depth).gt.(z(k_erosion + 1)+sumdz)) closest_k_depth = closest_k_depth - 1
                                    ! Interpolating values onto the layers with higher resolution
                                    interp_slope = (cc(closest_k_depth + 1, ip)- cc(closest_k_depth , ip))/ dz(closest_k_depth)                       
                                    cc(k, ip) = cc(k-1, ip) + dz(k-1) * interp_slope 
                                endif
                            else 
                            ! Repeat the values from the deepest sediment layer in the remaining bottom layers
                                cc(k, ip) = cc(k_max-1, ip)
                            endif
                        enddo
                        
                        ! Homogenizing concentrations of solids in the remaining sediment layer and mixing solutes with those in the water column                    
                        if (is_solid(ip).eq.1) then
                            mass_trawling = 0.0_rk
                            vol_trawling  = 0.0_rk
                            do k=k_bbl_sed+1,k_trawling
                                ! Calculate the total mass and volume of the remaining top sediment layer
                                mass_trawling = mass_trawling + (cc(k,ip) * (area_col*hz(k)))
                                vol_trawling = vol_trawling + (area_col * hz(k))
                            enddo

                            do k=k_bbl_sed+1,k_trawling
                                ! Homogenizing concentrations of solids in the top layer
                                cc(k,ip)= mass_trawling / vol_trawling                   
                            enddo
                        else 
                            !g If they are solutes, they are mixed with those dissolved in the resuspension layer
                            mass_trawling = 0.0_rk
                            vol_trawling  = 0.0_rk
                            do k=k_suspension,k_trawling
                                ! Calculate the total mass and volume of the remaining top sediment layer together with the suspension layer
                                mass_trawling = mass_trawling + (cc(k,ip) * (area_col*hz(k)))
                                vol_trawling = vol_trawling + (area_col * hz(k))
                            enddo
                            do k=k_suspension,k_trawling
                                ! Mixing solutes in the interface
                                cc(k,ip) = mass_trawling / vol_trawling
                            enddo
                        endif         
                    enddo                    
                endif
            endif

    
            !________Check for NaNs (stopping if any found)____________________!
            do ip=1,par_max                
                if (any(isnan(cc(1:k_max,ip)))) then
                    write(*,*) "Time step within day id = ", substep
                    write(*,*) "NaN detected in concentration array, ip = ", ip
                    write(*,*) "Variable name = ", model%interior_state_variables(ip)%name
                    if (show_nan.eq.1) write(*,*) cc(1:k_max,ip)
                    if (show_nan_kztCFL.gt.0) then
                        kztCFL(:,ip) = (kzti(2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
                        write(*,*) "maxval(kzti(:,ip)) = ", maxval(kzti(:,ip))
                        write(*,*) "maxval(kztCFL(:,ip)) = ", maxval(kztCFL(:,ip))
                        write(*,*) "fick = ", fick(:,ip)
                        if (show_nan_kztCFL.eq.2.and.maxval(kztCFL(:,ip)).gt.0.5_rk) then
                            write(*,*) "(z_L, kz, kztCFL) where kztCFL>0.5 = "
                            do k=1,k_max-1
                                if (kztCFL(k,ip).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kzti(k+1,ip), kztCFL(k,ip)
                            end do
                        end if
                    end if
                    if (show_nan_wCFL.gt.0) then
                        wCFL(:,ip) = (abs(wti(2:k_max,ip))*dt/freq_sed)/dz(1:k_max-1)
                        write(*,*) "maxval(wti(:,ip)) = ", maxval(wti(:,ip))
                        write(*,*) "maxval(wCFL(:,ip)) = ", maxval(wCFL(:,ip))
                        write(*,*) "wti = ", wti(:,ip)
                        if (show_nan_wCFL.eq.2.and.maxval(wCFL(:,ip)).gt.1.0_rk) then
                            write(*,*) "(z_L, wti, wCFL) where wCFL>1 = "
                            do k=1,k_max-1
                                if (wCFL(k,ip).gt.1.0_rk) write(*,*) z(k)+hz(k)/2, wti(k+1,ip), wCFL(k,ip)
                            end do
                        end if
                    end if
                    stop
                end if
            end do
    
            !--------------------------------------------------------------------------------------------------
            ! ---------------------------------------OUTPUT Netcdf---------------------------------------------
            if (sim_sec >= next_output_sec - 1.0e-6_rk) then
                fick_per_day = 86400.0_rk * fick
                sink_per_day = 86400.0_rk * sink
                ! here we save DIC (pCO2 in uM) air-sea flux
                air_sea_flux_CO2 = 86400.0_rk * surf_flux(9) !surf_flux(9)
                time_output = (real(next_output_sec,rk)/86400.0_rk) + real(first_day,rk) -1
                

                if (sediments_units_convert.eq.1) then
                    write(*,*) "Conversion to mass/pore water not supported at the moment. "
                endif

                call save_netcdf(k_max, cc, t(:,day_of_year), s(:,day_of_year), kz(:,day_of_year), &
                                 model, swradWm2(day_of_year), use_hice, hice(day_of_year), &
                                 fick_per_day, sink_per_day, air_sea_flux_CO2, time_output)

                next_output_sec = next_output_sec + output_step
            end if
            !--------------------------------------------------------------------------------------------------

            ! advance simulation time
            sim_sec = sim_sec + dt    

        end do ! end of Subloop over timesteps in the course of one day
    
    
            if (day_of_year == days_in_yr) then !Note: saving to ascii
                vv=1._rk
                call saving_state_variables(trim(outfile_name), model_year, day_of_year, &
                            k_max, par_max, par_name, z, hz, cc, vv, t, s, kz)
            endif
    
            !Write daily output to screen if required
        
            if (show_maxmin.eq.1) then
                write(*,*) "maxval(cc(:,:),1) = ", maxval(cc(:,:),1) !Good to keep an eye on max/min values for each parameter
                write(*,*) "minval(cc(:,:),1) = ", minval(cc(:,:),1)
            end if
            if (show_kztCFL.gt.0) then
                do ip=1,par_max
                    kztCFL(:,ip) = (kzti(2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
    
                end do
                if (show_kztCFL.eq.2) then
                    write(*,*) "maxval(kzti,2) = ", maxval(kzti,2)
                    write(*,*) "maxval(kztCFL,2) = ", maxval(kztCFL,2)
                else
                    write(*,*) "maxval(kzti) = ", maxval(kzti)
                    write(*,*) "maxval(kztCFL) = ", maxval(kztCFL)
                end if
            end if
            if (show_wCFL.gt.0) then
                do ip=1,par_max
                    wCFL(:,ip) = (abs(wti(2:k_max,ip))*dt/freq_sed)/dz(1:k_max-1)
                end do
                if (show_wCFL.eq.2) then
                    write(*,*) "maxval(wti(:,:),2) = ", maxval(wti(:,:),2)
                    write(*,*) "maxval(wCFL,2) = ", maxval(wCFL,2)
                else
                    write(*,*) "maxval(wti) = ", maxval(wti)
                    write(*,*) "maxval(wCFL) = ", maxval(wCFL)
                end if
            end if

    
        end do   ! end of BIG cycle

        call close_netcdf()
        close(8)
        end subroutine do_brom_transport
    !=======================================================================================================================
    
    
    !=======================================================================================================================
        subroutine clear_brom_transport()

            ! Grid arrays (water column, sediments, full column)
            if (allocated(z_w))        deallocate(z_w, dz_w, hz_w)
            if (allocated(t_w))        deallocate(t_w, s_w, kz_w)
            if (allocated(z))          deallocate(z, dz, hz, z1, z_s1)
        
            ! Forcing
            if (allocated(swradWm2))   deallocate(swradWm2)
            if (allocated(hice))       deallocate(hice)
            if (allocated(aice))       deallocate(aice)
        
            ! Physics (diffusivities, porosity, compaction, etc.)
            if (allocated(kz))         deallocate(kz)
            if (allocated(kzti))       deallocate(kzti)
            if (allocated(kztCFL))     deallocate(kztCFL)
            if (allocated(wCFL))       deallocate(wCFL)
            if (allocated(kz_bio))     deallocate(kz_bio)
            if (allocated(kz_mol))     deallocate(kz_mol)
            if (allocated(kz_molCFL))  deallocate(kz_molCFL)
            if (allocated(alpha))      deallocate(alpha)
            if (allocated(phi))        deallocate(phi)
            if (allocated(phi1))       deallocate(phi1)
            if (allocated(phi_inv))    deallocate(phi_inv)
            if (allocated(tortuosity)) deallocate(tortuosity)
            if (allocated(w_b))        deallocate(w_b)
            if (allocated(u_b))        deallocate(u_b)
            if (allocated(wti))        deallocate(wti)
            if (allocated(Izt))        deallocate(Izt)
            if (allocated(pressure))   deallocate(pressure)
            if (allocated(cell_thickness)) deallocate(cell_thickness)
        
            ! Tracers and fluxes
            if (allocated(cc))         deallocate(cc, cc_out, dcc, dcc_R)
            if (allocated(fick))       deallocate(fick, fick_per_day)
            if (allocated(sink))       deallocate(sink, sink_per_day)
            if (allocated(wbio))       deallocate(wbio)
            if (allocated(wbio_2d))    deallocate(wbio_2d)
            if (allocated(vv))         deallocate(vv, dVV)
        
            ! Biogeochem conversion factors
            if (allocated(pF1))        deallocate(pF1)
            if (allocated(pF2))        deallocate(pF2)
            if (allocated(pWC))        deallocate(pWC)
        
            ! Boundary conditions
            if (allocated(bc_top))     deallocate(bc_top, bc_bottom)
            if (allocated(bctype_top)) deallocate(bctype_top, bctype_bottom)
            if (allocated(bcpar_top))  deallocate(bcpar_top, bcpar_bottom)
            if (allocated(cc_top))     deallocate(cc_top, cc_bottom)
        
            ! Forcings at surface/bottom
            if (allocated(surf_flux))  deallocate(surf_flux, bott_flux, bott_source)
        
            ! Misc arrays
            if (allocated(cc_hmix))    deallocate(cc_hmix)
            if (allocated(par_name))   deallocate(par_name)
            if (allocated(is_solid))   deallocate(is_solid, is_gas, hmixtype)
            if (allocated(rho))        deallocate(rho)
        
            ! Index arrays
            if (allocated(k_wat))      deallocate(k_wat)
            if (allocated(k_sed))      deallocate(k_sed)
            if (allocated(k_sed1))     deallocate(k_sed1)
            if (allocated(k_bbl1))     deallocate(k_bbl1)
        
            ! Clean tridiagonal solver if needed
            if (diff_method.gt.0) call clean_tridiagonal()
        
        end subroutine clear_brom_transport
        
    !=======================================================================================================================
    
    
end module brom_transport
    