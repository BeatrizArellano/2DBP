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
        integer   :: i_day, start_year, freq_turb, freq_sed, freq_float, last_day, multiyears_physics  !time related ! ?? freq_sed, freq_turb
        integer   :: diff_method, bioturb_across_SWI  !vertical diffusivity related
        integer   :: h_relax, not_relax_centr  !horizontal transport  (relaxation) switches
        integer   :: use_swradWm2, use_hice ! use input for light, ice, calculate Kz
        integer   :: input_type, port_initial_state, ncoutfile_type !I/O related
        integer   :: bio_model ! basic ecosystem model: 0- for BROM_bio (default) 1- for OxyDep
        real(rk)  :: dt, water_layer_thickness
        real(rk)  :: K_O2s, gargett_a0, gargett_q, mult_Kz, Kz_storm
    
        ! Free timestep input and output
        integer   :: year_index, days_in_yr
        integer   :: ist, i_step, input_step, steps_in_yr, output_step ! ist - steps in the course of the day, i_step - step in the course of the year
    
        character(len=64) :: icfile_name, outfile_name, ncoutfile_name
        character :: hmix_file
    
        !Forcings to be provided to FABM: These must have the TARGET attribute
        real(rk), allocatable, target, dimension(:)   :: hice, aice, swradWm2
        real(rk), allocatable, target, dimension(:)   :: surf_flux , bott_flux, bott_source, Izt, pressure, depth, cell_thickness
        real(rk), allocatable, target, dimension(:,:) :: t, s
        real(rk), allocatable, target, dimension(:,:) :: vv, dVV, cc, cc_out, dcc, dcc_R, wbio, air_sea_flux ! add the description cc - all params, dcc - volumes of solids
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
        integer, allocatable, dimension(:,:)       :: kzCFL, kz_molCFL
        character(len=attribute_length), allocatable, dimension(:)    :: par_name
    
                ! variables for building the grid
        integer                                   :: sel, iday, istep
    !    integer                                   :: k_sed(k_max-k_bbl_sed), k_sed1(k_max+1-k_bbl_sed), k_bbl1(k_bbl_sed-k_wat_bbl)
        real(rk)                                  :: z_wat_bbl, z_bbl_sed, kz_gr !, z1(k_max+1), z_s1(k_max+1), phi1(k_max+1)
        real(rk)                                  :: hz_sed_min, dbl_thickness, kz_mol0
        real(rk)                                  :: a1_bioirr, a2_bioirr
        real(rk)                                  :: kz_bioturb_max, z_const_bioturb, z_decay_bioturb
        real(rk)                                  :: phi_0, phi_inf, z_decay_phi, w_binf, rho_def, wat_con_0, wat_con_inf
    !    real(rk)                                  :: kzCFL(k_bbl_sed-1,steps_in_yr), kz_molCFL(k_max-1,par_max)
    
        integer                                   :: inj_changing !for changing with time injection
        real(rk)                                  :: inj_square   ! square of the layer with injection
        !Constant forcings that can be read as parameters from brom.yaml
        real(rk) :: wind_speed, pco2_atm, mu0_musw, dphidz_SWI, area_col 
     
        ! Injection of something as a function or years
        real(rk)     :: inj_smth(400)          

        real(rk)     :: latitude, Io   !Variables used to calculate surface irradiance from latitude
           ! Environment
        real(rk),target :: decimal_yearday

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
    
    
        !Reading brom.yaml
        call init_common()
    
        !Get grid and numerical solution parameters from from brom.yaml
        dt = get_brom_par("dt")
        freq_turb = get_brom_par("freq_turb")
        freq_sed  = get_brom_par("freq_sed ")
        freq_float  = get_brom_par("freq_float ")
        last_day = get_brom_par("last_day")
        multiyears_physics = get_brom_par("multiyears_physics")
        water_layer_thickness = get_brom_par("water_layer_thickness")
        k_min = get_brom_par("k_min")
        k_storm = get_brom_par("k_storm")
        hz_sed_min = get_brom_par("hz_sed_min")
        hz_sed_min = get_brom_par("hz_sed_min")
        k_points_below_water = get_brom_par("k_points_below_water")
        area_col = get_brom_par("area_col")
        start_year = get_brom_par("start_year")
    
        ! for free length output (assumed to be a day fraction)
        input_step = get_brom_par("input_step")
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
        ncoutfile_name = get_brom_name("ncoutfile_name")
        ncoutfile_type = get_brom_par("ncoutfile_type")
        K_O2s = get_brom_par("K_O2s")
        h_relax =  get_brom_par("h_relax")
        not_relax_centr =  get_brom_par("not_relax_centr")
        ! light connected parameters (if not available in the forcing file)
        latitude = get_brom_par("latitude")
        Io = get_brom_par("Io")                    !W m-2 maximum surface downwelling irradiance at latitudes <= 23.5N,S
    
        ! vertical grid params    
        dbl_thickness = get_brom_par("dbl_thickness")
    
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
            ! Loading water column physics from netCDF forcing file
            !   - opens forcing file and figures out years and days per year
            !   - loads depth dimension (z_w)
            !-----------------------------------------------------------------
            call open_forcing_file(z_w)
            write(*,*) "NetCDF forcing successfully opened (depth axis and metadata)"
            !Note: This uses the netCDF file to set z_w = layer midpoints, dz_w = increments between layer midpoints, hz_w = layer thicknesses
        end if

        !Determine total number of vertical grid points (layers) now that k_wat_bbl is determined
        k_wat_bbl = size(z_w)
        k_max = k_wat_bbl + k_points_below_water

        !Determine number of days in the first year
        year_index = find_year_index(start_year)
        days_in_yr = days_in_year(year_index)
    
        !Allocate full grid variables now that k_max is knownk
        allocate(z(k_max))
        allocate(dz(k_max))
        allocate(hz(k_max))
        allocate(air_sea_flux(k_max,par_max))
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
        allocate(depth(k_max))
        allocate(kzti(k_max+1,par_max))
        allocate(kztCFL(k_max-1,par_max))
        allocate(wCFL(k_max-1,par_max))
        allocate(k_bbl1(k_bbl_sed-k_wat_bbl))
        allocate(z1(k_max+1))
        allocate(z_s1(k_max+1))
        allocate(kzCFL(k_bbl_sed-1,days_in_yr))
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

        ! Depth of interfaces relative to SWI (zero at sediment-water interface)
        z_s1 = z1 - z_bbl_sed
        !------------------------------------------------------------

        !------------------------------------------------------------
        ! Loading initial variables and building the forcing arrays for the full vertical grid
        call load_variable_year('temperature', start_year, t_w)
        call load_variable_year('salinity', start_year, s_w)
        call load_variable_year('Kz', start_year, kz_w)
        if (use_swradWm2 == 1) then
            call load_variable_year_1d('swradWm2', start_year, swradWm2)
        else
            allocate(swradWm2(days_in_yr))
            call build_swrad_year(Io, latitude, days_in_yr, swradWm2)
        end if
        if (use_hice.eq.1) then
            call load_variable_year_1d('hice', start_year, hice)
            call load_variable_year_1d('aice', start_year, aice)
        end if
        if (k_points_below_water>0) then
            ! Construct full-depth annual forcing arrays (T, S, Kz) for the model
            ! Below the water column, repeats the bottom value (constant T, S).
            call build_year_forcing(k_max, k_wat_bbl, k_bbl_sed, &
                                    z1, z_bbl_sed, dbl_thickness, &
                                    t_w, s_w, kz_w, t, s, kz)                               
        end if   


        !------------------------------------------------------------
!!!! Restart checking here 



    
        steps_in_yr = days_in_yr*24*3600/input_step ! determine how much timesteps in the course of the year
        if(multiyears_physics.gt.0) then
            steps_in_yr = last_day*24*3600/input_step ! determine how much timesteps in the course of the year
        endif 
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
            else if (bctype_top(ip).eq.2) then     !Model: bc_top = a1top + a2top*sin(omega*(julianday-a3top))
                bcpar_top(ip,1) = get_brom_par('a1top_' // trim(par_name(ip)))
                bcpar_top(ip,2) = get_brom_par('a2top_' // trim(par_name(ip)))
                bcpar_top(ip,3) = get_brom_par('a3top_' // trim(par_name(ip)))
                write(*,*) "Sinusoidal Dirichlet upper boundary condition for " // trim(par_name(ip))
                write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_top(ip,1), " + ", &
                      bcpar_top(ip,2), "*sin(omega*(julianday -", bcpar_top(ip,3), "))"
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
            else if (bctype_bottom(ip).eq.2) then  !Model: bc_bottom = a1bottom + a2bottom*sin(omega*julianday-a3bottom))
                bcpar_bottom(ip,1) = get_brom_par('a1bottom_' // trim(par_name(ip)))
                bcpar_bottom(ip,2) = get_brom_par('a2bottom_' // trim(par_name(ip)))
                bcpar_bottom(ip,3) = get_brom_par('a3bottom_' // trim(par_name(ip)))
                write(*,*) "Sinusoidal Dirichlet lower boundary condition for " // trim(par_name(ip))
                write(*,'(a, es10.3, a, es10.3, a, es10.3, a)') " = ", bcpar_bottom(ip,1), " + ", &
                bcpar_bottom(ip,2), "*sin(omega*(julianday -", bcpar_bottom(ip,3), "))"
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

    
        !Point FABM to array slices with biogeochemical state.
        do ip=1,par_max
            call model%link_interior_state_data(ip, cc(:,ip))
        end do
    
        !Link temperature and salinity data to FABM (needs to be redone every time julianday is updated below)
        call model%link_interior_data(fabm_standard_variables%temperature, t(:,1))
        call model%link_interior_data(fabm_standard_variables%practical_salinity, s(:,1))
    
        !Link other data needed by FABM
        call model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux, Izt)  !W m-2
        call model%link_interior_data(fabm_standard_variables%pressure, pressure)                              !dbar
        call model%link_interior_data(fabm_standard_variables%depth, depth)                            !dbar
        call model%link_interior_data(fabm_standard_variables%cell_thickness, cell_thickness)
        call model%link_horizontal_data(fabm_standard_variables%wind_speed, wind_speed)
        call model%link_horizontal_data(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air, pco2_atm)
        call model%link_horizontal_data(fabm_standard_variables%latitude, latitude)
        call model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, swradWm2(1))
        call model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, decimal_yearday)
        if (use_hice.eq.1) then
            call model%link_horizontal_data(type_horizontal_standard_variable(name='hice'), hice(1))
            call model%link_horizontal_data(type_horizontal_standard_variable(name='aice'), aice(1))
        endif    
        call model%link_interior_data(volume_of_cell, vv(:,1))
    
    
        !Check FABM is ready
        call model%start()    
    
        !Allow FABM models to use their default initialization (this sets cc)
        do k=1,k_max
            call model%initialize_interior_state(1, k)
        end do
    
        !Read initial values from ascii file if req'd
        if (port_initial_state.eq.1) call porting_initial_state_variables(trim(icfile_name), start_year, &
                                                    i_day, k_max, par_max, par_name, cc, vv)
    
        if (port_initial_state.eq.2) then
           call porting_initial_state_variables(trim(icfile_name), start_year, &
                                               i_day, k_max, par_max, par_name, cc, vv)
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
    
    
        !Initialize output
        call init_netcdf(trim(ncoutfile_name), k_max, model, use_swradWm2, use_hice, start_year)
    
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
    
    
        !Complete hydrophysical forcings
        cc_hmix=0.0_rk

        !Porosity (phi) (assumed constant in time)
        phi = 1.0_rk
        phi_inv = 1.0_rk/phi
        wat_content = 1.0_rk
        phi1 = 1.0_rk
        !Calculation below are for water column horizontal index

        phi(k_sed) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_sed)-z_bbl_sed)/z_decay_phi)
        dphidz_SWI = -1.0_rk * (phi_0-phi_inf) / z_decay_phi
        !water content (wat_content) (assumed constant in time)
        wat_content(k_sed) = wat_con_inf + (wat_con_0-wat_con_inf)*exp(-1.0_rk*(z(k_sed)-z_bbl_sed)/z_decay_phi)
        !Porosity on layer interfaces (phi1)
        phi1(k_sed) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_sed)-0.5_rk*hz(k_sed)-z_bbl_sed)/z_decay_phi)
        phi1(k_max+1) = phi_inf + (phi_0-phi_inf)*exp(-1.0_rk*(z(k_max)+0.5_rk*hz(k_max)-z_bbl_sed)/z_decay_phi)

        !Porosity factors used in diffusivity calculations (pF1, pF2)
        !(assumed constant in time but will vary between solutes vs. solids)
        !These allow us to use a single equation to model diffusivity updates in the water column and sediments, for both solutes and solids:
        ! dC/dt = d/dz(pF2*kzti*d/dz(pF1*C)) where C has units [mass per unit total volume (water+sediments)]
        pF1 = 1.0_rk
        pF2 = 1.0_rk
        pWC = 1.0_rk

        do ip=1,par_max
            if (is_solid(ip).eq.0) then !Factors for solutes
                pF1(k_sed,ip) = 1.0_rk/phi(k_sed)            !Factor to convert [mass per unit total volume] to [mass per unit volume pore water] for solutes in sediments
                pF2(k_sed1,ip) = phi1(k_sed1)                !Porosity-related area restriction factor for fluxes across layer interfaces
                pWC(k_sed,ip) = 1.0_rk/wat_content(k_sed)    !Factor to convert [mass per unit total volume] to [mass per unit volume pore water] for solutes in sediments (for analytical data)
                ! dC/dt = d/dz(kzti*dC/dz)                               in the water column
                ! dC/dt = d/dz(phi*kzti*d/dz(C/phi))                     in the sediments
            end if
            if (is_solid(ip).eq.1) then !Factors for solids
                pF1(k_sed,ip) = 1.0_rk/(1.0_rk - phi(k_sed)) !Factor to convert [mass per unit total volume] to [mass per unit volume solids] for solids in sediments
                pF2(k_sed1,ip) = 1.0_rk - phi1(k_sed1)       !Porosity-related area restriction factor for fluxes across layer interfaces
                pWC(k_sed,ip) = 1.0_rk/(1.0_rk - wat_content(k_sed)) !Factor to convert [mass per unit total volume] to [mass per unit volume solids] for solids in sediments (for analytical data)
                ! dC/dt = d/dz(kzti*dC/dz)                               in the water column
                ! dC/dt = d/dz((1-phi)*kzti*d/dz(C/(1-phi)))             in the sediments
            end if
        end do
        !Tortuosity on layer interfaces (following Boudreau 1996, equatoin 4.120)
        tortuosity(:) = sqrt(1.0_rk - 2.0_rk*log(phi1(:)))

    
        kz_mol = 0.0_rk
        kz_bio = 0.0_rk
        alpha = 0.0_rk
    


        !Vertical diffusivity due to effective molecular diffusivity of solutes in the sediments (kz_mol)
        !(assumed constant in time but in general may vary between solutes)
        do ip=1,par_max
            if (is_solid(ip).eq.0) then
                kz_mol(1:k_max+1,ip) = kz_mol0
                kz_mol(k_sed1,ip) = mu0_musw * kz_mol0/(tortuosity(k_sed1)**2)
    
                !Warning if the diffusive CFL condition is > 0.5
                kz_molCFL(:,ip) = (kz_mol(2:k_max,ip)*dt/freq_turb)/(dz(1:k_max-1)**2)
                if (diff_method.eq.0.and.maxval(kz_molCFL(:,ip)).gt.0.5_rk) then
                    write(*,*) "WARNING!!! CFL condition due to molecular diffusivity exceeds 0.5 for variable", trim(par_name(ip))
                    write(*,*) "z_L, kz_mol, kz_molCFL = "
                    do k=1,k_max-1
                        if (kz_molCFL(k,ip).gt.0.5_rk) write(*,*) z(k)+hz(k)/2, kz_mol(k+1,ip), kz_molCFL(k,ip)
                    end do
                end if
            end if
        end do
        !Vertical diffusivity due to bioturbation (kz_bio)
        !(assumed constant in time and between solutes/solids)
        do k=k_bbl_sed+(2-bioturb_across_SWI),k_max+1  !Note: bioturbation diffusivity is assumed to be non-zero on the SWI
            if (z_s1(k)<z_const_bioturb) then
                kz_bio(k) = kz_bioturb_max
            else
                kz_bio(k) = kz_bioturb_max*exp(-1.0_rk*(z_s1(k)-z_const_bioturb)/z_decay_bioturb)
            end if
        end do
    
        !Bioirrigation rate (alpha)
        !(assumed constant in time and between solutes)
        alpha(k_sed) = a1_bioirr*exp(-1.0_rk*a2_bioirr*(z(k_sed)-z_bbl_sed)) !Schluter et al. (2000), Eqn. 2


        !Background vertical advective velocities of particulates and solutes on layer interfaces in the sediments (w_b, u_b)
        !(these assume steady state compaction and neglect reaction terms)
        w_b = 0.0_rk
        u_b = 0.0_rk
        w_b(k_sed1) = ((1.0_rk - phi_inf)/(1.0_rk - phi1(k_sed1))) * w_binf !Boudreau (1997), Eqn 3.67; Holzbecher (2002) Eqn 3
        u_b(k_sed1) = (phi_inf/phi1(k_sed1)) * w_binf !Boudreau (1997), Eqn 3.68; Holzbecher (2002) Eqn 12

        !!Write to output file
!        open(12,FILE='Hydrophysics.dat')
!        write(12,'(6hiday  ,5hk    ,5hi    ,11hhz[m]      ,11hz[m]       ,11hz_H[m]     ,9ht[degC]  ,10hs[psu]    ,16hKz_H[m2/s]      ,18hKz_mol1_H[m2/s]   , 18hKz_bio_H[m2/s]    , 18hu[m/s]            , 18hstep              ,15halpha[/s]      ,5hphi  ,14htortuosity_H  ,17hw_bH [1E-10 m/s]  ,17hu_bH [1E-10 m/s]  )')
!            iday=1
!            do istep=1,steps_in_yr
!                if (mod(istep,int(86400/input_step)).eq.0) then
!    !            if (mod(istep,24).eq.0) then
!                    iday = iday + 1
!                end if
!                do k=1,k_max
!                    write (12,'(2(1x,i4))',advance='NO') istep, k!                    
!                        write(12,'(1x,i4,f10.4,1x,f10.4,1x,f10.4,2(1x,f8.4),1x,f15.11,1x,f17.13,1x,f17.13,12x,f7.4,12x,i5,1x,f15.11,f7.4,1x,f7.4,12x,f7.4,12x,f7.4)',advance='NO') &
!                            i, hz(k), z(k), (z(k)-0.5_rk*hz(k)),&
!                            t(k,istep), s(k,istep), kz(k,istep), kz_mol(k,1), &
!                            kz_bio(k), u_x(k,istep), istep, &
!                            alpha(k), phi(k), tortuosity(k), &
!                            1.E10*w_b(k), 1.E10*u_b(k)

!                    write(12,*)
!                end do
!            end do
!        close(12)
    
        write(*,*) "Made physics of BBL, sediments"
    
        !Get horizontal relaxation parameters from brom.yaml:
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

        ! read an array for injecting of something changing yearly
            
            inj_changing = get_brom_par("inj_changing")
            if (inj_changing.ne.0) then 
                open(21, file= get_brom_name("inj_smth_variable"))
                inj_smth=0.0_rk
                do iday=1,115
                    read(21, *) i_dummy, inj_smth(iday) ! 
                enddo               
                inj_square = get_brom_par("inj_square")
            endif

        open(8,FILE = 'burying_rate.dat')
        !Set constant forcings
        wind_speed = get_brom_par("wind_speed")    ! 10m wind speed [m s-1]
        pco2_atm   = get_brom_par("pco2_atm")      ! CO2 partical pressure [ppm]
        decimal_yearday=0.0_rk
    
        end subroutine init_brom_transport
    !=======================================================================================================================
    
    
    
    
    
    
    !=======================================================================================================================
        subroutine do_brom_transport()
    
        !Executes the offline vertical transport model BROM-transport
    
        use calculate, only:  calculate_phys, calculate_sed, calculate_sed_eya

        implicit none
    
        integer      :: id, idt, idf                     !time related
        integer      :: surf_flux_with_diff              !1 to include surface fluxes in diffusion update, 0 to include in bgc update
        integer      :: bott_flux_with_diff              !1 to include bottom fluxes in diffusion update, 0 to include in bgc update
        !integer      :: model_w_sed                      !1 to assume porosity effects for solutes and solids, 0 - sumplified approach
        integer      :: constant_w_sed                   !1 to assume constant burial (advection) velocities in the sediments
        integer      :: dynamic_w_sed                    !1 to assume dynamic burial (advection) velocities in the sediments depending on dVV(k_bbl_sed)
        integer      :: show_maxmin, show_kztCFL, show_wCFL, show_nan, show_nan_kztCFL, show_nan_wCFL     !options for runtime output to screen
        integer      :: bc_units_convert, sediments_units_convert !options for conversion of concentrations units in the sediment
        integer      :: julianday, model_year
    
        integer      :: ist, istep, istep_new ! AB

        integer      :: trawling_switch                  ! Switch to run a trawling experiment
        integer      :: k_trawling,k_erosion,i_trawling,start_trawling,k_suspension, closest_k_depth  ! Auxiliary integer variables for bottom trawling experiments
        integer      :: k_inj,i_inj,inj_switch,inj_num,start_inj,stop_inj    !#number of layer and column to inject into, start day, stop day number
        integer      :: start_inj2,stop_inj2,start_inj3,stop_inj3    !#start day, stop day number for several injections
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
        !model_w_sed = get_brom_par("model_w_sed")
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
        bc_units_convert = get_brom_par("bc_units_convert")
        sediments_units_convert = get_brom_par("sediments_units_convert")
        tau_relax = get_brom_par("tau_relax")
        injection_rate_ini = get_brom_par("injection_rate_ini")
        k_inj = get_brom_par("k_injection")
        i_inj = get_brom_par("i_injection")
        inj_switch = get_brom_par("injection_switch")
        start_inj = get_brom_par("start_inj")
        start_inj_part_day = get_brom_par("start_inj_part_day")
        stop_inj = get_brom_par("stop_inj")
        start_inj2 = get_brom_par("start_inj2")
        stop_inj2 = get_brom_par("stop_inj2")
        start_inj3 = get_brom_par("start_inj3")
        stop_inj3 = get_brom_par("stop_inj3")
        trawling_switch = get_brom_par("trawling_switch")
        depth_erosion = get_brom_par("depth_erosion")
        depth_trawling = get_brom_par("depth_trawling")
        thickness_suspension = get_brom_par("thickness_suspension")
        i_trawling = get_brom_par("i_trawling")
        start_inj = get_brom_par("start_inj")
        start_inj_part_day = get_brom_par("start_inj_part_day")
        stop_inj = get_brom_par("stop_inj")
        start_inj2 = get_brom_par("start_inj2")
        stop_inj2 = get_brom_par("stop_inj2")
        start_inj3 = get_brom_par("start_inj3")
        stop_inj3 = get_brom_par("stop_inj3")
        start_trawling = get_brom_par("start_trawling")


        idt = int(86400._rk/dt)                     !number of cycles per day
        kzti = 0.0_rk
        sink=0.0_rk
        wti=0.0_rk
        dVV = 0.0_rk
        model_year = 0

            !convert bottom boundary values from 'mass/pore water ml' for dissolved and 'mass/mass' for solids into 'mass/total volume'
            if (bc_units_convert.eq.1) then
                do ip=1,par_max
                    if (bctype_bottom(ip).eq.1) bc_bottom(ip) = bc_bottom(ip)/pF1(k_max,ip)
                enddo
            end if
    

        !Master time step loop over days
        write(*,*) "Starting time stepping"
        !_____Patch for STORM_______________________________!
        if(k_storm.gt.0) then
            kz(1:(k_storm),:)=Kz_storm
        endif
        !_______BIG Cycle ("i_day"=0,...,last_day-1)________!
        ist=1
        istep=1
      do i_day=0,(last_day-1)
    
            julianday = i_day - int(i_day/days_in_yr)*days_in_yr + 1    !"julianday" (1,2,...,days_in_yr)
            if (julianday==1) model_year = model_year + 1
    
        if (k_points_below_water.gt.0) then
            write (*,'(a, i4, a, i4, 3(a, f10.4))') " model year:", model_year, "; julianday:", julianday, &
                  "; w_sed 0 (cm/yr):", wti(k_bbl_sed,1)*365.*8640000., &
                  "; w_sed 1 (cm/yr):", wti(k_bbl_sed+1,1)*365.*8640000.        
        else
            write (*,'(a, i4, a, i4, a, f9.4)') " model year:", model_year, "; julianday:", julianday,"; w_sed (cm/yr):", wti(k_bbl_sed,1)*365.*8640000.
        endif
    
    !###########################################################################################
        do id=1,idt !in the course of 1 day
            !Note: The numerical approach here is Operator Splitting with tracer transport processes assumed to be
            !numerically more demanding than the biogeochemistry (hence freq_turb, freq_sed >= 1) (Butenschon et al., 2012)
    
    !                write (*,'(a, i4, a, i8)') "-julianday:", julianday,"; id:", id
    
            !Set time-varying Dirichlet boundary conditions for current julianday
            do ip=1,par_max                
                !Sinusoidal variations
                if (bctype_top(ip).eq.2) bc_top(ip) = bcpar_top(ip,1) + &
                    bcpar_top(ip,2)*sin(omega*(julianday-bcpar_top(ip,3)))
                if (bctype_bottom(ip).eq.2) bc_bottom(ip) = bcpar_bottom(ip,1) + &
                    bcpar_bottom(ip,2)*sin(omega*(julianday-bcpar_bottom(ip,3)))
    
                !Variations read from netcdf
                if (bctype_top(ip).eq.3) bc_top(ip) = cc_top(ip,julianday)
                if (bctype_bottom(ip).eq.3) bc_bottom(ip) = cc_bottom(ip,julianday)
    
                !Variations read from ascii file and/or calculated as a function of something
                if (bctype_top(ip).eq.4) then
                    bc_top(ip) = cc_hmix(ip,1,julianday)
                end if
    
                !SO4 in mmol/m3, SO4/Salt from Morris, A.W. and Riley, J.P.(1966) quoted in Dickson et al.(2007)
                if (bctype_top(ip).eq.5) bc_top(ip)=(0.1400_rk/96.062_rk)*(s(1,julianday)/1.80655_rk)*1.e6_rk !.and.ip.eq.id_SO4
                !if (bctype_top(ip).eq.4.and.ip.eq.id_Alk) bc_top(ip)=0.068*s(1,1,julianday)
                !         !Alk in mmol/m3, Alk/Salt from Murray, 2014
            enddo
    
            istep_new = max(1, 1 + nint((dble(i_day-1)*86400.d0 + dble(id-1)*dt) / dble(input_step)))
            if (istep_new.NE.istep) then
              ist = int((id*int(dt))/input_step)
              istep = istep_new

            if(multiyears_physics.gt.0) then
                istep = int((i_day)*(24*3600/input_step)) + ist
            endif
            
          ! Reload changes in t, s, light, ice (needed for FABM/biology)
            call model%link_interior_data(fabm_standard_variables%temperature, t(:,istep))
            call model%link_interior_data(fabm_standard_variables%practical_salinity, s(:,istep))

             decimal_yearday=real(julianday)
          !Calculate Izt = <PAR(z)>_24hr for this dayF

            if (use_swradWm2 == 0) then
                swradWm2(istep) = max(0.0_rk, Io*cos((latitude &
                        -23.5_rk*sin(2.0_rk*pi*(real(julianday,rk)-81.0_rk)/365.0_rk) & !Solar declination in degrees
                        )*pi/180.0_rk))
            endif
            call model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, swradWm2(istep))
    
            if (use_hice.eq.1) then
                call model%link_horizontal_data(type_horizontal_standard_variable(name='hice'), hice(istep))
                call model%link_horizontal_data(type_horizontal_standard_variable(name='aice'), aice(istep))
            endif
          endif
    
            !_______vertical diffusion________!

                call calculate_phys(k_max, par_max, model, cc, kzti, fick, &
                    dcc, bctype_top, bctype_bottom, bc_top, bc_bottom, &
                    surf_flux, bott_flux, bott_source, k_bbl_sed, dz, hz, kz, &
                    kz_mol, kz_bio, istep, julianday, id_O2, K_O2s, dt, freq_turb, &
                    diff_method, cnpar, surf_flux_with_diff,bott_flux_with_diff, &
                    bioturb_across_SWI, pF1, pF2, phi_inv, is_solid, cc0)

    
            !_______bioirrigation_____________!
                if (a1_bioirr.gt.0.0_rk) then
                    dcc = 0.0_rk                    
                    !Oxygen status of sediments set by O2 level just above sediment surface
                    O2stat = cc(k_bbl_sed,id_O2) / (cc(k_bbl_sed,id_O2) + K_O2s)
     !                   O2stat = cc(k_bbl_sed,id_Oxy) / (cc(k_bbl_sed,id_Oxy) + K_O2s)
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
    !!$OMP PARALLEL DO                
                !do k=1,k_max
                !   call model%get_interior_sources (1, k, dcc(k,:))
                !end do
    
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
                    bc_bottom, hz, dz, k_bbl_sed, wbio, w_b, u_b, julianday, &
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
                if (id.eq.1) write (8,'(a, i8,a, i4, a, i4,  a, f6.3,7(a, e9.3))') &
                    "i_day:", i_day, " year:", model_year, "; jday:", julianday, &
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
                                dcc(k,ip) = (cc_hmix(ip,k,julianday)-cc(k,ip))/tau_relax
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
              if (i_day.ge.start_inj.and.i_day.lt.stop_inj) then
                if (inj_switch.eq.1)  then
                    do ip = 1, par_max
                        if (par_name(ip).eq.get_brom_name("inj_var_name")) exit
                        inj_num = ip+1
                    end do
                    cc(k_inj,inj_num)=cc(k_inj,inj_num) &
                        +  dt*injection_rate_ini/(area_col*dz(k_inj))
                else
                  if (i_day.ne.start_inj.or.(real(id)/real(idt)).gt.start_inj_part_day) then
                    inj_switch=0 !!!!!! we do it only once
                    !print *, "injection num", inj_num
                    !"cc(i_inj,k_inj,inj_num)=cc(i_inj,k_inj,inj_num) & !+86400.0_rk*dt
                    !"         +86400.0_rk*dt/freq_float &
                    !"         *injection_rate/(dx(i_inj)*dy*dz(k_inj))
                    !cc(:,k_inj,inj_num)=cc(:,k_inj,inj_num)+86400.0_rk*dt*injection_rate/(dx(i_inj)*dx(i_inj)*dz(k_inj))
    !comment
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
            if (i_day+1.eq.start_trawling.and.id.eq.1) then
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
                    write(*,*) "Time step within day id = ", id
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

            vv=1._rk
    
        ! OUTPUT
        if (mod(int(id*dt),output_step).eq.0) then
    
        fick_per_day = 86400.0_rk * fick
        sink_per_day = 86400.0_rk * sink
    ! here we save DIC (pCO2 in uM) air-sea flux  in all the depth of the array air_sea_flux
        air_sea_flux = 0.0_rk
        do k=1, k_max          
            do ip=1,par_max
              air_sea_flux(k,ip) = 86400.0_rk * surf_flux(ip) !surf_flux(9)
            enddo
        enddo
        if (sediments_units_convert.eq.0) then
            if (ncoutfile_type == 1) then
                call save_netcdf(k_max, max(1,julianday), cc, t, s, kz, kzti, wti, model, z, hz, swradWm2, use_swradWm2, hice, use_hice, &
                fick_per_day, sink_per_day, ip_sol, ip_par, i_day+1, id, idt, output_step, input_step, air_sea_flux, multiyears_physics)
            else
                call save_netcdf(k_max, max(1,julianday), cc, t, s, kz, kzti, wti, model, z, hz, swradWm2, use_swradWm2, hice, use_hice, &
                fick_per_day, sink_per_day, ip_sol, ip_par, julianday, id, idt, output_step, input_step, air_sea_flux, multiyears_physics)
            endif
        endif
    
        endif ! END OF OUTPUT
    
        end do ! end of Subloop over timesteps in the course of one day (idt)
    
    
            if (julianday == 364) then !Note: saving to ascii
                vv=1._rk
                call saving_state_variables(trim(outfile_name), model_year, julianday, &
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
    
        end subroutine do_brom_transport
    !=======================================================================================================================
    
    
    
    
    
    
    
    !=======================================================================================================================
        subroutine clear_brom_transport()
    
        deallocate(z_w)
        deallocate(dz_w)
        deallocate(hz_w)
        deallocate(t_w)
        deallocate(s_w)
        deallocate(kz_w)
        deallocate(z)
        deallocate(dz)
        deallocate(hz)
        deallocate(t)
        deallocate(s)
        deallocate(kz)
        deallocate(air_sea_flux)
        deallocate(cc_hmix)
        deallocate(kz_bio)
        deallocate(kz_mol)
        deallocate(alpha)
        deallocate(phi)
        deallocate(phi1)
        deallocate(phi_inv)
        deallocate(tortuosity)
        deallocate(w_b)
        deallocate(u_b)
        deallocate(wti)
        deallocate(pF1)
        deallocate(pF2)
        deallocate(pWC)
        deallocate(cc)
        deallocate(cc_out)
        deallocate(dcc)
        deallocate(vv)
        deallocate(dVV)
        deallocate(fick)
        deallocate(fick_per_day)
        deallocate(sink)
        deallocate(sink_per_day)
        deallocate(wbio)
        deallocate(surf_flux)
        deallocate(bott_flux)
        deallocate(bc_top)
        deallocate(bc_bottom)
        deallocate(bctype_top)
        deallocate(bctype_bottom)
        deallocate(par_name)
        deallocate(Izt)
        deallocate(pressure)
        deallocate(cell_thickness)
        deallocate(hice)
        deallocate(swradWm2)
        deallocate(aice)
        deallocate(cc_top)
        deallocate(cc_bottom)
        deallocate(kzti)
        deallocate(kztCFL)
        deallocate(wCFL)
        deallocate(k_wat)
        deallocate(k_sed)
        deallocate(k_sed1)
        if (diff_method.gt.0) call clean_tridiagonal()
    
        end subroutine clear_brom_transport
    !=======================================================================================================================
    
    
        end module brom_transport
    