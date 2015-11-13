! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_network_init()

        use network

        call network_init()

      end subroutine ca_network_init


! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_extern_init(name,namlen)

        ! initialize the external runtime parameters in
        ! extern_probin_module

        integer :: namlen
        integer :: name(namlen)

        call runtime_init(name,namlen)

      end subroutine ca_extern_init

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_num_spec(nspec_out)

        use network, only : nspec

        implicit none 

        integer, intent(out) :: nspec_out

        nspec_out = nspec

      end subroutine get_num_spec

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_num_aux(naux_out)

        use network, only : naux

        implicit none 

        integer, intent(out) :: naux_out

        naux_out = naux

      end subroutine get_num_aux

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_spec_names(spec_names,ispec,len)

        use network, only : nspec, short_spec_names

        implicit none 

        integer, intent(in   ) :: ispec
        integer, intent(inout) :: len
        integer, intent(inout) :: spec_names(len)

        ! Local variables
        integer   :: i

        len = len_trim(short_spec_names(ispec+1))

        do i = 1,len
           spec_names(i) = ichar(short_spec_names(ispec+1)(i:i))
        end do

      end subroutine get_spec_names

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_spec_az(ispec,A,Z)

        use network, only : nspec, aion, zion

        implicit none 

        integer         , intent(in   ) :: ispec
        double precision, intent(inout) :: A, Z

        ! C++ is 0-based indexing, so increment
        A = aion(ispec+1)
        Z = zion(ispec+1)

      end subroutine get_spec_az

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_aux_names(aux_names,iaux,len)

        use network, only : naux, short_aux_names

        implicit none 

        integer, intent(in   ) :: iaux
        integer, intent(inout) :: len
        integer, intent(inout) :: aux_names(len)

        ! Local variables
        integer   :: i

        len = len_trim(short_aux_names(iaux+1))

        do i = 1,len
           aux_names(i) = ichar(short_aux_names(iaux+1)(i:i))
        end do

      end subroutine get_aux_names

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_method_params(nGrowHyp)

        ! Passing data from f90 back to C++

        use meth_params_module

        implicit none 

        integer, intent(out) :: ngrowHyp

        nGrowHyp = NHYP

      end subroutine get_method_params

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine allocate_outflow_data(np,nc)

        use meth_params_module, only: outflow_data_old, outflow_data_new, outflow_data_allocated

        implicit none

        integer, intent(in) :: np,nc

        if (.not. outflow_data_allocated) then
           allocate(outflow_data_old(nc,np))
           allocate(outflow_data_new(nc,np))
        end if

        outflow_data_allocated = .true.

      end subroutine allocate_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::
      subroutine set_old_outflow_data(radial,time,np,nc)

        ! Passing data from C++ to f90

        use meth_params_module, only: outflow_data_old, outflow_data_old_time

        implicit none

        double precision, intent(in) :: radial(nc,np)
        double precision, intent(in) :: time
        integer         , intent(in) :: np,nc

        ! Do this so the routine has the right size 
        deallocate(outflow_data_old)
          allocate(outflow_data_old(nc,np))

        outflow_data_old(1:nc,1:np) = radial(1:nc,1:np)

        outflow_data_old_time = time

      end subroutine set_old_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::
      subroutine set_new_outflow_data(radial,time,np,nc)

        ! Passing data from C++ to f90

        use meth_params_module, only: outflow_data_new, outflow_data_new_time

        implicit none

        double precision, intent(in) :: radial(nc,np)
        double precision, intent(in) :: time
        integer         , intent(in) :: np,nc

        ! Do this so the routine has the right size 
        deallocate(outflow_data_new)
          allocate(outflow_data_new(nc,np))

        outflow_data_new(1:nc,1:np) = radial(1:nc,1:np)

        outflow_data_new_time = time

      end subroutine set_new_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::
      subroutine swap_outflow_data()

        use meth_params_module, only: outflow_data_new, outflow_data_new_time, &
                                      outflow_data_old, outflow_data_old_time

        implicit none

        integer                       :: np,nc

        nc = size(outflow_data_new,dim=1)
        np = size(outflow_data_new,dim=2)

        if (size(outflow_data_old,dim=2) .ne. size(outflow_data_new,dim=2)) then
           ! Do this so the routine has the right size
           deallocate(outflow_data_old)
             allocate(outflow_data_old(nc,np))
        end if

        if (size(outflow_data_old,dim=2) .ne. size(outflow_data_new,dim=2)) then
           print *,'size of old and new dont match in swap_outflow_data '
            call bl_error("Error:: Castro_nd.f90 :: swap_outflow_data")
        end if

        outflow_data_old(1:nc,1:np) = outflow_data_new(1:nc,1:np)

        if (outflow_data_new_time .ge. 0.d0) &
           outflow_data_old_time = outflow_data_new_time
        outflow_data_new_time = -1.d0

      end subroutine swap_outflow_data

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_method_params(dm,Density,Xmom,Eden,Eint,Temp, &
                                   FirstAdv,FirstSpec,FirstAux,numadv, &
                                   difmag_in, small_dens_in, small_temp_in, small_pres_in, small_ener_in, &
                                   allow_negative_energy_in, &
                                   ppm_type_in,ppm_reference_in, &
                                   ppm_trace_grav_in, ppm_trace_rot_in, ppm_temp_fix_in, &
                                   ppm_tau_in_tracing_in, ppm_predict_gammae_in, &
                                   ppm_reference_edge_limit_in, &
                                   ppm_flatten_before_integrals_in, &
                                   ppm_reference_eigenvectors_in, &
                                   hybrid_riemann_in, use_colglaz_in, riemann_solver_in, use_flattening_in, &
                                   transverse_use_eos_in, transverse_reset_density_in, transverse_reset_rhoe_in, &
                                   cg_maxiter_in, cg_tol_in, &
                                   use_pslope_in, &
                                   do_grav_in, grav_source_type_in, &
                                   gravity_type_in, gravity_type_len, &
                                   do_sponge_in,normalize_species_in,fix_mass_flux_in,use_sgs, &
                                   burning_timestep_factor_in, &
                                   dual_energy_eta1_in,  dual_energy_eta2_in, dual_energy_eta3_in, dual_energy_update_E_from_E_in, &
                                   do_rotation_in, rot_source_type_in, rot_axis_in, &
                                   rot_period_in, rot_period_dot_in, &
                                   diffuse_cutoff_density_in, &
                                   const_grav_in, deterministic_in, do_acc_in)
!                                  phys_bc_lo,phys_bc_hi

        ! Passing data from C++ into f90

        use meth_params_module
        use network, only : nspec, naux
        use eos_module
        use parallel

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux
        integer, intent(in) :: numadv
        integer, intent(in) :: allow_negative_energy_in, ppm_type_in
        integer, intent(in) :: ppm_reference_in, ppm_trace_grav_in, ppm_trace_rot_in, ppm_temp_fix_in
        integer, intent(in) :: ppm_tau_in_tracing_in, ppm_predict_gammae_in
        integer, intent(in) :: ppm_reference_edge_limit_in
        integer, intent(in) :: ppm_flatten_before_integrals_in
        integer, intent(in) :: ppm_reference_eigenvectors_in
        integer, intent(in) :: hybrid_riemann_in, use_colglaz_in, riemann_solver_in, use_flattening_in
        integer, intent(in) :: transverse_use_eos_in, transverse_reset_density_in, transverse_reset_rhoe_in
        integer, intent(in) :: dual_energy_update_E_from_e_in
        double precision, intent(in) :: dual_energy_eta1_in, dual_energy_eta2_in, dual_energy_eta3_in
        integer, intent(in) :: use_pslope_in
        integer, intent(in) :: do_grav_in, grav_source_type_in, gravity_type_len
        integer, intent(in) :: gravity_type_in(gravity_type_len)
        integer, intent(in) :: cg_maxiter_in
        double precision, intent(in) :: cg_tol_in
        integer, intent(in) :: do_sponge_in
        double precision, intent(in) :: difmag_in
        double precision, intent(in) :: small_dens_in, small_temp_in, small_pres_in, small_ener_in
        integer, intent(in) :: normalize_species_in
        integer, intent(in) :: fix_mass_flux_in
        integer, intent(in) :: use_sgs
        double precision, intent(in) :: burning_timestep_factor_in
        double precision, intent(in) :: rot_period_in, rot_period_dot_in, const_grav_in, diffuse_cutoff_density_in
        integer, intent(in) :: do_rotation_in, rot_source_type_in, rot_axis_in
        integer, intent(in) :: deterministic_in, do_acc_in
        integer :: iadv, ispec

        integer             :: QLAST

        integer :: i
        
        call parallel_initialize()

        iorder = 2 

        difmag = difmag_in
        
        !---------------------------------------------------------------------
        ! conserved state components
        !---------------------------------------------------------------------

        ! NTHERM: number of thermodynamic variables (rho, 3 momenta, rho*e, rho*E, T)
        ! NVAR  : number of total variables in initial system
        NTHERM = 7
        if (use_sgs .eq. 1) NTHERM = NTHERM + 1
        NVAR = NTHERM + nspec + naux + numadv

        nadv = numadv

        ! We use these to index into the state "U"
        URHO  = Density   + 1
        UMX   = Xmom      + 1
        UMY   = Xmom      + 2
        UMZ   = Xmom      + 3
        UEDEN = Eden      + 1
        UEINT = Eint      + 1
        if (use_sgs .eq. 1) then
           UESGS = UEINT     + 1
        else
           UESGS = -1
        endif
        UTEMP = Temp      + 1

        if (numadv .ge. 1) then
          UFA   = FirstAdv  + 1
        else 
          UFA = 1
        end if

        UFS   = FirstSpec + 1

        if (naux .ge. 1) then
          UFX = FirstAux  + 1
        else 
          UFX = 1
        end if


        !---------------------------------------------------------------------
        ! primitive state components
        !---------------------------------------------------------------------

        ! QTHERM: number of primitive variables: rho, game, p, (rho e), T
        !         + 3 velocity components + 1 SGS components (if defined)
        ! QVAR  : number of total variables in primitive form
      
        QTHERM = NTHERM + 1  ! here the +1 is for QGAME always defined in primitive mode
                             ! the SGS component is accounted for already in NTHERM

        QVAR = QTHERM + nspec + naux + numadv

        ! We use these to index into the state "Q"
        QRHO  = 1

        QU    = 2
        QV    = 3
        QW    = 4

        ! we'll carry this around as an potential alternate to (rho e)
        QGAME   = 5
        QLAST   = QGAME

        QPRES   = QLAST + 1
        QREINT  = QLAST + 2

        if (use_sgs .eq. 1) then
           QESGS   = QLAST + 3
        else
           QESGS = -1
        endif

        QTEMP   = QTHERM 

        if (numadv .ge. 1) then
          QFA = QTHERM + 1
          QFS = QFA + numadv
        else 
          QFA = 1   ! density
          QFS = QTHERM + 1
        end if
        if (naux .ge. 1) then
          QFX = QFS + nspec
        else 
          QFX = 1
        end if

        ! easy indexing for the passively advected quantities.  This
        ! lets us loop over all four groups (SGS, advected, species, aux)
        ! in a single loop.
        allocate(qpass_map(QVAR))
        allocate(upass_map(NVAR))

        ! Transverse velocities

        if (dm .eq. 1) then
           upass_map(1) = UMY
           qpass_map(1) = QV

           upass_map(2) = UMZ
           qpass_map(2) = QW

           npassive = 2
        else if (dm .eq. 2) then
           upass_map(1) = UMZ
           qpass_map(1) = QW

           npassive = 1
        else
           npassive = 0
        endif

        if (QESGS > -1) then
           upass_map(npassive + 1) = UESGS
           qpass_map(npassive + 1) = QESGS
           npassive = npassive + 1
        endif
        do iadv = 1, nadv
           upass_map(npassive + iadv) = UFA + iadv - 1
           qpass_map(npassive + iadv) = QFA + iadv - 1
        enddo
        npassive = npassive + nadv
        if(QFS > -1) then
           do ispec = 1, nspec+naux
              upass_map(npassive + ispec) = UFS + ispec - 1
              qpass_map(npassive + ispec) = QFS + ispec - 1
           enddo
           npassive = npassive + nspec + naux
        endif
        

         
        !---------------------------------------------------------------------
        ! other initializations
        !---------------------------------------------------------------------

        if (small_dens_in > 0.d0) then
           small_dens = small_dens_in
        else
           call bl_warning("Warning:: small_dens has not been set, defaulting to 1.d-200.")           
           small_dens = 1.d-200
        endif

        if (small_temp_in > 0.d0) then
           small_temp = small_temp_in
        else
           call bl_warning("Warning:: small_temp has not been set, defaulting to 1.d-200.")
           small_temp = 1.d-200
        endif

        if (small_pres_in > 0.d0) then
           small_pres = small_pres_in
        else
           small_pres = 1.d-200
        endif

        if (small_ener_in > 0.d0) then
           small_ener = small_ener_in
        else
           small_ener = 1.d-200
        endif

        call eos_init(small_dens=small_dens, small_temp=small_temp)

        ! The EOS might have modified our choices because of its
        ! internal limitations, so let's get small_dens and small_temp
        ! again just to make sure we're consistent with the EOS.
        
        call eos_get_small_dens(small_dens)
        call eos_get_small_temp(small_temp)
        
        allow_negative_energy        = allow_negative_energy_in
        ppm_type                     = ppm_type_in
        ppm_reference                = ppm_reference_in
        ppm_trace_grav               = ppm_trace_grav_in
        ppm_temp_fix                 = ppm_temp_fix_in
        ppm_tau_in_tracing           = ppm_tau_in_tracing_in
        ppm_predict_gammae           = ppm_predict_gammae_in
        ppm_reference_edge_limit     = ppm_reference_edge_limit_in
        ppm_flatten_before_integrals = ppm_flatten_before_integrals_in

        ppm_reference_eigenvectors   = ppm_reference_eigenvectors_in
        hybrid_riemann               = hybrid_riemann_in
        use_colglaz                  = use_colglaz_in
        riemann_solver               = riemann_solver_in
        use_flattening               = use_flattening_in
        transverse_use_eos           = transverse_use_eos_in
        transverse_reset_density     = transverse_reset_density_in
        transverse_reset_rhoe        = transverse_reset_rhoe_in

        cg_tol                       = cg_tol_in
        cg_maxiter                   = cg_maxiter_in
        use_pslope                   = use_pslope_in
        do_grav                      = do_grav_in
        grav_source_type             = grav_source_type_in
        do_sponge                    = do_sponge_in
        normalize_species            = normalize_species_in
        fix_mass_flux                = fix_mass_flux_in
        burning_timestep_factor      = burning_timestep_factor_in
        do_rotation                  = do_rotation_in
        rot_period                   = rot_period_in
        rot_period_dot               = rot_period_dot_in
        rot_source_type              = rot_source_type_in
        rot_axis                     = rot_axis_in
        diffuse_cutoff_density       = diffuse_cutoff_density_in
        const_grav                   = const_grav_in
        deterministic                = deterministic_in .ne. 0
        do_acc                       = do_acc_in

        allocate(character(len=gravity_type_len) :: gravity_type)

        do i = 1, gravity_type_len
           gravity_type(i:i) = char(gravity_type_in(i))
        enddo
        
        dual_energy_eta1             = dual_energy_eta1_in
        dual_energy_eta2             = dual_energy_eta2_in
        dual_energy_eta3             = dual_energy_eta3_in
        dual_energy_update_E_from_e  = dual_energy_update_E_from_e_in .ne. 0

      end subroutine set_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_problem_params(dm,physbc_lo_in,physbc_hi_in,&
                                    Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in, &
                                    coord_type_in, &
                                    problo_in, probhi_in, center_in)

        ! Passing data from C++ into f90

        use bl_constants_module, only: ZERO
        use prob_params_module

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
        integer, intent(in) :: Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
        integer, intent(in) :: coord_type_in
        double precision, intent(in) :: problo_in(dm), probhi_in(dm), center_in(dm)

        dim = dm

        physbc_lo(1:dm) = physbc_lo_in(1:dm)
        physbc_hi(1:dm) = physbc_hi_in(1:dm)

        Outflow    = Outflow_in
        Symmetry   = Symmetry_in
        SlipWall   = SlipWall_in
        NoSlipWall = NoSlipWall_in

        coord_type = coord_type_in

        problo = ZERO
        probhi = ZERO
        center = ZERO

        problo(1:dm) = problo_in(1:dm)
        probhi(1:dm) = probhi_in(1:dm)
        center(1:dm) = center_in(1:dm)        

        dg(:) = 1

        if (dim .lt. 2) then
           dg(2) = 0
        endif

        if (dim .lt. 3) then
           dg(3) = 0
        endif

      end subroutine set_problem_params
      
! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_refinement_params(max_level_in, dx_level_in)

        use prob_params_module

        implicit none

        integer,          intent(in) :: max_level_in
        double precision, intent(in) :: dx_level_in(3*(max_level_in+1))

        integer :: lev, dir
        
        max_level = max_level_in

        allocate(dx_level(1:3, 0:max_level))

        do lev = 0, max_level
           do dir = 1, 3
              dx_level(dir,lev) = dx_level_in(3*lev + dir)
           enddo
        enddo

      end subroutine set_refinement_params      
      
! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine ca_set_special_tagging_flag(dummy,flag) 
      double precision :: dummy 
      integer          :: flag
      end subroutine ca_set_special_tagging_flag

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_tagging_params(name, namlen)

        use tagging_params_module

        ! Initialize the tagging parameters

        integer :: namlen
        integer :: name(namlen)
        
        integer :: un, i, status

        integer, parameter :: maxlen = 256
        character (len=maxlen) :: probin

        namelist /tagging/ &
           denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
           enterr,     entgrad,   max_enterr_lev,   max_entgrad_lev, &
           velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
           presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
           temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev, &
           raderr,     radgrad,   max_raderr_lev,   max_radgrad_lev

        ! Set namelist defaults
        denerr = 1.d20
        dengrad = 1.d20
        max_denerr_lev = 10
        max_dengrad_lev = 10

        enterr = 1.d20
        entgrad = 1.d20
        max_enterr_lev = -1
        max_entgrad_lev = -1

        presserr = 1.d20
        pressgrad = 1.d20
        max_presserr_lev = -1
        max_pressgrad_lev = -1

        velerr  = 1.d20
        velgrad = 1.d20
        max_velerr_lev = -1
        max_velgrad_lev = -1

        temperr  = 1.d20
        tempgrad = 1.d20
        max_temperr_lev = -1
        max_tempgrad_lev = -1

        raderr  = 1.d20
        radgrad = 1.d20
        max_raderr_lev = -1
        max_radgrad_lev = -1

        ! create the filename
        if (namlen > maxlen) then
           print *, 'probin file name too long'
           stop
        endif

        do i = 1, namlen
           probin(i:i) = char(name(i))
        end do

        ! read in the namelist
        un = 9
        open (unit=un, file=probin(1:namlen), form='formatted', status='old')
        read (unit=un, nml=tagging, iostat=status)

        if (status < 0) then
           ! the namelist does not exist, so we just go with the defaults
           continue

        else if (status > 0) then
           ! some problem in the namelist
           print *, 'ERROR: problem in the tagging namelist'
           stop
        endif

        close (unit=un)

      end subroutine get_tagging_params



! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_sponge_params(name, namlen)

        use sponge_params_module

        ! Initialize the sponge parameters

        integer :: namlen
        integer :: name(namlen)
        
        integer :: un, i, status

        integer, parameter :: maxlen = 256
        character (len=maxlen) :: probin

        namelist /sponge/ &
             sponge_lower_radius, sponge_upper_radius, &
             sponge_lower_density, sponge_upper_density, &
             sponge_timescale

        ! Set namelist defaults

        sponge_lower_radius = -1.d0
        sponge_upper_radius = -1.d0
        
        sponge_lower_density = -1.d0
        sponge_upper_density = -1.d0
        
        sponge_timescale    = -1.d0

        ! create the filename
        if (namlen > maxlen) then
           print *, 'probin file name too long'
           stop
        endif

        do i = 1, namlen
           probin(i:i) = char(name(i))
        end do

        ! read in the namelist
        un = 9
        open (unit=un, file=probin(1:namlen), form='formatted', status='old')
        read (unit=un, nml=sponge, iostat=status)

        if (status < 0) then
           ! the namelist does not exist, so we just go with the defaults
           continue

        else if (status > 0) then
           ! some problem in the namelist
           print *, 'ERROR: problem in the sponge namelist'
           stop
        endif

        close (unit=un)

      end subroutine get_sponge_params
