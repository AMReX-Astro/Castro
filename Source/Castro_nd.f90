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
                                   difmag_in, small_dens_in, small_temp_in, small_pres_in, &
                                   allow_negative_energy_in, &
                                   ppm_type_in,ppm_reference_in, &
                                   ppm_trace_grav_in, ppm_temp_fix_in, &
                                   ppm_tau_in_tracing_in, ppm_reference_edge_limit_in, &
                                   ppm_flatten_before_integrals_in, &
                                   ppm_reference_eigenvectors_in, &
                                   use_colglaz_in, &
                                   cg_maxiter_in, cg_tol_in, &
                                   use_pslope_in, &
                                   grav_source_type_in, &
                                   do_sponge_in,normalize_species_in,fix_mass_flux_in,use_sgs, &
                                   rot_period_in, const_grav_in)
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
        integer, intent(in) :: ppm_reference_in, ppm_trace_grav_in, ppm_temp_fix_in
        integer, intent(in) :: ppm_tau_in_tracing_in, ppm_reference_edge_limit_in
        integer, intent(in) :: ppm_flatten_before_integrals_in
        integer, intent(in) :: ppm_reference_eigenvectors_in
        integer, intent(in) :: use_colglaz_in, use_pslope_in, grav_source_type_in
        integer, intent(in) :: cg_maxiter_in
        double precision, intent(in) :: cg_tol_in
        integer, intent(in) :: do_sponge_in
        double precision, intent(in) :: difmag_in
        double precision, intent(in) :: small_dens_in, small_temp_in, small_pres_in
        integer, intent(in) :: normalize_species_in
        integer, intent(in) :: fix_mass_flux_in
        integer, intent(in) :: use_sgs
        double precision, intent(in) :: rot_period_in, const_grav_in

        integer             :: QLAST

        call parallel_initialize()

        iorder = 2 

!        difmag = 0.1d0
        difmag = difmag_in

        ! NTHERM: number of thermodynamic variables
        ! NVAR  : number of total variables in initial system
        ! dm refers to mometum components, '4' refers to rho, rhoE, rhoe and T
        NTHERM = dm + 4
        if (use_sgs .eq. 1) NTHERM = NTHERM + 1
        NVAR = NTHERM + nspec + naux + numadv

        nadv = numadv

        ! We use these to index into the state "U"
        URHO  = Density   + 1
        UMX   = Xmom      + 1
        if (dm .ge. 2) UMY = UMX + 1
        if (dm .eq. 3) UMZ = UMY + 1
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

        ! QTHERM: number of primitive variables, which includes pressure (+1), but
        !         not little e (-1)
        ! QVAR  : number of total variables in primitive form

        QTHERM = NTHERM
        if (use_colglaz_in == 1) QTHERM = QTHERM + 2

        QVAR = QTHERM + nspec + naux + numadv

        ! We use these to index into the state "Q"
        QRHO  = 1

        QU    = 2
        QLAST = 2

        if (dm .ge. 2) then
           QV    = 3
           QLAST = 3
        end if

        if (dm .eq. 3) then
           QW    = 4
           QLAST = 4
        end if

        if (use_colglaz_in == 1) then
           QGAME   = QLAST + 1
           QGAMC   = QLAST + 2
           QLAST   = QGAMC
        endif

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
          QFA = 1
          QFS = QTHERM + 1
        end if
        if (naux .ge. 1) then
          QFX = QFS + nspec
        else 
          QFX = 1
        end if

        if (small_pres_in > 0.d0) then
          small_pres = small_pres_in
        else
          small_pres = 1.d-8
        end if

        call eos_init(small_dens=small_dens_in, small_temp=small_temp_in)

        call eos_get_small_dens(small_dens)
        call eos_get_small_temp(small_temp)

        allow_negative_energy      = allow_negative_energy_in
        ppm_type                   = ppm_type_in
        ppm_reference              = ppm_reference_in
        ppm_trace_grav             = ppm_trace_grav_in
        ppm_temp_fix               = ppm_temp_fix_in
        ppm_tau_in_tracing         = ppm_tau_in_tracing_in
        ppm_reference_edge_limit   = ppm_reference_edge_limit_in
        ppm_flatten_before_integrals = ppm_flatten_before_integrals_in
        ppm_reference_eigenvectors = ppm_reference_eigenvectors_in
        use_colglaz                = use_colglaz_in
        cg_tol                     = cg_tol_in
        cg_maxiter                 = cg_maxiter_in
        use_pslope                 = use_pslope_in
        grav_source_type           = grav_source_type_in
        do_sponge                  = do_sponge_in
        normalize_species          = normalize_species_in
        fix_mass_flux              = fix_mass_flux_in
        rot_period                 = rot_period_in
        const_grav                 = const_grav_in
        

!       allocate(outflow_bc_lo(dm))
!       allocate(outflow_bc_hi(dm))

!       outflow_bc_lo(:) = phys_bc_lo(:)
!       outflow_bc_hi(:) = phys_bc_hi(:)

!       print *,'OUTFLOW_BC_LO ',outflow_bc_lo(:)
!       print *,'OUTFLOW_BC_HI ',outflow_bc_hi(:)

      end subroutine set_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_problem_params(dm,physbc_lo_in,physbc_hi_in,&
                                    Outflow_in,Symmetry_in,SlipWall_in,NoSlipWall_in, &
                                    coord_type_in)

        ! Passing data from C++ into f90

        use prob_params_module

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
        integer, intent(in) :: Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
        integer, intent(in) :: coord_type_in

        allocate(physbc_lo(dm))
        allocate(physbc_hi(dm))

        physbc_lo(:) = physbc_lo_in(:)
        physbc_hi(:) = physbc_hi_in(:)

        Outflow    = Outflow_in
        Symmetry   = Symmetry_in
        SlipWall   = SlipWall_in
        NoSlipWall = NoSlipWall_in

        coord_type = coord_type_in

      end subroutine set_problem_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine ca_set_special_tagging_flag(dummy,flag) 
      use probdata_module
      double precision :: dummy 
      integer          :: flag
      end subroutine ca_set_special_tagging_flag
