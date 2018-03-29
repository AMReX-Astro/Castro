  subroutine call_to_thornado(lo, hi, dt, &
                              S, dS, s_lo, s_hi, &
                              U_R_o, U_R_n, U_R_lo, U_R_hi, nr,
                              n_nodes, n_energy, n_species, n_moments) &
                              bind(C, name="call_to_thornado")

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UFX,NVAR

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) ::  U_R_lo(3),  U_R_hi(3)
    integer, intent(in) :: dU_F_lo(3), dU_F_hi(3)
    integer, intent(in) ::  nr
    real(rt), intent(in) :: dt !! KS: missing declaration

    ! Here we expect  nr = 20 x 16 x 6 x 4 (energy x nodes x species x moments)

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR) 
    real(rt), intent(inout) :: dS(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), nr)
    real(rt), intent(inout) ::  U_R_n(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), nr) 

    integer, parameter :: thor_density = 1
    integer, parameter :: thor_xmom    = 2
    integer, parameter :: thor_ymom    = 3
    integer, parameter :: thor_zmom    = 4
    integer, parameter :: thor_rhoe    = 5
    integer, parameter :: thor_ne      = 6
    integer, parameter :: thor_nfluid  = 6

    ! Temporary variables
    integer  :: ne,nn,ns,nm
    integer  :: iz_b0(4), iz_e0(4)
    integer  :: iz_b1(4), iz_e1(4)
    integer  :: i,j,k
    integer  :: ii,id,ie,im,is
    integer  :: ng !! KS: missing declaration

    real(rt), allocatable ::  U_F_thor(:,:,:,:,:)
    real(rt), allocatable :: dU_F_thor(:,:,:,:,:)

    real(rt), allocatable ::  U_R_thor(:,:,:,:,:,:,:)
    real(rt), allocatable :: dU_R_thor(:,:,:,:,:,:,:)

    ! Sanity check on size of arrays
    ! Note that we have set ngrow_thornado = ngrow_state in Castro_setup.cpp
    if (U_R_lo(1) .ne. S_lo(1) .or.  U_R_hi(1) .ne. S_hi(1) .or. &
        U_R_lo(2) .ne. S_lo(2) .or.  U_R_hi(2) .ne. S_hi(2) .or. &
        U_R_lo(3) .ne. S_lo(3) .or.  U_R_hi(3) .ne. S_hi(3)) then
        print *,'INCONSISTENT ARRAY BOUNDS ON FLUID AND RADIATION VARS' !! KS: added closing '
        stop
    endif
    ! End Sanity check 

    iz_b0(1) = lo(1)
    iz_b0(2) = lo(2)
    iz_b0(3) = lo(3)
    iz_b0(4) = 1

    iz_e0(1) = hi(1)
    iz_e0(2) = hi(2)
    iz_e0(3) = hi(3)
    iz_e0(4) = n_energy

    iz_b1(1) = lo(1)-1
    iz_b1(2) = lo(2)-1
    iz_b1(3) = lo(3)-1
    iz_b1(4) = 1

    iz_e1(1) = hi(1)+1
    iz_e1(2) = hi(2)+1
    iz_e1(3) = hi(3)+1
    iz_e1(4) = n_energy

    nn = n_nodes
    ne = n_energy
    ns = n_species
    nm = n_moments

    !! Should this be 1:nn instead of 1:4 given lines 113-118?
    !! Actually, why is this an array at all since only the first entry is used in lines 152-157?
    !! For clarity, should 1:6 be 1:thor_ne?

    !! ASA:I made up the idea that we only use the first component below -- should that be the
    !! ASA:    average of the nn components??? 
    !! ASA: I also made up the 1:4 -- I'm confused about what that should be

    ! ************************************************************************************
    ! ASA: Set ng for the temporaries here -- don't know what it should be
    ng = 2
    !! ASA: Another question -- should the definition of iz_b*, iz_e* above be consistent with the size 
    !!      of the arrays, i.e should iz_e1 be hi(1)+ng instead of hi(1)+1?
    ! ************************************************************************************
    allocate( U_F_thor(1:4, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,1:6))
    allocate(dU_F_thor(1:4, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,1:6))

    allocate( U_R_thor(1:nn, 1:ne, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1:nm, 1:ns))
    allocate(dU_R_thor(1:nn, 1:ne, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1:nm, 1:ns))

    ! ************************************************************************************
    ! Copy from the Castro arrays into temporary thornado arrays
    ! ************************************************************************************
    do k = lo(3)-ng,hi(3)+ng
    do j = lo(2)-ng,hi(2)+ng
    do i = lo(1)-ng,hi(1)+ng

         U_F_thor(1:nn,i,j,k,thor_density) = S_new(i,j,k,URHO)
         U_F_thor(1:nn,i,j,k,thor_xmom   ) = S_new(i,j,k,UMX)
         U_F_thor(1:nn,i,j,k,thor_ymom   ) = S_new(i,j,k,UMY)
         U_F_thor(1:nn,i,j,k,thor_zmom   ) = S_new(i,j,k,UMZ)
         U_F_thor(1:nn,i,j,k,thor_rhoe   ) = S_new(i,j,k,UEINT)
         U_F_thor(1:nn,i,j,k,thor_ne     ) = S_new(i,j,k,UFX)

         do is = 1, ns
         do im = 1, nm
         do ie = 1, ne
         do id = 1, nn
            ii = (is-1)*(nm*ne*nn) + (im-1)*(ne*nn) + (ie-1)*nn + (id-1)
            U_R_thor(id,ie,i,j,k,im,is) = U_R_o(i,j,k,ii) !! KS: I changed icomp to ii
         end do
         end do
         end do
         end do

    end do
    end do
    end do

    ! ************************************************************************************
    ! Call the Fortran interface that lives in the thornado repo
    ! ************************************************************************************
    ! call ComputeIncrement(iz_b0, iz_e0, iz_b1, iz_e1, U_F_thor, U_R_thor, dU_F_thor, dU_R_thor, nn, ns, nm)

    ! ************************************************************************************
    ! Copy back from the thornado arrays into Castro arrays
    ! ************************************************************************************
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

         ! We store dS as a source term which we can add to S_new outside of this routine

         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ASA: But why are we choosing 1 for the first component -- I made that up! 
         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         dS(i,j,k,URHO ) = dU_F_thor(1,i,j,k,thor_density) / dt !! KS: changed these to dS since dS_F doesn't exist
         dS(i,j,k,UMX  ) = dU_F_thor(1,i,j,k,thor_xmom)    / dt
         dS(i,j,k,UMY  ) = dU_F_thor(1,i,j,k,thor_ymom)    / dt
         dS(i,j,k,UMZ  ) = dU_F_thor(1,i,j,k,thor_zmom)    / dt
         dS(i,j,k,UEINT) = dU_F_thor(1,i,j,k,thor_rhoe)    / dt
         dS(i,j,k,UFX  ) = dU_F_thor(1,i,j,k,thor_ne)      / dt

         do is = 1, ns
         do im = 1, nm
         do ie = 1, ne
         do id = 1, nn
            ii = (is-1)*(nm*ne*nn) + (im-1)*(ne*nn) + (ie-1)*nn + (id-1)
            U_R_n(i,j,k,ii) = U_R_thor(id,ie,i,j,k,im,is) 
         end do
         end do
         end do
         end do

    end do
    end do
    end do

  end subroutine call_to_thornado
