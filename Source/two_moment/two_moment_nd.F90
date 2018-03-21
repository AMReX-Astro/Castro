  subroutine call_to_thornado(lo, hi, &
                                S_new, S_lo, S_hi, nf, &
                                U_R_o, &
                                U_R_n, U_R_lo,   U_R_hi, nr, &
                               dU_F, dU_F_lo, dU_F_hi, ndf) &
                               bind(C, name="call_to_thornado")

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,NE_COMP,NVAR

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  U_F_lo(3),  U_F_hi(3)
    integer, intent(in) ::  U_R_lo(3),  U_R_hi(3)
    integer, intent(in) :: dU_F_lo(3), dU_F_hi(3)
    integer, intent(in) ::  nf,nr, ndf

    ! Here we expect nf = ndf = 6
    !                nr = 20 x 16 x 6 x 4

    real(rt), intent(inout) ::  U_F  ( U_F_lo(1): U_F_hi(1),  U_F_lo(2): U_F_hi(2),  U_F_lo(3): U_F_hi(3), nf)
    real(rt), intent(inout) ::  U_R_o(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), nr) 
    real(rt), intent(inout) ::  U_R_n(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), nr) 
    real(rt), intent(inout) :: dU_F  (dU_F_lo(1):dU_F_hi(1), dU_F_lo(2):dU_F_hi(2), dU_F_lo(3):dU_F_hi(3), ndf)

    integer, parameter :: n_energy = 20
    integer, parameter :: n_nodes  = 16
    integer, parameter :: n_species = 6
    integer, parameter :: n_moments = 4

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
    real(rt) :: dU_R( U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),  U_R_lo(3): U_R_hi(3), nr)

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

    iz_b1(1) = hi(1)+1
    iz_b1(2) = hi(2)+1
    iz_b1(3) = hi(3)+1
    iz_b1(4) = n_energy

    ne = n_energy
    nn = n_nodes
    ns = n_species
    nm = n_moments

    allocate  U_F_thor(1:4, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,6)
    allocate dU_F_thor(1:4, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,6)

    allocate  U_R_thor(1:nn, 1:ne, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1:nm, 1:ns)
    allocate dU_R_thor(1:nn, 1:ne, lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1:nm, 1:ns)

    do k = lo(3)-2,hi(3)+2
    do j = lo(2)-2,hi(2)+2
    do i = lo(1)-2,hi(1)+2

         U_F_thor(1:nn,i,j,k,thor_density) = S_new(i,j,k,URHO)
         U_F_thor(1:nn,i,j,k,thor_xmom   ) = S_new(i,j,k,UMX)
         U_F_thor(1:nn,i,j,k,thor_ymom   ) = S_new(i,j,k,UMY)
         U_F_thor(1:nn,i,j,k,thor_zmom   ) = S_new(i,j,k,UMZ)
         U_F_thor(1:nn,i,j,k,thor_rhoe   ) = S_new(i,j,k,UEINT)
         U_F_thor(1:nn,i,j,k,thor_ne     ) = S_new(i,j,k,UFX)

         do is = 1, ns
            ii = (is-1)*(ncr*ndof*nenergy) 
            do im = 1, nm
               ii = ii + (im-1)*(ndof*nenergy)
               do ie = 1, ne
                  ii = ii + (ie-1)*ndof
                  do id = 1, nn
                     ii = ii + (id-1)
                     U_R_thor(id,ie,j,k,im,is) = U_R(i,j,k,ii)
                  end do
               end do
            end do
         end do

    end do
    end do
    end do

    ! This is calling the Fortran interface that lives in the thornado repo
    call ComputeIncrement(iz_b0, iz_e0, iz_b1, iz_e1, U_F_thor, U_R_thor, du_F_thor, du_R_thor, nn, ns, nm)

    do k = lo(3)-2,hi(3)+2
    do j = lo(2)-2,hi(2)+2
    do i = lo(1)-2,hi(1)+2
         S_new(i,j,k,URHO ) = S_new(i,j,k,URHO ) + dU_F_thor(1:dof,i,j,k,thor_density)
         S_new(i,j,k,UMX  ) = S_new(i,j,k,UMX  ) + dU_F_thor(1:dof,i,j,k,thor_xmom)
         S_new(i,j,k,UMY  ) = S_new(i,j,k,UMY  ) + dU_F_thor(1:dof,i,j,k,thor_ymom)
         S_new(i,j,k,UMZ  ) = S_new(i,j,k,UMZ  ) + dU_F_thor(1:dof,i,j,k,thor_zmom)
         S_new(i,j,k,UEINT) = S_new(i,j,k,UEINT) + dU_F_thor(1:dof,i,j,k,thor_rhoe)
         S_new(i,j,k,UFX  ) = S_new(i,j,k,UFX  ) + dU_F_thor(1:dof,i,j,k,thor_ne)

         do is = 1, ns
            ii = (is-1)*(ncr*ndof*nenergy) 
            do im = 1, nm
               ii = ii + (im-1)*(ndof*nenergy)
               do ie = 1, ne
                  ii = ii + (ie-1)*ndof
                  do id = 1, nn
                     ii = ii + (id-1)
                     U_R(i,j,k,nn) = U_R_thor(id,ie,j,k,im,is) 
                  end do
               end do
            end do
         end do

    end do
    end do
    end do


    ! Go ahead and update U_R here.
    U_r_n = U_r_o + dU_R

  end subroutine call_to_thornado

  subroutine ComputeIncrement(iz_b0, iz_e0, iz_b1, iz_e1, U_F, U_R, dU_F, dU_R, nn, ns, nm)

    integer, intent(in) :: iz_b0(4), iz_e0(4)
    integer, intent(in) :: iz_b1(4), iz_e1(4)
    integer, intent(in) :: nn, ns, nm

    double precision, intent(inout) ::  U_F(  iz_b0(1):iz_e0(1),  iz_b0(2):iz_e0(2), &
                                              iz_b0(3):iz_e0(3),  iz_b0(4):iz_e0(4), &
                                              1:nfluid_vars)
    double precision, intent(inout) ::  U_R(  iz_b1(1):iz_e1(1),  iz_b1(2):iz_e1(2), &
                                              iz_b1(3):iz_e1(3),  iz_b1(4):iz_e1(4), &
                                              1:nn, 1:ns, 1:nm)
    double precision, intent(inout) :: dU_F(  iz_b0(1):iz_e0(1),  iz_b0(2):iz_e0(2), &
                                              iz_b0(3):iz_e0(3),  iz_b0(4):iz_e0(4), &
                                              1:nfluid_vars)
    double precision, intent(inout) :: dU_R(  iz_b1(1):iz_e1(1),  iz_b1(2):iz_e1(2), &
                                              iz_b1(3):iz_e1(3),  iz_b1(4):iz_e1(4), &
                                              1:nn, 1:ns, 1:nm)

    ! get ordering of fluid fields in FluidFieldsModule.f90
    ! (DOF, i, j, k, 1:nfluid_vars)

    ! This will actually be what Eirik gives us

  end subroutine ComputeIncrement

  subroutine init_thornado(lo, hi, &
                           U_F, U_F_lo,   U_F_hi, nf, &
                           U_R, U_R_lo,   U_R_hi, nr) &
                           bind(C, name="init_thornado")

    use amrex_fort_module, only : rt => amrex_real

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  U_F_lo(3),  U_F_hi(3)
    integer, intent(in) ::  U_R_lo(3),  U_R_hi(3)
    integer, intent(in) ::  nf,nr

    ! Here we expect nf = ndf = 6
    !                nr = 20 x 16 x 6 x 4

    real(rt), intent(in   ) ::  U_F( U_F_lo(1): U_F_hi(1),  U_F_lo(2): U_F_hi(2),  U_F_lo(3): U_F_hi(3), nf)
    real(rt), intent(inout) ::  U_R(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), nr) 

    integer, parameter :: n_energy = 20
    integer, parameter :: n_nodes  = 16
    integer, parameter :: n_species = 6
    integer, parameter :: n_moments = 4

    integer  :: ne,nn,ns,nm
    integer  :: iz_b0(4), iz_e0(4)
    integer  :: iz_b1(4), iz_e1(4)

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

    iz_b1(1) = hi(1)+1
    iz_b1(2) = hi(2)+1
    iz_b1(3) = hi(3)+1
    iz_b1(4) = n_energy

    ne = n_energy
    nn = n_nodes
    ns = n_species
    nm = n_moments

    ! This is calling the Fortran interface that lives in the thornado repo
    call InitThornado(iz_b0, iz_e0, iz_b1, iz_e1, U_F, U_R, nn, ns, nm)

  end subroutine init_thornado

  subroutine InitThornado(iz_b0, iz_e0, iz_b1, iz_e1, U_F, U_R, nn, ns, nm)

    integer, intent(in) :: iz_b0(4), iz_e0(4)
    integer, intent(in) :: iz_b1(4), iz_e1(4)
    integer, intent(in) :: nn, ns, nm

    double precision, intent(inout) ::  U_F(  iz_b0(1):iz_e0(1),  iz_b0(2):iz_e0(2), &
                                              iz_b0(3):iz_e0(3),  iz_b0(4):iz_e0(4), &
                                              1:nn, 1:ns, 1:nm)
    double precision, intent(inout) ::  U_R(  iz_b1(1):iz_e1(1),  iz_b1(2):iz_e1(2), &
                                              iz_b1(3):iz_e1(3),  iz_b1(4):iz_e1(4), &
                                              1:nn, 1:ns, 1:nm)

  ! This will also be from Eirik

  end subroutine InitThornado
