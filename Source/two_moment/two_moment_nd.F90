  subroutine call_to_thornado(lo, hi, &
                                U_F, U_F_lo,   U_F_hi, nf, &
                                U_R_o, &
                                U_R_n, U_R_lo,   U_R_hi, nr, &
                               dU_F, dU_F_lo, dU_F_hi, ndf) &
                               bind(C, name="call_to_thornado")

    use amrex_fort_module, only : rt => amrex_real

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

    ! This is calling the Fortran interface that lives in the thornado repo
    call ComputeIncrement(iz_b0, iz_e0, iz_b1, iz_e1, U_F, U_R_o, du_F, du_R, nn, ns, nm)

    ! Go ahead and update U_R here.
    U_r_n = U_r_o + dU_R

  end subroutine call_to_thornado

  subroutine ComputeIncrement(iz_b0, iz_e0, iz_b1, iz_e1, U_F, U_R, dU_F, dU_R, nn, ns, nm)

    integer, intent(in) :: iz_b0(4), iz_e0(4)
    integer, intent(in) :: iz_b1(4), iz_e1(4)
    integer, intent(in) :: nn, ns, nm

    double precision, intent(inout) ::  U_F(  iz_b0(1):iz_e0(1),  iz_b0(2):iz_e0(2), &
                                              iz_b0(3):iz_e0(3),  iz_b0(4):iz_e0(4), &
                                              1:nn, 1:ns, 1:nm)
    double precision, intent(inout) ::  U_R(  iz_b1(1):iz_e1(1),  iz_b1(2):iz_e1(2), &
                                              iz_b1(3):iz_e1(3),  iz_b1(4):iz_e1(4), &
                                              1:nn, 1:ns, 1:nm)
    double precision, intent(inout) :: dU_F(  iz_b0(1):iz_e0(1),  iz_b0(2):iz_e0(2), &
                                              iz_b0(3):iz_e0(3),  iz_b0(4):iz_e0(4), &
                                              1:nn, 1:ns, 1:nm)
    double precision, intent(inout) :: dU_R(  iz_b1(1):iz_e1(1),  iz_b1(2):iz_e1(2), &
                                              iz_b1(3):iz_e1(3),  iz_b1(4):iz_e1(4), &
                                              1:nn, 1:ns, 1:nm)

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
