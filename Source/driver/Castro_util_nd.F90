module castro_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  function position(i, j, k, ccx, ccy, ccz)
    ! Given 3D indices (i,j,k), return the cell-centered spatial position.
    ! Optionally we can also be edge-centered in any of the directions.
    !

    use amrinfo_module, only: amr_level
    use prob_params_module, only: problo, probhi, physbc_lo, physbc_hi, dx_level, &
         domlo_level, domhi_level, Interior
    use amrex_constants_module, only: ZERO, HALF
    use amrex_fort_module, only: rt => amrex_real

    ! Input arguments

    integer :: i, j, k
    logical, optional :: ccx, ccy, ccz

    ! Local variables
    real(rt) :: position(3), dx(3), offset(3)
    integer  :: idx(3)
    logical  :: cc(3)
    integer  :: domlo(3), domhi(3)
    integer  :: dir

    !$gpu

    idx = [ i, j, k ]

    dx(:) = dx_level(:,amr_level)
    domlo = domlo_level(:,amr_level)
    domhi = domhi_level(:,amr_level)

    offset(:) = problo(:)

    cc(:) = .true.

    if (present(ccx)) then
       cc(1) = ccx
    endif

    if (present(ccy)) then
       cc(2) = ccy
    endif

    if (present(ccz)) then
       cc(3) = ccz
    endif

    do dir = 1, 3
       if (cc(dir)) then
          ! If we're cell-centered, we want to be in the middle of the zone.

          offset(dir) = offset(dir) + HALF * dx(dir)
       else
          ! Take care of the fact that for edge-centered indexing,
          ! we actually range from (domlo, domhi+1).

          domhi(dir) = domhi(dir) + 1
       endif
    enddo

    ! Be careful when using periodic boundary conditions. In that case,
    ! we need to loop around to the other side of the domain.

    do dir = 1, 3
       if      (physbc_lo(dir) .eq. Interior .and. idx(dir) .lt. domlo(dir)) then
          offset(dir) = offset(dir) + (probhi(dir) - problo(dir))
       else if (physbc_hi(dir) .eq. Interior .and. idx(dir) .gt. domhi(dir)) then
          offset(dir) = offset(dir) + (problo(dir) - probhi(dir))
       endif
    enddo

    position(:) = offset(:) + dble(idx(:)) * dx(:)

  end function position


  subroutine ca_find_center(data,new_center,icen,dx,problo) &
       bind(C, name="ca_find_center")
    !
    ! .. note::
    !    Binds to C function ``ca_find_center``

    use amrex_constants_module, only: ZERO, HALF, TWO
    use prob_params_module, only: dg, dim
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: data(-1:1,-1*dg(2):1*dg(2),-1*dg(3):1*dg(3))
    real(rt), intent(  out) :: new_center(3)
    real(rt), intent(in   ) :: dx(3),problo(3)

    real(rt) :: a,b,x,y,z,cen
    integer  :: icen(3)
    integer  :: i,j,k

    if (dim .eq. 1) then

       ! In 1-D it only make sense to have the center at the origin
       new_center = ZERO

    else if (dim .ge. 2) then

       ! We do this to take care of precision issues
       cen = data(0,0,0)
       do k = -1*dg(3),1*dg(3)
          do j = -1*dg(2),1*dg(2)
             do i = -1*dg(1),1*dg(1)
                data(i,j,k) = data(i,j,k) - cen
             end do
          end do
       end do

       ! This puts the "center" at the cell center
       new_center(1:dim) = problo(1:dim) +  (icen(1:dim)+HALF) * dx(1:dim)

       ! Fit parabola y = a x^2  + b x + c through three points
       ! a = 1/2 ( y_1 + y_-1)
       ! b = 1/2 ( y_1 - y_-1)
       ! x_vertex = -b / 2a

       ! ... in x-direction
       a = HALF * (data(1,0,0) + data(-1,0,0)) - data(0,0,0)
       b = HALF * (data(1,0,0) - data(-1,0,0)) - data(0,0,0)
       x = -b / (TWO*a)
       new_center(1) = new_center(1) +  x*dx(1)

       ! ... in y-direction
       a = HALF * (data(0,1,0) + data(0,-1,0)) - data(0,0,0)
       b = HALF * (data(0,1,0) - data(0,-1,0)) - data(0,0,0)
       y = -b / (TWO*a)
       new_center(2) = new_center(2) +  y*dx(2)

       if (dim .eq. 3) then

          ! ... in z-direction
          a = HALF * (data(0,0,1) + data(0,0,-1)) - data(0,0,0)
          b = HALF * (data(0,0,1) - data(0,0,-1)) - data(0,0,0)
          z = -b / (TWO*a)
          new_center(3) = new_center(3) +  z*dx(3)

       endif

    endif

  end subroutine ca_find_center


  subroutine ca_compute_avgstate(lo,hi,dx,dr,nc,&
       state,s_lo,s_hi,radial_state, &
       vol,v_lo,v_hi,radial_vol, &
       problo,numpts_1d) &
       bind(C, name="ca_compute_avgstate")
    !
    ! .. note::
    !    Binds to C function ``ca_compute_avgstate``

    use meth_params_module, only: URHO, UMX, UMY, UMZ
    use prob_params_module, only: dim
    use probdata_module, only: center
    use amrex_constants_module, only: HALF
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3),hi(3),nc
    real(rt), intent(in   ) :: dx(3),dr,problo(3)

    integer,  intent(in   ) :: numpts_1d
    real(rt), intent(inout) :: radial_state(nc,0:numpts_1d-1)
    real(rt), intent(inout) :: radial_vol(0:numpts_1d-1)

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nc)

    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer  :: i,j,k,n,index
    real(rt) :: x,y,z,r
    real(rt) :: x_mom,y_mom,z_mom,radial_mom

#ifndef AMREX_USE_CUDA
    if (dim .eq. 1) call castro_error("Error: cannot do ca_compute_avgstate in 1D.")
#endif

    !
    ! Do not OMP this.
    !
    do k = lo(3), hi(3)
       z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             r = sqrt(x**2 + y**2 + z**2)
             index = int(r/dr)
#ifndef AMREX_USE_CUDA
             if (index .gt. numpts_1d-1) then
                print *,'COMPUTE_AVGSTATE: INDEX TOO BIG ',index,' > ',numpts_1d-1
                print *,'AT (i,j,k) ',i,j,k
                print *,'R / DR ',r,dr
                call castro_error("Error:: Castro_util.F90 :: ca_compute_avgstate")
             end if
#endif
             radial_state(URHO,index) = radial_state(URHO,index) &
                  + vol(i,j,k)*state(i,j,k,URHO)
             !
             ! Store the radial component of the momentum in the
             ! UMX, UMY and UMZ components for now.
             !
             x_mom = state(i,j,k,UMX)
             y_mom = state(i,j,k,UMY)
             z_mom = state(i,j,k,UMZ)
             radial_mom = x_mom * (x/r) + y_mom * (y/r) + z_mom * (z/r)
             radial_state(UMX,index) = radial_state(UMX,index) + vol(i,j,k)*radial_mom
             radial_state(UMY,index) = radial_state(UMY,index) + vol(i,j,k)*radial_mom
             radial_state(UMZ,index) = radial_state(UMZ,index) + vol(i,j,k)*radial_mom

             do n = UMZ+1,nc
                radial_state(n,index) = radial_state(n,index) + vol(i,j,k)*state(i,j,k,n)
             end do
             radial_vol(index) = radial_vol(index) + vol(i,j,k)
          enddo
       enddo
    enddo

  end subroutine ca_compute_avgstate

end module castro_util_module
