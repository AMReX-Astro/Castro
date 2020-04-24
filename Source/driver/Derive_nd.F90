module derive_module
    ! All subroutines in this file must be threadsafe because they are called
    ! inside OpenMP parallel regions.

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains

#ifdef DIFFUSION
  subroutine dercond(lo, hi, &
                     cond, u_lo, u_hi, nd, &
                     state, d_lo, d_hi, nc, &
                     domlo, domhi, delta) bind(C, name="dercond")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO
    use meth_params_module, only: diffuse_cutoff_density, diffuse_cutoff_density_hi, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
    use conductivity_module, only: conductivity

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3)
    real(rt), intent(inout) :: cond(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in   ) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ), value :: nd, nc

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: multiplier

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%e      = state(i,j,k,UEINT)/state(i,j,k,URHO)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)/ state(i,j,k,URHO)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)/ state(i,j,k,URHO)
             call eos(eos_input_re,eos_state)

             if (eos_state%rho > diffuse_cutoff_density) then
                call conductivity(eos_state)

                if (eos_state%rho < diffuse_cutoff_density_hi) then
                    multiplier = (eos_state%rho - diffuse_cutoff_density) / &
                            (diffuse_cutoff_density_hi - diffuse_cutoff_density)
                    eos_state % conductivity = eos_state % conductivity * multiplier
                endif
             else
                eos_state % conductivity = ZERO
             endif

             cond(i,j,k,1) = eos_state % conductivity

          enddo
       enddo
    enddo

  end subroutine dercond


  subroutine derdiffcoeff(lo, hi, &
                          diff, u_lo, u_hi, nd, &
                          state, d_lo, d_hi, nc, &
                          domlo, domhi, delta) bind(C, name="derdiffcoeff")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO
    use meth_params_module, only: diffuse_cutoff_density, diffuse_cutoff_density_hi, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
    use conductivity_module, only: conductivity

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3)
    real(rt), intent(inout) :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in   ) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ), value :: nd, nc

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: multiplier

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%e      = state(i,j,k,UEINT)/state(i,j,k,URHO)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)/ state(i,j,k,URHO)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)/ state(i,j,k,URHO)
             call eos(eos_input_re,eos_state)

             if (eos_state%rho > diffuse_cutoff_density) then
                call conductivity(eos_state)

                if (eos_state%rho < diffuse_cutoff_density_hi) then
                    multiplier = (eos_state%rho - diffuse_cutoff_density) / &
                            (diffuse_cutoff_density_hi - diffuse_cutoff_density)
                    eos_state % conductivity = eos_state % conductivity * multiplier
                endif
             else
                eos_state % conductivity = ZERO
             endif

             diff(i,j,k,1) = eos_state % conductivity/(eos_state%rho * eos_state%cv)

          enddo
       enddo
    enddo

  end subroutine derdiffcoeff


  subroutine derdiffterm(lo, hi, &
                         diff, u_lo, u_hi, nd, &
                         state, d_lo, d_hi, nc, &
                         coeff_x, x_lo, x_hi, &
#if AMREX_SPACEDIM >= 2
                         coeff_y, y_lo, y_hi, &
#endif
#if AMREX_SPACEDIM == 3
                         coeff_z, z_lo, z_hi, &
#endif
                         domlo, domhi, delta) bind(C, name="derdiffterm")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: UTEMP
    use prob_params_module, only: problo, coord_type

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: x_lo(3), x_hi(3)
#if AMREX_SPACEDIM >= 2
    integer,  intent(in   ) :: y_lo(3), y_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer,  intent(in   ) :: z_lo(3), z_hi(3)
#endif
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3)
    real(rt), intent(inout) :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in   ) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    real(rt), intent(in   ) :: coeff_x(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(in   ) :: coeff_y(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(in   ) :: coeff_z(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
#endif
    integer,  intent(in   ), value :: nd, nc

    real(rt) :: diff_term
    real(rt) :: kgradT_xhi, kgradT_xlo, kgradT_yhi, kgradT_ylo, kgradT_zhi, kgradT_zlo
    integer  :: i, j, k

    real(rt) :: r, rp1, rm1

    !$gpu

    ! create the diff term
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             kgradT_xhi = coeff_x(i+1,j,k)*(state(i+1,j,k,UTEMP) - state(i,  j,k,UTEMP))/delta(1)
             kgradT_xlo = coeff_x(i  ,j,k)*(state(i  ,j,k,UTEMP) - state(i-1,j,k,UTEMP))/delta(1)
#if AMREX_SPACEDIM >= 2
             kgradT_yhi = coeff_y(i,j+1,k)*(state(i,j+1,k,UTEMP) - state(i,j  ,k,UTEMP))/delta(2)
             kgradT_ylo = coeff_y(i,j  ,k)*(state(i,j,  k,UTEMP) - state(i,j-1,k,UTEMP))/delta(2)
#endif
#if AMREX_SPACEDIM == 3
             kgradT_zhi = coeff_z(i,j,k+1)*(state(i,j,k+1,UTEMP) - state(i,j,k  ,UTEMP))/delta(3)
             kgradT_zlo = coeff_z(i,j,k  )*(state(i,j,k  ,UTEMP) - state(i,j,k-1,UTEMP))/delta(3)
#endif

             if (coord_type == 0) then
                diff_term = (kgradT_xhi - kgradT_xlo)/delta(1)
#if AMREX_SPACEDIM >= 2
                diff_term = diff_term + (kgradT_yhi - kgradT_ylo)/delta(2)
#endif
#if AMREX_SPACEDIM == 3
                diff_term = diff_term + (kgradT_zhi - kgradT_zlo)/delta(3)
#endif

             else if (coord_type == 1) then
                ! axisymmetric coords (2-d)
                r = dble(i + HALF)*delta(1) + problo(1)
                rm1 = dble(i - ONE + HALF)*delta(1) + problo(1)
                rp1 = dble(i + ONE + HALF)*delta(1) + problo(1)

                diff_term = (rp1*kgradT_xhi - rm1*kgradT_xlo)/(r*delta(1)) + &
                            (kgradT_yhi - kgradT_ylo)/delta(2)

             else if (coord_type == 2) then
                ! spherical coords (1-d)
                r = dble(i + HALF)*delta(1) + problo(1)
                rm1 = dble(i - ONE + HALF)*delta(1) + problo(1)
                rp1 = dble(i + ONE + HALF)*delta(1) + problo(1)

                diff_term = (rp1**2*kgradT_xhi - rm1**2*kgradT_xlo)/(r**2*delta(1))

             endif

             diff(i,j,k,1) = diff_term

          enddo
       enddo
    enddo

  end subroutine derdiffterm

#endif

#ifdef MHD  
   subroutine ca_dermagcenx(mag_cen_x,mag_x_l1,mag_x_l2,mag_x_l3,mag_x_h1,mag_x_h2,mag_x_h3, nbx, &
                         dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                            bind(C, name="ca_dermagcenx")

      !
      ! This routine will derive cell centered magnetic field in x direction
      ! mag_cen_x = 1/2 (mag_x(i,j,k) + mag_x(i+1,j,k))
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          mag_x_l1,mag_x_l2,mag_x_l3,mag_x_h1,mag_x_h2,mag_x_h3, nbx
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) mag_cen_x(mag_x_l1:mag_x_h1,mag_x_l2:mag_x_h2,mag_x_l3:mag_x_h3,nbx)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
      integer    i,j,k

      ! 
      ! Here dat contains (mag_x)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                mag_cen_x(i,j,k,1) = 0.5*(dat(i,j,k,1) + dat(i+1,j,k,1))
            end do
         end do
      end do

   end subroutine ca_dermagcenx

   subroutine ca_dermagceny(mag_cen_y,mag_y_l1,mag_y_l2,mag_y_l3,mag_y_h1,mag_y_h2,mag_y_h3, nby, &
                         dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                            bind(C, name="ca_dermagceny")

      !
      ! This routine will derive cell centered magnetic field in y direction
      ! mag_cen_y = 1/2 (mag_y(i,j,k) + mag_y(i,j+1,k))
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          mag_y_l1,mag_y_l2,mag_y_l3,mag_y_h1,mag_y_h2,mag_y_h3, nby
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) mag_cen_y(mag_y_l1:mag_y_h1,mag_y_l2:mag_y_h2,mag_y_l3:mag_y_h3,nby)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      ! 
      ! Here dat contains (mag_y)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                mag_cen_y(i,j,k,1) = 0.5*(dat(i,j,k,1) + dat(i,j+1,k,1))
            end do
         end do
      end do


   end subroutine ca_dermagceny

   subroutine ca_dermagcenz(mag_cen_z,mag_z_l1,mag_z_l2,mag_z_l3,mag_z_h1,mag_z_h2,mag_z_h3, nbz, &
                         dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                         lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                         bind(C, name="ca_dermagcenz")
      !
      ! This routine will derive cell centered magnetic field in z direction
      ! mag_cen_z = 1/2 (mag_z(i,j,k) + mag_z(i,j,k+1))
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          mag_z_l1,mag_z_l2,mag_z_l3,mag_z_h1,mag_z_h2,mag_z_h3, nbz
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) mag_cen_z(mag_z_l1:mag_z_h1,mag_z_l2:mag_z_h2,mag_z_l3:mag_z_h3,nbz)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      !
      ! Here dat contains (mag_z)
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                mag_cen_z(i,j,k,1) = 0.5*(dat(i,j,k,1) + dat(i,j,k+1,1))
            end do
         end do
      end do

  end subroutine ca_dermagcenz

  subroutine ca_derex(E_x,E_x_l1,E_x_l2,E_x_l3,E_x_h1,E_x_h2,E_x_h3, Ebx, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                         bind(C, name="ca_derex")

      !
      ! This routine will derive the x component of the electric field
      ! E_x = -(V X B)_x
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          E_x_l1,E_x_l2,E_x_l3,E_x_h1,E_x_h2,E_x_h3, Ebx
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) E_x(E_x_l1:E_x_h1,E_x_l2:E_x_h2,E_x_l3:E_x_h3,Ebx)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      real vy,vz
      !
      ! Here dat contains (mag_y,mag_z,density,ymom,zmom)
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                vy = dat(i,j,k,4)/dat(i,j,k,3)
                vz = dat(i,j,k,5)/dat(i,j,k,3)
                E_x(i,j,k,1) = -vy*dat(i,j,k,2) + vz*dat(i,j,k,1)
            end do
         end do
      end do

  end subroutine ca_derex

  subroutine ca_derey(E_y,E_y_l1,E_y_l2,E_y_l3,E_y_h1,E_y_h2,E_y_h3, Eby, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                       bind(C, name="ca_derey")
      !
      ! This routine will derive the y component of the electric field
      ! E_y = -(V X B)_y
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          E_y_l1,E_y_l2,E_y_l3,E_y_h1,E_y_h2,E_y_h3, Eby
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) E_y(E_y_l1:E_y_h1,E_y_l2:E_y_h2,E_y_l3:E_y_h3,Eby)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      real vx,vz
      !
      ! Here dat contains (mag_x,mag_z,density,xmom,zmom)
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                vx = dat(i,j,k,4)/dat(i,j,k,3)
                vz = dat(i,j,k,5)/dat(i,j,k,3)
                E_y(i,j,k,1) = -vz*dat(i,j,k,1) + vx*dat(i,j,k,2)
            end do
         end do
      end do

  end subroutine ca_derey

  subroutine ca_derez(E_z,E_z_l1,E_z_l2,E_z_l3,E_z_h1,E_z_h2,E_z_h3, Ebz, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                      bind(C, name="ca_derez") 
      !
      ! This routine will derive the z component of the electric field
      ! E_z = -(V X B)_z
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          lo(3), hi(3)
      integer          E_z_l1,E_z_l2,E_z_l3,E_z_h1,E_z_h2,E_z_h3, Ebz
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) E_z(E_z_l1:E_z_h1,E_z_l2:E_z_h2,E_z_l3:E_z_h3,Ebz)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      real vx,vy
      ! 
      ! Here dat contains (mag_x,mag_y,density,xmom,ymom)
      ! 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                        vx = dat(i,j,k,4)/dat(i,j,k,3)
                        vy = dat(i,j,k,5)/dat(i,j,k,3)
                E_z(i,j,k,1) = -vx*dat(i,j,k,2) + vy*dat(i,j,k,1)
            end do
         end do
      end do

  end subroutine ca_derez

  subroutine ca_derdivb(divb,divb_l1,divb_l2,divb_l3,divb_h1,divb_h2,divb_h3, Ebz, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)&
                         bind(C, name="ca_derdivb")
      !
      ! This routine will derive the divergence of the magnetic field.
      ! out = div(B)
      !
      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer          lo(3), hi(3)
      integer          divb_l1,divb_l2,divb_l3,divb_h1,divb_h2,divb_h3, Ebz
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      real(rt) delta(3), xlo(3), time, dt
      real(rt) divb(divb_l1:divb_h1,divb_l2:divb_h2,divb_l3:divb_h3,Ebz)
      real(rt)    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      real dx,dy,dz
      !
      ! Here dat contains (mag_x,mag_y,mag_z)
      !
      ! write(*,*) 'Computing Divergence!'
      divb = 0.0
      do k = lo(3) + 1, hi(3)
         do j = lo(2) + 1, hi(2)
            do i = lo(1) + 1, hi(1)
                        dx = dat(i,j,k,1) - dat(i-1,j,k,1)
                        dy = dat(i,j,k,2) - dat(i,j-1,k,2)
                        dz = dat(i,j,k,3) - dat(i,j,k-1,3)
                        divb(i,j,k,1) = dx/delta(1) + dy/delta(2) + dz/delta(3)
            end do
         end do
      end do

  end subroutine ca_derdivb



#endif

end module derive_module
