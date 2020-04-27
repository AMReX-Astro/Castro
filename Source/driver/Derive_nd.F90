module derive_module
    ! All subroutines in this file must be threadsafe because they are called
    ! inside OpenMP parallel regions.

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains


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
