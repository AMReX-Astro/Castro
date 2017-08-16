module interpolate_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  contains

      function interpolate(r, npts_model, model_r, model_var, iloc)
      
!     given the array of model coordinates (model_r), and variable (model_var),
!     find the value of model_var at point r using linear interpolation.
!     Eventually, we can do something fancier here.
      
      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , intent(in   ) :: r
      integer         , intent(in   ) :: npts_model
      real(rt)        , intent(in   ) :: model_r(npts_model), model_var(npts_model)
      integer, intent(in), optional   :: iloc
      real(rt)                        :: interpolate

      ! Local variables
      integer                         :: id
      real(rt)                        :: slope,minvar,maxvar
      
!     find the location in the coordinate array where we want to interpolate
      if (present(iloc)) then
         id = iloc
      else
         id = locate(r, npts_model, model_r)
      end if
      
      if (id .eq. 1) then

         slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
         interpolate = slope*(r - model_r(id)) + model_var(id)

         ! safety check to make sure interpolate lies within the bounding points
         minvar = min(model_var(id+1),model_var(id))
         maxvar = max(model_var(id+1),model_var(id))
         interpolate = max(interpolate,minvar)
         interpolate = min(interpolate,maxvar)
         
      else if (id .eq. npts_model) then

         slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
         interpolate = slope*(r - model_r(id)) + model_var(id)

         ! safety check to make sure interpolate lies within the bounding points
         minvar = min(model_var(id),model_var(id-1))
         maxvar = max(model_var(id),model_var(id-1))
         interpolate = max(interpolate,minvar)
         interpolate = min(interpolate,maxvar)

      else

         if (r .ge. model_r(id)) then

            slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
            interpolate = slope*(r - model_r(id)) + model_var(id)

            ! ! safety check to make sure interpolate lies within the bounding points
            ! minvar = min(model_var(id+1),model_var(id))
            ! maxvar = max(model_var(id+1),model_var(id))
            ! interpolate = max(interpolate,minvar)
            ! interpolate = min(interpolate,maxvar)

         else

            slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
            interpolate = slope*(r - model_r(id)) + model_var(id)

            ! ! safety check to make sure interpolate lies within the bounding points
            ! minvar = min(model_var(id),model_var(id-1))
            ! maxvar = max(model_var(id),model_var(id-1))
            ! interpolate = max(interpolate,minvar)
            ! interpolate = min(interpolate,maxvar)

          end if

      endif

    end function interpolate

    subroutine tri_interpolate(x, y, z, npts_x, npts_y, npts_z, &
                               model_x, model_y, model_z, model_var, &
                               interp_var, derivs, error)

      use bl_error_module
      use bl_constants_module, only: ONE

      ! tri-linear interpolation; useful for EOS tables
      ! this is stricly interpolation, so if the point (x,y,z) is outside
      ! the bounds of model_x,model_y,model_z, then we abort
      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , intent(in   ) :: x,y,z
      integer,          intent(in   ) :: npts_x, npts_y, npts_z
      real(rt)        , intent(in   ) :: model_x(npts_x), &
                                         model_y(npts_y), &
                                         model_z(npts_z), &
                                         model_var(npts_x,npts_y,npts_z)
      real(rt)        , intent(  out) :: interp_var, derivs(3)
      logical,          intent(  out) :: error

      integer :: ix,iy,iz
      real(rt)         :: deltax, deltay, deltaz
      real(rt)         :: c(8), delta(8)

      error = .false.

      ! find the indices below the point (x,y,z)
      do ix = 2, npts_x
         if (model_x(ix) .ge. x) exit
      enddo
      do iy = 2, npts_y
         if (model_y(iy) .ge. y) exit
      enddo
      do iz = 2, npts_z
         if (model_z(iz) .ge. z) exit
      enddo
      if ((ix>npts_x) .or. (iy>npts_y) .or. (iz>npts_z)) then !&
         error = .true.
         return
      endif

      ix=ix-1
      iy=iy-1
      iz=iz-1

      ! form the weights for each point
      deltax = (x - model_x(ix))/(model_x(ix+1)-model_x(ix))
      deltay = (y - model_y(iy))/(model_y(iy+1)-model_y(iy))
      deltaz = (z - model_z(iz))/(model_z(iz+1)-model_z(iz))
      delta = (/ ONE, deltax, deltay, deltaz, &
                 deltax*deltay, deltax*deltaz, deltay*deltaz, &
                 deltax*deltay*deltaz /)

      ! the model_var function values for Lagrange interpolation
      c(1) = model_var(ix  ,iy  ,iz  )
      c(2) = model_var(ix+1,iy  ,iz  ) - model_var(ix  ,iy  ,iz  )
      c(3) = model_var(ix  ,iy+1,iz  ) - model_var(ix  ,iy  ,iz  )
      c(4) = model_var(ix  ,iy  ,iz+1) - model_var(ix  ,iy  ,iz  )
      c(5) = model_var(ix+1,iy+1,iz  ) - model_var(ix  ,iy+1,iz  ) + &
             model_var(ix  ,iy  ,iz  ) - model_var(ix+1,iy  ,iz  )
      c(6) = model_var(ix+1,iy  ,iz+1) - model_var(ix  ,iy  ,iz+1) + &
             model_var(ix  ,iy  ,iz  ) - model_var(ix+1,iy  ,iz  )
      c(7) = model_var(ix  ,iy+1,iz+1) - model_var(ix  ,iy  ,iz+1) + &
             model_var(ix  ,iy  ,iz  ) - model_var(ix  ,iy+1,iz  )
      c(8) = model_var(ix+1,iy+1,iz+1) - model_var(ix  ,iy+1,iz+1) + &
             model_var(ix  ,iy  ,iz+1) - model_var(ix+1,iy  ,iz+1) + &
             model_var(ix  ,iy+1,iz  ) - model_var(ix+1,iy+1,iz  ) + &
             model_var(ix+1,iy  ,iz  ) - model_var(ix  ,iy  ,iz  )

      ! interpolated value
      interp_var = sum(c*delta)
      ! derivs: dvar/dx, dvar/dy, dvar/dz
      derivs(1) = c(2) / (model_x(ix+1)-model_x(ix))
      derivs(2) = c(3) / (model_y(iy+1)-model_y(iy))
      derivs(3) = c(4) / (model_z(iz+1)-model_z(iz))
      
    end subroutine tri_interpolate


    function locate(x, n, xs)
      use amrex_fort_module, only : rt => amrex_real
      integer, intent(in) :: n
      real(rt)        , intent(in) :: x, xs(n)
      integer :: locate      

      integer :: ilo, ihi, imid

      if (x .le. xs(1)) then
         locate = 1
      else if (x .gt. xs(n-1)) then
         locate = n
      else

         ilo = 1
         ihi = n-1

         do while (ilo+1 .ne. ihi)
            imid = (ilo+ihi)/2
            if (x .le. xs(imid)) then
               ihi = imid
            else 
               ilo = imid
            end if
         end do
         
         locate = ihi

      end if

    end function locate

end module interpolate_module
