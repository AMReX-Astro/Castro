      subroutine compute_temp(lo,hi,state,state_l1,state_l2,state_l3, &
                              state_h1,state_h2,state_h3)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, UMX, UMY, UMZ, allow_negative_energy
      use bl_constants_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: state_l1,state_l2,state_l3
      integer         , intent(in   ) :: state_h1,state_h2,state_h3
      double precision, intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,&
                                               state_l3:state_h3,NVAR)

      integer          :: i,j,k,n
      double precision :: rhoInv
      integer          :: pt_index(3)

      type (eos_t), allocatable :: eos_state(:)
      integer :: eos_state_len(1), nx(3), n
      
      nx = hi - lo + 1
      eos_state_len(1) = nx(1) * nx(2) * nx(3)
      allocate(eos_state(eos_state_len(1)))

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        if (state(i,j,k,URHO) <= ZERO) then
           print *,'   '
           print *,'>>> Error: Castro_3d::compute_temp ',i,j,k
           print *,'>>> ... negative density ',state(i,j,k,URHO)
           print *,'    '
           call bl_error("Error:: Castro_3d.f90 :: compute_temp")
        end if
      enddo
      enddo
      enddo

      if (allow_negative_energy.eq.0) then
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (state(i,j,k,UEINT) <= ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_3d::compute_temp ',i,j,k
                print *,'>>> ... (rho e) is negative '
                call bl_error("Error:: Castro_3d.f90 :: compute_temp")
            end if
         enddo
         enddo
         enddo
      end if

      eos_state(:) % rho = reshape(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),URHO ), eos_state_len)
      eos_state(:) % T   = reshape(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UTEMP), eos_state_len) ! Initial guess for the EOS
      eos_state(:) % e   = reshape(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UEINT), eos_state_len) / eos_state(:) % rho
      do n = 1, nspec
         eos_state(:) % xn(n)  = reshape(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UFS+n-1), eos_state_len) / eos_state(:) % rho
      enddo
      do n = 1, naux
         eos_state(:) % aux(n) = reshape(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UFX+n-1), eos_state_len) / eos_state(:) % rho
      enddo

      call eos(eos_input_re, eos_state, state_len = eos_state_len(1))

      !$OMP PARALLEL WORKSHARE

      state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UTEMP) = reshape(eos_state(:) % T, nx)

      ! Reset energy in case we floored
      state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UEINT) = state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),URHO) * &
                                                         reshape(eos_state(:) % e, nx)
      state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UEDEN) = state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UEINT) + HALF * &
                                                         (state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UMX)**2 + &
                                                          state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UMY)**2 + &
                                                          state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),UMZ)**2) / &
                                                          state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),URHO)
      
      !$OMP END PARALLEL WORKSHARE

      end subroutine compute_temp
