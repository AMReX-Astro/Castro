
      subroutine ca_react_state(lo,hi, &
                                s_in ,so_l1,so_l2,so_l3,so_h1,so_h2,so_h3, &
                                s_out,sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                reaction_terms,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3, &
                                time,dt_react)

      use eos_module
      use network           , only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, small_dens, small_temp, allow_negative_energy
      use castro_burner_module
      use bl_constants_module

      implicit none

      integer lo(3),hi(3)
      integer so_l1,so_h1,so_l2,so_h2,so_l3,so_h3
      integer sn_l1,sn_h1,sn_l2,sn_h2,sn_l3,sn_h3
      integer  r_l1, r_h1, r_l2, r_h2, r_l3, r_h3
      double precision s_in (so_l1:so_h1,so_l2:so_h2,so_l3:so_h3,NVAR)
      double precision s_out(sn_l1:sn_h1,sn_l2:sn_h2,sn_l3:sn_h3,NVAR)
      double precision reaction_terms(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3,nspec+2)
      double precision time, dt_react

      integer          :: i,j,k
      integer          :: pt_index(3)
      double precision :: rho, rhoInv, u, v, w, ke, e_in, e_out, T
      double precision :: x_in(nspec+naux), x_out(nspec+naux)

      type (eos_t) :: eos_state

      !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index,rho,rhoInv,u,v,w,ke,T,x_in,e_in,x_out,e_out)
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

        if (s_in(i,j,k,URHO) .gt. small_dens) then

           ! Define T from e
           rho           = s_in(i,j,k,URHO)
           rhoInv        = ONE / rho
           u             = s_in(i,j,k,UMX)*rhoInv
           v             = s_in(i,j,k,UMY)*rhoInv
           w             = s_in(i,j,k,UMZ)*rhoInv
           ke            = HALF * (u**2 + v**2 + w**2)
           T             = s_in(i,j,k,UTEMP)
           x_in(1:nspec) = s_in(i,j,k,UFS:UFS+nspec-1) * rhoInv
           if (naux > 0) &
             x_in(nspec+1:nspec+naux)  = s_in(i,j,k,UFX:UFX+naux-1) * rhoInv

           e_in          = s_in(i,j,k,UEINT) * rhoInv

           pt_index(1) = i
           pt_index(2) = j
           pt_index(3) = k

           eos_state % T   = T
           eos_state % rho = rho
           eos_state % e   = e_in
           eos_state % xn  = x_in(1:nspec)
           eos_state % aux = s_in(nspec+1:nspec+naux)

           if (allow_negative_energy .eq. 0 .and. e_in .le. ZERO) then
              print *,'... e negative in react_state: ', i, j, k, e_in
              T = max(T, small_temp)
              eos_state % T = T
              call eos(eos_input_rt, eos_state, pt_index = pt_index)
              e_in = eos_state % e
              if (e_in .lt. ZERO) then
                 print *,'... call to eos (input_rt) with small_temp still gives neg e ', e_in
                 call bl_error("Error:: React_3d.f90 :: ca_react_state")
              else
                 print *,'... able to re-set using eos (input_rt) with small_temp ', e_in
              end if
           end if

           ! Use this call to define T

           call eos(eos_input_re, eos_state, pt_index = pt_index)

           T = eos_state % T
           e = eos_state % e

           call burner(rho, T, x_in, e_in, dt_react, time, x_out, e_out)

           ! Make sure that species emerge in the proper range: [0,1]
           do n = 1, nspec
             x_out(n) = max(min(x_out(n),ONE),ZERO)
           end do

           if (i.ge.r_l1 .and. i.le.r_h1 .and. j.ge.r_l2 .and. j.le.r_h2 .and. &
               k.ge.r_l3 .and. k.le.r_h3) then
               reaction_terms(i,j,k,1:nspec) = reaction_terms(i,j,k,1:nspec) &
                + (x_out(1:nspec) - x_in(1:nspec))
               reaction_terms(i,j,k,nspec+1) = reaction_terms(i,j,k,nspec+1) &
               +     (e_out - e_in)
               reaction_terms(i,j,k,nspec+2) = reaction_terms(i,j,k,nspec+2) &
               +      rho*(e_out - e_in)
           end if

           s_out(i,j,k,URHO)            = rho
           s_out(i,j,k,UEINT)           = rho * e_out
           s_out(i,j,k,UEDEN)           = rho * (e_out + ke)
           s_out(i,j,k,UFS:UFS+nspec-1) = rho * x_out(1:nspec)
           s_out(i,j,k,UFX:UFX+naux -1) = s_in(i,j,k,UFX:UFX+naux-1)
   
           ! Now update the temperature to match the new internal energy

           eos_state % rho = rho
           eos_state % e   = e_out
           eos_state % xn  = x_out(1:nspec)
           eos_state % aux = x_out(nspec+1:nspec+naux)

           call eos(eos_input_re, eos_state, pt_index = pt_index)

           s_out(i,j,k,UTEMP)           = eos_state % T

        else

           s_out(i,j,k,1:NVAR)          = s_in(i,j,k,1:NVAR)

        end if

      end do
      end do
      end do
      !$OMP END PARALLEL DO

      end subroutine ca_react_state
