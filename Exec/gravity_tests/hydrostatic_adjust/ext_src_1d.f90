
  subroutine ca_ext_src(lo,hi, &
                        old_state,old_state_l1,old_state_h1, &
                        new_state,new_state_l1,new_state_h1, &
                        src,src_l1,src_h1,problo,dx,time,dt)

    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UFS
    use probdata_module   , only : xmin, &
                                   heating_time, heating_rad, &
                                   heating_peak, heating_sigma, prob_type
    use prob_params_module, only: center
    use network
    
    use amrex_fort_module, only : rt => c_real
    implicit none
    integer         ,intent(in   ) :: lo(1),hi(1)
    integer         ,intent(in   ) :: old_state_l1,old_state_h1
    integer         ,intent(in   ) :: new_state_l1,new_state_h1
    integer         ,intent(in   ) :: src_l1,src_h1
    real(rt)        ,intent(in   ) :: old_state(old_state_l1:old_state_h1,NVAR)
    real(rt)        ,intent(in   ) :: new_state(new_state_l1:new_state_h1,NVAR)
    real(rt)        ,intent(  out) :: src(src_l1:src_h1,NVAR)
    real(rt)        ,intent(in   ) :: problo(1),dx(1),time,dt
    
    integer          :: i
    real(rt)         :: x, r_0, H_0, W_0, Hext, t_stop

    integer :: ihe4_p

    r_0 = heating_rad

    H_0 = heating_peak
    W_0 = heating_sigma

    t_stop = heating_time

    if (prob_type == 1) then

       ! For heating at the center
       if (time .le. t_stop) then
          if (lo(1) .eq. 0) print *,'TIME vs TSTOP ',time, t_stop 
          do i = lo(1), hi(1)

             x = problo(1) + ((dble(i)+0.5e0_rt)*dx(1) + xmin) - center(1)

             Hext = H_0*exp(-((x-r_0)**2)/W_0**2)

             src(i,UEINT) = old_state(i,URHO)*Hext
             src(i,UEDEN) = old_state(i,URHO)*Hext
             
          enddo
       else

          src(lo(1):hi(1),UEINT) = 0.e0_rt
          src(lo(1):hi(1),UEDEN) = 0.e0_rt

       end if

    else if (prob_type == 3) then

       ihe4_p = network_species_index("helium-4")

       ! sub-chandra heating -- modulate by He
       if (time .le. t_stop) then
          if (lo(1) .eq. 0) print *,'TIME vs TSTOP ',time, t_stop 
          do i = lo(1), hi(1)

             x = problo(1) + ((dble(i)+0.5e0_rt)*dx(1) + xmin) - center(1)

             Hext = H_0*exp(-((x-r_0)**2)/W_0**2)*old_state(i,UFS-1+ihe4_p)/old_state(i,URHO)

             src(i,UEINT) = old_state(i,URHO)*Hext
             src(i,UEDEN) = old_state(i,URHO)*Hext
             
          enddo
       else

          src(lo(1):hi(1),UEINT) = 0.e0_rt
          src(lo(1):hi(1),UEDEN) = 0.e0_rt

       end if


    endif


  end subroutine ca_ext_src

