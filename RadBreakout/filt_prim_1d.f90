subroutine ca_filt_prim(lo, hi, &
     Stmp, Stmp_l1, Stmp_h1, &
     Snew, Snew_l1, Snew_h1, &
     mask, mask_l1, mask_h1, &
     filt_T, S, domlo,domhi, &
     delta,xlo,problo,time,level)

  use probdata_module, only : filter_rhomax, filter_timemax
  use network, only : naux, nspec
  use eos_module, only : eos_given_RTX
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
       small_temp, small_dens, nadv
  use filter_module

  implicit none

  integer, intent(in) :: lo(1), hi(1), domlo(1), domhi(1), level
  integer, intent(in) :: filt_T, S
  double precision, intent(in) :: delta(1), xlo(1), problo(1), time
  integer, intent(in) ::   Stmp_l1,Stmp_h1
  integer, intent(in) ::   Snew_l1,Snew_h1
  integer, intent(in) ::   mask_l1,mask_h1
  double precision :: Stmp(Stmp_l1:Stmp_h1,NVAR)
  double precision :: Snew(Snew_l1:Snew_h1,NVAR)
  double precision :: mask(mask_l1:mask_h1)
  ! mask has three possible values: -1.d0, 0.d0, and 1.d0.
  ! -1.d0 appears only in cells that are covered by neither this level nor the finer level.
  !       It can only appear in ghost cells. 
  !  0.d0 appears only in cells that are covered by only this level, but not the finer level.
  !  1.d0 appears only in cells that are covered by the finer level.
  !       It can appear in either valid cells or ghost cells. 

  integer :: i
  double precision :: rhotmpinv, p, e, X(nspec+naux)
  logical, allocatable :: filtered(:)

  allocate(filtered(lo(1):hi(1)))
  filtered = .false.

  if (filt_T .eq. 1) then

     do i=lo(1), hi(1)
        if (time < filter_timemax .or. Snew(i,URHO) < filter_rhomax) then
           filtered(i) = .true.
           Snew(i,URHO ) = ff1(0) *  Stmp(i  ,URHO ) &
                &        + ff1(1) * (Stmp(i-1,URHO ) + Stmp(i+1,URHO ))           
           Snew(i,UMX  ) = ff1(0) *  Stmp(i  ,UMX  ) &
                &        + ff1(1) * (Stmp(i-1,UMX  ) + Stmp(i+1,UMX  ))
           Snew(i,UTEMP) = ff1(0) *  Stmp(i  ,UTEMP) &
                &        + ff1(1) * (Stmp(i-1,UTEMP) + Stmp(i+1,UTEMP))
        end if
     end do

  else if (filt_T .eq. 2) then

     do i=lo(1), hi(1)
        if (time < filter_timemax .or. Snew(i,URHO) < filter_rhomax) then
           filtered(i) = .true.
           Snew(i,URHO ) = ff2(0,S) *  Stmp(i  ,URHO ) &
                &        + ff2(1,S) * (Stmp(i-1,URHO ) + Stmp(i+1,URHO )) &
                &        + ff2(2,S) * (Stmp(i-2,URHO ) + Stmp(i+2,URHO ))
           Snew(i,UMX  ) = ff2(0,S) *  Stmp(i  ,UMX  ) &
                &        + ff2(1,S) * (Stmp(i-1,UMX  ) + Stmp(i+1,UMX  )) &
                &        + ff2(2,S) * (Stmp(i-2,UMX  ) + Stmp(i+2,UMX  ))
           Snew(i,UTEMP) = ff2(0,S) *  Stmp(i  ,UTEMP) &
                &        + ff2(1,S) * (Stmp(i-1,UTEMP) + Stmp(i+1,UTEMP)) &
                &        + ff2(2,S) * (Stmp(i-2,UTEMP) + Stmp(i+2,UTEMP))
        end if
     end do

  else if (filt_T .eq. 3) then

     do i=lo(1), hi(1)
        if (time < filter_timemax .or. Snew(i,URHO) < filter_rhomax) then
           filtered(i) = .true.
           Snew(i,URHO ) = ff3(0,S) *  Stmp(i  ,URHO ) &
                &        + ff3(1,S) * (Stmp(i-1,URHO ) + Stmp(i+1,URHO )) &
                &        + ff3(2,S) * (Stmp(i-2,URHO ) + Stmp(i+2,URHO )) &
                &        + ff3(3,S) * (Stmp(i-3,URHO ) + Stmp(i+3,URHO ))
           Snew(i,UMX  ) = ff3(0,S) *  Stmp(i  ,UMX  ) &
                &        + ff3(1,S) * (Stmp(i-1,UMX  ) + Stmp(i+1,UMX  )) &
                &        + ff3(2,S) * (Stmp(i-2,UMX  ) + Stmp(i+2,UMX  )) &
                &        + ff3(3,S) * (Stmp(i-3,UMX  ) + Stmp(i+3,UMX  ))
           Snew(i,UTEMP) = ff3(0,S) *  Stmp(i  ,UTEMP) &
                &        + ff3(1,S) * (Stmp(i-1,UTEMP) + Stmp(i+1,UTEMP)) &
                &        + ff3(2,S) * (Stmp(i-2,UTEMP) + Stmp(i+2,UTEMP)) &
                &        + ff3(3,S) * (Stmp(i-3,UTEMP) + Stmp(i+3,UTEMP))
        end if
     end do

  else if (filt_T .eq. 4) then

     do i=lo(1), hi(1)
        if (time < filter_timemax .or. Snew(i,URHO) < filter_rhomax) then
           filtered(i) = .true.
           Snew(i,URHO ) = ff4(0,S) *  Stmp(i  ,URHO ) &
                &        + ff4(1,S) * (Stmp(i-1,URHO ) + Stmp(i+1,URHO )) &
                &        + ff4(2,S) * (Stmp(i-2,URHO ) + Stmp(i+2,URHO )) &
                &        + ff4(3,S) * (Stmp(i-3,URHO ) + Stmp(i+3,URHO )) &
                &        + ff4(4,S) * (Stmp(i-4,URHO ) + Stmp(i+4,URHO ))
           Snew(i,UMX  ) = ff4(0,S) *  Stmp(i  ,UMX  ) &
                &        + ff4(1,S) * (Stmp(i-1,UMX  ) + Stmp(i+1,UMX  )) &
                &        + ff4(2,S) * (Stmp(i-2,UMX  ) + Stmp(i+2,UMX  )) &
                &        + ff4(3,S) * (Stmp(i-3,UMX  ) + Stmp(i+3,UMX  )) &
                &        + ff4(4,S) * (Stmp(i-4,UMX  ) + Stmp(i+4,UMX  ))
           Snew(i,UTEMP) = ff4(0,S) *  Stmp(i  ,UTEMP) &
                &        + ff4(1,S) * (Stmp(i-1,UTEMP) + Stmp(i+1,UTEMP)) &
                &        + ff4(2,S) * (Stmp(i-2,UTEMP) + Stmp(i+2,UTEMP)) &
                &        + ff4(3,S) * (Stmp(i-3,UTEMP) + Stmp(i+3,UTEMP)) &
                &        + ff4(4,S) * (Stmp(i-4,UTEMP) + Stmp(i+4,UTEMP))
        end if
     end do

  end if


  do i=lo(1), hi(1)
     if (filtered(i)) then

        Snew(i,URHO ) = max(Snew(i,URHO ), small_dens)
        Snew(i,UTEMP) = max(Snew(i,UTEMP), small_temp)
        
        rhotmpinv = 1.d0 / Stmp(i,URHO)
        
        X(1:nspec) = Stmp(i,UFS:UFS-1+nspec) * rhotmpinv
        
        if (naux > 0) &
             X(nspec+1:nspec+naux)  = Stmp(i,UFX:UFX+naux-1) * rhotmpinv
        
        call eos_given_RTX(e,p,Snew(i,URHO),Snew(i,UTEMP),X)
        
        Snew(i,UEINT) = Snew(i,URHO) * e
        Snew(i,UEDEN) = Snew(i,UEINT) + &
             0.5*(Snew(i,UMX)**2)/Snew(i,URHO) 
        if (nadv > 0) then
           Snew(i,UFA:UFA-1+nadv) = Stmp(i,UFA:UFA-1+nadv) * rhotmpinv * Snew(i,URHO)
        end if
        Snew(i,UFS:UFS-1+nspec) = X(1:nspec) * Snew(i,URHO)
        if (naux > 0) then
           Snew(i,UFX:UFX-1+naux) = X(nspec+1:) * Snew(i,URHO)
        end if
        
     end if
  end do

  deallocate(filtered)

end subroutine ca_filt_prim
