module dg_amr_helper_module

  implicit none

contains

  subroutine dg_refine(fblo, fbhi, &
                       fine, fine_lo, fine_hi, &
                       crse, crse_lo, crse_hi, &
                       nFineX, &
                       nvar, n_dof, n_species, n_moments, n_energy_bins) bind(c,name='ca_dg_refine')

    use amrex_fort_module, only: rt => amrex_real
    use TwoMoment_MeshRefinementModule, only: Refine_TwoMoment

    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: fblo(3), fbhi(3)
    integer :: nFineX(3)
    integer :: nvar, n_dof, n_species, n_moments, n_energy_bins
    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: idof, ispec, imom, iebin, ic, jc, kc, icomp, ifx, ify, ifz, nfx, nfy, nfz
    real(rt) :: U_Coarse(1:n_dof, 1:n_energy_bins)
    real(rt) :: U_Fine(1:n_dof, 1:n_energy_bins, 1:nFineX(1), 1:nFineX(2), 1:nFineX(3))

    ! loop over coarse grid ...
    do ic = crse_lo(1), crse_hi(1)
       do jc = crse_lo(2), crse_hi(2)
          do kc = crse_lo(3), crse_hi(3)

             ! loop over species and moments ...
             do ispec = 1, n_species
                do imom = 1, n_moments

                   ! pack crse into U_Coarse
                   do iebin = 1, n_energy_bins
                      do idof = 1, n_dof
                         icomp = (ispec-1)*(n_moments*n_energy_bins*n_dof) + (imom-1)*(n_energy_bins*n_dof) + (iebin-1)*n_dof + (idof-1)
                         U_Coarse(idof, iebin) = crse(ic, jc, kc, icomp)
                      enddo
                   enddo

                   ! Derive U_Fine from U_Coarse
                   call Refine_TwoMoment(n_energy_bins, nFineX, U_Coarse, U_Fine)

                   ! for (ifine, jfine) in the U_Fine we obtained ...
                   ! unpack U_Fine into fine if (ifine, jfine) lies within (fblo, fbhi)
                   do iebin = 1, n_energy_bins
                      do idof = 1, n_dof
                         do nfx = 1, nFineX(1)
                            do nfy = 1, nFineX(2)
                               do nfz = 1, nFineX(3)
                                  ifx = fine_lo(1) + (nfx - 1) + nFineX(1) * (ic - crse_lo(1))
                                  ify = fine_lo(2) + (nfy - 1) + nFineX(2) * (jc - crse_lo(2))
                                  ifz = fine_lo(3) + (nfz - 1) + nFineX(3) * (kc - crse_lo(3))
                                  if (fblo(1) .le. ifx .and. ifx .le. fbhi(1) .and. &
                                      fblo(2) .le. ify .and. ify .le. fbhi(2) .and. &
                                      fblo(3) .le. ifz .and. ifz .le. fbhi(3)) then
                                     icomp = (ispec-1)*(n_moments*n_energy_bins*n_dof) + (imom-1)*(n_energy_bins*n_dof) + (iebin-1)*n_dof + (idof-1)
                                     fine(ifx, ify, ifz, icomp) = U_Fine(idof, iebin, nfx, nfy, nfz)
                                  endif
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                   ! ... end fine unpacking loop

                enddo
             enddo
             ! ... end species/moments loop

          enddo
       enddo
    enddo
    ! ... end coarse grid loop

  end subroutine dg_refine


  subroutine dg_coarsen(cblo, cbhi, &
                        fine, fine_lo, fine_hi, &
                        crse, crse_lo, crse_hi, &
                        nFineX, &
                        nvar, n_dof, n_species, n_moments, n_energy_bins) bind(c,name='ca_dg_coarsen')

    use amrex_fort_module, only: rt => amrex_real
    use TwoMoment_MeshRefinementModule, only: Coarsen_TwoMoment    
    
    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: cblo(3), cbhi(3)
    integer :: nFineX(3)
    integer :: nvar, n_dof, n_species, n_moments, n_energy_bins
    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: idof, ispec, imom, iebin, ic, jc, kc, icomp, ifx, ify, ifz, nfx, nfy, nfz
    real(rt) :: U_Coarse(1:n_dof, 1:n_energy_bins)
    real(rt) :: U_Fine(1:n_dof, 1:n_energy_bins, 1:nFineX(1), 1:nFineX(2), 1:nFineX(3))

    ! loop over coarse grid ...
    do ic = crse_lo(1), crse_hi(1)
       do jc = crse_lo(2), crse_hi(2)
          do kc = crse_lo(3), crse_hi(3)

             ! loop over species and moments ...
             do ispec = 1, n_species
                do imom = 1, n_moments

                   ! pack fine into U_Fine
                   do iebin = 1, n_energy_bins
                      do idof = 1, n_dof
                         do nfx = 1, nFineX(1)
                            do nfy = 1, nFineX(2)
                               do nfz = 1, nFineX(3)
                                  ifx = fine_lo(1) + (nfx - 1) + nFineX(1) * (ic - crse_lo(1))
                                  ify = fine_lo(2) + (nfy - 1) + nFineX(2) * (jc - crse_lo(2))
                                  ifz = fine_lo(3) + (nfz - 1) + nFineX(3) * (kc - crse_lo(3))
                                  icomp = (ispec-1)*(n_moments*n_energy_bins*n_dof) + (imom-1)*(n_energy_bins*n_dof) + (iebin-1)*n_dof + (idof-1)
                                  U_Fine(idof, iebin, nfx, nfy, nfz) = fine(ifx, ify, ifz, icomp)
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                   ! ... end fine packing loop
                   
                   ! Derive U_Coarse from U_Fine
                   call Coarsen_TwoMoment(n_energy_bins, nFineX, U_Coarse, U_Fine)

                   ! pack U_Coarse into crse
                   do iebin = 1, n_energy_bins
                      do idof = 1, n_dof
                         icomp = (ispec-1)*(n_moments*n_energy_bins*n_dof) + (imom-1)*(n_energy_bins*n_dof) + (iebin-1)*n_dof + (idof-1)
                         crse(ic, jc, kc, icomp) = U_Coarse(idof, iebin)
                      enddo
                   enddo
                   ! ... end crse packing loop

                enddo
             enddo
             ! ... end species/moments loop

          enddo
       enddo
    enddo
    ! ... end coarse grid loop

  end subroutine dg_coarsen
  
end module dg_amr_helper_module
