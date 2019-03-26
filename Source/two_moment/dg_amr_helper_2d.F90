module dg_amr_helper_module

  implicit none

contains

  subroutine dg_refine(fblo, fbhi, &
                       fine, fine_lo, fine_hi, &
                       crse, crse_lo, crse_hi, &
                       nFineX, &
                       nvar, n_moments, n_energy_bins) bind(c,name='ca_dg_refine')

    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only: rt => amrex_real
    use ProgramHeaderModule           , only: nNodesE, nE, nDOF
    use RadiationFieldsModule         , only: nSpecies
    use TwoMoment_MeshRefinementModule, only: Refine_TwoMoment

    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: fblo(3), fbhi(3)
    integer :: nFineX(3)
    integer :: nvar, n_moments, n_energy_bins
    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: idof, ispec, imom, ie, ic, jc, kc, icomp, ifx, ify, ifz, nfx, nfy, nfz
    real(rt) :: U_Coarse(1:nDOF, 1:n_energy_bins)

    real(rt) :: U_Fine(1:nDOF, 1:n_energy_bins, 1:nFineX(1), 1:nFineX(2), 1:1)

    ! loop over coarse grid ...
    do ic = crse_lo(1), crse_hi(1)
       do jc = crse_lo(2), crse_hi(2)
          do kc = crse_lo(3), crse_hi(3)

             ! loop over species and moments ...
             do ispec = 1, nSpecies
                do imom = 1, n_moments

                   ! pack crse into U_Coarse
                   do ie = 1, n_energy_bins
                      do idof = 1, nDOF
                         icomp = (ispec-1)*(n_moments*n_energy_bins*nDOF) + (imom-1)*(n_energy_bins*nDOF) + (ie-1)*nDOF + (idof-1)
                         U_Coarse(idof, ie) = crse(ic, jc, kc, icomp)
                      enddo
                   enddo

                   ! Derive U_Fine from U_Coarse
                   call Refine_TwoMoment(n_energy_bins, nFineX, U_Coarse, U_Fine)

                   ! for (ifine, jfine) in the U_Fine we obtained ...
                   ! unpack U_Fine into fine if (ifine, jfine) lies within (fblo, fbhi)
                   do ie = 1, n_energy_bins
                      do idof = 1, nDOF
                         do nfx = 1, nFineX(1)
                            do nfy = 1, nFineX(2)
                               do nfz = 1, 1
                                  ifx = fblo(1) + (nfx-1) + nFineX(1) * (ic - crse_lo(1))
                                  ify = fblo(2) + (nfy-1) + nFineX(2) * (jc - crse_lo(2))
                                  ifz = fblo(3) + (nfz-1) + nFineX(3) * (kc - crse_lo(3))
                                  if (fblo(1) .le. ifx .and. ifx .le. fbhi(1) .and. &
                                      fblo(2) .le. ify .and. ify .le. fbhi(2) .and. &
                                      fblo(3) .le. ifz .and. ifz .le. fbhi(3)) then
                                     icomp = (ispec-1)*(n_moments*n_energy_bins*nDOF) + (imom-1)*(n_energy_bins*nDOF) + (ie-1)*nDOF + (idof-1)
                                     fine(ifx, ify, ifz, icomp) = U_Fine(idof, ie, nfx, nfy, nfz)
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
                        nvar) bind(c,name='ca_dg_coarsen')

    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only: rt => amrex_real
    use ProgramHeaderModule           , only: nNodesX, nNodesE, nE, nDOF
    use RadiationFieldsModule         , only: nSpecies, nCR
    use TwoMoment_MeshRefinementModule, only: Coarsen_TwoMoment    
    
    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: cblo(3), cbhi(3)
    integer :: nFineX(3)
    integer :: nvar

    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: ispec, imom, ic, jc
    integer :: icc, jcc, ifc, jfc
    integer :: id, ie, ii
    integer :: nFineX_local(3)

    real(rt) :: U_FineCoarse(1:nDOF, 1:nE, 1:nFineX(1), 1:nFineX(2), 1:1)
    real(rt) :: U_CrseCoarse(1:nDOF, 1:nE)

    ! The C++ sends in zero but Coarsen_TwoMoment expects 1
    nFineX_local(1) = nFineX(1)
    nFineX_local(2) = nFineX(2)
    nFineX_local(3) = 1

    if ( nvar .ne. nNodesE * nE * nCR ) &
       call amrex_error("Mismatch of nvar with nNodesE*nE*nCR in dg_coarsen")
    if ( nDOF .ne. nNodesE * nNodesX(1) * nNodesX(2) ) &
       call amrex_error("Mismatch of nDOF with nNodesE*nNodesX(1)*nNodesX(2) in dg_coarsen")

    ! loop over coarse grid -- note that we must skip here because we need to construct an intermediate
    !      array at half of this resolution
    do ic = cblo(1), cbhi(1), 2
       do jc = cblo(2), cbhi(2), 2

          ! loop over species and moments ...
          do ispec = 1, nSpecies
             do imom = 1, nCR
                do ie = 1, nE

                   do jcc = 1,2
                   do icc = 1,2

                       ifc = 2*(ic+icc-1)
                       jfc = 2*(jc+jcc-1)

                       do id = 1, nNodesE
                          ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-1)
                          U_FineCoarse(id,ie,icc,jcc,1) = fine(ifc,jfc,0,ii)
                       end do

                       do id = nNodesE+1, 2*nNodesE
                          ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,1) = fine(ifc+1,jfc,0,ii)
                       end do

                       do id = 2*nNodesE+1, 3*nNodesE
                          ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-2*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,1) = fine(ifc,jfc+1,0,ii)
                       end do

                       do id = 3*nNodesE+1, 4*nNodesE
                          ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-3*nNodesE-1)
                          U_FineCoarse(id,ie,icc,jcc,1) = fine(ifc+1,jfc+1,0,ii)
                       end do

                   end do ! icc
                   end do ! jcc
                enddo ! ie
                   
                ! Coarsen U_FineCoarse --> U_CrseCoarse

                call Coarsen_TwoMoment(nE, nFineX_local, U_CrseCoarse, U_FineCoarse)

                do ie = 1, nE
                   do id = 1, nNodesE
                      ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-1)
                      crse(ic,jc,0,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = nNodesE+1, 2*nNodesE
                      ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-nNodesE-1)
                      crse(ic+1,jc,0,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 2*nNodesE+1, 3*nNodesE
                      ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-2*nNodesE-1)
                      crse(ic,jc+1,0,ii) = U_CrseCoarse(id,ie)
                   end do
   
                   do id = 3*nNodesE+1, 4*nNodesE
                      ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-3*nNodesE-1)
                      crse(ic+1,jc+1,0,ii) = U_CrseCoarse(id,ie)
                   end do
                enddo ! ie

             enddo ! imom
          enddo ! ispec

       enddo ! jc
    enddo ! ic
    ! ... end coarse grid loop

  end subroutine dg_coarsen
  
end module dg_amr_helper_module
