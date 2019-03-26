module dg_amr_helper_module

  implicit none

contains

  subroutine dg_refine(fblo, fbhi, &
                       fine, fine_lo, fine_hi, &
                       crse, crse_lo, crse_hi, &
                       nFineX, nvar) &
                       bind(c,name='ca_dg_refine')

    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only: rt => amrex_real
    use ProgramHeaderModule           , only: nNodesX, nNodesE, nE, nDOF
    use RadiationFieldsModule         , only: nSpecies, nCR
    use TwoMoment_MeshRefinementModule, only: Refine_TwoMoment

    implicit none

    ! arguments
    integer :: crse_lo(3), crse_hi(3)
    integer :: fine_lo(3), fine_hi(3)
    integer :: fblo(3), fbhi(3)
    integer :: nFineX(3)
    integer :: nvar

    real(rt) :: fine(fine_lo(1):fine_hi(1), fine_lo(2):fine_hi(2), fine_lo(3):fine_hi(3), 0:nvar-1)
    real(rt) :: crse(crse_lo(1):crse_hi(1), crse_lo(2):crse_hi(2), crse_lo(3):crse_hi(3), 0:nvar-1)

    ! local variables
    integer :: ispec, imom, ic, jc
    integer :: icc, jcc, ifc, jfc
    integer :: id, ie, ii
    integer :: nFineX_local(3)
    integer :: cblo(2), cbhi(2)

    real(rt) :: U_FineCoarse(1:nDOF, 1:nE, 1:nFineX(1), 1:nFineX(2), 1:1)
    real(rt) :: U_CrseCoarse(1:nDOF, 1:nE)

    ! The C++ sends in zero but Refine_TwoMoment expects 1
    nFineX_local(1) = nFineX(1)
    nFineX_local(2) = nFineX(2)
    nFineX_local(3) = 1

    if ( nvar .ne. nNodesE * nE * nCR ) &
       call amrex_error("Mismatch of nvar with nNodesE*nE*nCR in dg_refine")
    if ( nDOF .ne. nNodesE * nNodesX(1) * nNodesX(2) ) &
       call amrex_error("Mismatch of nDOF with nNodesE*nNodesX(1)*nNodesX(2) in dg_refine")

    ! print *,'FINE DIMS ', fine_lo(:), fine_hi(:)
    ! print *,'CRSE DIMS ', crse_lo(:), crse_hi(:)
    ! print *,'FB LO/HI  ',fblo(:), fbhi(:)
 
    cblo(1) = fblo(1) / nFineX(1)
    cblo(2) = fblo(2) / nFineX(2)
    cbhi(1) = fbhi(1) / nFineX(1)
    cbhi(2) = fbhi(2) / nFineX(2)

    ! loop over coarse grid under fine region...
    do ic = cblo(1), cbhi(1), 2
       do jc = cblo(2), cbhi(2), 2

             do ispec = 1, nSpecies
                do imom = 1, nCR

                   do ie = 1, nE

                      do jcc = 1,2
                      do icc = 1,2

                          ifc = 2*(ic+icc-1)
                          jfc = 2*(jc+jcc-1)
   
                          do id = 1, nNodesE
                             ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-1)
                             U_CrseCoarse(id,ie) = crse(ic,jc,0,ii)
                          end do
   
                          do id = nNodesE+1, 2*nNodesE
                             ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic+1,jc,0,ii)
                          end do
   
                          do id = 2*nNodesE+1, 3*nNodesE
                             ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-2*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic,jc+1,0,ii)
                          end do
   
                          do id = 3*nNodesE+1, 4*nNodesE
                          ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-3*nNodesE-1)
                             U_CrseCoarse(id,ie) = crse(ic+1,jc+1,0,ii)
                          end do
   
                      end do ! icc
                      end do ! jcc
                   enddo ! ie

                   ! Coarsen U_FineCoarse --> U_CrseCoarse
                   call Refine_TwoMoment(nE, nFineX, U_CrseCoarse, U_FineCoarse)

                   do ie = 1, nE

                      do jcc = 1,2
                      do icc = 1,2

                         ifc = 2*(ic+icc-1)
                         jfc = 2*(jc+jcc-1)

                         do id = 1, nNodesE
                            ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-1)
                            fine(ifc,jfc,0,ii) = U_FineCoarse(id,ie,icc,jcc,1)
                         end do
         
                         do id = nNodesE+1, 2*nNodesE
                            ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-nNodesE-1)
                            fine(ifc+1,jfc,0,ii) = U_FineCoarse(id,ie,icc,jcc,1)
                         end do
      
                         do id = 2*nNodesE+1, 3*nNodesE
                            ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-2*nNodesE-1)
                            fine(ifc,jfc+1,0,ii) = U_FineCoarse(id,ie,icc,jcc,1)
                         end do
         
                         do id = 3*nNodesE+1, 4*nNodesE
                            ii   = (ispec-1)*(nCR*nE*nNodesE) + (imom-1)*(nE*nNodesE) + (ie-1)*nNodesE + (id-3*nNodesE-1)
                            fine(ifc+1,jfc+1,0,ii) = U_FineCoarse(id,ie,icc,jcc,1)
                         end do
                      end do
                      end do
                   enddo ! ie

             enddo ! imom
          enddo ! ispec

       enddo ! jc
    enddo ! ic

  end subroutine dg_refine


  subroutine dg_coarsen(cblo, cbhi, &
                        fine, fine_lo, fine_hi, &
                        crse, crse_lo, crse_hi, &
                        nFineX, nvar) &
                        bind(c,name='ca_dg_coarsen')

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

  end subroutine dg_coarsen
  
end module dg_amr_helper_module
