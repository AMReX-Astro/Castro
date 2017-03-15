module probdata_module

!     These determine the refinement criteria
      use bl_fort_module, only : rt => c_real
      character (len=80), save  :: model_name
      
      real(rt)    , allocatable, save :: rho_0(:)
      real(rt)    , allocatable, save :: T_0(:)
      integer       , save :: npts
end module probdata_module
