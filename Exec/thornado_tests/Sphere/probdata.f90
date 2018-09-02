module probdata_module

      use amrex_fort_module, only : rt => amrex_real
      use network, only: nspec

      real(rt), save :: xn(nspec)
      real(rt), save :: rho_i, T_i, Ye_i, rhoe_i, p_i
      real(rt), save :: centx, centy

end module probdata_module
