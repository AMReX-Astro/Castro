module probdata_module

    use amrex_fort_module, only : rt => amrex_real

    character (len=80), save  :: model_name

    real(rt), save :: min_density
    real(rt), save :: min_temperature
    real(rt), save :: fluff_ye

end module probdata_module
