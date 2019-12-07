module probdata_module

    use amrex_fort_module, only : rt => amrex_real

    character (len=80), save  :: model_name

    ! Initial model mapping parameters
    real(rt), save :: min_density
    real(rt), save :: min_temperature
    real(rt), save :: fluff_ye

    ! Tagging parameters
    integer,  save :: max_base_tagging_level
    real(rt), save :: tag_density
    real(rt), save :: tag_max_density_fraction

end module probdata_module
