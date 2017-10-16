#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void find_ignited_zones(const int* lo, const int* hi,
                            BL_FORT_FAB_ARG_3D(ignition_radius),
                            BL_FORT_FAB_ARG_3D(state),
                            const amrex::Real* dx, int* num_zones, const int* level);
    
    void set_num_zones_ignited(const int* num_zones, const int* level);

    void get_num_zones_ignited(int* num_zones, int* level);
    
#ifdef __cplusplus
}
#endif

