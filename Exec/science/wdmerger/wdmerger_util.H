#ifndef wdmerger_util_H
#define wdmerger_util_H

#include <math.H>

// If we are in a rotating reference frame, then rotate a vector
// by an amount corresponding to the time that has passed
// since the beginning of the simulation.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
GpuArray<Real, 3> inertial_rotation(const GpuArray<Real, 3>& vec, const GpuArray<Real, 3>& omega, Real time)
{
    GpuArray<Real, 3> vec_i{};

    Array1D<Real, 0, 2> theta{};
    Array2D<Real, 0, 2, 0, 2> rot_matrix{};

    if (castro::do_rotation == 1) {
        for (int n = 0; n < 3; ++n) {
            theta(n) = omega[n] * time;
        }
    }

    // This is the 3D rotation matrix for converting between reference frames.
    // It is the composition of rotations along the x, y, and z axes. Therefore 
    // it allows for the case where we are rotating about multiple axes. Normally 
    // we use the right-hand convention for constructing the usual rotation matrix, 
    // but this is the transpose of that rotation matrix to account for the fact 
    // that we are rotating *back* to the inertial frame, rather than from the 
    // inertial frame to the rotating frame.

    rot_matrix(0,0) =  cos(theta(1)) * cos(theta(2));
    rot_matrix(0,1) = -cos(theta(1)) * sin(theta(2));
    rot_matrix(0,2) =  sin(theta(1));
    rot_matrix(1,0) =  cos(theta(0)) * sin(theta(2)) + sin(theta(0)) * sin(theta(1)) * cos(theta(2));
    rot_matrix(1,1) =  cos(theta(0)) * cos(theta(2)) - sin(theta(0)) * sin(theta(1)) * sin(theta(2));
    rot_matrix(1,2) = -sin(theta(0)) * cos(theta(1));
    rot_matrix(2,0) =  sin(theta(0)) * sin(theta(2)) - cos(theta(0)) * sin(theta(1)) * cos(theta(2));
    rot_matrix(2,1) =  sin(theta(0)) * cos(theta(2)) + cos(theta(0)) * sin(theta(1)) * sin(theta(2));
    rot_matrix(2,2) =  cos(theta(0)) * cos(theta(1));

    for (int l = 0; l < 3; ++l) {
        for (int m = 0; m < 3; ++m) {
            vec_i[l] += rot_matrix(l,m) * vec[m];
        }
    }

    return vec_i;
}

// Given a rotating frame velocity, get the inertial frame velocity.
// Note that we simply return the original velocity if we're
// already in the inertial frame.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
GpuArray<Real, 3> inertial_velocity(const GpuArray<Real, 3>& loc, const GpuArray<Real, 3>& vel, const GpuArray<Real, 3>& omega)
{
    GpuArray<Real, 3> vel_i{};

    if (castro::do_rotation == 1 && castro::state_in_rotating_frame == 1) {
        cross_product(omega, loc, vel_i);
    }

    for (int n = 0; n < 3; ++n) {
        vel_i[n] += vel[n];
    }

    return vel_i;
}

#endif
