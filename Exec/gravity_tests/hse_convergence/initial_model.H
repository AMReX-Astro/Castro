// create a simple hydrostatic initial model that is 4th order accurate
// at cell-centers by using RK4 to integrate:
//
//   dp/dr = rho(p, s) g
//   ds/dr = 0


struct model_t {
    Real dens_base = -1;
    Real T_base = -1;
    Real xn[NumSpec] = {0.0};
}


///
/// construct an initial model in HSE.  Note: this does not return
/// anything, but rather updates the model_parser globals with the
/// model information.
///
AMREX_INLINE
void
generate_initial_model(const int nx, const Real xmin, const Real xmax,
                       const model_t model_params, const int nbuf) {

    constexpr int MAX_ITER = 100;
    constexpr Real TOL = 1.e-10_rt;

    // allocate the storage in the model_parser_module
    int npts_model = nx + 2*nbuf;

    allocate_model_parser_data(npts_model);

    int ibase = nbuf;
    int itop = ibase + nx - 1;

    // get the base conditions

    eos_t eos_state;
    eos_state.rho = model_params.dens_base;
    eos_state.T = model_params.T_base;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model_params.xn[n];
    }

    eos(eos_input_rt, eos_state);

    Real entropy_fixed = eos_state.s;

    // create the grid -- cell centers

    Real dx = (xmax - xmin) / nx;

    for (int i = 0; i < npts_model; i++) {
        model::r[i] = xmin + (static_cast<Real>(i - nbuf) + 0.5_rt) * dx;
    }

    // note, those conditions are the lower boundary.  This means we
    // need to integrate dx/2 in the first step to get to the first
    // zone's cell-center.  After that, we integrate from
    // center-to-center.

    Real p = eos_state.p;

    if (time_integration_method == 2 && sdc_order == 4) {

        // a fourth order accurate method

        // do RK 4 integration up from the lower boundary

        for (int i = ibase; i <= itop + nbuf; i++) {

            // rho and T here are guesses for the EOS call

            Real rho;
            Real T;

            if (i == ibase) {
                rho = eos_state.rho;
                T = eos_state.T;
            } else {
                rho = model::state(MODELINDEX(i-1, model::idens));
                T = model::state(MODELINDEX(i-1, model::itemp));
            }

            Real xn[NumSpec];
            for (int n = 0; n < NumSpec; n++) {
                xn[n] = model_params.xn[n];
            }

            // step size

            Real h;
            if (i == ibase) {
                h = 0.5_rt * dx;
            } else {
                h = dx;
            }

            // entropy never changes in this model

            Real s = entropy_fixed;

            Real k1 = f(p, s, const_grav, rho, T, xn);
            Real k2 = f(p + 0.5_rt * h * k1, s, const_grav, rho, T, xn);
            Real k3 = f(p + 0.5_rt * h * k2, s, const_grav, rho, T, xn);
            Real k4 = f(p + h * k3, s, const_grav, rho, T, xn);

            Real pnew = p + (1.0_rt/6.0_rt) * h * (k1 + 2.0_rt * k2 + 2.0_rt * k3 + k4);

            // call the EOS to get the remainder of the thermodynamics

            eos_t eos_state;
            eos_state.T = T; // initial guess
            eos_state.rho = rho; // initial guess
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = model_params.xn[n];
            }
            eos_state.p = pnew;
            eos_state.s = s;

            eos(eos_input_ps, eos_state);

            // update the thermodynamics in this zone

            model::state(MODELINDEX(i, idens_model)) = eos_state.rho;
            model::state(MODELINDEX(i, itemp_model)) = eos_state.T;
            model::state(MODELINDEX(i, ipres_model)) = eos_state.p;
            for (int n = 0; n < NumSpec; n++) {
                model::state(MODELINDEX(i, ispec_model+n)) = eos_state.xn[n];
            }

            // reset for the next iteration

            p = pnew;

        }

        // now integrate down

        p = model::state(MODELINDEX(ibase, ipres_model));

        for (int i = ibase-1; i >= 0, --i) {

            Real rho = model::state(MODELINDEX(i+1, idens_model));
            Real T = model::state(MODELINDEX(i+1, itemp_model));
            for (int n = 0; n < NumSpec; n++) {
                xn[n] = model_params.xn[n];
            }

            Real h = -dx;

            // entropy never changes in this model

            Real s = entropy_fixed;

            Real k1 = f(p, s, const_grav, rho, T, xn);
            Real k2 = f(p + 0.5_rt * h * k1, s, const_grav, rho, T, xn);
            Real k3 = f(p + 0.5_rt * h * k2, s, const_grav, rho, T, xn);
            Real k4 = f(p + h * k3, s, const_grav, rho, T, xn);

            Real pnew = p + (1.0_rt/6.0_rt) * h * (k1 + 2.0_rt*k2 + 2.0_rt*k3 + k4);

            //! call the EOS to get the remainder of the thermodynamics
            eos_t eos_state;
            eos_state.T = T; // initial guess
            eos_state.rho = rho; // initial guess
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = model_params.xn[n];
            }
            eos_state.p = pnew;
            eos_state.s = s;

            eos(eos_input_ps, eos_state);

            // update the thermodynamics in this zone

            model::state(MODELINDEX(i, idens_model)) = eos_state.rho;
            model::state(MODELINDEX(i, itemp_model)) = eos_state.T;
            model::state(MODELINDEX(i, ipres_model)) = eos_state.p;
            for (int n = 0; n < NumSpec; n++) {
                model::state(MODELINDEX(i, ispec_model+n)) = eos_state.xn[n];
            }

            // reset for the next iteration

            p = pnew;

        }

    } else {

        //  a second order accurate scheme

        for (int i = ibase; i <= itop + nbuf; i++) {

            if (i == ibase) {
                model_state(MODELINDEX(i, idens_model)) = model_params.dens_base; // initial guess
                model_state(MODELINDEX(i, ipres_model)) =
                    p + 0.5_rt * dx * model_params.dens_base * const_grav;

                eos_t eos_state;
                eos_state.rho = model_state(MODELINDEX(i, idens_model));
                eos_state.T = model_params.dens_base;
                eos_state.p = model_state(MODELINDEX(i, ipres_model));
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = model_params.xn[n];
                }
                eos_state.s = entropy_fixed;

                eos(eos_input_ps, eos_state);

                model_state(MODELINDEX(i, idens_model)) = eos_state.rho;
                model_state(MODELINDEX(i, itemp_model)) = eos_state.T;
                for (int n = 0; n < NumSpec; n++) {
                    model::state(MODELINDEX(i, ispec_model+n)) = model_params.xn[n];
                }

            } else {

                Real dens_zone = model_state(MODELINDEX(i-1, idens_model));

                bool converged_hse = false;

                for (int iter = 0; iter < MAX_ITER; iter++) {

                    Real p_want = model_state(MODELINDEX(i-1, ipres_model)) +
                        dx * 0.5_rt * (dens_zone + model_state(MODELINDEX(i-1, idens_model))) * const_grav;

                    // use the EOS with constant entropy to find corrected state

                    eos_t eos_state;
                    eos_state.p = p_want;
                    eos_state.s = entropy_fixed;
                    for (int n = 0; n < NumSpec; n++) {
                        eos_state.xn[n] = model_params.xn[n];
                    }

                    eos(eos_input_ps, eos_state);

                    Real drho = eos_state.rho - dens_zone;
                    dens_zone = eos_state.rho;
                    if (std::abs(drho) < TOL * dens_zone) {
                        converged_hse = true;
                        break;
                    }

                }

                if (!converged_hse) {
                    amrex::Error("ERROR: failed to convergence in initial model generation");
                }

            }

            // initialze zone

            model_state(MODELINDEX(i, model::idens)) = dens_zone;
            model_state(MODELINDEX(i, model::ipres)) = p_want;
            model_state(MODELINDEX(i, model::itemp)) = eos_state.T;
            for (int n = 0; n < NumSpec; n++) {
                model_state(MODELINDEX(i, model::ispec+n)) = model_params.xn[n];
            }

        }

        amrex::Print() << i << " " model_state(MODELINDEX(i, idens_model)) << " " << model_state(MODELINDEX(i, ipres_model)) << std::endl;

        //  now integrate down

        for (int i = ibase-1; i >= 0; --i) {

            Real dens_zone = model::state(MODELINDEX(i+1, model::idens));

            bool converged_hse = false;

            for (int iter = 0; iter < MAX_ITER; iter++) {

                Real p_want = model_state(MODELINDEX(i+1, model::ipres)) -
                    dx * 0.5_rt * (dens_zone + model::state(MODELINDEX(i+1, model::idens))) * const_grav;

                // use the EOS with constant entropy to find corrected state

                eos_t eos_state;
                eos_state.p = p_want;
                eos_state.s = entropy_fixed;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = model_params.xn[n];
                }

                eos(eos_input_ps, eos_state);

                Real drho = eos_state.rho - dens_zone;
                Real dens_zone = eos_state.rho;
                if (std::abs(drho) < TOL * dens_zone) {
                    converged_hse = true;
                    break;
                }
            }

            if (!converged_hse) {
                amrex::Error("failed to convergence in initial model generation");
            }

            // initialze zone

            model::state(MODELINDEX(i, model::idens)) = dens_zone;
            model::state(MODELINDEX(i, model::ipres)) = p_want;
            model::state(MODELINDEX(i, model::itemp)) = eos_state.T;
            for (int n = 0; n < NumSpec; n++) {
                model_state(MODELINDEX(i, model::ispec+n)) = model_params.xn[n];
            }

        }

    }
}


  function f(p, s, g, rho, T, xn) result (rhs)

    use eos_type_module, only : eos_t, eos_input_ps
    use eos_module, only : eos

    real(rt), intent(in) :: p, s, g, rho, T, xn(nspec)

    real(rt) :: rhs
    type(eos_t) :: eos_state

    eos_state % T = T ! initial guess
    eos_state % rho = rho ! initial guess
    eos_state % xn(:) = xn(:)
    eos_state % p = p
    eos_state % s = s

    call eos(eos_input_ps, eos_state)

    rhs = eos_state % rho * g

  end function f

end module initial_model_module
