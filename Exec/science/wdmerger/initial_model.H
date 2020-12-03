#ifndef initial_model_H
#define initial_model_H

#include <AMReX.H>

#include <network_properties.H>
#include <eos.H>
#include <model_parser_data.H>
#include <model_parser.H>
#include <prob_parameters.H>

void establish_hse (model::initial_model_t& model,
                    Real& mass, Real& central_density,
                    Real envelope_mass, Real& radius,
                    const Real core_comp[NumSpec], const Real envelope_comp[NumSpec]);

#endif
