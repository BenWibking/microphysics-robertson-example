#ifndef EXTERN_PARAMETERS_H
#define EXTERN_PARAMETERS_H

#include <AMReX_REAL.H>

#include <string>

struct extern_t {};

inline extern_t init_extern_parameters() { return extern_t{}; }

namespace integrator_rp {
    inline amrex::Real ode_max_dt = 1.0e99;
    inline int ode_max_steps = 100000;
    inline int linalg_do_pivoting = 1;
    inline int use_jacobian_caching = 0;
}

#endif
