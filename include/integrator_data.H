#ifndef INTEGRATOR_DATA_H
#define INTEGRATOR_DATA_H

#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <AMReX_ANSIEscCode.H>

#include <cstdint>

#include <ArrayUtilities.H>
#include <network.H>
#include <extern_parameters.H>

#ifndef NET_LOOP_UNROLL_LEN
#define NET_LOOP_UNROLL_LEN 1
#endif

#ifndef MICROPHYSICS_CUSTOM_JAC_T
#ifdef SINGLE_PRECISION_JACOBIAN
using jac_t = float;
#else
using jac_t = amrex::Real;
#endif
#define MICROPHYSICS_CUSTOM_JAC_T 1
#endif

using namespace amrex::literals;

constexpr int INT_NEQS = NumSpec + 1;
constexpr amrex::Real species_failure_tolerance = 1.0e-2_rt;

enum integrator_errors : std::int8_t {
    IERR_SUCCESS = 1,
    IERR_BAD_INPUTS = -1,
    IERR_DT_UNDERFLOW = -2,
    IERR_SPRAD_CONVERGENCE = -3,
    IERR_TOO_MANY_STEPS = -4,
    IERR_TOO_MUCH_ACCURACY_REQUESTED = -5,
    IERR_CORRECTOR_CONVERGENCE = -6,
    IERR_LU_DECOMPOSITION_ERROR = -7
};

template <typename BurnT>
constexpr int integrator_neqs()
{
    (void)sizeof(BurnT);
    return INT_NEQS;
}

using IArray1D = amrex::Array1D<short, 1, INT_NEQS>;
using RArray1D = amrex::Array1D<amrex::Real, 1, INT_NEQS>;
using RArray2D = ArrayUtil::MathArray2D<jac_t, 1, INT_NEQS, 1, INT_NEQS>;

#endif
