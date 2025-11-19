#include <AMReX_REAL.H>

#include <iomanip>
#include <iostream>
#include <vector>

#include <integrator_data.H>

struct RobertsonState {};

template <typename BurnT, typename IntegratorT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs(const amrex::Real time, BurnT& state, IntegratorT& int_state,
         RArray1D& ydot, [[maybe_unused]] const bool in_jacobian=false);

template <typename BurnT, typename IntegratorT, class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac(const amrex::Real time, BurnT& state, IntegratorT& int_state, MatrixType& jac);

#include <vode_type.H>
#include <vode_dvode.H>

using namespace amrex::literals;

template <typename BurnT, typename IntegratorT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs([[maybe_unused]] const amrex::Real time, [[maybe_unused]] BurnT& state,
         IntegratorT& int_state, RArray1D& ydot, [[maybe_unused]] const bool in_jacobian)
{
    const amrex::Real y1 = int_state.y(1);
    const amrex::Real y2 = int_state.y(2);
    const amrex::Real y3 = int_state.y(3);

    ydot(1) = -0.04_rt * y1 + 1.0e4_rt * y2 * y3;
    ydot(2) = 0.04_rt * y1 - 1.0e4_rt * y2 * y3 - 3.0e7_rt * y2 * y2;
    ydot(3) = 3.0e7_rt * y2 * y2;
    ydot(net_ienuc) = 0.0_rt;
}

template <typename BurnT, typename IntegratorT, class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac([[maybe_unused]] const amrex::Real time, [[maybe_unused]] BurnT& state,
         IntegratorT& int_state, MatrixType& jac)
{
    jac.zero();

    const amrex::Real y2 = int_state.y(2);
    const amrex::Real y3 = int_state.y(3);

    jac(1,1) = -0.04_rt;
    jac(1,2) = 1.0e4_rt * y3;
    jac(1,3) = 1.0e4_rt * y2;

    jac(2,1) = 0.04_rt;
    jac(2,2) = -1.0e4_rt * y3 - 6.0e7_rt * y2;
    jac(2,3) = -1.0e4_rt * y2;

    jac(3,1) = 0.0_rt;
    jac(3,2) = 6.0e7_rt * y2;
    jac(3,3) = 0.0_rt;
}

int main()
{
    RobertsonState state{};

    constexpr int int_neqs = integrator_neqs<RobertsonState>();
    dvode_t<int_neqs> vode_state{};

    vode_state.atol_spec = 1.0e-12_rt;
    vode_state.rtol_spec = 1.0e-12_rt;
    vode_state.atol_enuc = 1.0e-16_rt;
    vode_state.rtol_enuc = 1.0e-12_rt;
    vode_state.jacobian_type = 1;

    vode_state.t = 0.0_rt;
    vode_state.tout = 0.0_rt;

    vode_state.y(1) = 1.0_rt;
    vode_state.y(2) = 0.0_rt;
    vode_state.y(3) = 0.0_rt;
    vode_state.y(net_ienuc) = 0.0_rt;

    std::vector<amrex::Real> output_times = {
        4.0e-1_rt, 4.0e0_rt, 4.0e1_rt, 4.0e2_rt, 4.0e3_rt,
        4.0e4_rt, 4.0e5_rt, 4.0e6_rt, 4.0e7_rt, 4.0e8_rt,
        4.0e9_rt, 4.0e10_rt
    };

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Robertson problem solved with Microphysics VODE\n";
    std::cout << std::setw(14) << "time"
              << std::setw(14) << "y1"
              << std::setw(14) << "y2"
              << std::setw(14) << "y3" << '\n';
    std::cout << std::setw(14) << 0.0
              << std::setw(14) << vode_state.y(1)
              << std::setw(14) << vode_state.y(2)
              << std::setw(14) << vode_state.y(3) << '\n';

    for (amrex::Real target_time : output_times) {
        vode_state.tout = target_time;
        const int istate = dvode(state, vode_state);
        if (istate != IERR_SUCCESS) {
            std::cerr << "dvode integration failed with code " << istate << '\n';
            return istate;
        }

        std::cout << std::setw(14) << vode_state.t
                  << std::setw(14) << vode_state.y(1)
                  << std::setw(14) << vode_state.y(2)
                  << std::setw(14) << vode_state.y(3) << '\n';
    }

    return 0;
}
