#ifndef NETWORK_H
#define NETWORK_H

#include <AMReX_REAL.H>

#include <string>
#include <vector>

constexpr int NumSpec = 3;
constexpr int NumAux = 0;
constexpr int net_ienuc = NumSpec + 1;

namespace Species {
    enum NetworkSpecies {
        Y1 = 1,
        Y2,
        Y3,
        NumberSpecies = Y3
    };
}

inline void network_init() {}

inline const std::vector<std::string> short_spec_names_cxx = {"Y1", "Y2", "Y3"};
inline const std::vector<std::string> spec_names_cxx = {"y1", "y2", "y3"};

#endif
