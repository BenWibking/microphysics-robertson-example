cmake_minimum_required(VERSION 3.20)

project(RobertsonVODE
    DESCRIPTION "Self-contained Robertson problem using Microphysics VODE"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

set(AMREX_TAG "development" CACHE STRING "AMReX tag to fetch")
set(AMReX_SPACEDIM 1 CACHE STRING "Spatial dimension for AMReX" FORCE)
set(AMReX_MPI OFF CACHE BOOL "Enable AMReX MPI support" FORCE)
set(AMReX_OMP OFF CACHE BOOL "Enable AMReX OpenMP support" FORCE)
set(AMReX_GPU_BACKEND NONE CACHE STRING "Disable GPU backends" FORCE)
FetchContent_Declare(
    amrex
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex
    GIT_TAG        ${AMREX_TAG}
)
FetchContent_MakeAvailable(amrex)

add_executable(robertson_vode
    src/main.cpp
)

# Locate the Microphysics tree. Prefer an in-tree checkout (../interfaces) but fall back to
# ../Microphysics/* so this example can build standalone.
set(MICROPHYSICS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
if(EXISTS "${MICROPHYSICS_ROOT}/interfaces")
    set(MICROPHYSICS_TREE "${MICROPHYSICS_ROOT}")
elseif(EXISTS "${MICROPHYSICS_ROOT}/Microphysics/interfaces")
    set(MICROPHYSICS_TREE "${MICROPHYSICS_ROOT}/Microphysics")
else()
    message(FATAL_ERROR "Could not locate Microphysics interfaces. Expected either '../interfaces' "
        "or '../Microphysics/interfaces'.")
endif()

target_include_directories(robertson_vode PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MICROPHYSICS_TREE}/interfaces
    ${MICROPHYSICS_TREE}/integration
    ${MICROPHYSICS_TREE}/integration/VODE
    ${MICROPHYSICS_TREE}/util
)

target_link_libraries(robertson_vode PRIVATE AMReX::amrex)

target_compile_features(robertson_vode PRIVATE cxx_std_17)
