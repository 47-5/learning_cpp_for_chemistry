cmake_minimum_required(VERSION 3.20)

project(
    learning_cpp_for_chemistry
    VERSION 0.1.0
    DESCRIPTION "Small C++ exercises for computational chemistry"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CTest)

add_library(project_warnings INTERFACE)
if(MSVC)
    target_compile_options(project_warnings INTERFACE /W4 /permissive-)
else()
    target_compile_options(project_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_subdirectory(01_basic)

option(LCC_BUILD_EIGEN_EXAMPLE "Build the Eigen matrix example" ON)
if(LCC_BUILD_EIGEN_EXAMPLE)
    find_package(Eigen3 3.3 QUIET NO_MODULE)
    if(Eigen3_FOUND)
        add_subdirectory(02_eigen)
    else()
        message(STATUS "Eigen3 was not found; skipping 02_eigen")
        message(STATUS "On Ubuntu/WSL, install it with: sudo apt install libeigen3-dev")
    endif()
endif()

add_subdirectory(03_miniHF)
