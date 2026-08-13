add_library(
        project01_molecule
        STATIC
        atomic_masses.h
        molecule.cpp
        molecule.h
)

target_include_directories(
        project01_molecule
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
        project01_molecule
        PUBLIC Eigen3::Eigen
        PRIVATE project_warnings
)

add_executable(
        project01
        main.cpp
        physical_constants.h
)

target_link_libraries(
        project01
        PRIVATE
        project01_molecule
)