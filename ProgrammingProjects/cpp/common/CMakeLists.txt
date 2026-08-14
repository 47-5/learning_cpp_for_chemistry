add_library(
    chemistry_core
    STATIC
    src/geometry_analyzer.cpp
    src/molecule.cpp
    src/rotational_analyzer.cpp
        include/chemistry/vibrational_analysis.h
    src/vibrational_analysis.cpp
)

target_include_directories(
    chemistry_core
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
    chemistry_core
    PUBLIC
    Eigen3::Eigen
    PRIVATE
    project_warnings
)
