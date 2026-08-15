add_library(
    project03_core
    STATIC
    src/integrals.cpp
    src/rhf_calculator.cpp
)

target_include_directories(
    project03_core
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(
    project03_core
    PUBLIC
    Eigen3::Eigen
    PRIVATE
    project_warnings
)

add_executable(project03 main.cpp)
target_link_libraries(
    project03
    PRIVATE
    project03_core
    project_warnings
)

if(BUILD_TESTING)
    add_executable(
        project03_scaffold_test
        tests/scaffold_test.cpp
    )
    target_link_libraries(
        project03_scaffold_test
        PRIVATE
        project03_core
        project_warnings
    )

    add_test(
        NAME project03_scaffold
        COMMAND
        project03_scaffold_test
        ${CMAKE_SOURCE_DIR}/ProgrammingProjects/resources/Project_03/input/h2o/STO-3G
    )
endif()
