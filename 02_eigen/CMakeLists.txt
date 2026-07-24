add_executable(eigen_matrix_demo src/main.cpp)
target_link_libraries(eigen_matrix_demo PRIVATE Eigen3::Eigen project_warnings)

if(BUILD_TESTING)
    add_test(NAME eigen_matrix_demo_runs COMMAND eigen_matrix_demo)
endif()
