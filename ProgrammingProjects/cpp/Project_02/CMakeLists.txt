add_executable(
    project02
    main.cpp
)

target_link_libraries(
    project02
    PRIVATE
    chemistry_core
    project_warnings
)
if(BUILD_TESTING)
    add_test(
        NAME project02_h2o_runs
        COMMAND
        project02
        ${CMAKE_SOURCE_DIR}/ProgrammingProjects/resources/Project_02/input/h2o_geom.txt
        ${CMAKE_SOURCE_DIR}/ProgrammingProjects/resources/Project_02/input/h2o_hessian.txt
    )
endif()