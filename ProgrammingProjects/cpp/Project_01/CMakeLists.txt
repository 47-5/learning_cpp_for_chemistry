add_executable(
    project01
    main.cpp
)

target_link_libraries(
    project01
    PRIVATE
    chemistry_core
    project_warnings
)

if(BUILD_TESTING)
    add_test(
        NAME project01_allene_runs
        COMMAND
        project01
        ${CMAKE_SOURCE_DIR}/ProgrammingProjects/resources/Project_01/input/allene.dat
    )
endif()
