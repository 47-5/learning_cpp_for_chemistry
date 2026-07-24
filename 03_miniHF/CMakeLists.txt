add_executable(minihf_placeholder src/main.cpp)
target_link_libraries(minihf_placeholder PRIVATE project_warnings)

if(BUILD_TESTING)
    add_test(NAME minihf_placeholder_runs COMMAND minihf_placeholder)
endif()
