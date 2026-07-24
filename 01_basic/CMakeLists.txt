add_executable(basic_hello src/main.cpp)
target_link_libraries(basic_hello PRIVATE project_warnings)

if(BUILD_TESTING)
    add_test(NAME basic_hello_runs COMMAND basic_hello)
endif()
