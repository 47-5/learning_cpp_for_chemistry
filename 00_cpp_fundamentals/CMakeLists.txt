# 第一节练习：创建一个名为 study_hello_world 的可执行程序。
add_executable(
    study_hello_world
    01_types_and_variables/hello_world.cpp
)

add_executable(
        fizzbuzz
        02_functions_conditionals_loops/FizzBuzz.cpp
)


# 让这个练习继承顶层 CMakeLists.txt 中统一定义的编译警告选项。
target_link_libraries(study_hello_world PRIVATE project_warnings)

# 如果测试功能已开启，就把“程序能够正常运行”注册为一个测试。
if(BUILD_TESTING)
    add_test(NAME study_hello_world_runs COMMAND study_hello_world)
    add_test(NAME study_fizzbuzz_runs COMMAND fizzbuzz)
endif()
