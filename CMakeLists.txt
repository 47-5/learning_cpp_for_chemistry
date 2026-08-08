# 这是整个仓库最上层的 CMake 配置文件。
# CMake 会先读取这里的配置，生成适合当前平台的构建系统，
# 然后再由构建系统调用编译器，把 .cpp 源文件编译成可执行程序。

# 要求使用 CMake 3.20 或更高版本。
# 如果版本太旧，CMake 会在配置阶段停止并提示升级。
cmake_minimum_required(VERSION 3.20)

# 定义项目的基本信息。
project(
    # 项目名称。它不等于最终生成的可执行文件名。
    learning_cpp_for_chemistry

    # 当前项目版本：主版本号.次版本号.修订号。
    VERSION 0.1.0

    # 项目描述，主要作为项目的元信息。
    DESCRIPTION "Small C++ exercises for computational chemistry"

    # 本项目只使用 C++；CXX 是 CMake 中对 C++ 的称呼。
    LANGUAGES CXX
)

# 默认使用 C++17 标准编译后面创建的所有 C++ 目标。
set(CMAKE_CXX_STANDARD 17)

# 必须支持 C++17；不允许编译器自动退回到较旧的 C++ 标准。
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 关闭编译器专有扩展，尽量只使用标准 C++17，提高代码可移植性。
# 例如 GCC/Clang 会倾向于使用 -std=c++17，而不是 -std=gnu++17。
set(CMAKE_CXX_EXTENSIONS OFF)

# 生成 compile_commands.json，记录每个源文件的完整编译命令。
# clangd、clang-tidy 和一些编辑器插件可以用它理解项目。
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 加载 CTest 测试支持，同时创建 BUILD_TESTING 选项（默认开启）。
# 子目录可以通过 if(BUILD_TESTING) 决定是否注册测试。
include(CTest)

# 创建一个“接口库”。它没有自己的源文件，也不会生成真正的库文件；
# 它只用来集中保存编译警告选项，再把这些选项传给使用它的程序。
add_library(project_warnings INTERFACE)

# MSVC 是 Microsoft Visual C++ 编译器。
if(MSVC)
    # /W4：开启较严格的警告。
    # /permissive-：让 MSVC 更严格地遵守 C++ 标准。
    target_compile_options(project_warnings INTERFACE /W4 /permissive-)
else()
    # 对 GCC、Clang 等编译器开启常用、额外和标准一致性警告。
    # 这些选项只产生警告，不会把警告自动当成编译错误。
    target_compile_options(project_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()

# 读取 00_cpp_fundamentals/CMakeLists.txt，并加入 StudyPlan.dev 基础练习。
add_subdirectory(00_cpp_fundamentals)

# 读取 01_basic/CMakeLists.txt，并把基础示例加入当前构建。
add_subdirectory(01_basic)

# 创建一个可开关的选项，决定是否构建 Eigen 矩阵示例；默认开启。
# 配置时可以使用下面的命令关闭它：
# cmake -S . -B build -DLCC_BUILD_EIGEN_EXAMPLE=OFF
option(LCC_BUILD_EIGEN_EXAMPLE "Build the Eigen matrix example" ON)

# 只有选项开启时才查找 Eigen。
if(LCC_BUILD_EIGEN_EXAMPLE)
    # 查找 3.3 或更高版本的 Eigen3。
    # QUIET 表示没找到时不让 find_package 自己报出警告；
    # NO_MODULE 表示使用 Eigen 提供的配置文件来查找它。
    find_package(Eigen3 3.3 QUIET NO_MODULE)

    # Eigen3_FOUND 是 find_package 设置的变量，表示是否成功找到 Eigen。
    if(Eigen3_FOUND)
        # 找到 Eigen 后，读取 02_eigen/CMakeLists.txt 并构建该示例。
        add_subdirectory(02_eigen)

        # Project_01 的 Molecule 类型使用 Eigen::Vector3d，因此也在找到 Eigen 后构建。
        add_subdirectory(ProgrammingProjects/Project_01)
    else()
        # 没找到时只显示提示并跳过示例，不让整个项目配置失败。
        message(STATUS "Eigen3 was not found; skipping 02_eigen")
        message(STATUS "On Ubuntu/WSL, install it with: sudo apt install libeigen3-dev")
    endif()
endif()

# 读取 03_miniHF/CMakeLists.txt，并把 miniHF 练习加入当前构建。
add_subdirectory(03_miniHF)
