# 将 Molecule 和 Atom 编译成静态库，以便后续项目或测试复用。
add_library(
    project01_molecule
    STATIC
    molecule.cpp
    molecule.h
)

# molecule.h 是这个库的公开接口，使用者需要能够通过 #include "molecule.h" 找到它。
target_include_directories(
    project01_molecule
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# molecule.h 公开使用 Eigen::Vector3d，因此 Eigen 依赖也必须传递给库的使用者。
# project_warnings 只负责这个库自身的编译警告，不需要继续传播。
target_link_libraries(
    project01_molecule
    PUBLIC Eigen3::Eigen
    PRIVATE project_warnings
)
