# CMake 速记

CMake 的核心不是逐条描述编译命令，而是声明“目标”及其依赖。

```cmake
add_executable(my_program src/main.cpp)
target_link_libraries(my_program PRIVATE some_library)
```

常用流程：

```bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

- `-S .`：源码目录。
- `-B build`：独立构建目录，避免生成文件污染源码。
- `target_*`：把头文件路径、编译选项和链接库绑定到具体目标。
- `PRIVATE`：依赖仅用于当前目标。
- `PUBLIC`：当前目标及其使用者都需要该依赖。
- `INTERFACE`：仅目标的使用者需要该依赖。
