# CMakeLists.txt 入门教程

这篇教程面向第一次接触 CMake 的 C++ 学习者，并以本仓库为实际例子。

读完后，你应该能够回答下面几个问题：

1. CMake、编译器和构建工具分别负责什么？
2. `CMakeLists.txt` 描述的“目标”是什么？
3. `add_executable()`、`add_library()` 和 `target_*()` 如何配合？
4. `PRIVATE`、`PUBLIC`、`INTERFACE` 分别表示什么？
5. `project_warnings` 是怎样把警告选项传给程序的？
6. 怎样配置、构建、运行和测试这个仓库？

---

## 1. CMake 到底是什么？

C++ 源代码不能直接运行。一个简单程序通常要经历：

```text
main.cpp
   │
   │ 编译（compile）
   ▼
main.o 或 main.obj
   │
   │ 链接（link）
   ▼
可执行程序
```

真正负责把 C++ 源代码变成机器代码的是编译器，例如：

- GCC 的 `g++`
- Clang 的 `clang++`
- Visual Studio 的 `cl.exe`

当项目只有一个文件时，可以手动调用编译器：

```bash
g++ -std=c++17 -Wall main.cpp -o hello
```

但是项目变大以后，我们还需要处理：

- 多个 `.cpp` 文件；
- 头文件搜索路径；
- 第三方库；
- 不同操作系统和编译器；
- Debug 和 Release 配置；
- 测试；
- 文件之间的依赖关系。

CMake 的工作是读取 `CMakeLists.txt`，然后生成真正的构建系统：

```text
CMakeLists.txt
      │
      │ CMake 配置和生成
      ▼
Ninja / Makefile / Visual Studio 工程
      │
      │ 调用编译器和链接器
      ▼
可执行程序或库
```

因此要记住：

> CMake 通常不是编译器，也不直接替代编译器。
>
> CMake 是“构建系统生成器”和“项目构建规则的管理工具”。

---

## 2. 第一个最小项目

假设目录结构是：

```text
hello_project/
├── CMakeLists.txt
└── main.cpp
```

`main.cpp`：

```cpp
#include <iostream>

int main()
{
    std::cout << "Hello, CMake!\n";
    return 0;
}
```

最小的 `CMakeLists.txt` 可以写成：

```cmake
cmake_minimum_required(VERSION 3.20)

project(hello_project LANGUAGES CXX)

add_executable(hello main.cpp)
```

逐句理解：

```cmake
cmake_minimum_required(VERSION 3.20)
```

要求 CMake 版本至少为 3.20。

```cmake
project(hello_project LANGUAGES CXX)
```

定义项目：

- 项目名称是 `hello_project`；
- 使用的语言是 C++，在 CMake 中写作 `CXX`。

```cmake
add_executable(hello main.cpp)
```

定义一个可执行程序目标：

- 目标名称：`hello`
- 源文件：`main.cpp`

这里的 `hello` 是一个 **CMake 目标（target）**。它既是 CMake 内部引用这个程序时使用的名字，也通常是最终可执行文件的名字。

---

## 3. 配置、构建和运行是三个不同动作

在项目根目录执行：

```bash
cmake -S . -B build -G Ninja
```

这一步称为“配置和生成”：

- `-S .`：源码目录是当前目录；
- `-B build`：把生成文件放进 `build` 目录；
- `-G Ninja`：生成 Ninja 构建文件。

目录会变成：

```text
hello_project/
├── CMakeLists.txt
├── main.cpp
└── build/               # CMake 生成的文件
```

接下来执行：

```bash
cmake --build build
```

这一步才真正调用编译器和链接器。

最后运行程序。在 Linux/WSL + Ninja 中通常是：

```bash
./build/hello
```

三个动作可以概括为：

```text
cmake -S . -B build    配置项目并生成构建文件
cmake --build build    编译和链接
./build/hello          运行程序
```

把生成文件放在独立的 `build` 目录称为 **out-of-source build**。它可以避免大量临时文件污染源码目录。

---

## 4. CMake 最重要的概念：目标 target

现代 CMake 的核心思想不是“到处设置编译参数”，而是：

> 先定义目标，再说明每个目标需要什么。

常见目标有三类。

### 4.1 可执行程序

```cmake
add_executable(chemistry_app
    src/main.cpp
    src/molecule.cpp
)
```

这表示创建一个名为 `chemistry_app` 的可执行程序，它由两个 `.cpp` 文件组成。

### 4.2 真正的库

```cmake
add_library(chemistry_math STATIC
    src/matrix.cpp
    src/integral.cpp
)
```

`STATIC` 表示静态库。它会编译源文件，并生成类似下面的文件：

```text
libchemistry_math.a    # Linux
chemistry_math.lib     # Windows
```

也可以使用：

```cmake
add_library(chemistry_math SHARED ...)
```

`SHARED` 表示动态库，通常生成 `.so`、`.dll` 或 `.dylib`。

### 4.3 接口库

```cmake
add_library(project_warnings INTERFACE)
```

接口库比较特殊：

- 没有需要编译的 `.cpp` 文件；
- 不会生成 `.a`、`.lib` 或 `.dll`；
- 它只负责保存并传递使用要求。

“使用要求”可以包括：

- 编译选项；
- 预处理宏；
- 头文件搜索目录；
- 其他依赖目标。

所以对初学者来说，可以暂时把接口库理解为：

> 一个有名字、可以复用的“编译配置包”。

---

## 5. `project_warnings` 到底怎样工作？

本仓库的顶层 [`CMakeLists.txt`](../CMakeLists.txt) 中有：

```cmake
add_library(project_warnings INTERFACE)
```

这一步创建一个空的接口目标：

```text
project_warnings
└── 暂时没有任何配置
```

然后根据当前编译器，把警告选项保存到这个目标中：

```cmake
if(MSVC)
    target_compile_options(
        project_warnings
        INTERFACE
        /W4
        /permissive-
    )
else()
    target_compile_options(
        project_warnings
        INTERFACE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()
```

`target_compile_options()` 的基本形式是：

```cmake
target_compile_options(目标名 传播范围 编译选项...)
```

因此：

```cmake
target_compile_options(project_warnings INTERFACE -Wall -Wextra)
```

可以读成：

> 把 `-Wall` 和 `-Wextra` 保存为 `project_warnings` 提供给使用者的编译选项。

如果当前使用 GCC 或 Clang，此时可以想象成：

```text
project_warnings
└── 提供给使用者的编译选项
    ├── -Wall
    ├── -Wextra
    └── -Wpedantic
```

接下来，[`01_basic/CMakeLists.txt`](../01_basic/CMakeLists.txt) 中写了：

```cmake
add_executable(basic_hello src/main.cpp)

target_link_libraries(
    basic_hello
    PRIVATE
    project_warnings
)
```

第二条命令建立了两个目标之间的关系：

```text
project_warnings
        │
        │ 提供编译警告选项
        ▼
basic_hello
```

CMake 看到 `basic_hello` 使用了 `project_warnings`，就会读取后者的 `INTERFACE` 设置，并在编译 `basic_hello` 时加入这些选项。

如果使用 GCC，最终编译命令可能近似于：

```bash
g++ -std=c++17 -Wall -Wextra -Wpedantic src/main.cpp -o basic_hello
```

这里有一个容易困惑的地方：

```cmake
target_link_libraries(basic_hello PRIVATE project_warnings)
```

虽然命令名称中有 `link_libraries`，但它的作用不只是传递真正的库文件。在现代 CMake 中，它也用来建立目标之间的依赖关系，并传递目标的使用要求。

因为 `project_warnings` 没有真正的库文件，所以链接器不会收到一个叫 `project_warnings.lib` 的文件；CMake 只会传递它携带的编译设置。

如果不用接口库，就要对每个程序重复编写：

```cmake
if(MSVC)
    target_compile_options(basic_hello PRIVATE /W4 /permissive-)
    target_compile_options(eigen_matrix_demo PRIVATE /W4 /permissive-)
    target_compile_options(minihf_placeholder PRIVATE /W4 /permissive-)
else()
    target_compile_options(basic_hello PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(eigen_matrix_demo PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(minihf_placeholder PRIVATE -Wall -Wextra -Wpedantic)
endif()
```

使用接口库后，只需要定义一次：

```cmake
target_compile_options(project_warnings INTERFACE ...)
```

各个程序再选择使用它：

```cmake
target_link_libraries(basic_hello PRIVATE project_warnings)
target_link_libraries(eigen_matrix_demo PRIVATE project_warnings)
target_link_libraries(minihf_placeholder PRIVATE project_warnings)
```

如果用接近 C++ 的伪代码表示，整个过程类似：

```cpp
CompileSettings project_warnings;

if (compiler == MSVC) {
    project_warnings.options = {"/W4", "/permissive-"};
} else {
    project_warnings.options = {"-Wall", "-Wextra", "-Wpedantic"};
}

basic_hello.use(project_warnings);
eigen_matrix_demo.use(project_warnings);
minihf_placeholder.use(project_warnings);
```

接口库与其说像“函数”，不如说更像一个保存了配置的对象。

---

## 6. `PRIVATE`、`PUBLIC` 和 `INTERFACE`

很多 `target_*()` 命令都要求写下面三个关键词之一：

- `PRIVATE`
- `PUBLIC`
- `INTERFACE`

它们回答的是：

> 这项要求是当前目标自己需要，还是使用当前目标的其他目标也需要？

| 关键词 | 当前目标自己使用 | 传递给使用者 |
|---|---:|---:|
| `PRIVATE` | 是 | 否 |
| `PUBLIC` | 是 | 是 |
| `INTERFACE` | 否 | 是 |

### 6.1 `PRIVATE`

```cmake
target_compile_options(my_app PRIVATE -Wall)
```

含义：

> `my_app` 自己编译时使用 `-Wall`，但不把它继续传给其他目标。

对于可执行程序，绝大多数要求都可以先从 `PRIVATE` 开始，因为通常不会有其他目标链接一个可执行程序。

### 6.2 `PUBLIC`

假设库的公开头文件中直接使用了 Eigen：

```cpp
// include/chemistry/matrix.hpp
#include <Eigen/Dense>

Eigen::MatrixXd make_matrix();
```

那么这个库自己需要 Eigen，使用这个库的程序也需要 Eigen：

```cmake
target_link_libraries(chemistry_math PUBLIC Eigen3::Eigen)
```

### 6.3 `INTERFACE`

```cmake
target_compile_options(project_warnings INTERFACE -Wall)
```

含义：

> `project_warnings` 自己不需要编译，但使用它的目标需要 `-Wall`。

注意，下面两个 `INTERFACE` 出现在不同命令中：

```cmake
add_library(project_warnings INTERFACE)
```

这里的 `INTERFACE` 表示目标类型是接口库。

```cmake
target_compile_options(project_warnings INTERFACE -Wall)
```

这里的 `INTERFACE` 表示 `-Wall` 只提供给目标的使用者。

### 6.4 一个简单判断方法

遇到作用域时可以依次问：

1. 当前目标自己编译时需要它吗？
2. 使用当前目标的其他目标也需要它吗？

答案分别是：

```text
自己需要，使用者不需要  -> PRIVATE
自己需要，使用者也需要  -> PUBLIC
自己不需要，使用者需要  -> INTERFACE
```

---

## 7. 常见的 `target_*()` 命令

现代 CMake 推荐把设置绑定到具体目标。

### 7.1 添加编译选项

```cmake
target_compile_options(my_app PRIVATE -Wall -Wextra)
```

### 7.2 添加预处理宏

```cmake
target_compile_definitions(my_app PRIVATE ENABLE_DEBUG_OUTPUT)
```

C++ 中就可以使用：

```cpp
#ifdef ENABLE_DEBUG_OUTPUT
std::cout << "debug information\n";
#endif
```

### 7.3 添加头文件搜索目录

```cmake
target_include_directories(chemistry_math
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
```

这样源代码可以写：

```cpp
#include <chemistry/matrix.hpp>
```

### 7.4 链接或使用其他目标

```cmake
target_link_libraries(my_app
    PRIVATE
    chemistry_math
    Eigen3::Eigen
    project_warnings
)
```

这些依赖可能是真正的库，也可能是只传递设置的接口目标。

### 7.5 指定目标使用的 C++ 特性

除了设置全局默认标准，还可以写：

```cmake
target_compile_features(my_app PRIVATE cxx_std_17)
```

这表示 `my_app` 至少需要 C++17。对于可复用库，目标级写法通常比全局变量更精确。

---

## 8. 阅读本仓库的顶层 `CMakeLists.txt`

下面按功能块理解本仓库的 [`CMakeLists.txt`](../CMakeLists.txt)。

### 8.1 CMake 版本和项目信息

```cmake
cmake_minimum_required(VERSION 3.20)

project(
    learning_cpp_for_chemistry
    VERSION 0.1.0
    DESCRIPTION "Small C++ exercises for computational chemistry"
    LANGUAGES CXX
)
```

这部分定义：

- 最低 CMake 版本；
- 项目名称；
- 项目版本；
- 项目描述；
- 使用 C++。

### 8.2 C++ 标准

```cmake
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
```

含义分别是：

- 默认使用 C++17；
- 必须支持 C++17，不能退回旧标准；
- 尽量关闭编译器专有语言扩展。

这些变量会初始化后面创建的 C++ 目标，包括通过 `add_subdirectory()` 创建的目标。

### 8.3 生成编译命令数据库

```cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
```

这要求 CMake 生成 `compile_commands.json`。clangd、clang-tidy 和编辑器插件可以通过它知道每个 `.cpp` 文件的编译参数。

它主要适用于 Ninja 和 Makefile 生成器。

### 8.4 加载测试功能

```cmake
include(CTest)
```

它会提供 `BUILD_TESTING` 选项，默认开启，并允许项目使用 `add_test()` 注册测试。

### 8.5 统一警告选项

```cmake
add_library(project_warnings INTERFACE)

if(MSVC)
    target_compile_options(project_warnings INTERFACE /W4 /permissive-)
else()
    target_compile_options(project_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()
```

这里使用接口库统一保存不同编译器的警告选项。

### 8.6 加载基础练习子目录

```cmake
add_subdirectory(01_basic)
```

CMake 会进入 `01_basic`，继续读取其中的 `CMakeLists.txt`。

这不是启动另一个独立项目，而是把该目录定义的目标加入当前构建。

### 8.7 创建 Eigen 示例开关

```cmake
option(
    LCC_BUILD_EIGEN_EXAMPLE
    "Build the Eigen matrix example"
    ON
)
```

它创建一个布尔选项，默认值为 `ON`。

可以在配置时关闭：

```bash
cmake -S . -B build -DLCC_BUILD_EIGEN_EXAMPLE=OFF
```

选项会保存在 `build/CMakeCache.txt` 中。再次运行 CMake 时，如果没有明确修改，它通常会继续使用缓存中的值。

### 8.8 查找 Eigen

```cmake
if(LCC_BUILD_EIGEN_EXAMPLE)
    find_package(Eigen3 3.3 QUIET NO_MODULE)

    if(Eigen3_FOUND)
        add_subdirectory(02_eigen)
    else()
        message(STATUS "Eigen3 was not found; skipping 02_eigen")
    endif()
endif()
```

流程是：

```text
是否启用 Eigen 示例？
        │
        ├── 否：不查找，也不构建
        │
        └── 是：查找 Eigen 3.3+
                    │
                    ├── 找到：加载 02_eigen
                    └── 没找到：显示提示并跳过
```

`QUIET` 表示找不到时不让 `find_package()` 自己输出警告；项目随后用自己的 `message()` 给出更友好的提示。

`NO_MODULE` 表示使用 Eigen 安装时提供的 CMake 配置文件。

### 8.9 加载 miniHF 子目录

```cmake
add_subdirectory(03_miniHF)
```

它没有放在 `if()` 中，因此总会被加载。

---

## 9. 子目录中的 `CMakeLists.txt`

[`01_basic/CMakeLists.txt`](../01_basic/CMakeLists.txt) 内容很短：

```cmake
add_executable(basic_hello src/main.cpp)
target_link_libraries(basic_hello PRIVATE project_warnings)

if(BUILD_TESTING)
    add_test(NAME basic_hello_runs COMMAND basic_hello)
endif()
```

第一句：

```cmake
add_executable(basic_hello src/main.cpp)
```

创建可执行程序目标 `basic_hello`。

第二句：

```cmake
target_link_libraries(basic_hello PRIVATE project_warnings)
```

让 `basic_hello` 使用统一的警告配置。

第三部分：

```cmake
if(BUILD_TESTING)
    add_test(NAME basic_hello_runs COMMAND basic_hello)
endif()
```

如果测试功能开启，就注册一个名为 `basic_hello_runs` 的测试。CTest 运行这个测试时，会启动 `basic_hello`；如果程序返回 0，测试通常视为通过。

[`02_eigen/CMakeLists.txt`](../02_eigen/CMakeLists.txt) 多了一个依赖：

```cmake
target_link_libraries(
    eigen_matrix_demo
    PRIVATE
    Eigen3::Eigen
    project_warnings
)
```

`Eigen3::Eigen` 是 `find_package(Eigen3 ...)` 成功后提供的导入目标。它主要携带 Eigen 头文件路径等使用要求。

程序不需要自己猜测 Eigen 安装在 `/usr/include/eigen3` 还是其他位置；只要使用 `Eigen3::Eigen`，CMake 就会传递正确的信息。

---

## 10. 测试：CTest 和 `add_test()`

顶层文件加载：

```cmake
include(CTest)
```

子目录注册：

```cmake
add_test(NAME basic_hello_runs COMMAND basic_hello)
```

构建后运行所有测试：

```bash
ctest --test-dir build --output-on-failure
```

参数含义：

- `--test-dir build`：测试构建目录中的项目；
- `--output-on-failure`：测试失败时显示程序输出。

只运行名字匹配的测试：

```bash
ctest --test-dir build -R basic_hello --output-on-failure
```

关闭测试：

```bash
cmake -S . -B build -DBUILD_TESTING=OFF
```

注意：`add_test()` 只是注册如何运行测试，不负责构建被测试的程序。因此一般要先执行：

```bash
cmake --build build
```

---

## 11. 本仓库的常用构建命令

推荐在 WSL 的仓库根目录执行。

### 11.1 第一次配置

```bash
cmake -S . -B build -G Ninja
```

### 11.2 构建全部目标

```bash
cmake --build build
```

### 11.3 只构建一个目标

```bash
cmake --build build --target basic_hello
```

### 11.4 运行测试

```bash
ctest --test-dir build --output-on-failure
```

### 11.5 运行程序

```bash
./build/01_basic/basic_hello
./build/02_eigen/eigen_matrix_demo
./build/03_miniHF/minihf_placeholder
```

如果没有安装 Eigen，第二个程序不会生成。

### 11.6 关闭 Eigen 示例

```bash
cmake -S . -B build -DLCC_BUILD_EIGEN_EXAMPLE=OFF
cmake --build build
```

### 11.7 重新开启 Eigen 示例

```bash
cmake -S . -B build -DLCC_BUILD_EIGEN_EXAMPLE=ON
cmake --build build
```

### 11.8 使用 Debug 配置

对于 Ninja 这样的单配置生成器：

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

Debug 配置通常包含调试信息，便于 GDB、CLion 或 VS Code 调试。

---

## 12. CMake 变量、选项和目标不要混淆

### 12.1 普通变量

```cmake
set(MY_NAME "chemistry")
message(STATUS "Project topic: ${MY_NAME}")
```

使用 `${变量名}` 读取变量。

### 12.2 布尔选项

```cmake
option(BUILD_DEMO "Build the demo program" ON)

if(BUILD_DEMO)
    add_subdirectory(demo)
endif()
```

选项通常可以用 `-D` 修改：

```bash
cmake -S . -B build -DBUILD_DEMO=OFF
```

### 12.3 目标

```cmake
add_executable(my_app main.cpp)
add_library(my_library STATIC library.cpp)
add_library(my_settings INTERFACE)
```

`my_app`、`my_library` 和 `my_settings` 都是目标名。引用目标时一般直接写名字，不需要 `${}`：

```cmake
target_link_libraries(my_app PRIVATE my_library my_settings)
```

错误的思维是：

```cmake
target_link_libraries(${my_app} ...)
```

除非 `my_app` 本身是一个变量，否则不应该写 `${}`。

---

## 13. CMake 中的条件判断

基本形式是：

```cmake
if(条件)
    # 条件为真时执行
else()
    # 条件为假时执行
endif()
```

本项目使用：

```cmake
if(MSVC)
```

`MSVC` 是 CMake 提供的条件，表示当前是否使用 Microsoft Visual C++ 风格的编译器。

还使用：

```cmake
if(Eigen3_FOUND)
```

这是 `find_package(Eigen3 ...)` 设置的结果变量。

需要注意：

> CMake 的 `if()` 在“配置阶段”执行，不是在 C++ 程序运行时执行。

因此 CMake 会先决定使用哪组编译选项，再生成 Ninja 或其他构建文件。

---

## 14. 常见错误及理解方法

### 14.1 `cmake: command not found`

说明当前环境没有安装 CMake，或者 `cmake` 不在 `PATH` 中。

在 Ubuntu/WSL 中可以安装：

```bash
sudo apt update
sudo apt install cmake ninja-build build-essential
```

### 14.2 找不到 Eigen

本项目会显示：

```text
Eigen3 was not found; skipping 02_eigen
```

在 Ubuntu/WSL 中安装：

```bash
sudo apt install libeigen3-dev
```

然后重新配置：

```bash
cmake -S . -B build
```

### 14.3 修改了选项，但结果没有变化

CMake 会把配置保存在：

```text
build/CMakeCache.txt
```

请明确传递新的值：

```bash
cmake -S . -B build -DLCC_BUILD_EIGEN_EXAMPLE=ON
```

如果构建目录确实已经混乱，也可以删除**明确的项目构建目录**后重新配置。删除之前要确认路径，避免误删源码。

### 14.4 修改 `.cpp` 后是否要重新运行 CMake？

通常不需要。只要文件列表和 CMake 配置没有变化，直接执行：

```bash
cmake --build build
```

如果添加了新的 `.cpp` 文件，并且需要把它写进 `add_executable()` 或 `add_library()`，就要修改 `CMakeLists.txt`。构建工具通常会自动重新运行 CMake，也可以手动再次执行配置命令。

### 14.5 为什么不推荐到处使用全局编译选项？

例如：

```cmake
add_compile_options(-Wall)
include_directories(some/path)
```

这些目录级命令可能影响后面许多无关目标。项目变大后，很难判断某个设置来自哪里。

现代 CMake 更推荐：

```cmake
target_compile_options(my_app PRIVATE -Wall)
target_include_directories(my_app PRIVATE some/path)
```

这样每个设置都明确属于某个目标。

---

## 15. 一个带库的完整小例子

目录结构：

```text
demo/
├── CMakeLists.txt
├── include/
│   └── add.hpp
└── src/
    ├── add.cpp
    └── main.cpp
```

`include/add.hpp`：

```cpp
#pragma once

double add(double a, double b);
```

`src/add.cpp`：

```cpp
#include "add.hpp"

double add(double a, double b)
{
    return a + b;
}
```

`src/main.cpp`：

```cpp
#include "add.hpp"

#include <iostream>

int main()
{
    std::cout << add(1.0, 2.0) << '\n';
    return 0;
}
```

`CMakeLists.txt`：

```cmake
cmake_minimum_required(VERSION 3.20)
project(demo LANGUAGES CXX)

add_library(math STATIC
    src/add.cpp
)

target_include_directories(math
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(calculator
    src/main.cpp
)

target_link_libraries(calculator
    PRIVATE
    math
)
```

依赖关系是：

```text
math 静态库
├── 编译 src/add.cpp
└── PUBLIC 提供 include/ 头文件路径
             │
             ▼
calculator 可执行程序
└── 编译 src/main.cpp，并链接 math
```

这里头文件目录使用 `PUBLIC`，是因为：

1. `math` 自己编译 `add.cpp` 时需要找到 `add.hpp`；
2. 使用 `math` 的 `calculator` 编译 `main.cpp` 时也需要找到 `add.hpp`。

---

## 16. 建议的学习练习

### 练习一：观察警告配置的传播

在 `01_basic/src/main.cpp` 中临时定义一个未使用变量：

```cpp
int unused_value = 42;
```

重新构建：

```bash
cmake --build build
```

观察编译器是否给出未使用变量警告。这可以帮助理解 `project_warnings` 确实影响了 `basic_hello`。

完成观察后再删除这个临时代码。

### 练习二：创建自己的接口配置目标

尝试在顶层添加：

```cmake
add_library(project_debug_options INTERFACE)
target_compile_definitions(
    project_debug_options
    INTERFACE
    LEARNING_CMAKE
)
```

然后让 `basic_hello` 使用它：

```cmake
target_link_libraries(
    basic_hello
    PRIVATE
    project_warnings
    project_debug_options
)
```

C++ 中可以检查：

```cpp
#ifdef LEARNING_CMAKE
std::cout << "CMake definition received\n";
#endif
```

### 练习三：增加一个新程序

创建：

```text
04_practice/
├── CMakeLists.txt
└── src/
    └── main.cpp
```

子目录的 `CMakeLists.txt`：

```cmake
add_executable(cmake_practice src/main.cpp)
target_link_libraries(cmake_practice PRIVATE project_warnings)
```

顶层加入：

```cmake
add_subdirectory(04_practice)
```

然后配置、构建并运行它。

---

## 17. 初学阶段最值得记住的内容

不必一次记住所有 CMake 命令。先掌握下面这些：

```cmake
# 定义项目
cmake_minimum_required(VERSION 3.20)
project(my_project LANGUAGES CXX)

# 创建程序和库
add_executable(my_app main.cpp)
add_library(my_library STATIC library.cpp)
add_library(my_settings INTERFACE)

# 给目标添加设置或依赖
target_compile_options(my_app PRIVATE -Wall)
target_include_directories(my_library PUBLIC include)
target_link_libraries(my_app PRIVATE my_library my_settings)

# 加载子目录
add_subdirectory(subdirectory)

# 查找第三方库
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(my_app PRIVATE Eigen3::Eigen)
```

以及三个命令：

```bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

最后记住现代 CMake 的基本阅读方法：

```text
1. 找到 add_executable() 和 add_library()：项目有哪些目标？
2. 查看每个 target_*()：每个目标需要什么？
3. 查看 target_link_libraries()：目标之间怎样连接？
4. 查看 PRIVATE/PUBLIC/INTERFACE：设置怎样传播？
5. 查看 add_subdirectory()：项目被拆成了哪些部分？
```

只要始终围绕“目标及其依赖”阅读，`CMakeLists.txt` 就不会再像一组彼此无关的命令。
