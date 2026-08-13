# Learning C++ for Computational Chemistry

这个仓库用于系统学习面向计算化学与科学计算的现代 C++。当前重点不是实现完整量化程序，而是逐步打通下面这条链路：

```text
C++ 基础
  -> CMake 与调试
  -> Eigen 数值线性代数
  -> RHF/SCF 玩具程序
  -> 可测试、可维护的计算化学软件
```

## 仓库结构

```text
.
├── notes/                  # C++、CMake、Eigen 与 HF 实现笔记
├── 00_cpp_fundamentals/    # StudyPlan.dev C++ 基础练习
├── 01_basic/               # 不依赖第三方库的 C++ 基础示例
├── 02_eigen/               # Eigen 矩阵与对称本征值问题示例
├── 03_miniHF/              # miniHF 的代码骨架与分阶段路线
├── ProgrammingProjects/    # 按语言组织的项目、公共库与共享输入
└── examples/               # 后续独立小实验和输入示例
```

## 推荐环境

建议使用 WSL2 Ubuntu、CLion、GCC/G++、CMake、Ninja、GDB 和 Eigen：

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build gdb git libeigen3-dev
```

在 CLion 中将 Toolchain 设为 WSL，并确认 CMake、C++ Compiler 和 Debugger 均指向 WSL 内的工具。

换电脑时的完整安装、Git/SSH、CLion Profile、构建和故障恢复步骤见
[`ENVIRONMENT_SETUP.md`](ENVIRONMENT_SETUP.md)。

## 构建与运行

在仓库根目录执行：

```bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

运行示例：

```bash
./build/01_basic/basic_hello
./build/02_eigen/eigen_matrix_demo
./build/03_miniHF/minihf_placeholder
./build/ProgrammingProjects/cpp/Project_01/project01 \
    ProgrammingProjects/resources/Project_01/input/allene.dat
```

如果尚未安装 Eigen，CMake 会跳过 `eigen_matrix_demo`，其余示例仍然可以构建。

## 学习里程碑

1. 配置 CLion + WSL + CMake，并运行 `basic_hello`。
2. 安装 Eigen，完成一次实对称矩阵对角化。
3. 熟悉函数、引用、类、`std::vector` 和基本文件组织。
4. 用给定积分实现闭壳层 RHF 的 SCF 循环。
5. 加入正交化、密度矩阵、能量收敛判断和 DIIS。
6. 与 PySCF 的小分子/STO-3G 结果交叉验证。

## 当前边界

`03_miniHF` 目前只是工程骨架。第一版应保持小而完整：闭壳层、有限基组、给定一电子与二电子积分、单点 RHF，不急于加入 DFT、UHF、解析积分或并行化。
