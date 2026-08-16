# Project 03：Restricted Hartree–Fock（C++ / Eigen）

这是 Project 3 的 C++17/Eigen 实现。程序读取预先计算好的 AO 积分，
对闭壳层体系执行 Restricted Hartree–Fock（RHF）自洽场计算，并打印每一步
SCF 迭代的能量与密度变化。

当前实现包括：

- 从目录读取几何、核排斥能、重叠积分、单电子积分和双电子积分；
- 构造对称正交化矩阵 \(X=S^{-1/2}\)；
- 使用核 Hamiltonian 的本征轨道生成初始密度；
- 构造 Fock 矩阵并求解 Roothaan–Hall 方程；
- 同时使用能量变化和密度矩阵变化判断收敛；
- 保存并打印完整的 SCF 迭代历史。

## 目录结构

```text
Project_03/
├── include/minihf/
│   ├── integrals.h        # AOIntegrals 和 ERI 访问接口
│   └── rhf_calculator.h   # SCF 数据结构与 RHFCalculator
├── src/
│   ├── integrals.cpp      # Project 3 输入文件读取
│   └── rhf_calculator.cpp # RHF 数值步骤与 SCF 循环
├── tests/
│   └── scaffold_test.cpp  # 输入读取和对象构造的冒烟测试
├── main.cpp               # 命令行入口与 SCF 结果输出
└── CMakeLists.txt
```

## 数值约定

本项目的密度矩阵不包含闭壳层占据因子 2：

\[
D_{\mu\nu}=\sum_{i\in\mathrm{occ}}C_{\mu i}C_{\nu i}.
\]

与此对应，Fock 矩阵和电子能量采用

\[
F_{\mu\nu}=H^{\mathrm{core}}_{\mu\nu}
+\sum_{\lambda\sigma}D_{\lambda\sigma}
\left[2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)\right],
\]

\[
E_\mathrm{elec}=\sum_{\mu\nu}D_{\mu\nu}
\left(H^{\mathrm{core}}_{\mu\nu}+F_{\mu\nu}\right).
\]

双电子积分以一维 `std::vector<double>` 稠密存储，并通过

```cpp
integrals.eri(mu, nu, lambda, sigma)
```

访问化学家记号下的 \((\mu\nu|\lambda\sigma)\)。稠密存储便于将第一版
C++ 实现与 Python 四维数组逐项比较；八重对称压缩可以作为后续优化。

## 构建

需要 CMake 3.20、支持 C++17 的编译器和 Eigen3。以下命令均从仓库根目录执行。

使用 Ninja：

```bash
cmake -S . -B build -G Ninja
cmake --build build --target project03
```

如果没有安装 Ninja，可以省略 `-G Ninja`，让 CMake 使用系统默认生成器：

```bash
cmake -S . -B build
cmake --build build --target project03
```

## 运行

程序接受一个参数：包含 Project 3 输入文件的目录。

```bash
./build/ProgrammingProjects/cpp/Project_03/project03 \
  ProgrammingProjects/resources/Project_03/input/h2o/STO-3G
```

如果使用 CLion 的 WSL Debug 构建目录，等价的命令是：

```bash
./cmake-build-wsl-debug/ProgrammingProjects/cpp/Project_03/project03 \
  ProgrammingProjects/resources/Project_03/input/h2o/STO-3G
```

建议从仓库根目录运行上述命令，这样示例中的相对输入路径才能正确解析。

输入目录应包含：

```text
geom.dat  # 原子数、核电荷和坐标
enuc.dat  # 核间排斥能
s.dat     # AO 重叠积分
t.dat     # 电子动能积分
v.dat     # 核对电子吸引积分
eri.dat   # 双电子排斥积分
```

当前入口默认按中性分子计算电子数，并要求电子数为偶数，以满足闭壳层 RHF
假设。

## 输出

程序会打印每次迭代的电子能量、总能量、能量变化和密度矩阵的 Frobenius
范数变化。H₂O/STO-3G 示例的最终结果应为：

```text
Converged: True
SCF iterations: 25
Final electronic energy: -82.944446990002 Eh
Nuclear repulsion: 8.002367061810 Eh
Final total energy: -74.942079928192 Eh
```

不同线性代数库、编译器或平台可能使最后几位和收敛步数出现很小差异。

默认 SCF 参数定义在 `SCFOptions` 中：

```cpp
double energy_tolerance{1.0e-10};
double density_tolerance{1.0e-8};
int max_iterations{100};
double overlap_eigenvalue_tolerance{1.0e-10};
```

## 测试

当前测试验证输入读取、矩阵尺寸、电子数和 `RHFCalculator` 的构造：

```bash
cmake --build build --target project03_scaffold_test
ctest --test-dir build -R project03_scaffold --output-on-failure
```

目前尚未加入完整 SCF 数值结果的自动化回归测试。
