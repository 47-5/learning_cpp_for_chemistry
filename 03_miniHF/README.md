# miniHF

这里将逐步实现一个教学用途的闭壳层 Hartree–Fock 程序。

建议的最终结构：

```text
03_miniHF/
├── include/minihf/
│   ├── molecule.hpp
│   ├── integrals.hpp
│   ├── scf.hpp
│   └── diis.hpp
├── src/
│   ├── main.cpp
│   ├── molecule.cpp
│   ├── scf.cpp
│   └── diis.cpp
└── tests/
```

当前 `minihf_placeholder` 只验证工程骨架可以被 CMake 编译。正式实现前，先在测试中固定一组来自 PySCF 的小尺寸 AO 积分和参考能量。
