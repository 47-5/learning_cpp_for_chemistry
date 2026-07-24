# miniHF 实现路线

第一版目标：对一个闭壳层小分子，使用给定 AO 积分完成 RHF 单点能计算，并与 PySCF 对照。

## 数据流

```text
geometry and basis metadata
  -> S, T, V, ERI
  -> H_core = T + V
  -> initial density P
  -> Fock matrix F(P)
  -> solve F C = S C epsilon
  -> new density P
  -> electronic and total energy
  -> convergence check
```

## 推荐实现顺序

1. 从文本或测试夹具读取小尺寸积分。
2. 计算核间排斥能和核心 Hamiltonian。
3. 用 `S^{-1/2}` 进行对称正交化。
4. 从密度矩阵构造 Coulomb、exchange 和 Fock 矩阵。
5. 生成占据轨道密度矩阵。
6. 同时检查能量变化和密度矩阵 RMS 变化。
7. 添加单元测试，再加入 DIIS。

## 暂不加入

- Gaussian 积分递推
- UHF/ROHF
- DFT 数值积分
- MP2/CC
- 周期性边界、GPU 或 MPI
