# Eigen 速记

Eigen 是头文件为主的 C++ 线性代数库。对于 miniHF，最先需要掌握：

- `Eigen::MatrixXd`：动态尺寸实矩阵。
- `Eigen::VectorXd`：动态尺寸实向量。
- `Eigen::SelfAdjointEigenSolver`：实对称矩阵本征值问题。
- 矩阵转置、乘法、分块与逐元素操作。

`02_eigen` 中的示例求解普通本征值问题

```math
F C = C \varepsilon .
```

RHF 随后需要处理非正交 AO 基下的广义本征值问题

```math
F C = S C \varepsilon ,
```

通常先由重叠矩阵构造对称正交化矩阵 `X = S^{-1/2}`，再对 `X^T F X` 对角化。
