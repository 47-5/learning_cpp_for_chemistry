# Project 3：从线性变分法到 Roothaan–Hall 方程与 RHF/SCF

## 1. 这份文档要回答什么

Project 3 的目标不是从 Gaussian 基函数计算积分，而是读取已经给出的 AO 积分，完成一个闭壳层限制性 Hartree–Fock（restricted Hartree–Fock, RHF）自洽场程序。

完整的逻辑链条是：

```text
线性变分法
    ↓
在有限基组中展开未知波函数
    ↓
广义本征值问题 Hc = ESc
    ↓
用 Slater 行列式保证多电子波函数反对称
    ↓
对 Slater 行列式中的轨道作变分
    ↓
Hartree–Fock 方程 fφᵢ = εᵢφᵢ
    ↓
用 AO 基函数展开分子轨道
    ↓
Roothaan–Hall 方程 FC = SCε
    ↓
F 依赖由 C 产生的密度矩阵 D
    ↓
通过 SCF 迭代寻找自洽解
```

本文的重点不只是列出公式，而是说明每个新对象为什么出现、它解决了什么问题，以及它怎样对应到 Project 3 的代码。

---

## 2. 符号、单位和研究范围

### 2.1 原子单位

除非特别说明，Hartree–Fock 推导通常使用原子单位：

\[
\hbar=m_e=e=4\pi\varepsilon_0=1.
\]

因此电子动能算符写成

\[
-\frac12\nabla^2,
\]

电子间 Coulomb 排斥写成

\[
\frac{1}{r_{12}}.
\]

能量单位是 hartree，长度单位是 bohr。

这里的“等于 1”指这些物理量在原子单位制中的**数值**为 1，并不改变物理关系。例如

\[
\hbar=\frac{h}{2\pi}
\]

始终成立；由于原子单位制取 \(\hbar=1\)，所以在这套单位制中 \(h=2\pi\)。类似地，\(e\) 表示元电荷的正值大小，因此电子的电荷为 \(-e=-1\)。

#### bohr 和 hartree 是怎样得到的

bohr 长度 \(a_0\) 和 hartree 能量 \(E_h\) 不是凭空规定的；它们可以从氢原子的 Schrödinger Hamiltonian 中自然地读出。SI 单位下，氢原子的电子 Hamiltonian 为

\[
\hat H
=-\frac{\hbar^2}{2m_e}\nabla_{\mathbf r}^2
-\frac{e^2}{4\pi\varepsilon_0r}.
\]

现在希望用某个典型原子长度 \(a_0\) 衡量坐标，定义无量纲坐标

\[
\boldsymbol\rho=\frac{\mathbf r}{a_0},
\qquad
\mathbf r=a_0\boldsymbol\rho.
\]

由于每次对坐标求导都会产生一个 \(1/a_0\)，有

\[
\nabla_{\mathbf r}^2
=\frac{1}{a_0^2}\nabla_{\boldsymbol\rho}^2.
\]

代回 Hamiltonian：

\[
\hat H
=-\frac{\hbar^2}{2m_ea_0^2}\nabla_{\boldsymbol\rho}^2
-\frac{e^2}{4\pi\varepsilon_0a_0}\frac{1}{\rho}.
\]

动能和 Coulomb 势能分别给出两个自然的能量尺度

\[
\frac{\hbar^2}{m_ea_0^2},
\qquad
\frac{e^2}{4\pi\varepsilon_0a_0}.
\]

原子束缚态中二者处在同一数量级。选择 \(a_0\) 使这两个尺度完全相等：

\[
\frac{\hbar^2}{m_ea_0^2}
=\frac{e^2}{4\pi\varepsilon_0a_0}.
\]

解出

\[
\boxed{
a_0
=\frac{4\pi\varepsilon_0\hbar^2}{m_ee^2}
}.
\]

这就是 bohr。上述匹配没有把动能算符中的 \(1/2\) 吸收到长度定义中，因此无量纲 Hamiltonian 的动能项最终仍是 \(-\tfrac12\nabla^2\)。也可以从量级估计

\[
E(r)\sim
\frac{\hbar^2}{2m_er^2}
-\frac{e^2}{4\pi\varepsilon_0r}
\]

出发，对 \(r\) 求极小值；同样会得到 \(r=a_0\)。这反映了电子的量子动能阻止它无限靠近原子核，而 Coulomb 吸引又试图把它拉近，两者平衡决定了典型原子尺度。

把 \(a_0\) 代入任一个能量尺度，定义 hartree 能量

\[
E_h
=\frac{e^2}{4\pi\varepsilon_0a_0}
=\frac{\hbar^2}{m_ea_0^2}.
\]

利用

\[
\frac{1}{a_0}
=\frac{m_ee^2}{4\pi\varepsilon_0\hbar^2},
\]

可得

\[
\boxed{
E_h
=\frac{m_ee^4}
{(4\pi\varepsilon_0)^2\hbar^2}
}.
\]

因此，bohr 是氢原子中动能与 Coulomb 势能平衡产生的自然长度，hartree 则是该长度对应的自然能量。数值上

\[
1\ a_0\approx 0.529177\ \text{\AA},
\qquad
1\ E_h\approx 27.2114\ \mathrm{eV}.
\]

若再令

\[
\mathbf r=a_0\boldsymbol\rho,
\qquad
\hat H=E_h\hat H_{\mathrm{au}},
\]

上面的氢原子 Hamiltonian 就变成

\[
\hat H_{\mathrm{au}}
=-\frac12\nabla_{\boldsymbol\rho}^2-\frac1\rho.
\]

量子化学中通常省略下标 \(\mathrm{au}\) 和无量纲坐标符号，直接写成 \(-\tfrac12\nabla^2-1/r\)。单位并没有消失，而是已经包含在 bohr 和 hartree 的定义中。

### 2.2 下标约定

| 下标 | 含义 |
|---|---|
| \(A,B\) | 原子核 |
| \(i,j\) | 占据分子轨道 |
| \(a,b\) | 未占据（虚）分子轨道 |
| \(p,q\) | 任意分子轨道 |
| \(\mu,\nu,\lambda,\sigma\) | AO 基函数 |
| \(x=(\mathbf r,\omega)\) | 电子的空间和自旋坐标 |

本文假设：

- Born–Oppenheimer 近似，原子核固定；
- 实 AO 基函数和实轨道系数；
- 偶数电子、闭壳层 RHF；
- 每个占据空间轨道中有一对自旋相反的电子。

---

## 3. 线性变分法

### 3.1 变分原理

设 Hamiltonian 为 \(\hat H\)，其精确基态满足

\[
\hat H\Psi_0=E_0\Psi_0.
\]

对任意满足归一化条件的试探波函数 \(\Phi\)，都有

\[
E[\Phi]
=\frac{\langle\Phi|\hat H|\Phi\rangle}
       {\langle\Phi|\Phi\rangle}
\ge E_0.
\]

这意味着：只要选择一族允许的试探波函数，并在其中最小化能量，就能得到基态能量的上界。

变分法并不保证试探空间中含有精确波函数；它保证的是，在选定试探空间内找到能量最低的那个波函数。

### 3.2 在有限基组中展开试探波函数

选取 \(K\) 个已知基函数 \(\{\chi_\mu\}\)，把未知波函数写成

\[
\Phi=\sum_{\mu=1}^{K}c_\mu\chi_\mu.
\]

未知函数问题于是变成未知系数 \(c_\mu\) 的问题。

定义 Hamiltonian 矩阵和重叠矩阵：

\[
H_{\mu\nu}
=\langle\chi_\mu|\hat H|\chi_\nu\rangle,
\]

\[
S_{\mu\nu}
=\langle\chi_\mu|\chi_\nu\rangle.
\]

如果基函数彼此不正交，\(S\ne I\)。下面从变分原理中的能量期望值出发，推导能量泛函的矩阵形式：

\[
E[\Phi]
=\frac{\langle\Phi|\hat H|\Phi\rangle}
       {\langle\Phi|\Phi\rangle}.
\]

分母之所以保留，是因为任意选择的系数 \(\mathbf c\) 不一定已经使试探波函数归一化。把

\[
|\Phi\rangle=\sum_\nu c_\nu|\chi_\nu\rangle,
\qquad
\langle\Phi|=\sum_\mu c_\mu^*\langle\chi_\mu|
\]

代入分子，得到

\[
\begin{aligned}
\langle\Phi|\hat H|\Phi\rangle
&=\left(\sum_\mu c_\mu^*\langle\chi_\mu|\right)
  \hat H
  \left(\sum_\nu c_\nu|\chi_\nu\rangle\right)\\
&=\sum_{\mu\nu}c_\mu^*c_\nu
  \langle\chi_\mu|\hat H|\chi_\nu\rangle\\
&=\sum_{\mu\nu}c_\mu^*H_{\mu\nu}c_\nu\\
&=\mathbf c^\dagger\mathbf H\mathbf c.
\end{aligned}
\]

同样地，分母为

\[
\begin{aligned}
\langle\Phi|\Phi\rangle
&=\sum_{\mu\nu}c_\mu^*c_\nu
  \langle\chi_\mu|\chi_\nu\rangle\\
&=\sum_{\mu\nu}c_\mu^*S_{\mu\nu}c_\nu\\
&=\mathbf c^\dagger\mathbf S\mathbf c.
\end{aligned}
\]

因此一般的复数形式是

\[
E(\mathbf c)
=\frac{\mathbf c^\dagger\mathbf H\mathbf c}
       {\mathbf c^\dagger\mathbf S\mathbf c}.
\]

本文假设基函数和展开系数都是实数，此时共轭转置 \(\dagger\) 就退化为普通转置 \({\mathrm T}\)，于是

\[
E(\mathbf c)
=\frac{\mathbf c^{\mathrm T}\mathbf H\mathbf c}
       {\mathbf c^{\mathrm T}\mathbf S\mathbf c}.
\]

这里 \(\mathbf c^{\mathrm T}\mathbf H\mathbf c\) 是 Hamiltonian 期望值的矩阵写法，而 \(\mathbf c^{\mathrm T}\mathbf S\mathbf c\) 是试探波函数的范数。重叠矩阵 \(S\) 的出现仅仅是因为基函数之间一般不正交；若基函数正交，则 \(S=I\)，分母变成 \(\mathbf c^{\mathrm T}\mathbf c\)。这个比值也称为 Rayleigh quotient（Rayleigh 商）。

### 3.3 用 Lagrange 乘子保持归一化

归一化条件是

\[
\mathbf c^{\mathrm T}\mathbf S\mathbf c=1.
\]

构造 Lagrange 函数

\[
\mathcal L
=\mathbf c^{\mathrm T}\mathbf H\mathbf c
-E\left(\mathbf c^{\mathrm T}\mathbf S\mathbf c-1\right).
\]

对 \(\mathbf c\) 求驻值：

\[
\frac{\partial\mathcal L}{\partial\mathbf c}=0.
\]

由于 \(H\) 和 \(S\) 都是对称矩阵，得到

\[
2\mathbf H\mathbf c-2E\mathbf S\mathbf c=0,
\]

即

\[
\boxed{\mathbf H\mathbf c=E\mathbf S\mathbf c}.
\]

这就是线性变分法产生的广义本征值问题。

如果基组正交，\(S=I\)，它退化为普通本征值问题

\[
\mathbf H\mathbf c=E\mathbf c.
\]

### 3.4 这与 Hartree–Fock 有什么关系

在线性变分法中，\(H\) 是给定的矩阵，因此求解一次本征值问题即可。

Hartree–Fock 中对应的矩阵是 Fock 矩阵 \(F\)，但 \(F\) 又依赖它自己的本征向量所生成的电子密度。因此每一次对角化仍然是线性的，但整个问题是非线性的，必须迭代到自洽。

---

## 4. 固定核多电子 Hamiltonian

### 4.1 Born–Oppenheimer 分离

固定原子核后，电子 Hamiltonian 为

\[
\hat H_\mathrm{elec}
=\sum_{i=1}^{N_e}\hat h(i)
+\sum_{i<j}^{N_e}\frac{1}{r_{ij}},
\]

其中一电子算符

\[
\hat h(i)
=-\frac12\nabla_i^2
-\sum_A\frac{Z_A}{r_{iA}}.
\]

它包含：

- 电子动能；
- 电子与所有原子核之间的吸引。

二电子项

\[
\frac{1}{r_{ij}}
\]

表示电子之间的 Coulomb 排斥。

原子核之间的排斥能是常数：

\[
E_\mathrm{nuc}
=\sum_{A<B}\frac{Z_AZ_B}{R_{AB}}.
\]

最终 Born–Oppenheimer 总能量为

\[
E_\mathrm{total}
=E_\mathrm{elec}+E_\mathrm{nuc}.
\]

Project 3 直接从 `enuc.dat` 读取 \(E_\mathrm{nuc}\)。

### 4.2 困难来自电子相互作用

如果没有 \(1/r_{ij}\)，Hamiltonian 是一电子算符之和，每个电子可以独立求解。

加入电子排斥后，一个电子的运动取决于所有其他电子的位置，多电子 Schrödinger 方程一般不能分离。Hartree–Fock 的核心近似是：用一个 Slater 行列式描述多电子波函数，让每个电子在其他电子产生的平均场中运动。

---

## 5. 从 Hartree 乘积到 Slater 行列式

### 5.1 Hartree 乘积的问题

最简单的独立电子波函数是

\[
\Psi_H(x_1,\ldots,x_{N_e})
=\chi_1(x_1)\chi_2(x_2)\cdots\chi_{N_e}(x_{N_e}).
\]

但电子是不可区分的费米子，交换两个电子时总波函数必须变号：

\[
\Psi(\ldots,x_i,\ldots,x_j,\ldots)
=-\Psi(\ldots,x_j,\ldots,x_i,\ldots).
\]

Hartree 乘积不满足这个要求。

### 5.2 Slater 行列式

用自旋轨道 \(\chi_p(x)\) 构造

\[
\Psi(x_1,\ldots,x_{N_e})
=\frac{1}{\sqrt{N_e!}}
\begin{vmatrix}
\chi_1(x_1)&\chi_2(x_1)&\cdots&\chi_{N_e}(x_1)\\
\chi_1(x_2)&\chi_2(x_2)&\cdots&\chi_{N_e}(x_2)\\
\vdots&\vdots&\ddots&\vdots\\
\chi_1(x_{N_e})&\chi_2(x_{N_e})&\cdots&\chi_{N_e}(x_{N_e})
\end{vmatrix}.
\]

交换两个电子等价于交换行列式的两行，因此自动变号。

如果两个电子占据完全相同的自旋轨道，行列式出现两列相同，波函数为零。这就是 Pauli 不相容原理在行列式中的体现。

---

## 6. Slater 行列式的能量

假设占据自旋轨道彼此正交。根据 Slater–Condon 规则，单行列式能量为

\[
E_\mathrm{elec}
=\sum_i\langle i|\hat h|i\rangle
+\frac12\sum_{i,j}
\left[(ii|jj)-(ij|ji)\right].
\]

定义二电子积分

\[
(pq|rs)
=\iint
\chi_p^*(x_1)\chi_q(x_1)
\frac{1}{r_{12}}
\chi_r^*(x_2)\chi_s(x_2)
\,dx_1dx_2.
\]

其中：

- \((ii|jj)\) 是 Coulomb 项，描述轨道电荷分布之间的经典排斥；
- \((ij|ji)\) 是 exchange 项，来自波函数的反对称性，没有简单的经典对应物。

当 \(i=j\) 时 Coulomb 和 exchange 项相消，因此电子不会与自己产生虚假的自相互作用。

---

## 7. Hartree–Fock 方程的变分来源

### 7.1 为什么需要一整组 Lagrange 乘子

Hartree–Fock 要在“所有允许的单个 Slater 行列式”中寻找能量最低者，因此需要改变占据轨道，同时保持它们正交：

\[
\langle\chi_i|\chi_j\rangle=\delta_{ij}.
\]

这里不只有每条轨道自身的归一化条件 \(\langle\chi_i|\chi_i\rangle=1\)，还有任意两条不同轨道之间的正交条件 \(\langle\chi_i|\chi_j\rangle=0\)。因此不能只为每条轨道引入一个标量乘子，而要为每一对 \((i,j)\) 引入一个乘子。把它们排列起来便得到 Lagrange 乘子矩阵 \(\boldsymbol\varepsilon\)：

\[
\mathcal L
=E_\mathrm{HF}
-\sum_{ij}\varepsilon_{ji}
\left(\langle\chi_i|\chi_j\rangle-\delta_{ij}\right).
\]

这里把乘子的指标写成 \(\varepsilon_{ji}\)，是为了使后面的矩阵方程具有自然的列向量形式；由于 \(i,j\) 都是求和哑指标，这只是记号约定。

### 7.2 对轨道作函数变分

第 6 节的单行列式能量是

\[
E_\mathrm{HF}
=\sum_i\langle i|\hat h|i\rangle
+\frac12\sum_{ij}\left[(ii|jj)-(ij|ji)\right].
\]

现在把 \(\chi_i(x)\) 和它的复共轭 \(\chi_i^*(x)\) 暂时当作独立变量，对 \(\chi_i^*(1)\) 作函数变分。能量泛函的导数可以整理为

\[
\frac{\delta E_\mathrm{HF}}
     {\delta\chi_i^*(1)}
=\hat f(1)\chi_i(1),
\]

而正交约束项的导数为

\[
\frac{\delta}{\delta\chi_i^*(1)}
\sum_{pq}\varepsilon_{qp}
\langle\chi_p|\chi_q\rangle
=\sum_j\varepsilon_{ji}\chi_j(1).
\]

令 \(\delta\mathcal L/\delta\chi_i^*(1)=0\)，得到

\[
\boxed{
\hat f\chi_i
=\sum_j\varepsilon_{ji}\chi_j
}.
\]

这个方程还不是通常的本征值方程。右侧可能含有所有占据轨道，表示在保持轨道彼此正交时，对第 \(i\) 条轨道的优化会与其他占据轨道耦合。此时一般的 \(\varepsilon_{ji}\) 只是正交约束的 Lagrange 乘子，不能立刻把每个对角元解释成轨道能量。

### 7.3 Fock 算符从哪里出现

对自旋轨道形式，Fock 算符为

\[
\hat f(1)
=\hat h(1)
+\sum_{j\in\mathrm{occ}}
\left[\hat J_j(1)-\hat K_j(1)\right].
\]

其中 \(\hat h\) 是一电子 Hamiltonian。变分一电子能量产生 \(\hat h\chi_i\)，变分二电子 Coulomb 项产生 \(\hat J_j\chi_i\)，变分 exchange 项产生 \(-\hat K_j\chi_i\)。这就是为什么对整个 HF 能量泛函求导后，可以把结果合并成一个 Fock 算符。

Coulomb 算符的作用是

\[
\hat J_j(1)\chi_i(1)
=\left[
\int\frac{|\chi_j(2)|^2}{r_{12}}dx_2
\right]\chi_i(1),
\]

exchange 算符的作用是

\[
\hat K_j(1)\chi_i(1)
=\left[
\int\frac{\chi_j^*(2)\chi_i(2)}{r_{12}}dx_2
\right]\chi_j(1).
\]

对于一组给定的占据轨道，\(\hat h\)、\(\hat J_j\) 和 \(\hat K_j\) 都是 Hermitian 算符，因此 \(\hat f\) 也是 Hermitian 算符。不过 \(\hat f\) 本身依赖占据轨道，所以 HF 方程仍然是非线性的。

### 7.4 为什么 \(\boldsymbol\varepsilon\) 是 Hermitian 矩阵

从变分方程

\[
\hat f\chi_i=\sum_j\varepsilon_{ji}\chi_j
\]

出发，在左侧乘 \(\langle\chi_k|\)。利用占据轨道正交性，得到

\[
\langle\chi_k|\hat f|\chi_i\rangle
=\sum_j\varepsilon_{ji}
\langle\chi_k|\chi_j\rangle
=\varepsilon_{ki}.
\]

因此 Lagrange 乘子就是 Fock 算符在占据轨道空间中的矩阵元：

\[
\boxed{
\varepsilon_{ki}
=\langle\chi_k|\hat f|\chi_i\rangle
}.
\]

由于 \(\hat f\) 是 Hermitian 算符，

\[
\begin{aligned}
\varepsilon_{ki}
&=\langle\chi_k|\hat f|\chi_i\rangle\\
&=\langle\chi_i|\hat f|\chi_k\rangle^*\\
&=\varepsilon_{ik}^*.
\end{aligned}
\]

这正是 Hermitian 矩阵的定义：

\[
\boxed{
\boldsymbol\varepsilon
=\boldsymbol\varepsilon^\dagger
}.
\]

对于本文采用的实轨道，Hermitian 就退化为实对称，即 \(\boldsymbol\varepsilon=\boldsymbol\varepsilon^{\mathrm T}\)。

### 7.5 为什么允许旋转占据轨道

设有 \(N_\mathrm{occ}\) 条占据轨道，并在这个占据空间内部作幺正变换：

\[
\chi'_p=\sum_i\chi_iU_{ip},
\qquad
U^\dagger U=I.
\]

新轨道仍然正交，因为

\[
\langle\chi'_p|\chi'_q\rangle
=(U^\dagger U)_{pq}
=\delta_{pq}.
\]

更重要的是，这个变换没有改变被占据的轨道子空间。相应的占据空间投影算符保持不变：

\[
\begin{aligned}
\hat P'
&=\sum_p|\chi'_p\rangle\langle\chi'_p|\\
&=\sum_{ij}|\chi_i\rangle
   (UU^\dagger)_{ij}
   \langle\chi_j|\\
&=\sum_i|\chi_i\rangle\langle\chi_i|
=\hat P.
\end{aligned}
\]

所以电子密度、Fock 算符和 HF 能量都不变。从 Slater 行列式本身来看，占据轨道的幺正变换只会让行列式乘上 \(\det U\)。由于 \(|\det U|=1\)，这至多是一个不影响任何可观测量的整体相位。

因此，HF 解真正确定的是“哪个 \(N_\mathrm{occ}\) 维轨道子空间被占据”，而不是这个子空间中唯一的一组轨道。我们可以在其中自由选择更方便的轨道表示。

### 7.6 从耦合方程到 canonical HF 方程

把占据轨道看成组成矩阵 \(\boldsymbol\chi\) 的列，变分方程可以写成

\[
\hat f\boldsymbol\chi
=\boldsymbol\chi\boldsymbol\varepsilon.
\]

作占据轨道旋转

\[
\boldsymbol\chi'
=\boldsymbol\chi U
\]

以后，Lagrange 乘子矩阵变为

\[
\boldsymbol\varepsilon'
=U^\dagger\boldsymbol\varepsilon U.
\]

谱定理保证任何 Hermitian 矩阵都存在一个幺正矩阵 \(U\)，使得

\[
U^\dagger\boldsymbol\varepsilon U
=\operatorname{diag}(\varepsilon_1,\varepsilon_2,\ldots).
\]

选择这个特殊的 \(U\) 后，原来相互耦合的轨道方程就逐列分离：

\[
\boxed{
\hat f\chi'_i
=\varepsilon_i\chi'_i
}.
\]

这组使 \(\boldsymbol\varepsilon\) 对角化的轨道称为 **canonical Hartree–Fock orbitals**。通常重新去掉撇号，写成

\[
\boxed{
\hat f\chi_i
=\varepsilon_i\chi_i
}.
\]

现在 \(\varepsilon_i\) 才是 canonical 轨道能量。canonical 轨道不是另一种 HF 近似，而只是同一个最优占据空间的一种特别方便的表示：在这个表示中，Fock 算符在轨道基中是对角的。

在有限 AO 基组的实际程序中，我们通常直接对整个 Fock 矩阵求解本征值问题。这样不仅得到占据 canonical 轨道，也得到一组未占据的 virtual canonical 轨道。占据轨道用于构造密度，而虚轨道不进入 HF 基态的 Slater 行列式。

最后要注意这个循环依赖：\(\hat f\) 由所有占据轨道构造，而 canonical 占据轨道又必须是 \(\hat f\) 的最低本征函数。也就是

```text
占据轨道 → 密度与 Fock 算符 → Fock 本征轨道 → 新的占据轨道
```

这正是必须进行 self-consistent field 迭代的根本原因。

---

## 8. 闭壳层 RHF

### 8.1 空间轨道双占据

RHF 假设每个占据空间轨道 \(\phi_i(\mathbf r)\) 同时容纳一个 \(\alpha\) 电子和一个 \(\beta\) 电子：

\[
\chi_{i\alpha}(x)=\phi_i(\mathbf r)\alpha(\omega),
\qquad
\chi_{i\beta}(x)=\phi_i(\mathbf r)\beta(\omega).
\]

若电子数为 \(N_e\)，占据空间轨道数为

\[
N_\mathrm{occ}=\frac{N_e}{2}.
\]

### 8.2 RHF 能量

对自旋积分后，RHF 电子能量为

\[
E_\mathrm{RHF}
=2\sum_i h_{ii}
+\sum_{i,j}
\left(2J_{ij}-K_{ij}\right),
\]

其中

\[
h_{ii}=\langle\phi_i|\hat h|\phi_i\rangle,
\]

\[
J_{ij}=(ii|jj),
\qquad
K_{ij}=(ij|ji).
\]

空间轨道 Fock 算符为

\[
\hat f
=\hat h
+\sum_{j\in\mathrm{occ}}
\left(2\hat J_j-\hat K_j\right).
\]

注意：轨道能量之和不是总电子能量，因为每一条轨道能量已经包含它与其他占据轨道的平均相互作用，直接相加会重复计算电子–电子作用。

---

## 9. LCAO：用 AO 基组展开分子轨道

### 9.1 分子轨道展开

实际计算中不能直接在无限维函数空间中优化轨道。选择有限组 AO 基函数 \(\{\chi_\mu\}\)，写成

\[
\phi_i=\sum_{\mu=1}^{K}C_{\mu i}\chi_\mu.
\]

\(C_{\mu i}\) 表示第 \(i\) 个分子轨道中第 \(\mu\) 个 AO 基函数的系数。矩阵 \(C\) 的每一列是一条分子轨道。

AO 基函数通常不是正交的，因此定义

\[
S_{\mu\nu}
=\langle\chi_\mu|\chi_\nu\rangle.
\]

分子轨道正交条件变成

\[
\boxed{C^{\mathrm T}SC=I}.
\]

它不是 \(C^{\mathrm T}C=I\)，因为系数向量所处的 AO 坐标系具有度量矩阵 \(S\)。

### 9.2 投影 Hartree–Fock 方程

把 LCAO 展开代入

\[
\hat f\phi_i=\varepsilon_i\phi_i,
\]

得到

\[
\hat f\sum_\nu C_{\nu i}\chi_\nu
=\varepsilon_i\sum_\nu C_{\nu i}\chi_\nu.
\]

左乘 \(\langle\chi_\mu|\)，得到

\[
\sum_\nu
\langle\chi_\mu|\hat f|\chi_\nu\rangle C_{\nu i}
=\varepsilon_i
\sum_\nu
\langle\chi_\mu|\chi_\nu\rangle C_{\nu i}.
\]

定义

\[
F_{\mu\nu}
=\langle\chi_\mu|\hat f|\chi_\nu\rangle,
\]

就得到

\[
\sum_\nu F_{\mu\nu}C_{\nu i}
=\varepsilon_i\sum_\nu S_{\mu\nu}C_{\nu i}.
\]

矩阵形式是

\[
\boxed{FC=SC\varepsilon}.
\]

这就是 Roothaan–Hall 方程。

它本质上与第 3 节的线性变分方程相同：非正交基组必然带来重叠矩阵 \(S\)。不同之处是 Fock 矩阵依赖未知系数 \(C\)，所以必须自洽迭代。

---

## 10. AO 矩阵和双电子积分

### 10.1 一电子积分

重叠矩阵：

\[
S_{\mu\nu}
=\int\chi_\mu(\mathbf r)\chi_\nu(\mathbf r)d\mathbf r.
\]

动能矩阵：

\[
T_{\mu\nu}
=\int\chi_\mu(\mathbf r)
\left(-\frac12\nabla^2\right)
\chi_\nu(\mathbf r)d\mathbf r.
\]

核吸引矩阵：

\[
V_{\mu\nu}
=\int\chi_\mu(\mathbf r)
\left(-\sum_A\frac{Z_A}{r_A}\right)
\chi_\nu(\mathbf r)d\mathbf r.
\]

把电子动能与电子–原子核吸引合在一起，得到**一电子核心 Hamiltonian 矩阵**（one-electron core Hamiltonian matrix）：

\[
\boxed{H^\mathrm{core}=T+V}.
\]

这里的 “core Hamiltonian” 是量子化学中的固定术语，但很容易被误解。它不是“原子核的 Hamiltonian”，也不描述原子核的运动；更明确的中文说法是“一电子核心 Hamiltonian”或“电子的一电子 Hamiltonian”。它来自一电子算符

\[
\hat h(i)
=-\frac12\nabla_i^2
-\sum_A\frac{Z_A}{r_{iA}},
\]

其中两部分都作用在电子坐标 \(\mathbf r_i\) 上：

- \(-\tfrac12\nabla_i^2\) 是第 \(i\) 个电子的动能算符；
- \(-\sum_A Z_A/r_{iA}\) 是该电子受到所有固定原子核吸引的势能算符。

相应的 AO 矩阵元为

\[
H^\mathrm{core}_{\mu\nu}
=\langle\chi_\mu|\hat h|\chi_\nu\rangle
=T_{\mu\nu}+V_{\mu\nu}.
\]

名称中的 “core” 可以理解为构造 Fock 矩阵时预先存在的“一电子核心部分”：

\[
F=H^\mathrm{core}+G[D],
\]

其中 \(H^\mathrm{core}\) 不依赖电子密度，而 \(G[D]\) 才是由电子–电子 Coulomb 和 exchange 作用产生的平均场部分。

真正与原子核自由度有关的项需要另行区分。在完整的分子 Hamiltonian 中有原子核动能；采用 Born–Oppenheimer 近似并固定原子核后，核动能不再出现在电子 Hamiltonian 中，而核–核排斥能

\[
E_\mathrm{nuc}
=\sum_{A<B}\frac{Z_AZ_B}{R_{AB}}
\]

作为一个与电子坐标无关的常数单独加入总能量。因此不能把 \(H^\mathrm{core}=T+V\) 称为“核运动的 Hamiltonian”或简称为容易产生歧义的“核 Hamiltonian”。

### 10.2 双电子排斥积分

Project 3 使用 Mulliken/chemist 记号：

\[
(\mu\nu|\lambda\sigma)
=\iint
\chi_\mu(1)\chi_\nu(1)
\frac{1}{r_{12}}
\chi_\lambda(2)\chi_\sigma(2)
\,d\mathbf r_1d\mathbf r_2.
\]

实 AO 基函数下具有八重置换对称性：

\[
(\mu\nu|\lambda\sigma)
=(\nu\mu|\lambda\sigma)
=(\mu\nu|\sigma\lambda)
=(\lambda\sigma|\mu\nu)
=\cdots.
\]

因此文件只需给出置换唯一的积分。

### 10.3 Project 3 的密度矩阵约定

本项目定义

\[
\boxed{
D_{\mu\nu}
=\sum_{i=1}^{N_\mathrm{occ}}C_{\mu i}C_{\nu i}
}.
\]

这里没有显式的占据数 2。为了看清 Fock 矩阵公式的来源，先回到第 8 节的闭壳层空间轨道 Fock 算符：

\[
\hat f
=\hat h
+\sum_{i\in\mathrm{occ}}
\left(2\hat J_i-\hat K_i\right).
\]

其中 \(i\) 遍历占据的**空间轨道**，而不是逐个遍历电子。每条占据空间轨道中有一个 \(\alpha\) 电子和一个 \(\beta\) 电子，因此 Coulomb 作用出现两次，给出 \(2\hat J_i\)；exchange 只发生在相同自旋的电子之间，对正在求解的某一自旋轨道只出现一次，给出 \(-\hat K_i\)。这就是 \(2J-K\) 的来源。

Fock 矩阵元是 Fock 算符在 AO 基组中的投影：

\[
F_{\mu\nu}
=\langle\chi_\mu|\hat f|\chi_\nu\rangle.
\]

把算符的三部分分别投影，得到

\[
F_{\mu\nu}
=H^\mathrm{core}_{\mu\nu}
+\sum_{i\in\mathrm{occ}}
\left(2J^i_{\mu\nu}-K^i_{\mu\nu}\right),
\]

其中

\[
H^\mathrm{core}_{\mu\nu}
=\langle\chi_\mu|\hat h|\chi_\nu\rangle.
\]

下面分别把 \(J^i_{\mu\nu}\) 和 \(K^i_{\mu\nu}\) 写成 AO 双电子积分。

#### Coulomb 部分

第 \(i\) 条占据空间轨道展开为

\[
\phi_i=\sum_\lambda C_{\lambda i}\chi_\lambda.
\]

先回忆 Coulomb 算符的作用定义。对任意待作用函数 \(g(1)\)，

\[
\boxed{
\hat J_i(1)g(1)
=\left[
\int\frac{|\phi_i(2)|^2}{r_{12}}d\mathbf r_2
\right]g(1)
}.
\]

方括号中的积分是占据轨道 \(\phi_i\) 的电子密度在电子 1 所在位置产生的 Coulomb 势。\(\hat J_i\) 作用在 \(g(1)\) 上，就是先计算这个势，再把它乘到 \(g(1)\) 上。

现在令 \(g(1)=\chi_\nu(1)\)，得到

\[
\hat J_i(1)\chi_\nu(1)
=\left[
\int\frac{|\phi_i(2)|^2}{r_{12}}d\mathbf r_2
\right]\chi_\nu(1).
\]

矩阵元 \(J^i_{\mu\nu}\) 的定义是左乘 \(\langle\chi_\mu|\)，也就是再对电子 1 的坐标积分：

\[
\begin{aligned}
J^i_{\mu\nu}
&=\langle\chi_\mu|\hat J_i|\chi_\nu\rangle\\
&=\int
\chi_\mu^*(1)
\left[
\hat J_i(1)\chi_\nu(1)
\right]d\mathbf r_1\\
&=\int
\chi_\mu^*(1)
\left[
\int\frac{|\phi_i(2)|^2}{r_{12}}d\mathbf r_2
\right]
\chi_\nu(1)d\mathbf r_1\\
&=\iint
\chi_\mu^*(1)\chi_\nu(1)
\frac{|\phi_i(2)|^2}{r_{12}}
\,d\mathbf r_1d\mathbf r_2.
\end{aligned}
\]

这解释了第一个等号：它只是把 bra–ket 矩阵元和 \(\hat J_i\) 的作用定义写成显式积分。

接着处理积分中的 \(|\phi_i(2)|^2\)。由于

\[
\phi_i(2)
=\sum_\sigma C_{\sigma i}\chi_\sigma(2),
\qquad
\phi_i^*(2)
=\sum_\lambda C_{\lambda i}^*\chi_\lambda^*(2),
\]

所以

\[
\begin{aligned}
|\phi_i(2)|^2
&=\phi_i^*(2)\phi_i(2)\\
&=\left(\sum_\lambda
C_{\lambda i}^*\chi_\lambda^*(2)\right)
\left(\sum_\sigma
C_{\sigma i}\chi_\sigma(2)\right)\\
&=\sum_{\lambda\sigma}
C_{\lambda i}^*C_{\sigma i}
\chi_\lambda^*(2)\chi_\sigma(2).
\end{aligned}
\]

这里出现两个下标 \(\lambda,\sigma\)，是因为 \(|\phi_i|^2\) 是两个 AO 线性组合的乘积；这与 \((a+b)^2=a^2+ab+ba+b^2\) 会产生所有两两组合是同一件事。

将它代回 Coulomb 矩阵元，并把与积分坐标无关的系数移到积分号外：

\[
\begin{aligned}
J^i_{\mu\nu}
&=\sum_{\lambda\sigma}
C_{\lambda i}^*C_{\sigma i}
\iint
\chi_\mu^*(1)\chi_\nu(1)
\frac{1}{r_{12}}
\chi_\lambda^*(2)\chi_\sigma(2)
\,d\mathbf r_1d\mathbf r_2\\
&=\sum_{\lambda\sigma}
C_{\lambda i}^*C_{\sigma i}
(\mu\nu|\lambda\sigma).
\end{aligned}
\]

最后一个等号直接使用了 AO 双电子积分的定义：

\[
(\mu\nu|\lambda\sigma)
=\iint
\chi_\mu^*(1)\chi_\nu(1)
\frac{1}{r_{12}}
\chi_\lambda^*(2)\chi_\sigma(2)
\,d\mathbf r_1d\mathbf r_2.
\]

本文和 Project 3 假设 AO、分子轨道和系数都是实数，因此星号可以去掉，最终得到

\[
\boxed{
J^i_{\mu\nu}
=\sum_{\lambda\sigma}
C_{\lambda i}C_{\sigma i}
(\mu\nu|\lambda\sigma)
}.
\]

这里 \(\mu,\nu\) 同时出现在电子 1 的坐标上，而 \(\lambda,\sigma\) 同时出现在电子 2 的坐标上，所以得到 Coulomb 排列 \((\mu\nu|\lambda\sigma)\)。

#### Exchange 部分

继续采用 Project 3 的实轨道假设。对任意待作用函数 \(g(1)\)，exchange 算符定义为

\[
\boxed{
\hat K_i(1)g(1)
=\left[
\int
\frac{\phi_i(2)g(2)}{r_{12}}
d\mathbf r_2
\right]\phi_i(1)
}.
\]

它与 Coulomb 算符有一个重要区别。Coulomb 算符使用密度 \(|\phi_i(2)|^2\) 产生一个普通的乘法势；exchange 算符则先把待作用函数 \(g(2)\) 与占据轨道 \(\phi_i(2)\) 配对积分，再让结果乘上电子 1 处的 \(\phi_i(1)\)。因此 \(\hat K_i\) 的输出不仅取决于位置 1，也取决于 \(g\) 在所有位置 2 的取值；它是一个非局域算符。

现在令 \(g(1)=\chi_\nu(1)\)：

\[
\hat K_i(1)\chi_\nu(1)
=\left[
\int
\frac{\phi_i(2)\chi_\nu(2)}{r_{12}}
d\mathbf r_2
\right]\phi_i(1).
\]

矩阵元 \(K^i_{\mu\nu}\) 的定义是左乘 \(\langle\chi_\mu|\)，即再对电子 1 的坐标积分：

\[
\begin{aligned}
K^i_{\mu\nu}
&=\langle\chi_\mu|\hat K_i|\chi_\nu\rangle\\
&=\int
\chi_\mu(1)
\left[\hat K_i(1)\chi_\nu(1)\right]
d\mathbf r_1\\
&=\int
\chi_\mu(1)
\left[
\int
\frac{\phi_i(2)\chi_\nu(2)}{r_{12}}
d\mathbf r_2
\right]
\phi_i(1)d\mathbf r_1\\
&=\iint
\chi_\mu(1)\phi_i(1)
\frac{1}{r_{12}}
\chi_\nu(2)\phi_i(2)
\,d\mathbf r_1d\mathbf r_2.
\end{aligned}
\]

这一步只是把 bra–ket 矩阵元和 \(\hat K_i\) 的作用定义写成显式的二重积分。

接下来，积分中出现了两次占据分子轨道，但它们位于不同的电子坐标上：

\[
\phi_i(1)
=\sum_\lambda C_{\lambda i}\chi_\lambda(1),
\qquad
\phi_i(2)
=\sum_\sigma C_{\sigma i}\chi_\sigma(2).
\]

将这两个展开式代入：

\[
\begin{aligned}
K^i_{\mu\nu}
&=\iint
\chi_\mu(1)
\left(\sum_\lambda
C_{\lambda i}\chi_\lambda(1)\right)
\frac{1}{r_{12}}
\chi_\nu(2)
\left(\sum_\sigma
C_{\sigma i}\chi_\sigma(2)\right)
\,d\mathbf r_1d\mathbf r_2\\
&=\sum_{\lambda\sigma}
C_{\lambda i}C_{\sigma i}
\iint
\chi_\mu(1)\chi_\lambda(1)
\frac{1}{r_{12}}
\chi_\nu(2)\chi_\sigma(2)
\,d\mathbf r_1d\mathbf r_2.
\end{aligned}
\]

按照 AO 双电子积分的定义，积分部分正是

\[
(\mu\lambda|\nu\sigma)
=\iint
\chi_\mu(1)\chi_\lambda(1)
\frac{1}{r_{12}}
\chi_\nu(2)\chi_\sigma(2)
\,d\mathbf r_1d\mathbf r_2.
\]

所以

\[
\boxed{
K^i_{\mu\nu}
=\sum_{\lambda\sigma}
C_{\lambda i}C_{\sigma i}
(\mu\lambda|\nu\sigma)
}.
\]

现在可以直接从被积函数读出指标排列：

- 在电子 1 的坐标上相乘的是 \(\chi_\mu(1)\chi_\lambda(1)\)，形成第一对指标 \(\mu\lambda\)；
- 在电子 2 的坐标上相乘的是 \(\chi_\nu(2)\chi_\sigma(2)\)，形成第二对指标 \(\nu\sigma\)。

因此 exchange 积分是 \((\mu\lambda|\nu\sigma)\)，而不是 Coulomb 项中的 \((\mu\nu|\lambda\sigma)\)。从直观上说，Coulomb 项把两个 AO 矩阵指标 \(\mu,\nu\) 留在同一个电子坐标上，而 exchange 项把 \(\nu\) 与一个占据轨道指标交换到了电子 2 的配对中。这正是实现 Fock 构造时最需要注意的指标交换。

把两部分代回：

\[
\begin{aligned}
F_{\mu\nu}
&=H^\mathrm{core}_{\mu\nu}
+\sum_{i\in\mathrm{occ}}
 \sum_{\lambda\sigma}
C_{\lambda i}C_{\sigma i}
\left[
2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)
\right]\\
&=H^\mathrm{core}_{\mu\nu}
+\sum_{\lambda\sigma}
\left(\sum_{i\in\mathrm{occ}}
C_{\lambda i}C_{\sigma i}\right)
\left[
2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)
\right].
\end{aligned}
\]

括号里的占据轨道求和正是本项目定义的密度矩阵元：

\[
\sum_{i\in\mathrm{occ}}
C_{\lambda i}C_{\sigma i}
=D_{\lambda\sigma}.
\]

因此最终得到

\[
\boxed{
F_{\mu\nu}
=H^\mathrm{core}_{\mu\nu}
+\sum_{\lambda\sigma}
D_{\lambda\sigma}
\left[
2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)
\right]
}.
\]

可以把这个推导压缩记成

```text
空间轨道 Fock 算符  f = h + Σocc(2J - K)
          ↓ 投影到 AO 基组
Fμν = Hcoreμν + Σocc Σλσ Cλi Cσi [2(μν|λσ) - (μλ|νσ)]
          ↓ 用 Dλσ = Σocc Cλi Cσi 收拢占据轨道求和
Fμν = Hcoreμν + Σλσ Dλσ [2(μν|λσ) - (μλ|νσ)]
```

#### 电子能量公式的来源

不能把 Fock 矩阵简单地当作多电子 Hamiltonian，然后直接把占据轨道的 Fock 期望值相加，因为每条轨道的 Fock 算符中已经包含它与其他占据轨道的平均相互作用。直接相加会把每一对电子的相互作用重复计算。

应从第 8 节的闭壳层 RHF 电子能量出发：

\[
E_\mathrm{RHF}
=2\sum_{i\in\mathrm{occ}}h_{ii}
+\sum_{i,j\in\mathrm{occ}}
\left(2J_{ij}-K_{ij}\right).
\]

先处理一电子部分。因为

\[
h_{ii}
=\langle\phi_i|\hat h|\phi_i\rangle
=\sum_{\mu\nu}
C_{\mu i}C_{\nu i}
H^\mathrm{core}_{\mu\nu},
\]

对占据空间轨道求和后得到

\[
\begin{aligned}
2\sum_i h_{ii}
&=2\sum_{\mu\nu}
\left(\sum_iC_{\mu i}C_{\nu i}\right)
H^\mathrm{core}_{\mu\nu}\\
&=2\sum_{\mu\nu}
D_{\mu\nu}H^\mathrm{core}_{\mu\nu}.
\end{aligned}
\]

这里的因子 2 来自每条占据空间轨道容纳两个电子，而本项目定义的 \(D\) 本身不包含这个占据数 2。

再定义 Fock 矩阵中的二电子平均场部分

\[
\boxed{
G_{\mu\nu}[D]
=\sum_{\lambda\sigma}D_{\lambda\sigma}
\left[
2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)
\right]
},
\]

于是刚才推导的 Fock 矩阵可以简写为

\[
\boxed{
F=H^\mathrm{core}+G[D]
}.
\]

为了看清二电子部分，先分别展开轨道间 Coulomb 积分和 exchange 积分：

\[
\begin{aligned}
J_{ij}
&=(ii|jj)\\
&=\sum_{\mu\nu\lambda\sigma}
C_{\mu i}C_{\nu i}
C_{\lambda j}C_{\sigma j}
(\mu\nu|\lambda\sigma),
\end{aligned}
\]

\[
\begin{aligned}
K_{ij}
&=(ij|ji)\\
&=\sum_{\mu\nu\lambda\sigma}
C_{\mu i}C_{\lambda i}
C_{\nu j}C_{\sigma j}
(\mu\lambda|\nu\sigma).
\end{aligned}
\]

对所有占据轨道 \(i,j\) 求和时，属于轨道 \(i\) 的两个系数收拢成一个密度矩阵元，属于轨道 \(j\) 的两个系数收拢成另一个密度矩阵元。例如

\[
\sum_{ij}
C_{\mu i}C_{\nu i}
C_{\lambda j}C_{\sigma j}
=\left(\sum_iC_{\mu i}C_{\nu i}\right)
 \left(\sum_jC_{\lambda j}C_{\sigma j}\right)
=D_{\mu\nu}D_{\lambda\sigma}.
\]

因此 RHF 能量中的整个二电子部分为

\[
\begin{aligned}
\sum_{ij}(2J_{ij}-K_{ij})
&=\sum_{\mu\nu\lambda\sigma}
D_{\mu\nu}D_{\lambda\sigma}
\left[
2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)
\right]\\
&=\sum_{\mu\nu}
D_{\mu\nu}G_{\mu\nu}[D].
\end{aligned}
\]

因此总电子能量是

\[
\begin{aligned}
E_\mathrm{elec}
&=2\sum_{\mu\nu}
D_{\mu\nu}H^\mathrm{core}_{\mu\nu}
+\sum_{\mu\nu}
D_{\mu\nu}G_{\mu\nu}[D]\\
&=\sum_{\mu\nu}D_{\mu\nu}
\left(2H^\mathrm{core}_{\mu\nu}
+G_{\mu\nu}[D]\right).
\end{aligned}
\]

最后利用

\[
F_{\mu\nu}
=H^\mathrm{core}_{\mu\nu}+G_{\mu\nu}[D],
\]

可得

\[
2H^\mathrm{core}_{\mu\nu}+G_{\mu\nu}[D]
=H^\mathrm{core}_{\mu\nu}+F_{\mu\nu}.
\]

所以

\[
\boxed{
E_\mathrm{elec}
=\sum_{\mu\nu}
D_{\mu\nu}
\left(
H^\mathrm{core}_{\mu\nu}+F_{\mu\nu}
\right)
}.
\]

这个公式的结构很值得记住：

\[
\underbrace{H^\mathrm{core}+F}_{H^\mathrm{core}+(H^\mathrm{core}+G)}
=2H^\mathrm{core}+G.
\]

也就是说，一电子能量按双占据计算两次，而电子–电子相互作用只计算一次。若错误地使用 \(2DF\)，就会得到 \(2DH+2DG\)，把二电子相互作用重复计算；若只使用 \(DF\)，又会得到 \(DH+DG\)，漏掉一半的一电子贡献。

另一种教材和程序常定义

\[
P_{\mu\nu}=2D_{\mu\nu}.
\]

此时公式变为

\[
F_{\mu\nu}
=H^\mathrm{core}_{\mu\nu}
+\sum_{\lambda\sigma}P_{\lambda\sigma}
\left[
(\mu\nu|\lambda\sigma)
-\frac12(\mu\lambda|\nu\sigma)
\right],
\]

\[
E_\mathrm{elec}
=\frac12\sum_{\mu\nu}P_{\mu\nu}
\left(H^\mathrm{core}_{\mu\nu}+F_{\mu\nu}\right).
\]

两种写法完全等价，但绝对不能把一种密度定义和另一种能量/Fock 公式混合。

---

## 11. 怎样求解 Roothaan–Hall 方程

### 11.1 为什么不能直接对角化 F

我们要求解

\[
FC=SC\varepsilon,
\]

而 Eigen 的 `SelfAdjointEigenSolver` 求解的是

\[
AC=C\varepsilon.
\]

只有在 \(S=I\) 时才能直接对角化 \(F\)。AO 基组通常非正交，所以要先构造一个正交化变换。

### 11.2 对称正交化

对重叠矩阵对角化：

\[
S=UsU^{\mathrm T},
\]

其中 \(s\) 是正的重叠本征值组成的对角矩阵。

定义

\[
\boxed{X=S^{-1/2}=Us^{-1/2}U^{\mathrm T}}.
\]

可以验证

\[
X^{\mathrm T}SX=I.
\]

把 Fock 矩阵变换到正交 AO 基：

\[
F'=X^{\mathrm T}FX.
\]

求解普通对称本征值问题

\[
F'C'=C'\varepsilon.
\]

再变换回原 AO 基：

\[
\boxed{C=XC'}.
\]

于是

\[
C^{\mathrm T}SC
=C'^{\mathrm T}X^{\mathrm T}SXC'
=C'^{\mathrm T}C'
=I.
\]

### 11.3 近线性相关

如果 \(S\) 的某些本征值非常小，\(s^{-1/2}\) 会非常大，放大数值误差。这说明 AO 基组接近线性相关。

教学项目的测试数据通常不会出现严重问题，但完整程序应设置阈值，删除或特殊处理过小的重叠本征值。

---

## 12. 为什么需要 SCF

Fock 矩阵依赖密度：

\[
F=F[D].
\]

密度依赖占据分子轨道：

\[
D=D[C_\mathrm{occ}].
\]

分子轨道又来自 Fock 矩阵的本征向量：

\[
C=C[F].
\]

所以我们要求的是不动点

\[
D^*=D[C[F[D^*]]].
\]

SCF 就是通过迭代寻找这个不动点。

---

## 13. RHF/SCF 的完整流程

### 13.1 输入预处理

Project 3 读取：

```text
enuc.dat  → E_nuc
s.dat     → S
t.dat     → T
v.dat     → V
eri.dat   → (μν|λσ)
```

构造

\[
H^\mathrm{core}=T+V.
\]

一电子积分文件只包含对称矩阵的唯一元素，文件下标从 1 开始；放进 C++/Eigen 前要减 1，并同时填充上下三角。

ERI 文件只包含八重对称下的唯一积分，可以压缩存入一维数组。

### 13.2 一电子核心 Hamiltonian 初猜

构造正交化矩阵 \(X=S^{-1/2}\)，然后使用

\[
F'_0=X^{\mathrm T}H^\mathrm{core}X
\]

作为第一次轨道求解的矩阵。

对角化：

\[
F'_0C'_0=C'_0\varepsilon_0,
\]

回到 AO 基：

\[
C_0=XC'_0.
\]

取最低的 \(N_\mathrm{occ}\) 列构造初始密度：

\[
D^0_{\mu\nu}
=\sum_{i=1}^{N_\mathrm{occ}}
C^0_{\mu i}C^0_{\nu i}.
\]

### 13.3 一次 SCF 迭代

给定旧密度 \(D^{n-1}\)：

1. 构造 Fock 矩阵

\[
F^n_{\mu\nu}
=H^\mathrm{core}_{\mu\nu}
+\sum_{\lambda\sigma}D^{n-1}_{\lambda\sigma}
\left[2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)\right].
\]

2. 变换到正交基

\[
F'^n=X^{\mathrm T}F^nX.
\]

3. 对角化

\[
F'^nC'^n=C'^n\varepsilon^n.
\]

Eigen 通常按特征值从小到大返回，因此最低的轨道首先出现；代码仍应明确把前 \(N_\mathrm{occ}\) 列视为占据轨道。

4. 回到原 AO 基

\[
C^n=XC'^n.
\]

5. 构造新密度

\[
D^n_{\mu\nu}
=\sum_{i=1}^{N_\mathrm{occ}}C^n_{\mu i}C^n_{\nu i}.
\]

6. 计算电子能和总能量

\[
E^n_\mathrm{elec}
=\sum_{\mu\nu}D^n_{\mu\nu}
\left(H^\mathrm{core}_{\mu\nu}+F^n_{\mu\nu}\right),
\]

\[
E^n_\mathrm{total}
=E^n_\mathrm{elec}+E_\mathrm{nuc}.
\]

这里必须使用与能量所对应的 Fock 和密度。实现时应清楚规定一次迭代中“旧密度构造 F、新轨道构造新密度、何时计算能量”的顺序，避免把不同迭代的数据混在一起。

### 13.4 收敛判断

能量变化：

\[
\Delta E^n
=E^n_\mathrm{elec}-E^{n-1}_\mathrm{elec}.
\]

密度变化可使用 Frobenius 范数：

\[
\Delta D_\mathrm{F}
=\left[
\sum_{\mu\nu}
\left(D^n_{\mu\nu}-D^{n-1}_{\mu\nu}\right)^2
\right]^{1/2}.
\]

有些程序将它再除以矩阵元素数的平方根并称为 RMS。Project 3 的参考公式没有这个除数，因此比较参考输出时要沿用项目定义。

收敛条件通常同时要求

\[
|\Delta E^n|<\delta_E,
\qquad
\Delta D_\mathrm{F}<\delta_D.
\]

还应设置最大迭代次数，避免不收敛时无限循环。

### 13.5 伪代码

```text
read E_nuc, S, T, V, ERI
H_core = T + V
X = symmetric_orthogonalizer(S)

C = diagonalize(Xᵀ H_core X)
C = X C
D = density_from_occupied_orbitals(C)
E_old = initial_energy

for iteration = 1 ... max_iterations:
    F = build_fock(H_core, ERI, D)
    F_orth = Xᵀ F X
    epsilon, C_orth = diagonalize(F_orth)
    C_new = X C_orth
    D_new = density_from_occupied_orbitals(C_new)
    E_elec = electronic_energy(D_new, H_core, F)
    E_total = E_elec + E_nuc

    delta_E = E_elec - E_old
    delta_D = norm(D_new - D)

    print iteration, E_elec, E_total, delta_E, delta_D

    if abs(delta_E) < energy_tolerance
       and delta_D < density_tolerance:
        converged = true
        break

    D = D_new
    E_old = E_elec
```

---

## 14. 与 Project 3 十个步骤逐项对应

### Step 1：核排斥能

读取 `enuc.dat`，得到 \(E_\mathrm{nuc}\)。虽然它也能由几何和核电荷计算，但作业要求先读取给定结果。

### Step 2：一电子积分

读取 \(S,T,V\)，恢复完整对称矩阵并构造

\[
H^\mathrm{core}=T+V.
\]

### Step 3：双电子积分

读取 `eri.dat`，利用八重置换对称性进行压缩存储。

对从 0 开始的 AO 下标，可定义

\[
\operatorname{pair}(i,j)
=\frac{p(p+1)}{2}+q,
\qquad
p=\max(i,j),\ q=\min(i,j).
\]

然后令

\[
P=\operatorname{pair}(i,j),
\qquad
Q=\operatorname{pair}(k,l),
\]

最终 ERI 下标为

\[
\operatorname{eri}(i,j,k,l)
=\operatorname{pair}(P,Q).
\]

这样所有八个等价排列都映射到同一个位置。

### Step 4：正交化矩阵

对角化 \(S\)，构造

\[
X=S^{-1/2}.
\]

验证

\[
X^{\mathrm T}SX\approx I.
\]

### Step 5：初始密度

用一电子核心 Hamiltonian 作为 Fock 初猜，在正交基中对角化，回变换系数，并从占据轨道构造 \(D^0\)。

### Step 6：初始 SCF 能量

用项目指定的密度约定计算初始电子能与总能量，打印第 0 次迭代。

### Step 7：新 Fock 矩阵

用旧密度和 ERI 构造包含 Coulomb 与 exchange 贡献的新 \(F\)。这是四重 AO 下标最集中的部分。

### Step 8：新密度矩阵

计算

\[
F'=X^{\mathrm T}FX,
\]

对角化，回变换 \(C=XC'\)，再由占据轨道形成新密度。

### Step 9：新 SCF 能量

计算 \(E_\mathrm{elec}\) 和 \(E_\mathrm{total}\)。

### Step 10：收敛判断

同时检查能量和密度变化；未收敛则返回 Step 7。

---

## 15. 关键数值与物理检查

### 15.1 矩阵对称性

在实 RHF 中应满足

\[
S=S^{\mathrm T},\quad
H^\mathrm{core}=H^{\mathrm{core,T}},\quad
F=F^{\mathrm T},\quad
D=D^{\mathrm T}.
\]

明显不对称通常说明文件恢复、ERI 索引或 Fock 构造有错误。

### 15.2 正交化检查

\[
X^{\mathrm T}SX\approx I.
\]

### 15.3 分子轨道正交检查

\[
C^{\mathrm T}SC\approx I.
\]

### 15.4 电子数检查

由于 Project 3 的 \(D\) 不含因子 2：

\[
\boxed{2\operatorname{Tr}(DS)=N_e}.
\]

若采用含因子 2 的 \(P\)，则是

\[
\operatorname{Tr}(PS)=N_e.
\]

### 15.5 收敛后的 canonical MO 检查

把 AO Fock 矩阵变换到 MO 基：

\[
F_\mathrm{MO}=C^{\mathrm T}FC.
\]

对 canonical 轨道，收敛后它应接近对角矩阵，对角元是轨道能量 \(\varepsilon_i\)。

### 15.6 H₂O/STO-3G 参考能量

Project 3 参考输出最终给出

\[
E_\mathrm{elec}
\approx-82.944446990003\ E_h,
\]

\[
E_\mathrm{nuc}
\approx 8.0023670618\ E_h,
\]

\[
\boxed{
E_\mathrm{total}
\approx-74.942079928192\ E_h
}.
\]

这是第一版程序最重要的端到端验证目标。

---

## 16. 最容易混淆的问题

### 16.1 AO 基函数不是分子轨道

AO 基函数 \(\chi_\mu\) 是预先选定的展开函数；分子轨道 \(\phi_i\) 是它们的线性组合。SCF 优化的是系数 \(C_{\mu i}\)，不是改变输入基函数本身。

### 16.2 为什么有 S

不是因为 Hartree–Fock 特殊，而是因为 AO 基函数非正交。任何在非正交基组中的线性变分问题都会产生

\[
HC=SC\varepsilon.
\]

### 16.3 为什么一次对角化不够

因为

\[
F\rightarrow C\rightarrow D\rightarrow F.
\]

第一次一电子核心 Hamiltonian 对角化只是初猜，不是最终 Hartree–Fock 解。

### 16.4 轨道能量不等于总能量分解

\(\sum_i\varepsilon_i\) 重复包含电子相互作用。必须使用 RHF 能量表达式计算总电子能量。

### 16.5 密度矩阵的因子 2

这是 Project 3 最危险的约定问题。务必把下面三者当成一整套：

```text
D = Σocc C_i C_iᵀ
F = Hcore + ΣD(2J - K)
Eelec = ΣD(Hcore + F)
```

不要混入使用 \(P=2D\) 的公式。

### 16.6 Coulomb 与 exchange 的 ERI 下标不同

Fock 元素中：

\[
2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma).
\]

exchange 项不是简单地把 Coulomb 项乘一个系数；它交换了指标。下标写错时矩阵仍可能看起来对称，但最终能量会错误。

### 16.7 本征向量按列存储

Eigen 的本征向量矩阵中，每一列对应一个本征值。构造密度时要选取最低 \(N_\mathrm{occ}\) 个列向量。

### 16.8 收敛要使用绝对能量变化

参考输出可以保留带符号的 \(\Delta E\)，但判断时应使用

\[
|\Delta E|<\delta_E.
\]

---

## 17. Hartree–Fock 做了哪些近似

Hartree–Fock 并不是精确求解多电子 Schrödinger 方程。它限制波函数只能是一个 Slater 行列式，因此只完整包含由反对称性带来的 exchange，缺少通常所说的电子相关能。

在固定 Hamiltonian 和给定有限基组中：

- SCF 收敛给出该基组中某个 HF 驻点；
- 稳定的最低 RHF 解是单行列式空间中的变分最优解；
- 有限基组 RHF 能量通常高于完备基组 HF 极限；
- HF 能量又高于同一非相对论 Born–Oppenheimer Hamiltonian 的精确基态能量。

Project 3 的目标不是消除这些近似，而是亲手实现产生 HF 解的自洽场数值过程。

---

## 18. 最核心公式汇总

线性变分：

\[
Hc=ESc.
\]

固定核电子 Hamiltonian：

\[
\hat H_\mathrm{elec}
=\sum_i\hat h(i)+\sum_{i<j}\frac1{r_{ij}}.
\]

RHF 空间轨道方程：

\[
\hat f\phi_i=\varepsilon_i\phi_i.
\]

LCAO 展开：

\[
\phi_i=\sum_\mu C_{\mu i}\chi_\mu.
\]

Roothaan–Hall：

\[
\boxed{FC=SC\varepsilon}.
\]

对称正交化：

\[
X=S^{-1/2},
\qquad
F'=X^{\mathrm T}FX,
\qquad
C=XC'.
\]

Project 3 密度：

\[
D_{\mu\nu}
=\sum_{i\in\mathrm{occ}}C_{\mu i}C_{\nu i}.
\]

Fock 矩阵：

\[
F_{\mu\nu}
=H^\mathrm{core}_{\mu\nu}
+\sum_{\lambda\sigma}D_{\lambda\sigma}
\left[2(\mu\nu|\lambda\sigma)
-(\mu\lambda|\nu\sigma)\right].
\]

电子能量：

\[
E_\mathrm{elec}
=\sum_{\mu\nu}D_{\mu\nu}
\left(H^\mathrm{core}_{\mu\nu}+F_{\mu\nu}\right).
\]

总能量：

\[
E_\mathrm{total}=E_\mathrm{elec}+E_\mathrm{nuc}.
\]

整个 Project 3 可以浓缩成一句话：

> 在非正交 AO 基组中，用当前密度构造 Fock 矩阵，求解 Roothaan–Hall 方程得到新轨道和新密度，直到输入密度与输出密度自洽。

---

## 19. 推荐阅读和对照材料

- Project 3 正式说明：`project3-instructions.pdf`。
- Project 3 README、hints 和 H₂O/STO-3G 参考输出。
- Szabo and Ostlund, *Modern Quantum Chemistry*, Chapter 3。
- 编码时优先逐项对照 Project 3 的参考矩阵和每次迭代能量，不要只比较最终总能量。
