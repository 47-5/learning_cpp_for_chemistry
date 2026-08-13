# Project 2：谐振动分析的数学基础

## 1. 这个项目究竟在求什么

Project 2 的目标是从分子在某一几何构型附近的势能曲率出发，求出分子的简正振动模式及其谐振动频率。

计算链条可以概括为：

```text
平衡几何 + Cartesian Hessian + 原子质量
                    ↓
              质量加权 Hessian
                    ↓ 对角化
          特征值 λ 和特征向量 L
                    ↓
       角频率 ω = √λ，波数 ν̃ = ω/(2πc)
```

Hessian 告诉我们势能面在各个方向上“弯得有多厉害”，质量告诉我们原子受到同样的回复力后运动得有多快。两者缺一不可。

本文先建立经典小振动理论，再把它逐项对应到 Project 2。

## 2. 从一维弹簧开始

### 2.1 Hooke 定律

质量为 \(m\) 的小球连接在弹簧上，令 \(x=0\) 为平衡位置。平衡点附近的势能为

\[
V(x)=V(0)+\frac12kx^2,
\]

其中 \(k\) 是力常数。力等于势能的负梯度：

\[
F=-\frac{dV}{dx}=-kx.
\]

Newton 方程给出

\[
m\ddot{x}=-kx,
\]

也就是

\[
m\ddot{x}+kx=0.
\]

它的解可以写成

\[
x(t)=A\cos(\omega t+\phi),
\]

这里的 \(A\) 是**振幅**：它表示小球相对平衡位置能够达到的最大位移，单位与 \(x\) 相同。\(A\) 不是由弹簧常数 \(k\) 和质量 \(m\) 唯一决定的，而是由初始条件决定的。例如，若在 \(t=0\) 时把小球拉得更远再释放，\(A\) 就更大；但只要仍处于谐近似范围内，振动频率 \(\omega\) 并不会因此改变。

\(\phi\) 是初相位，它说明我们把计时起点选在振动周期的什么位置。于是这个解可以拆成两个部分来看：

\[
\underbrace{A}_{\text{振动有多大}}
\underbrace{\cos(\omega t+\phi)}_{\text{位移怎样随时间往复变化}}.
\]

这个区分对后面的多原子体系非常重要：在一维中，振动“有多大”只需要一个数 \(A\)；在 \(3N\) 维中，各个 Cartesian 坐标的振幅一般不同，因此要把 \(A\) 推广成一个振幅向量 \(\mathbf a\)。

代回运动方程得到

\[
\omega^2=\frac{k}{m},
\qquad
\omega=\sqrt{\frac{k}{m}}.
\]

这里要区分三个相关量：

| 名称 | 符号 | 单位 | 关系 |
|---|---:|---:|---:|
| 角频率 | \(\omega\) | rad s\(^{-1}\) | \(\omega=\sqrt{k/m}\) |
| 普通频率 | \(\nu\) | Hz | \(\nu=\omega/(2\pi)\) |
| 光谱波数 | \(\tilde\nu\) | cm\(^{-1}\) | \(\tilde\nu=\nu/c=\omega/(2\pi c)\) |

Project 2 最终打印的是光谱波数 \(\tilde\nu\)，而不是角频率 \(\omega\)。

### 2.2 一维模型已经揭示了核心

\[
\omega=\sqrt{\frac{\text{势能曲率}}{\text{质量}}}.
\]

- 曲率越大，势能井越陡，振动越快；
- 质量越大，相同回复力造成的加速度越小，振动越慢。

多原子分子的数学只是把一个坐标、一个质量和一个力常数推广成许多坐标、质量矩阵和 Hessian 矩阵。

## 3. 分子的 Cartesian 坐标

一个含 \(N\) 个原子的分子有 \(3N\) 个 Cartesian 坐标。把它们按顺序排成一个列向量：

\[
\mathbf{x}=
\begin{pmatrix}
x_1&y_1&z_1&x_2&y_2&z_2&\cdots&x_N&y_N&z_N
\end{pmatrix}^{\mathrm T}.
\]

设参考构型为 \(\mathbf{x}_0\)，定义位移

\[
\Delta\mathbf{x}=\mathbf{x}-\mathbf{x}_0.
\]

为了简化下标，后面用单个复合下标 \(p=1,\ldots,3N\) 表示“某个原子的某个 Cartesian 分量”。例如：

```text
p = 0 → 原子 0 的 x
p = 1 → 原子 0 的 y
p = 2 → 原子 0 的 z
p = 3 → 原子 1 的 x
```

因此坐标下标 \(p\) 所属的原子编号是整数除法 \(p/3\)。这正是后面质量加权时需要的映射。

## 4. 势能面的 Taylor 展开

分子的电子能量在 Born-Oppenheimer 近似下是核坐标的函数：

\[
V=V(\mathbf{x}).
\]

在参考构型 \(\mathbf{x}_0\) 附近作 Taylor 展开：

\[
V(\mathbf{x}_0+\Delta\mathbf{x})
=V(\mathbf{x}_0)
+\sum_p
\left(\frac{\partial V}{\partial x_p}\right)_0\Delta x_p
+\frac12\sum_{p,q}
\left(\frac{\partial^2V}{\partial x_p\partial x_q}\right)_0
\Delta x_p\Delta x_q
+\cdots.
\]

三个层次分别是：

- 零阶项：参考点的能量；
- 一阶项：梯度，决定该点是否还受净力；
- 二阶项：曲率，决定小位移后的回复力；
- 三阶及更高项：非谐性。

### 4.1 驻点为什么重要

若 \(\mathbf{x}_0\) 是势能面的驻点，则

\[
\left(\frac{\partial V}{\partial x_p}\right)_0=0
\]

对所有 \(p\) 都成立。势能展开于是简化为

\[
V\approx V_0+\frac12\sum_{p,q}H_{pq}\Delta x_p\Delta x_q,
\]

其中

\[
H_{pq}=\left(\frac{\partial^2V}{\partial x_p\partial x_q}\right)_0
\]

就是 Cartesian Hessian。

矩阵形式为

\[
V\approx V_0+\frac12\Delta\mathbf{x}^{\mathrm T}
\mathbf H\Delta\mathbf{x}.
\]

### 4.2 “谐近似”是什么意思

谐近似就是只保留到二阶项，丢弃三阶及更高阶项。因此势能被近似为一个多维抛物面。

它适用于：

- 位移较小；
- 几何接近平衡结构；
- 振动量子数不太高；
- 化学键没有接近断裂。

真实分子势能不是完美抛物线，所以谐频率与实验基频通常不完全一致。实际量化化学中常使用经验缩放因子修正系统误差，但这不属于本项目的核心算法。

## 5. Hessian 的含义与性质

### 5.1 对角元素和非对角元素

\[
H_{pp}=\frac{\partial^2V}{\partial x_p^2}
\]

描述沿单个坐标方向的曲率。

\[
H_{pq}=\frac{\partial^2V}{\partial x_p\partial x_q},\qquad p\ne q
\]

描述两个坐标之间的耦合。例如移动一个氢原子可能同时改变两条键和一个键角，因此不同坐标通常不是独立的。

### 5.2 Hessian 为什么对称

只要势能足够光滑，混合偏导可交换：

\[
\frac{\partial^2V}{\partial x_p\partial x_q}
=
\frac{\partial^2V}{\partial x_q\partial x_p}.
\]

所以

\[
\mathbf H^{\mathrm T}=\mathbf H.
\]

实对称矩阵有两个关键性质：

1. 所有特征值都为实数；
2. 可以选择一组彼此正交的特征向量。

这也是 C++ 中应使用 `Eigen::SelfAdjointEigenSolver`，而不是一般非对称矩阵求解器的原因。

### 5.3 Hessian 的单位

Project 2 的能量使用 hartree \(E_h\)，坐标使用 bohr \(a_0\)，因此

\[
[H_{pq}]=\frac{E_h}{a_0^2}.
\]

它是“能量除以长度平方”，也就是力常数的单位。

## 6. 从势能得到耦合运动方程

谐近似下，力向量是

\[
\mathbf f=-\nabla V=-\mathbf H\Delta\mathbf{x}.
\]

### 6.1 第二个等号是怎样得到的

谐近似下的势能为

\[
V(\Delta\mathbf x)
=V_0+
\frac12\Delta\mathbf x^{\mathrm T}
\mathbf H\Delta\mathbf x.
\]

为了看清楚对它求梯度后为什么得到 \(\mathbf H\Delta\mathbf x\)，先考虑只有两个位移坐标的情况：

\[
\Delta\mathbf x=
\begin{pmatrix}
\Delta x_1\\
\Delta x_2
\end{pmatrix},
\qquad
\mathbf H=
\begin{pmatrix}
H_{11}&H_{12}\\
H_{21}&H_{22}
\end{pmatrix}.
\]

将二次型完全展开：

\[
\begin{aligned}
V
=V_0+\frac12\bigl(&
H_{11}\Delta x_1^2
+H_{12}\Delta x_1\Delta x_2\\
&+H_{21}\Delta x_2\Delta x_1
+H_{22}\Delta x_2^2
\bigr).
\end{aligned}
\]

由于 Hessian 是对称矩阵，\(H_{12}=H_{21}\)，所以

\[
V=V_0+
\frac12H_{11}\Delta x_1^2
+H_{12}\Delta x_1\Delta x_2
+\frac12H_{22}\Delta x_2^2.
\]

分别对两个位移坐标求偏导：

\[
\frac{\partial V}{\partial\Delta x_1}
=H_{11}\Delta x_1+H_{12}\Delta x_2,
\]

\[
\frac{\partial V}{\partial\Delta x_2}
=H_{21}\Delta x_1+H_{22}\Delta x_2.
\]

把这两个偏导数组合成梯度向量：

\[
\begin{aligned}
\nabla V
&=
\begin{pmatrix}
H_{11}\Delta x_1+H_{12}\Delta x_2\\
H_{21}\Delta x_1+H_{22}\Delta x_2
\end{pmatrix}\\
&=
\begin{pmatrix}
H_{11}&H_{12}\\
H_{21}&H_{22}
\end{pmatrix}
\begin{pmatrix}
\Delta x_1\\
\Delta x_2
\end{pmatrix}\\
&=\mathbf H\Delta\mathbf x.
\end{aligned}
\]

因此，根据力是势能的负梯度，便有

\[
\boxed{
\mathbf f=-\nabla V=-\mathbf H\Delta\mathbf x
}.
\]

对一般的 \(3N\) 维分子，势能可以写成求和形式：

\[
V=V_0+
\frac12
\sum_{p=1}^{3N}\sum_{q=1}^{3N}
H_{pq}\Delta x_p\Delta x_q.
\]

对第 \(r\) 个位移坐标求偏导，利用 Hessian 的对称性可得

\[
\frac{\partial V}{\partial\Delta x_r}
=\sum_{q=1}^{3N}H_{rq}\Delta x_q.
\]

所以第 \(r\) 个力分量为

\[
\boxed{
f_r=-\sum_{q=1}^{3N}H_{rq}\Delta x_q
}.
\]

把全部 \(3N\) 个分量组合起来，就是矩阵形式

\[
\mathbf f=-\mathbf H\Delta\mathbf x.
\]

它正是一维 Hooke 定律

\[
F=-kx
\]

在多维、坐标之间可以相互耦合时的推广。负号表示势能升高的方向与力的方向相反；在稳定平衡点附近，它表现为把分子拉回平衡构型的回复力。

Cartesian 坐标的质量矩阵为

\[
\mathbf M=\operatorname{diag}(
m_1,m_1,m_1,
m_2,m_2,m_2,
\ldots,
m_N,m_N,m_N).
\]

每个原子的质量连续出现三次，因为该原子的 \(x,y,z\) 三个方向具有同一个质量。

Newton 方程写成

\[
\mathbf M\Delta\ddot{\mathbf{x}}
+\mathbf H\Delta\mathbf{x}=0.
\]

### 6.2 振幅向量 \(\mathbf a\) 是从哪里来的

一维弹簧的解写成

\[
x(t)=A\cos(\omega t+\phi).
\]

对含 \(N\) 个原子的分子，位移 \(\Delta\mathbf x(t)\) 有 \(3N\) 个分量。我们先寻找一种特别重要的集体运动：所有坐标都以**同一个频率** \(\omega\) 和同一个相位 \(\phi\) 往复变化，但每个坐标可以具有不同的振幅和方向。逐分量写就是

\[
\begin{aligned}
\Delta x_1(t)&=a_1\cos(\omega t+\phi),\\
\Delta x_2(t)&=a_2\cos(\omega t+\phi),\\
&\ \vdots\\
\Delta x_{3N}(t)&=a_{3N}\cos(\omega t+\phi).
\end{aligned}
\]

因为每一行都有相同的时间因子，可以把它提出来：

\[
\Delta\mathbf{x}(t)
=
\begin{pmatrix}
a_1\\a_2\\\vdots\\a_{3N}
\end{pmatrix}
\cos(\omega t+\phi)
=\mathbf a\cos(\omega t+\phi).
\]

所以 \(\mathbf a\) 并不是突然出现的新假设；它正是一维振幅 \(A\) 的多维推广：

- \(|a_p|\) 表示第 \(p\) 个 Cartesian 坐标在该集体运动中的相对振幅；
- \(a_p\) 的正负表示该坐标相对于所选方向和其他坐标的运动方向；
- 各分量之间的比例决定分子振动时的“形状”，也就是哪些原子一起动、各自沿什么方向动；
- 整个 \(\mathbf a\) 同时乘一个常数只会改变总振幅，不会改变模式的形状或频率。

例如，若某个瞬间两个坐标满足 \(a_1>0\)、\(a_2<0\)，它们便沿各自坐标轴的相反方向运动；若 \(a_3=0\)，该坐标在这个模式中不运动。可以把 \(\mathbf a\) 想成一张固定的“集体动作图”，而 \(\cos(\omega t+\phi)\) 负责让整张动作图随时间在正向和反向之间来回缩放。

我们之所以先尝试这种形式，是因为线性、常系数的耦合振动方程允许把运动分解成若干这样的简正模式。这里并不是预先知道任意运动都只有一个频率；而是先寻找每一种能够保持固定动作比例、以单一频率振动的基本模式。一般运动随后可以由这些基本模式线性叠加得到。

现在将这个试探形式写为

\[
\Delta\mathbf{x}(t)=\mathbf a\cos(\omega t+\phi),
\]

由于 \(\mathbf a\) 不随时间变化，对时间求两次导数得到

\[
\Delta\ddot{\mathbf{x}}(t)
=-\omega^2\mathbf a\cos(\omega t+\phi)
=-\omega^2\Delta\mathbf{x}(t).
\]

代回运动方程得到

\[
\left(\mathbf H\mathbf a-\omega^2\mathbf M\mathbf a\right)
\cos(\omega t+\phi)=0.
\]

这个等式必须在整个振动过程中都成立，而余弦函数并不恒为零，所以它前面的向量必须为零：

\[
\mathbf H\mathbf a=\omega^2\mathbf M\mathbf a.
\]

这一步也解释了为什么要求特征值：并不是我们为了使用矩阵算法而凭空构造了特征值问题，而是“所有坐标保持固定比例 \(\mathbf a\)，并以单一频率共同振动”这一物理要求，自动产生了这个方程。只有某些特殊的方向 \(\mathbf a\) 能在运动中保持形状不变；这些方向就是简正模式，而相应的 \(\omega^2\) 就是广义特征值。

这叫作广义特征值问题。它与普通的

\[
\mathbf A\mathbf v=\lambda\mathbf v
\]

不同，因为右边多了质量矩阵 \(\mathbf M\)。

## 7. 为什么必须质量加权

### 7.1 直接对角化 Hessian 为什么不够

直接求

\[
\mathbf H\mathbf v=\lambda\mathbf v
\]

只分析了势能曲率，却没有正确包含不同原子的惯性。

例如 O-H 与 O-D 在近似相同的电子势能面上具有近似相同的力常数，但氘更重，所以 O-D 振动频率明显更低。只对角化原始 Hessian 无法给出这种正确的同位素效应。

### 7.2 定义质量加权坐标

为了避免混淆，先区分两类量：

- \(\Delta\mathbf x(t)\) 和 \(\mathbf q(t)\) 是随时间变化的坐标；
- \(\mathbf a\) 和 \(\mathbf l\) 是某一个简正模式中固定不变的模式向量，也就是相应坐标空间里的振幅方向。

定义质量加权坐标

\[
\mathbf q(t)=\mathbf M^{1/2}\Delta\mathbf{x}(t),
\qquad
\Delta\mathbf{x}(t)=\mathbf M^{-1/2}\mathbf q(t).
\]

质量矩阵及其平方根为

\[
\mathbf M=\operatorname{diag}
(m_1,m_1,m_1,m_2,m_2,m_2,\ldots),
\]

\[
\mathbf M^{1/2}=\operatorname{diag}
\left(
\sqrt{m_1},\sqrt{m_1},\sqrt{m_1},\ldots
\right),
\]

\[
\mathbf M^{-1/2}=\operatorname{diag}
\left(
\frac1{\sqrt{m_1}},\frac1{\sqrt{m_1}},\frac1{\sqrt{m_1}},\ldots
\right).
\]

因为所有原子质量都为正，所以 \(\mathbf M^{1/2}\) 和 \(\mathbf M^{-1/2}\) 都存在，并满足

\[
\mathbf M^{-1/2}\mathbf M^{1/2}
=\mathbf I.
\]

第 6 节在 Cartesian 坐标中把一个简正模式写成

\[
\Delta\mathbf x(t)
=\mathbf a\cos(\omega t+\phi).
\]

把它变换到质量加权坐标：

\[
\begin{aligned}
\mathbf q(t)
&=\mathbf M^{1/2}\Delta\mathbf x(t)\\
&=\mathbf M^{1/2}\mathbf a
  \cos(\omega t+\phi)\\
&=\mathbf l\cos(\omega t+\phi),
\end{aligned}
\]

这里定义

\[
\boxed{
\mathbf l\equiv\mathbf M^{1/2}\mathbf a
},
\qquad
\boxed{
\mathbf a=\mathbf M^{-1/2}\mathbf l
}.
\]

因此，\(\mathbf a\) 是 Cartesian 坐标空间中的模式方向，\(\mathbf l\) 是同一个物理模式在质量加权坐标空间中的方向。\(\mathbf l\) 不是突然出现的新模式，也不等于随时间变化的 \(\mathbf q(t)\)；它是 \(\mathbf q(t)\) 的固定振幅向量。

现在从第 6 节的广义特征值方程出发：

\[
\mathbf H\mathbf a
=\omega^2\mathbf M\mathbf a.
\]

代入 \(\mathbf a=\mathbf M^{-1/2}\mathbf l\)：

\[
\mathbf H\mathbf M^{-1/2}\mathbf l
=\omega^2\mathbf M\mathbf M^{-1/2}\mathbf l
=\omega^2\mathbf M^{1/2}\mathbf l.
\]

在等式两边从左侧乘以 \(\mathbf M^{-1/2}\)，得到

\[
\mathbf M^{-1/2}\mathbf H\mathbf M^{-1/2}\mathbf l
=\omega^2\mathbf l.
\]

定义质量加权 Hessian

\[
\mathbf F=\mathbf M^{-1/2}\mathbf H\mathbf M^{-1/2},
\]

于是问题变成普通特征值问题：

\[
\mathbf F\mathbf l=\lambda\mathbf l,
\qquad
\lambda=\omega^2.
\]

所以符号之间的完整关系是

\[
\boxed{
\Delta\mathbf x(t)=\mathbf a\cos(\omega t+\phi),
\quad
\mathbf q(t)=\mathbf l\cos(\omega t+\phi),
\quad
\mathbf l=\mathbf M^{1/2}\mathbf a
}.
\]

特征向量的整体倍数是任意的。若数值求解器把 \(\mathbf l\) 归一化为 \(\mathbf l^{\mathrm T}\mathbf l=1\)，对应的 Cartesian 模式向量便满足 \(\mathbf a^{\mathrm T}\mathbf M\mathbf a=1\)。这种归一化只规定模式向量的尺度，不决定分子在实际运动中的振幅。

这就是 Project 2 Step 3 和 Step 4 的理论来源。

### 7.3 元素形式

若坐标 \(p\) 属于原子 \(A(p)\)，坐标 \(q\) 属于原子 \(A(q)\)，则

\[
F_{pq}=\frac{H_{pq}}
{\sqrt{m_{A(p)}m_{A(q)}}}.
\]

程序中的原子映射为

```text
atom_of_coordinate_p = p / 3
atom_of_coordinate_q = q / 3
```

这里是整数除法。

### 7.4 质量加权坐标的另一个好处

Cartesian 动能为

\[
T=\frac12\Delta\dot{\mathbf{x}}^{\mathrm T}
\mathbf M\Delta\dot{\mathbf{x}}.
\]

代入 \(\mathbf q=\mathbf M^{1/2}\Delta\mathbf{x}\) 后：

\[
T=\frac12\dot{\mathbf q}^{\mathrm T}\dot{\mathbf q}.
\]

质量矩阵消失了。质量加权实质上把不同质量的坐标重新缩放，使动能在新坐标中具有统一的欧氏形式。

## 8. 对角化为什么得到简正模式

质量加权坐标中的运动方程为

\[
\ddot{\mathbf q}(t)+\mathbf F\mathbf q(t)=0,
\]

其中 \(\mathbf F=\mathbf M^{-1/2}\mathbf H\mathbf M^{-1/2}\) 是质量加权 Hessian。这个方程仍然是耦合的：\(\mathbf F\mathbf q\) 的第 \(p\) 个分量通常依赖 \(\mathbf q\) 的许多分量。

### 8.1 从特征向量基底得到简正坐标

由于 \(\mathbf F\) 是实对称矩阵，它有 \(3N\) 个彼此正交且可以归一化的特征向量：

\[
\mathbf F\mathbf l_k
=\lambda_k\mathbf l_k,
\qquad
\mathbf l_i^{\mathrm T}\mathbf l_j=\delta_{ij}.
\]

把这些特征向量作为列组成矩阵

\[
\mathbf L
=
\begin{pmatrix}
\mathbf l_1&\mathbf l_2&\cdots&\mathbf l_{3N}
\end{pmatrix}.
\]

于是

\[
\mathbf L^{\mathrm T}\mathbf L
=\mathbf L\mathbf L^{\mathrm T}
=\mathbf I,
\]

并且

\[
\mathbf F\mathbf L
=\mathbf L\boldsymbol\Lambda,
\qquad
\mathbf L^{\mathrm T}\mathbf F\mathbf L
=\boldsymbol\Lambda,
\]

其中

\[
\boldsymbol\Lambda
=\operatorname{diag}(\lambda_1,\ldots,\lambda_{3N}).
\]

因为 \(\{\mathbf l_k\}\) 构成质量加权坐标空间的一组完备正交基，任意质量加权位移 \(\mathbf q(t)\) 都可以在这组基底中展开：

\[
\boxed{
\mathbf q(t)
=\sum_{k=1}^{3N}\mathbf l_k Q_k(t)
=\mathbf L\mathbf Q(t)
}.
\]

这里

\[
\mathbf Q(t)
=
\begin{pmatrix}
Q_1(t)&Q_2(t)&\cdots&Q_{3N}(t)
\end{pmatrix}^{\mathrm T}
\]

称为**简正坐标向量**，而 \(Q_k(t)\) 是第 \(k\) 个简正坐标。它们不是新的原子 Cartesian 坐标，而是 \(\mathbf q(t)\) 沿各个特征向量方向的展开系数。由于 \(\mathbf L\) 是正交矩阵，

\[
\boxed{
\mathbf Q(t)=\mathbf L^{\mathrm T}\mathbf q(t)
},
\qquad
\boxed{
Q_k(t)=\mathbf l_k^{\mathrm T}\mathbf q(t)
}.
\]

所以 \(Q_k(t)\) 可以理解为：在时刻 \(t\)，整个分子位移中“第 \(k\) 个简正模式含有多少”。若只有 \(Q_k(t)\) 非零，分子就在做纯粹的第 \(k\) 个简正振动；若多个 \(Q_k(t)\) 非零，实际运动就是多个简正模式的叠加。

第 7.2 节写出的 \(\mathbf q(t)=\mathbf l\cos(\omega t+\phi)\)，正是这里只保留一个模式时的特殊情况。更明确地写，对第 \(k\) 个模式有

\[
\mathbf q(t)=\mathbf l_k Q_k(t),
\qquad
Q_k(t)=A_k\cos(\omega_k t+\phi_k),
\]

其中 \(A_k\) 和 \(\phi_k\) 由初始条件决定。

### 8.2 为什么换成 \(\mathbf Q\) 后运动会解耦

把 \(\mathbf q(t)=\mathbf L\mathbf Q(t)\) 代入质量加权运动方程。因为 \(\mathbf L\) 不随时间变化，

\[
\mathbf L\ddot{\mathbf Q}
+\mathbf F\mathbf L\mathbf Q=0.
\]

利用 \(\mathbf F\mathbf L=\mathbf L\boldsymbol\Lambda\)，得到

\[
\mathbf L\ddot{\mathbf Q}
+\mathbf L\boldsymbol\Lambda\mathbf Q=0.
\]

在等式两边从左侧乘以 \(\mathbf L^{\mathrm T}\)，并使用 \(\mathbf L^{\mathrm T}\mathbf L=\mathbf I\)：

\[
\boxed{
\ddot{\mathbf Q}
+\boldsymbol\Lambda\mathbf Q=0
}.
\]

逐分量写就是

\[
\ddot Q_k(t)+\lambda_kQ_k(t)=0,
\qquad k=1,\ldots,3N.
\]

这是一组互不耦合的一维振子方程。原来 \(\mathbf q\) 的各分量会通过非对角矩阵 \(\mathbf F\) 相互影响；换到特征向量基底后，矩阵变成对角矩阵 \(\boldsymbol\Lambda\)，每个 \(Q_k\) 只出现在自己的方程中。

因此每个模式满足

\[
\omega_k=\sqrt{\lambda_k}.
\]

Cartesian 位移可以由简正坐标恢复：

\[
\boxed{
\Delta\mathbf x(t)
=\mathbf M^{-1/2}\mathbf L\mathbf Q(t)
=\sum_{k=1}^{3N}
\mathbf M^{-1/2}\mathbf l_k Q_k(t)
}.
\]

这条式子完整地区分了三种对象：

- \(\mathbf M^{-1/2}\mathbf l_k\)：第 \(k\) 个模式的 Cartesian 运动形状；
- \(Q_k(t)\)：该模式随时间变化的幅度；
- \(\Delta\mathbf x(t)\)：所有模式叠加后，原子在时刻 \(t\) 的实际 Cartesian 位移。

### 8.3 特征值代表什么

特征值表示相应简正方向上的质量加权势能曲率：

- \(\lambda_k>0\)：稳定的谐振动；
- \(\lambda_k=0\)：势能不随该运动改变，理想情况下对应整体平移或转动；
- \(\lambda_k<0\)：该方向曲率向下，参考构型不是极小值。

### 8.4 特征向量代表什么

特征向量 \(\mathbf l_k\) 给出质量加权坐标中的运动方向。要恢复 Cartesian 位移方向：

\[
\Delta\mathbf{x}_k
\propto\mathbf M^{-1/2}\mathbf l_k.
\]

所以 Eigen 返回的特征向量不能在不去质量加权的情况下直接当作原子位移。若将来要制作振动动画，这一步很重要。

特征向量整体乘以 \(-1\) 仍表示同一个模式，因此两个程序给出的模式箭头可能方向完全相反，但物理上没有区别。

## 9. 为什么有 3N-6 个振动模式

分子共有 \(3N\) 个 Cartesian 自由度，但并非所有自由度都是内部振动。

### 9.1 整体平移

整个分子一起沿 \(x,y,z\) 移动，不改变任何内部距离和势能，因此有 3 个平移自由度。

### 9.2 整体转动

对于非线性分子，绕三条独立轴整体旋转也不改变内部结构，因此有 3 个转动自由度。

所以非线性分子的振动自由度为

\[
3N-3-3=3N-6.
\]

例如水分子 \(N=3\)：

\[
3N-6=3,
\]

对应对称伸缩、弯曲和反对称伸缩三个内部振动。

### 9.3 线性分子为什么是 3N-5

线性分子绕分子轴旋转时，所有位于轴上的点都不发生位移；在经典点核模型中这不是独立的有限转动自由度。线性分子只有 2 个有效整体转动自由度，因此

\[
3N-3-2=3N-5.
\]

### 9.4 Hessian 中的零特征值

如果满足以下理想条件：

- 几何是真正的驻点；
- 势能严格满足平移和转动不变性；
- Hessian 数值完全精确；

那么非线性分子应有 6 个零特征值，线性分子应有 5 个。

但 Project 2 的 H₂O 参考说明指出，使用的结构不是该计算势能面上的严格驻点，因此只有三个平移频率接近零，三个“转动方向”出现明显非零值。这不是 \(3N-6\) 理论失效，而是参考几何和 Hessian 没有处在理想驻点条件下。

实际量化化学程序通常还会显式投影掉整体平移与转动，再报告内部振动频率；Project 2 的基础算法没有做这一步。

## 10. 负特征值与虚频

若

\[
\lambda_k<0,
\]

则 \(\sqrt{\lambda_k}\) 不是实数。量化化学输出通常把它报告为“虚频”，常见写法是负的 cm\(^{-1}\) 数值。

它的几何含义是：沿这个模式稍微移动，势能会降低，因而参考点不是局部极小值。

- 0 个显著虚频：通常是局部极小值；
- 1 个显著虚频：通常是一阶鞍点，即过渡态候选；
- 多个显著虚频：更高阶鞍点或几何尚未充分优化。

非常小的负特征值也可能只是数值误差，不能不加阈值地解释为真实虚频。

Project 2 的参考数据主要展示非负特征值，但健壮实现仍应明确决定如何处理负数，而不是直接对负数调用 `sqrt()` 得到 `NaN`。

## 11. 从特征值换算到 cm⁻¹

### 11.1 质量加权 Hessian 的单位

作业中：

\[
[H]=\frac{E_h}{a_0^2},
\qquad [m]=\mathrm{amu}.
\]

所以质量加权 Hessian 特征值的单位为

\[
[\lambda]
=\frac{E_h}{a_0^2\,\mathrm{amu}}.
\]

这在量纲上等于 \(\mathrm{s^{-2}}\)，因为

\[
\frac{\text{energy}}{\text{mass}\times\text{length}^2}
=
\frac{\mathrm{kg\,m^2\,s^{-2}}}
{\mathrm{kg\,m^2}}
=\mathrm{s^{-2}}.
\]

因此平方根对应角频率。

### 11.2 完整换算公式

若程序中的数值特征值为 \(\lambda_k\)，则 SI 角频率是

\[
\omega_k=
\sqrt{
\lambda_k
\frac{E_h\;[\mathrm J]}
{a_0^2\;[\mathrm{m^2}]\,m_u\;[\mathrm{kg}]}
}.
\]

转换为光谱波数：

\[
\tilde\nu_k
=\frac{\omega_k}{2\pi c},
\]

其中若 \(c\) 使用 cm s\(^{-1}\)，结果就是 cm\(^{-1}\)。合并后：

\[
\boxed{
\tilde\nu_k
=
\frac{1}{2\pi c}
\sqrt{
\lambda_k
\frac{E_h}{a_0^2m_u}
}
}
\]

使用现代物理常数，常写成方便的数值形式：

\[
\boxed{
\tilde\nu_k\;[\mathrm{cm^{-1}}]
\approx5140.49\sqrt{\lambda_k}
}
\]

这里的前提必须同时满足：

1. 原始 Hessian 数值的单位是 \(E_h/a_0^2\)；
2. 质量加权时使用的数值单位是 amu；
3. \(\lambda_k\) 是该质量加权 Hessian 的数值特征值。

若任何输入单位改变，\(5140.49\) 就不能直接使用。

### 11.3 用 H₂O 的一个特征值检查

参考答案最大的特征值约为

\[
\lambda=0.2351542439.
\]

于是

\[
\tilde\nu
\approx5140.49\sqrt{0.2351542439}
\approx2492.76\ \mathrm{cm^{-1}},
\]

与作业参考输出一致。

### 11.4 为什么一定有 \(2\pi\)

对运动方程求平方根得到的是角频率 \(\omega\)，而不是普通频率 \(\nu\)：

\[
\omega=2\pi\nu.
\]

光谱波数则是

\[
\tilde\nu=\frac{\nu}{c}
=\frac{\omega}{2\pi c}.
\]

忘记 \(2\pi\) 会让最终结果大约大 6.283 倍。

## 12. 量子力学中的谐振子

前面的 Hessian 和简正模式来自经典小振动理论。量子力学并不改变简正频率，而是把每个独立模式的能量量子化。

对于第 \(k\) 个简正模式：

\[
E_{v_k}=\hbar\omega_k\left(v_k+\frac12\right),
\qquad v_k=0,1,2,\ldots
\]

整个谐振动能量为

\[
E_{\mathrm{vib}}
=\sum_k\hbar\omega_k
\left(v_k+\frac12\right).
\]

即使所有模式都处于 \(v_k=0\)，仍有零点振动能：

\[
E_{\mathrm{ZPVE}}=\frac12\sum_k\hbar\omega_k.
\]

因此原子核在绝对零度也不会静止在势能最低点。Project 2 只计算频率，没有要求累加零点能，但频率正是进一步计算零点能、振动热容和振动熵的基础。

## 13. Project 2 五个步骤与理论的对应

### Step 1：读取几何

几何文件提供：

- 原子数 \(N\)；
- 每个原子的原子序数；
- Cartesian 坐标。

在本项目中，几何的首要用途是确定原子种类和质量。基础频率计算使用现成 Hessian，因此不会再次从坐标数值计算二阶导数。

### Step 2：读取 Cartesian Hessian

创建

\[
(3N)\times(3N)
\]

的动态矩阵。H₂O 有三个原子，所以矩阵是 \(9\times9\)。

作业文件每行按三个数分组，但内存中仍是普通方阵。矩形的文本排版只是文件格式，不代表数学上的 Hessian 是矩形矩阵。

还应检查 Hessian 文件首行的原子数是否与几何文件一致。

### Step 3：质量加权

对每个元素执行

\[
F_{pq}=\frac{H_{pq}}
{\sqrt{m_{p/3}m_{q/3}}}.
\]

结果仍是 \(3N\times3N\) 实对称矩阵。

### Step 4：对角化

求解

\[
\mathbf F\mathbf l_k=\lambda_k\mathbf l_k.
\]

特征值给出 \(\omega_k^2\)，特征向量给出简正模式。Eigen 的自伴随求解器通常把特征值按从小到大排列；具体代码不应在没有确认文档时假定其他顺序。

### Step 5：计算频率

对正特征值：

\[
\tilde\nu_k\approx5140.49\sqrt{\lambda_k}.
\]

对接近零的特征值，应使用容差将数值噪声视为零；对显著负特征值，则应按虚频处理。

## 14. 数值计算中容易忽略的问题

### 14.1 浮点数中的“零”

理论上的零特征值可能计算成

```text
 2.0e-12
-7.0e-11
```

因此判断零时应使用容差，而不是只写 `value == 0.0`。

### 14.2 Hessian 的轻微不对称

理论上 \(\mathbf H\) 对称，但文本舍入或上游程序可能造成

\[
H_{pq}\ne H_{qp}
\]

的极小误差。可以先检查最大不对称量。必要时使用

\[
\mathbf H_{\mathrm{sym}}
=\frac12(\mathbf H+\mathbf H^{\mathrm T})
\]

进行对称化，但不应在存在巨大差异时悄悄掩盖输入错误。

### 14.3 简并模式

如果两个特征值相同或非常接近，相应特征子空间中的任意正交线性组合都是合法特征向量。因此不同软件可能给出看起来不同的两个模式，但它们张成同一个简并子空间，频率仍相同。

### 14.4 频率排序

参考输出可能按降序打印，而 Eigen 可能按升序返回。排序不同不代表数值错误。应先比较特征值集合，再比较打印顺序。

### 14.5 原子质量的选择

频率依赖质量。使用最丰富同位素质量、平均原子量或指定同位素质量会得到略有不同的结果。Project 2 要求沿用 Project 1 的 amu 质量表，因此参考答案比较时应保持相同选择。

## 15. Hessian 与惯量张量虽然都要对角化，但不是同一个问题

Project 1 对角化的是 \(3\times3\) 惯量张量：

\[
\mathbf I\mathbf u=I\mathbf u,
\]

得到主转动轴和主惯量。

Project 2 对角化的是 \(3N\times3N\) 质量加权 Hessian：

\[
\mathbf F\mathbf l=\omega^2\mathbf l,
\]

得到简正振动方向和频率平方。

共同点只是“对称矩阵的特征值问题”；它们的物理含义、维度和单位完全不同。

| 项目 | 矩阵 | 大小 | 特征值含义 |
|---|---|---:|---|
| Project 1 | 惯量张量 \(\mathbf I\) | \(3\times3\) | 主转动惯量 |
| Project 2 | 质量加权 Hessian \(\mathbf F\) | \(3N\times3N\) | \(\omega^2\) |

## 16. 一套适合实现时使用的检查清单

1. 几何文件和 Hessian 文件的原子数一致。
2. Hessian 读入了恰好 \((3N)^2\) 个元素。
3. Hessian 接近对称。
4. 每个 Cartesian 坐标使用了所属原子的质量。
5. 质量加权时分母是 \(\sqrt{m_i m_j}\)，不是 \(m_i m_j\)。
6. 使用实对称矩阵特征值求解器。
7. 检查求解器是否成功。
8. 平方根之前区分正值、数值零和显著负值。
9. 从角频率转波数时包含 \(2\pi c\)。
10. 比较参考结果时确认排序、质量表和物理常数一致。

## 17. 最核心的公式

如果只保留整篇文档的主线，就是下面四步。

势能的谐近似：

\[
V\approx V_0+\frac12
\Delta\mathbf{x}^{\mathrm T}\mathbf H\Delta\mathbf{x}.
\]

耦合运动方程：

\[
\mathbf M\Delta\ddot{\mathbf{x}}
+\mathbf H\Delta\mathbf{x}=0.
\]

质量加权特征值问题：

\[
\left(\mathbf M^{-1/2}\mathbf H\mathbf M^{-1/2}\right)
\mathbf l_k
=\lambda_k\mathbf l_k,
\qquad \lambda_k=\omega_k^2.
\]

转换为光谱波数：

\[
\tilde\nu_k
=\frac{\omega_k}{2\pi c}
\approx5140.49\sqrt{\lambda_k}
\quad(\mathrm{cm^{-1}}),
\]

其中最后一个数值公式只适用于 Hessian 使用 \(E_h/a_0^2\)、质量使用 amu 的本项目单位约定。

## 参考资料

- Project 2 作业说明：`project2-instructions.pdf`
- E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations*, McGraw-Hill, 1955.
- Project 1 的单位说明：`../Project_01/UNIT_CONVERSIONS.md`
