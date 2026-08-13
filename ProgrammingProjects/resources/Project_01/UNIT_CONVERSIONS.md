# Project 1 单位体系与换算

## 1. 本项目的单位流

Project 1 中各阶段使用的单位并不完全相同：

| 物理量 | 当前程序中的单位 |
|---|---|
| Cartesian 坐标 | bohr，记作 \(a_0\) |
| 原子间距离 | bohr |
| 键角、面外角、二面角 | 内部先用弧度，打印时用 degree |
| 原子质量 | unified atomic mass unit，\(u\) 或 amu |
| 质心坐标 | bohr |
| 惯量张量 | \(\mathrm{amu\,bohr^2}\) |
| 主转动惯量 | \(\mathrm{amu\,bohr^2}\)，随后转换为其他单位 |
| 转动常数 | MHz 或 \(\mathrm{cm^{-1}}\) |

最重要的一条链是：

```text
坐标 (bohr)
    + 原子质量 (amu)
          ↓
惯量 (amu·bohr²)
          ↓
惯量 (amu·Å²) 或 (g·cm²)
          ↓
转动常数 (Hz → MHz) 或 (cm⁻¹)
```

## 2. 输入坐标为什么是 bohr

输入文件中的坐标，例如：

```text
6  0.000000000000  0.000000000000  1.889725988579
```

最后三个数的单位是 bohr。bohr 是 Hartree 原子单位制中的长度单位，即玻尔半径 \(a_0\)。

常用换算为：

\[
1\ \mathrm{bohr}=0.5291772105\ \mathrm{\mathring{A}}
\]

以及：

\[
1\ \mathrm{bohr}=0.5291772105\times10^{-8}\ \mathrm{cm}
\]

所以程序直接用输入坐标计算距离时，结果自然仍是 bohr。两个 bohr 坐标相减得到 bohr，求向量长度后也还是 bohr。

## 3. 角度为什么没有长度单位

键角公式中：

\[
\cos\theta=
\frac{\mathbf u\cdot\mathbf v}{|\mathbf u||\mathbf v|}
\]

分子的单位为 \(\mathrm{bohr^2}\)，分母也是 \(\mathrm{bohr^2}\)，因此单位约掉，\(\cos\theta\) 是无量纲数。

`np.arccos()`、`np.arcsin()` 返回弧度。弧度本身也是无量纲的；为了便于阅读，程序再用 `np.degrees()` 转为角度制。

## 4. amu 到底是什么

amu 在这里指 unified atomic mass unit，也写作 \(u\)。它定义为一个碳-12 原子静止质量的十二分之一：

\[
1\ \mathrm{u}\approx1.66053906892\times10^{-27}\ \mathrm{kg}
\]

或：

\[
1\ \mathrm{u}\approx1.66053906892\times10^{-24}\ \mathrm{g}
\]

本项目质量表使用最丰富同位素的质量，例如：

```python
H = 1.00782503223
C = 12.0
O = 15.9949146196
```

这些数值的单位是 amu。

“相对原子质量”严格来说是相对于 \(1\ \mathrm{u}\) 的无量纲比值，但其数值与用 amu 表示的单原子质量相同。因此在程序中通常直接把表中的数值当作 amu 使用。

## 5. 质心为什么仍以 bohr 表示

质心为：

\[
\mathbf R_{\mathrm{COM}}
=
\frac{\sum_i m_i\mathbf r_i}{\sum_i m_i}
\]

单位分析为：

\[
\frac{\mathrm{amu}\times\mathrm{bohr}}{\mathrm{amu}}
=\mathrm{bohr}
\]

所以 `RotationalAnalyzer.center_of_mass()` 返回的坐标仍是 bohr。

计算惯量前使用：

```python
centered_positions = molecule.positions - analyzer.center_of_mass()
```

平移只改变坐标原点，不改变坐标单位。

## 6. 惯量为什么是 amu·bohr²

点质量绕某转轴的惯量形式为：

\[
I=mr_\perp^2
\]

本项目中：

```text
m 的单位：amu
r 的单位：bohr
```

因此：

\[
[I]=\mathrm{amu\,bohr^2}
\]

例如：

\[
I_{xx}=\sum_i m_i(y_i^2+z_i^2)
\]

其中 \(y_i^2+z_i^2\) 的单位为 \(\mathrm{bohr^2}\)，乘质量后就是 \(\mathrm{amu\,bohr^2}\)。惯量张量对角化只是在旋转坐标系，不改变单位，因此三个特征值仍是 \(\mathrm{amu\,bohr^2}\)。

## 7. 惯量单位的两种转换

### 7.1 转换为 amu·Å²

质量单位不变，只转换长度：

\[
I(\mathrm{amu\,\mathring{A}^2})
=I(\mathrm{amu\,bohr^2})
\left(0.5291772105\right)^2
\]

代码为：

```python
moments_amu_angstrom2 = (
    moments_amu_bohr2 * BOHR_TO_ANGSTROM ** 2
)
```

注意换算因子必须平方，因为惯量包含长度的平方。

### 7.2 转换为 g·cm²

质量和长度都要转换：

\[
I(\mathrm{g\,cm^2})
=I(\mathrm{amu\,bohr^2})
\times
\frac{\mathrm{g}}{\mathrm{amu}}
\times
\left(\frac{\mathrm{cm}}{\mathrm{bohr}}\right)^2
\]

代码为：

```python
moments_g_cm2 = (
    moments_amu_bohr2
    * AMU_TO_GRAM
    * BOHR_TO_CM ** 2
)
```

## 8. 转动常数是什么，以及公式从哪里来（详细版）

### 8.1 从经典转动动能开始

刚体绕某条主惯性轴旋转时：

\[
E_{\mathrm{rot}}=\frac{1}{2}I\omega^2
\]

角动量为 \(L=I\omega\)，所以：

\[
E_{\mathrm{rot}}=\frac{L^2}{2I}
\]

这说明在角动量相同时，惯量越小，转动能量越大。

### 8.2 量子化之后

量子力学中：

\[
L^2=\hbar^2J(J+1),\qquad J=0,1,2,\ldots
\]

因此线性刚性转子的能级为：

\[
E_J=\frac{\hbar^2}{2I}J(J+1)
\]

把与 \(J\) 无关的系数定义为能量形式的转动常数：

\[
B_E=\frac{\hbar^2}{2I}
\]

于是 \(E_J=B_EJ(J+1)\)。此处 \(B_E\) 的单位是能量，例如 J 或 erg。

### 8.3 为什么公式变成 h/(8π²I)

光谱学常用 \(E=h\nu\) 把能量表示成频率，因此：

\[
B_{\mathrm{Hz}}=\frac{B_E}{h}
=\frac{\hbar^2}{2Ih}
\]

利用：

\[
\hbar=\frac{h}{2\pi}
\]

得到：

\[
B_{\mathrm{Hz}}
=\frac{(h/2\pi)^2}{2Ih}
=\frac{h}{8\pi^2I}
\]

所以作业公式最终使用 \(h\)，但它来自含 \(\hbar^2\) 的量子转动能级。不能直接把公式中的 \(h\) 换成 \(\hbar\)。

### 8.4 转动常数的物理意义

因为：

\[
B\propto\frac{1}{I}
\]

所以：

```text
质量靠近转轴 → I 小 → B 大 → 转动能级间隔较大
质量远离转轴 → I 大 → B 小 → 转动能级间隔较小
```

线性刚性转子的频率形式能级为：

\[
\frac{E_J}{h}=B_{\mathrm{Hz}}J(J+1)
\]

对于允许的 \(\Delta J=+1\) 跃迁：

\[
\nu_{J\rightarrow J+1}=2B_{\mathrm{Hz}}(J+1)
\]

因此，\(B\) 通常不是某一条谱线的频率，而是决定整套转动谱线位置和间隔的参数。最低的 \(J=0\rightarrow1\) 跃迁位于 \(2B\)，不是 \(B\)。

### 8.5 三个主惯量对应 A、B、C

一般非线性分子有：

\[
I_A\le I_B\le I_C
\]

分别定义：

\[
A=\frac{h}{8\pi^2I_A},\qquad
B=\frac{h}{8\pi^2I_B},\qquad
C=\frac{h}{8\pi^2I_C}
\]

因此：

\[
A\ge B\ge C
\]

一般刚性转子的 Hamiltonian 为：

\[
\hat H_{\mathrm{rot}}
=\frac{\hat J_A^2}{2I_A}
+\frac{\hat J_B^2}{2I_B}
+\frac{\hat J_C^2}{2I_C}
\]

| 转子类型 | 惯量关系 | 转动常数关系 |
|---|---|---|
| spherical top | \(I_A=I_B=I_C\) | \(A=B=C\) |
| prolate symmetric top | \(I_A<I_B=I_C\) | \(A>B=C\) |
| oblate symmetric top | \(I_A=I_B<I_C\) | \(A=B>C\) |
| asymmetric top | 三者不同 | 三者不同 |

### 8.6 单位为什么能得到 Hz

若使用 CGS 单位：

```text
h：erg·s
I：g·cm²
```

由于：

\[
1\ \mathrm{erg}=1\ \mathrm{g\,cm^2\,s^{-2}}
\]

所以：

\[
\frac{\mathrm{erg\,s}}{\mathrm{g\,cm^2}}
=\mathrm{s^{-1}}=\mathrm{Hz}
\]

代码为：

```python
constants_hz = PLANCK_CONSTANT_ERG_S / (
    8 * np.pi ** 2 * moments_g_cm2
)

constants_mhz = constants_hz / 1e6
```

也可以完全使用 SI：\(h\) 用 J·s，\(I\) 用 kg·m²，结果同样是 Hz。关键是所有量必须使用同一套单位。

### 8.7 Hz 与 cm⁻¹ 表示的是同一个量

频率与波数满足：

\[
\nu=c\tilde\nu
\]

所以：

\[
\tilde B=\frac{B_{\mathrm{Hz}}}{c}
\]

如果 \(c\) 使用 cm/s，则结果自然是 cm⁻¹：

```python
constants_cm_inverse = (
    constants_hz / SPEED_OF_LIGHT_CM_S
)
```

MHz 和 cm⁻¹ 不是两种不同的转动常数，只是同一能量尺度的两种表示方法。

### 8.8 allene 的 A 常数完整演算

allene 的最小主转动惯量约为：

\[
I_A=10.7970237\ \mathrm{amu\,bohr^2}
\]

先转换为 g·cm²：

\[
I_A=10.7970237
\times1.66053906892\times10^{-24}
\times(0.5291772105\times10^{-8})^2
\]

得到：

\[
I_A\approx5.02059764\times10^{-40}\ \mathrm{g\,cm^2}
\]

代入转动常数公式：

\[
A=\frac{6.62607015\times10^{-27}}
{8\pi^2(5.02059764\times10^{-40})}
\]

得到：

\[
A\approx1.67151731\times10^{11}\ \mathrm{Hz}
\]

即：

\[
A\approx167151.731\ \mathrm{MHz}
\]

再转换成波数：

\[
\tilde A=\frac{1.67151731\times10^{11}}
{2.99792458\times10^{10}}
\approx5.57558\ \mathrm{cm^{-1}}
\]

课程答案尾数略有差异，是因为课程使用较旧且经过舍入的物理常数。

### 8.9 一句话总结

```text
几何和质量 → 转动惯量 I
角动量量子化 → 转动能级
B = h/(8π²I) → 能级的频率尺度
除以 10⁶ → MHz
除以光速 → cm⁻¹
```

## 参考资料

- [NIST：2022 CODATA fundamental constants](https://physics.nist.gov/cuu/Constants/)
- [NIST：Introduction to fundamental constants](https://physics.nist.gov/cuu/Constants/introduction.html)

## 9. 为什么除以光速得到 cm⁻¹

频率和波数满足：

\[
\nu=c\tilde\nu
\]

所以：

\[
\tilde\nu=\frac{\nu}{c}
\]

如果频率用 \(\mathrm{s^{-1}}\)，光速使用 \(\mathrm{cm\,s^{-1}}\)：

\[
\frac{\mathrm{s^{-1}}}{\mathrm{cm\,s^{-1}}}
=\mathrm{cm^{-1}}
\]

代码为：

```python
constants_cm_inverse = (
    constants_hz / SPEED_OF_LIGHT_CM
)
```

## 10. 为什么 A ≥ B ≥ C

三个主转动惯量通常按升序排列：

\[
I_A\le I_B\le I_C
\]

而转动常数与惯量成反比：

\[
A=\frac{h}{8\pi^2I_A},\qquad
B=\frac{h}{8\pi^2I_B},\qquad
C=\frac{h}{8\pi^2I_C}
\]

所以：

\[
A\ge B\ge C
\]

小惯量对应大转动常数，大惯量对应小转动常数。

## 11. 这里与 QM 原子单位制有什么不同

Hartree 原子单位制通常令：

\[
\hbar=m_e=e=4\pi\epsilon_0=1
\]

其基本单位包括：

| 物理量 | 原子单位 |
|---|---|
| 长度 | bohr，\(a_0\) |
| 质量 | 电子质量，\(m_e\) |
| 电荷 | 基本电荷，\(e\) |
| 能量 | hartree，\(E_h\) |
| 作用量 | \(\hbar\) |

关键区别是：

```text
原子单位制的质量单位：电子质量 me
本项目惯量的质量单位：amu
```

它们不是同一个单位：

\[
1\ \mathrm{u}\approx1822.888486\,m_e
\]

因此，若真正要把本项目惯量写成纯原子单位，应把核质量从 amu 转换为电子质量：

\[
I(\mathrm{a.u.})
=I(\mathrm{amu\,bohr^2})
\times1822.888486
\]

所得惯量单位可以理解为 \(m_ea_0^2\)。但本作业没有这样做，而是混合使用：

```text
长度：原子单位 bohr
核质量：化学中常用的 amu
惯量：amu·bohr²
```

这种混合单位在分子光谱与量化化学后处理中很常见，只要在代入转动常数公式前完整转换到一致的 SI 或 CGS 单位即可。

## 12. 实际 QM 程序常见的单位

不同程序和输入格式可能不同，但常见情况是：

| 输出量 | 常见单位 |
|---|---|
| 分子轨道与总能量 | hartree |
| 几何坐标 | bohr 或 Å |
| 梯度 | hartree/bohr |
| Hessian | hartree/bohr² |
| 偶极矩 | 原子单位或 debye |
| 原子质量 | amu，或在严格原子单位计算中转换为 \(m_e\) |
| 振动频率 | cm⁻¹ |
| 转动常数 | MHz、GHz 或 cm⁻¹ |

因此看到“程序使用原子单位”时，不能自动认为每一项输出都没有单位。需要检查具体物理量和程序文档，尤其要确认：

1. 坐标是 bohr 还是 Å；
2. 核质量是 amu 还是电子质量；
3. 能量是 hartree、eV 还是其他单位；
4. 频率是 Hz、MHz、GHz 还是 cm⁻¹。

## 13. 本项目建议使用的常数

```python
BOHR_TO_ANGSTROM = 0.5291772105
BOHR_TO_CM = BOHR_TO_ANGSTROM * 1e-8

AMU_TO_GRAM = 1.66053906892e-24

PLANCK_CONSTANT_ERG_S = 6.62607015e-27
SPEED_OF_LIGHT_CM_S = 2.99792458e10
```

其中 Planck 常数和真空光速在现行 SI 中是精确定义值；bohr 和原子质量常数取自 CODATA 推荐值。课程参考答案使用了较旧且经过舍入的常数，因此最后几位可能与使用现行 CODATA 常数的结果略有不同，这不代表公式错误。

## 14. 对应的 Python 实现

```python
def get_principal_moments_in_different_units(self):
    moments_amu_bohr2 = self.get_principal_moments_of_inertia()

    moments_amu_angstrom2 = (
        moments_amu_bohr2 * BOHR_TO_ANGSTROM ** 2
    )

    moments_g_cm2 = (
        moments_amu_bohr2
        * AMU_TO_GRAM
        * BOHR_TO_CM ** 2
    )

    return (
        moments_amu_bohr2,
        moments_amu_angstrom2,
        moments_g_cm2,
    )


def get_rotational_constants(self):
    _, _, moments_g_cm2 = (
        self.get_principal_moments_in_different_units()
    )

    constants_hz = PLANCK_CONSTANT_ERG_S / (
        8 * np.pi ** 2 * moments_g_cm2
    )

    constants_mhz = constants_hz / 1e6
    constants_cm_inverse = constants_hz / SPEED_OF_LIGHT_CM_S

    return constants_mhz, constants_cm_inverse
```

## 15. 常见错误清单

1. 把 bohr 当成 Å，造成长度差约 1.8897 倍。
2. 惯量中的长度换算忘记平方。
3. 把 amu 直接当成 g 或 kg。
4. 把原子单位质量 \(m_e\) 与 amu 混为一谈。
5. 在 \(h/(8\pi^2I)\) 中误用 \(\hbar\)。
6. 计算 \(\mathrm{cm^{-1}}\) 时使用 \(c\) 的 m/s 数值，忘记额外换算 100。
7. 用未平移到质心的坐标计算惯量。
8. 课程答案与现行常数略有尾数差异时，误以为算法错误。

## 参考资料

- [NIST：2022 CODATA fundamental constants](https://physics.nist.gov/cuu/Constants/)
- [NIST：Atomic Weights and Isotopic Compositions](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii2&ele=&isotype=some)
- [IUPAC：Periodic Table of Elements](https://iupac.org/what-we-do/periodic-table-of-elements/)
