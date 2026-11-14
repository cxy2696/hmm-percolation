# 数学笔记

这份笔记只处理代码里实际用到的结论。我尽量把“证明了的”“引用文献的”“只是数值观察的”分开写，避免把有限尺寸实验说成定理。

约定：状态数 \(N\)，观测字母表大小 \(M\)，序列长度 \(T\)。\(\Delta_N=\{x\in[0,1]^N:\sum_i x_i=1\}\)。概率矩阵按行随机。取 \(0\log 0=0\)，\(\log 0=-\infty\)。所有随机变量定义在同一个概率空间上。

---

## 1. 隐马尔可夫模型

### 1.1 生成假设与联合分布

**假设 1.1.** 随机变量 \(Q_t\in\{1,\dots,N\}\)、\(O_t\in\{1,\dots,M\}\) 满足：

1. \(P(Q_1=i)=\pi_i\)；
2. 对 \(t<T\)，\(P(Q_{t+1}=j\mid Q_t=i,\ Q_{1:t-1},O_{1:t})=a_{ij}\)；
3. 给定 \(Q_{1:T}\)，\(O_1,\dots,O_T\) 条件独立，且
   \[
   P(O_t=k\mid Q_{1:T})=b_{Q_t k}.
   \]

条件 2 比“只依赖 \(Q_t\)”稍强：它同时排除了过去观测对下一状态的影响。条件 3 是通常的发射独立性。

**引理 1.2（联合分布因子化）.** 对任意 \(q\in\{1,\dots,N\}^T\)、\(o\in\{1,\dots,M\}^T\)，
\[
P(Q=q,O=o)
=\pi_{q_1}b_{q_1}(o_1)\prod_{t=2}^{T}a_{q_{t-1}q_t}b_{q_t}(o_t).
\]

**证明.** 用链式法则按 \(Q_1,O_1,Q_2,O_2,\dots\) 展开：
\[
P(Q=q,O=o)
=P(Q_1=q_1)\prod_{t=1}^{T}P(O_t=o_t\mid Q_{1:t}=q_{1:t},O_{1:t-1})
\prod_{t=2}^{T}P(Q_t=q_t\mid Q_{1:t-1}=q_{1:t-1},O_{1:t-1}).
\]
由条件 2，第二个乘积的每一项是 \(a_{q_{t-1}q_t}\)；由条件 3，第一个乘积在给定 \(Q_{1:T}\) 时是 \(\prod_t b_{q_t}(o_t)\)。合并即得。\(\square\)

### 1.2 前向递推

定义
\[
\alpha_t(j)=P(O_1=o_1,\dots,O_t=o_t,\ Q_t=j).
\]

**定理 1.3.** \(\alpha_1(j)=\pi_j b_j(o_1)\)，且对 \(t<T\)，
\[
\alpha_{t+1}(j)=b_j(o_{t+1})\sum_{i=1}^{N}\alpha_t(i)a_{ij}.
\]
从而
\[
P(O=o)=\sum_{j=1}^{N}\alpha_T(j).
\]

**证明.** 对 \(t+1\) 用全概率公式，按 \(Q_t=i\) 分类：
\[
\begin{aligned}
\alpha_{t+1}(j)
&=P(O_{1:t+1},Q_{t+1}=j)\\
&=\sum_i P(O_{1:t},Q_t=i)\,
P(Q_{t+1}=j\mid Q_t=i,O_{1:t})\,
P(O_{t+1}=o_{t+1}\mid Q_t=i,Q_{t+1}=j,O_{1:t}).
\end{aligned}
\]
第一项是 \(\alpha_t(i)\)，第二项由假设 1.1(2) 等于 \(a_{ij}\)，第三项由假设 1.1(3) 等于 \(b_j(o_{t+1})\)。代回即得递推。最后对 \(Q_T\) 求和得到边缘分布。\(\square\)

### 1.3 后向递推、平滑与转移期望

定义
\[
\beta_t(i)=P(O_{t+1:T}=o_{t+1:T}\mid Q_t=i),\qquad \beta_T(i)=1.
\]

**定理 1.4.**
\[
\beta_t(i)=\sum_{j=1}^{N}a_{ij}b_j(o_{t+1})\beta_{t+1}(j),\qquad t<T.
\]

**证明.** 对 \(Q_{t+1}=j\) 用全概率公式，并利用给定 \(Q_{t+1}\) 后 \(O_{t+2:T}\) 与 \(Q_t,O_{t+1}\) 条件独立（这由假设 1.1 的 Markov 性与发射独立性推出）：
\[
\begin{aligned}
\beta_t(i)
&=\sum_j P(Q_{t+1}=j\mid Q_t=i)\,
P(O_{t+1}=o_{t+1}\mid Q_{t+1}=j)\\
&\qquad\times P(O_{t+2:T}\mid Q_{t+1}=j)\\
&=\sum_j a_{ij}b_j(o_{t+1})\beta_{t+1}(j).\qquad\square
\end{aligned}
\]

**命题 1.5（平滑）.** 当 \(P(O=o)>0\) 时，
\[
\gamma_t(i):=P(Q_t=i\mid O=o)=\frac{\alpha_t(i)\beta_t(i)}{P(O=o)}.
\]
对 \(t<T\)，
\[
\xi_t(i,j):=P(Q_t=i,Q_{t+1}=j\mid O=o)
=\frac{\alpha_t(i)a_{ij}b_j(o_{t+1})\beta_{t+1}(j)}{P(O=o)}.
\]

**证明.** 由条件独立性，
\[
P(O=o,Q_t=i)
=P(O_{1:t},Q_t=i)P(O_{t+1:T}\mid Q_t=i,O_{1:t})
=\alpha_t(i)\beta_t(i).
\]
两边除以 \(P(O=o)\) 得第一式。第二式把 \(O_{1:t}\)、边 \(i\to j\)、\(O_{t+1}\)、\(O_{t+2:T}\) 四段乘起来，再除以 \(P(O=o)\) 即可。\(\square\)

### 1.4 Viterbi：路径最优与逐点最优不是一回事

定义
\[
\delta_t(j)=\max_{q_1,\dots,q_{t-1}}
P(Q_1=q_1,\dots,Q_{t-1}=q_{t-1},Q_t=j,O_{1:t}=o_{1:t}).
\]

**定理 1.6.**
\[
\delta_1(j)=\pi_j b_j(o_1),\qquad
\delta_{t+1}(j)=b_j(o_{t+1})\max_i\bigl[\delta_t(i)a_{ij}\bigr].
\]
若
\[
\psi_{t+1}(j)\in\arg\max_i \delta_t(i)a_{ij},
\]
则从 \(\arg\max_j\delta_T(j)\) 沿 \(\psi\) 回溯得到一条使 \(P(q,O)\) 最大的路径。

**证明.** 任一在 \(t+1\) 时刻结束于 \(j\) 的路径，在 \(t\) 时刻必结束于某个 \(i\)。其联合概率等于“到 \(t\) 为止的前缀概率”乘 \(a_{ij}b_j(o_{t+1})\)。先对前缀取最大，再对 \(i\) 取最大；有限集合上的最大值可以交换次序，于是
\[
\delta_{t+1}(j)=b_j(o_{t+1})\max_i\bigl(\delta_t(i)a_{ij}\bigr).
\]
归纳得第一式。回溯指针记录每个 \((t+1,j)\) 的最优前驱；从终点的最大值出发，逐层取前驱，得到全局最优路径。\(\square\)

**注 1.7（与逐点后验的区别）.** 逐点取 \(\arg\max_i\gamma_t(i)\) 不是 Viterbi。取 \(N=2\)、\(T=2\)、\(\pi=(1/2,1/2)\)、\(A=\begin{pmatrix}0&1\\1&0\end{pmatrix}\)、\(B\) 的两行相同。此时每个时刻的后验都是均匀分布；若平局规则在 \(t=1,2\) 都选状态 1，得到的路径 \((1,1)\) 的概率为 0，因为 \(a_{11}=0\)。所以“逐点最优”甚至不一定是一条可行路径；Viterbi 优化的是整条路径的联合概率，两者不能互相替代。

### 1.5 Baum-Welch 与 EM 的单调性

固定观测 \(o\)，记 \(\lambda=(\pi,A,B)\)，\(P_\lambda(q,o)\) 为联合概率，\(P_\lambda(q\mid o)\) 为后验。定义完全数据对数似然
\[
L_c(\lambda;q,o)=\log\pi_{q_1}
+\sum_{t=1}^{T-1}\log a_{q_t q_{t+1}}
+\sum_{t=1}^{T}\log b_{q_t}(o_t).
\]
给定当前参数 \(\lambda'\)，定义辅助函数
\[
Q(\lambda\mid\lambda')
=\mathbb E_{\lambda'}\bigl[L_c(\lambda;Q,o)\mid O=o\bigr].
\]
其中期望是对 \(P_{\lambda'}(Q=\cdot\mid O=o)\) 取的。以下先假设 \(P_{\lambda'}(q\mid o)>0\) 对所有 \(q\) 成立。代码的概率下限让参数为正，但后验仍可能数值为 0，所以实现与定理有差距。

**引理 1.8.**
\[
Q(\lambda\mid\lambda')
=\sum_i\gamma_1(i)\log\pi_i
+\sum_{t=1}^{T-1}\sum_{i,j}\xi_t(i,j)\log a_{ij}
+\sum_{t=1}^{T}\sum_j\gamma_t(j)\log b_j(o_t).
\]

**证明.** 对完全数据对数似然逐项取条件期望，用 \(\gamma_t(i)=P_{\lambda'}(Q_t=i\mid O=o)\) 与 \(\xi_t(i,j)=P_{\lambda'}(Q_t=i,Q_{t+1}=j\mid O=o)\) 的定义即可。\(\square\)

**定理 1.9（M 步）.** 在概率约束下最大化 \(Q(\lambda\mid\lambda')\)，若所有分母为正，则
\[
\hat\pi_i=\gamma_1(i),\qquad
\hat a_{ij}=\frac{\sum_{t=1}^{T-1}\xi_t(i,j)}
{\sum_{t=1}^{T-1}\gamma_t(i)},\qquad
\hat b_j(k)=\frac{\sum_{t:o_t=k}\gamma_t(j)}
{\sum_{t=1}^{T}\gamma_t(j)}.
\]

**证明.** 三个参数块在约束下可分离，逐块用 Lagrange 乘子。

- 对 \(\pi\)：约束 \(\sum_i\pi_i=1\)。拉格朗日函数对 \(\pi_i\) 求导得 \(\gamma_1(i)/\pi_i-\mu=0\)，故 \(\pi_i=\gamma_1(i)/\mu\)。对 \(i\) 求和，左边为 1，右边为 \(1/\mu\)，所以 \(\mu=1\)，即 \(\hat\pi_i=\gamma_1(i)\)。
- 对 \(A\) 的第 \(i\) 行：约束 \(\sum_j a_{ij}=1\)。令 \(C_{ij}=\sum_{t=1}^{T-1}\xi_t(i,j)\)。求导得 \(C_{ij}/a_{ij}-\mu_i=0\)，于是 \(a_{ij}=C_{ij}/\mu_i\)。对 \(j\) 求和得 \(\mu_i=\sum_j C_{ij}=\sum_{t<T}\gamma_t(i)\)，即分母。
- 对 \(B\) 的第 \(j\) 行：约束 \(\sum_k b_j(k)=1\)。令 \(D_{jk}=\sum_{t:o_t=k}\gamma_t(j)\)。求导得 \(D_{jk}/b_j(k)-\nu_j=0\)，求和得 \(\nu_j=\sum_k D_{jk}=\sum_t\gamma_t(j)\)。

目标函数是 \(\sum_x c_x\log x\) 的和，在单纯形上是凹函数；上述驻点就是全局最大点（若某系数为 0，则该坐标的目标为常数，最优不唯一，代码保留旧值）。\(\square\)

**定理 1.10（EM 单调性）.** 令 \(\lambda_{n+1}\) 为上述 M 步的解。则
\[
\log P_{\lambda_{n+1}}(o)\ge \log P_{\lambda_n}(o).
\]
若后验不同且 \(Q\) 严格增加，则不等式严格。

**证明.** 对任意 \(\lambda\)，逐点有
\[
\log P_\lambda(o)=\log P_\lambda(Q,o)-\log P_\lambda(Q\mid o).
\]
对 \(P_{\lambda'}(Q\mid o)\) 取期望：
\[
\log P_\lambda(o)=Q(\lambda\mid\lambda')+H(\lambda\mid\lambda'),
\]
其中
\[
H(\lambda\mid\lambda')=-\mathbb E_{\lambda'}\bigl[\log P_\lambda(Q\mid o)\bigr].
\]
于是
\[
\begin{aligned}
\log P_\lambda(o)-\log P_{\lambda'}(o)
&=\bigl[Q(\lambda\mid\lambda')-Q(\lambda'\mid\lambda')\bigr]\\
&\quad+\mathbb E_{\lambda'}\left[\log\frac{P_{\lambda'}(Q\mid o)}{P_\lambda(Q\mid o)}\right].
\end{aligned}
\]
第二项是 \(\mathrm{KL}(P_{\lambda'}(\cdot\mid o)\|P_\lambda(\cdot\mid o))\ge 0\)。取 \(\lambda=\lambda_{n+1}\)，由 M 步的定义 \(Q(\lambda_{n+1}\mid\lambda_n)\ge Q(\lambda_n\mid\lambda_n)\)，故对数似然不降。等号需要 KL 为 0（后验相同）且 \(Q\) 取等。\(\square\)

**推论 1.11.** \(\log P_{\lambda_n}(o)\) 单调有界，故收敛；在正则条件下，极限点满足 \(Q\) 的驻点方程，即 Baum-Welch 不动点方程。它不保证全局最优，也不保证参数收敛（只保证似然收敛）。

**GMM-HMM 的扩展.** 在完全数据里再加分量标签 \(Z_t\)，完全数据对数似然多出
\[
\sum_t \log w_{Q_t Z_t}
+\log \mathcal N(x_t;\mu_{Q_t Z_t},\Sigma_{Q_t Z_t}).
\]
E 步算分量责任 \(r_t(j,k)=P(Q_t=j,Z_t=k\mid X,\lambda')\)；M 步对 \(w,\mu,\Sigma\) 分别是加权频率、加权均值和加权协方差。EM 单调性证明不变。代码只实现对角 \(\Sigma\)，并在 M 步加方差下限。

**实现上的偏差.** 代码对方差、权重和转移概率加了下限，并在数值回退时保留历史最优参数。因此定理 1.10 的严格形式只适用于“无下限、精确 M 步”的算法；工程实现只能保证“通常不降”，测试里用 `max_iter` 对照做经验检查。这一点不应被模糊掉。

### 1.6 可辨识性：标签与退化

**命题 1.12（标签置换）.** 对任意置换 \(\sigma\)，定义
\[
\pi'_i=\pi_{\sigma(i)},\quad
a'_{ij}=a_{\sigma(i)\sigma(j)},\quad
b'_i(k)=b_{\sigma(i)}(k).
\]
则对所有观测序列 \(o\)，\(P_{\lambda'}(o)=P_\lambda(o)\)。

**证明.** 在联合分布的求和里令 \(q'_t=\sigma^{-1}(q_t)\)，逐项对应；联合概率不变，对 \(q\) 求和得到相同的观测边缘。\(\square\)

所以状态编号没有内在含义；代码按状态均值排序只是固定一个报告约定。

**例 1.13（发射退化）.** 若所有状态的发射分布相同，即 \(b_j(k)=b(k)\)，则
\[
P(O=o)=\prod_{t=1}^{T}b(o_t),
\]
与 \(A,\pi\) 无关。此时转移矩阵不可辨识。这个例子说明：没有“发射分布有差异”这类条件，参数恢复是不可能的。

一般 HMM 在什么条件下可辨识，是比这里更深的课题；本仓库不依赖一般可辨识性定理，也不声称能恢复“真实状态”。

---

## 2. 渗流

### 2.1 耦合与单调性

设图 \(G=(V,E)\)，每条边 \(e\) 有独立均匀变量 \(U_e\sim\mathrm{Unif}[0,1]\)。给定 \(p\)，令 \(e\) 开放当且仅当 \(U_e\le p\)。

**命题 2.1（单调耦合）.** 若 \(p\le q\)，则存在耦合使 \(p\)-开放边集包含于 \(q\)-开放边集。因此对任意顶点 \(v\)：
\[
\theta_v(p):=P_p(v\text{ 属于无限开放团簇})\le \theta_v(q),
\]
有限图上的最大团簇尺寸、spanning 事件也单调不减。

**证明.** 用同一组 \(U_e\) 同时构造两个构型。\(U_e\le p\Rightarrow U_e\le q\)，故边集包含。连通性是边集的单调性质：增加边不会拆散已有路径。最大团簇、spanning 同理。\(\square\)

在 \(\mathbb Z^d\) 这类可迁图上，\(p_c=\inf\{p:\theta(p)>0\}\) 与顶点选择无关，并且 \(\theta\) 非降。代码里的 `hp_lattice` 与相关图扫描都依赖这条单调性。

### 2.2 一维：\(p_c=1\)

**命题 2.2.** 在 \(\mathbb Z\) 上的键渗流或点渗流，若 \(p<1\)，则几乎必然没有无限团簇；\(p=1\) 时整条线连通。因此 \(p_c=1\)。

**证明（键渗流）.** 若原点属于无限开放团簇，则对每个 \(n\ge1\)，要么边 \(0,1,\dots,n-1\) 全开，要么边 \(-1,\dots,-n\) 全开。由并集界，
\[
P(\text{原点连通到距离 }n)\le 2p^n.
\]
令 \(n\to\infty\)，概率为 0。可数并（对所有顶点）仍为 0。\(p=1\) 时所有边开，显然连通。点渗流把“边全开”换成“连续 \(n+1\) 个点全开”，概率 \(p^{n+1}\)，同样趋于 0。\(\square\)

### 2.3 Galton-Watson、\(R_0\) 与配置模型

设后代分布 \((p_k)_{k\ge0}\)，概率母函数
\[
f(s)=\sum_{k\ge0}p_k s^k,\qquad s\in[0,1].
\]
令 \(Z_n\) 为第 \(n\) 代个体数，\(Z_0=1\)，灭绝概率 \(q=P(\exists n:Z_n=0)\)。

**引理 2.3.** \(q\) 是 \(f\) 在 \([0,1]\) 上的最小不动点。

**证明.** 令 \(q_n=P(Z_n=0)=f^n(0)\)。由 \(f\) 单调且 \(q_1=f(0)\ge0\)，归纳得 \(q_n\) 单调递增；事件 \(\{Z_n=0\}\) 也单调递增，故 \(q_n\uparrow q\)。由连续性，\(q=f(q)\)。若 \(r\in[0,1]\) 满足 \(f(r)=r\)，则 \(q_0=0\le r\)，且由 \(f\) 单调，\(q_n\le r\Rightarrow q_{n+1}=f(q_n)\le f(r)=r\)。取极限得 \(q\le r\)。所以 \(q\) 是最小不动点。\(\square\)

**定理 2.4（均值判据）.** 设 \(m=f'(1)\in(0,\infty)\)，且后代分布非退化。则
\[
q=1\iff m\le1,\qquad q<1\iff m>1.
\]

**证明.** 若 \(m\le1\)：\(f\) 是凸函数，在 \(s=1\) 处的切线为 \(1+m(s-1)\)，凸性给出
\[
f(s)\ge 1+m(s-1)\ge s,\qquad s\in[0,1].
\]
于是 \(q_n=f(q_{n-1})\ge q_{n-1}\)，且 \(q_n\le1\)；极限必为 1。若 \(m=1\) 且分布非退化，则 \(f\) 严格凸，\(s<1\) 时 \(f(s)>s\)，仍得 \(q=1\)。

若 \(m>1\)：令 \(g(s)=f(s)-s\)。则 \(g(1)=0\)、\(g'(1)=m-1>0\)，故存在 \(s_0<1\) 使 \(g(s_0)<0\)。由 \(q_0=0\le s_0\) 与 \(f\) 单调，归纳得 \(q_n\le s_0\)：若 \(q_n\le s_0\)，则 \(q_{n+1}=f(q_n)\le f(s_0)<s_0\)。所以 \(q\le s_0<1\)。\(\square\)

**注 2.5（临界情形）.** \(m=1\) 且退化（每个个体恰有一个后代）时 \(f(s)=s\)，\(q=0\)；这对应一条无限直线，而不是通常的临界分支过程。定理中的“非退化”不能省。

#### 2.3.1 带传播率的扩展

每条边以概率 \(T\in[0,1]\) 保留。对“邻居”的剩余度分支过程，其母函数为
\[
h(s)=\sum_{k\ge1}\frac{k p_k}{\mu}(1-T+Ts)^{k-1},
\qquad \mu=\sum_k k p_k.
\]
其中 \(k p_k/\mu\) 是“随机一条边通向的顶点度数为 \(k\)”的尺寸偏置分布，\(k-1\) 是除来边外的子代数。于是
\[
h'(1)=T\,\frac{\sum_k k(k-1)p_k}{\mu}=:T\nu.
\]
由定理 2.4，灭绝概率 \(q\)（\(h\) 的最小不动点）小于 1 当且仅当 \(T\nu>1\)。这就是代码里的
\[
R_0=T\cdot\frac{\langle k(k-1)\rangle}{\langle k\rangle}.
\]
在树近似下，度为 \(k\) 的根不与巨分量相连的概率是 \((1-T+Tq)^k\)，所以
\[
S=1-\sum_k p_k(1-T+Tq)^k.
\]
对 Poisson 度分布 \(p_k=e^{-c}c^k/k!\)，有 \(h(s)=e^{-cT(1-s)}\)，\(q=1-S\)，两式化为
\[
S=1-e^{-cTS}.
\]

**定理 2.6（配置模型，引用）.** 在适当的有限二阶矩与正则条件下，配置模型的巨分量比例依概率收敛到上式给出的 \(S\)；当 \(T\nu>1\) 时 \(S>0\)，当 \(T\nu<1\) 时 \(S=0\)。这是 Molloy–Reed 型结果，本笔记不证明。代码计算的是树近似值；在 \(T=1\) 且存在大量二度链时，树近似会忽略环，可能高估 \(S\)，此时应把结果当作平均场诊断而不是有限图的精确值。

**定理 2.7（ER 随机图，引用）.** 对 \(G(n,p)\)，若 \(np\to c\)，则 \(c>1\) 时最大团簇占比趋于 \(S(c)\)，其中 \(S=1-e^{-cS}\)；\(c<1\) 时最大团簇为 \(o(n)\)。这是 Erdős–Rényi 定理，代码用 `hp_perc_er_giant_theory` 计算其树近似解并用模拟核对。

### 2.4 二维方格：\(p_c=1/2\) 是定理，不是本文的证明

2D 键渗流的 \(p_c=1/2\) 来自自对偶性与临界点唯一性（Kesten）。完整证明需要 RSW 估计与对偶格点论证，篇幅远超这份笔记。代码只做数值标定：对 \(L=32,64,128,256\) 用 spanning 概率 0.5 的交叉点外推，得到 0.4997；这是验证，不是证明。点渗流的 \(p_c\approx0.5927\) 同样引用已知结果。

### 2.5 相关矩阵阈值图

给定相关矩阵 \(\rho\) 与阈值 \(\theta\)，令
\[
E_\theta=\{\{i,j\}:|\rho_{ij}|\ge\theta\}.
\]
若 \(\theta_1\le\theta_2\)，则 \(E_{\theta_2}\subseteq E_{\theta_1}\)。由命题 2.1 的同一论证，巨分量比例 \(G(\theta)\) 关于 \(\theta\) 非增。因此若存在 \(\theta\) 使 \(G(\theta)\ge1/2>G(\theta')\)，则“首次跌破 1/2”的位置是良定义的；代码用线性插值给出 `theta_half`。

注意 susceptibility \(\chi(\theta)=\sum_{s\ne s_{\max}}s^2/n\) 并不单调，且在小 \(n\) 上取值很粗。`theta_peak=argmax χ` 只作为诊断，不作为主指标。

---

## 3. 因果滤波与回测时序

令 \(\mathcal F_t=\sigma(O_1,\dots,O_t)\)。

**命题 3.1.** 固定模型参数时，由前向递推，\(\alpha_t(j)\) 只依赖 \(O_1,\dots,O_t\)。因此滤波概率
\[
P(Q_t=j\mid O_1,\dots,O_t)=\frac{\alpha_t(j)}{\sum_i\alpha_t(i)}
\]
是 \(\mathcal F_t\)-可测的。

**证明.** 对 \(t\) 归纳：\(\alpha_1\) 只含 \(O_1\)；\(\alpha_{t+1}\) 由 \(\alpha_t\) 与 \(O_{t+1}\) 线性组合而成。\(\square\)

**命题 3.2（回测时序）.** 设 \(s_t\) 是 \(\mathcal F_t\)-可测的信号，\(r_{t+1}\) 是 \(t+1\) 期收益。策略在 \(t+1\) 的收益取 \(s_t r_{t+1}\)（再扣换手成本）时，只使用了 \(t\) 时刻及以前的信息。

**证明.** \(s_t\) 的构造只涉及 \(\mathcal F_t\)；\(r_{t+1}\) 在 \(t+1\) 收盘后才已知。乘积中没有任何 \(t+1\) 时刻之后的信息进入 \(s_t\)。\(\square\)

实现上有两个容易破坏这个性质的细节：

1. 平滑后验 \(\gamma_t\) 用了全序列，不能用于 \(s_t\)；代码的仓位只用 `filtered`。
2. 标准化 scaler 若在全样本上估计，会把未来信息带入早期特征；代码只在训练段估计并固定。

模型参数的重估也必须是“只用重估点之前的数据”。代码在样本外每个 block 开始时重估，block 内的信号用重估前已经固定的参数和截至当日的观测。报告里的 `oos_start` 是第一个有样本外收益的日期；第一个决策在 `train_end` 收盘做出，收益从 `train_end+1` 开始。

---

## 4. 哪些结论没有证明

为避免误读，明确列出：

- 2D 键渗流 \(p_c=1/2\)、点渗流 \(p_c\approx0.5927\)：引用经典结果，代码只做数值验证。
- ER/配置模型巨分量的渐近定理：引用 Erdős–Rényi、Molloy–Reed；代码计算树近似，不证明有限图收敛。
- 一般 HMM 的可辨识性：只证明了标签置换与一个退化反例；没有证明“满足某条件即可恢复参数”。
- EM 全局收敛、参数收敛速率：只证明似然单调；加下限后的实现不在定理的严格假设内。
- 有限尺寸标度 \(p_c(L)=p_c+aL^{-1/\nu}\) 与 \(\nu=4/3\)：作为拟合模型使用，未证明误差界。
- 回测收益、Sharpe、最大回撤：单一合成样本上的数值结果，不是期望业绩，也不是统计显著性的证明。
- 状态图的中位阈值、渗流闸门的 80% 分位：工程启发式，没有最优性定理。

这份清单会随代码一起维护；新增结论时，要么补证明，要么补引用，要么明确标成经验观察。

---

## 参考文献

- H. Kesten. The critical probability of bond percolation on the square lattice equals 1/2. *Communications in Mathematical Physics* 74 (1980), 41–59.
- G. Grimmett. *Percolation*. 2nd ed., Springer, 1999.（耦合、一维、临界现象的标准参考）
- P. Erdős, A. Rényi. On the evolution of random graphs. *Publications of the Mathematical Institute of the Hungarian Academy of Sciences* 5 (1960), 17–60.
- M. Molloy, B. Reed. A critical point for random graphs with a given degree sequence. *Random Structures & Algorithms* 6 (1995), 161–180.
- A. P. Dempster, N. M. Laird, D. B. Rubin. Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society B* 39 (1977), 1–38.
