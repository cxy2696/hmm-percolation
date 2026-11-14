# 用法

## 1. 构建

```bash
make -j                 # 静态库、共享库、CLI、examples
make check              # C + C++ 测试
make shared             # 只构建 libhp.dylib / libhp.so
make python-test        # 需要共享库
make figures            # 需要 build/venv 里的 matplotlib
```

CMake：

```bash
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake -j
ctest --test-dir build/cmake --output-on-failure
```

## 2. C API

总入口：

```c
#include "hp/hp.h"

printf("%s\n", hp_version());   /* HP_VERSION_STRING */
```

### 2.1 离散 HMM

```c
hp_hmm *hmm = hp_hmm_new(3, 4);       /* 3 states, 4 symbols */
hp_hmm_train_opts opt;
hp_hmm_train_opts_default(&opt);
opt.n_restarts = 4;
opt.seed = 42;

hp_hmm_report rep;
hp_status st = hp_hmm_train(obs, T, 3, 4, &opt, hmm, &rep);
if (st != HP_OK) { /* 处理错误 */ }

uint32_t *path = malloc((size_t)T * sizeof *path);
hp_hmm_viterbi(hmm, obs, T, path, NULL);

double *logalpha = malloc((size_t)T * 3 * sizeof *logalpha);
double logL = 0;
hp_hmm_forward_log(hmm, obs, T, logalpha, &logL);

hp_hmm_free(hmm);
free(path);
free(logalpha);
```

- `obs[t]` 必须在 `[0, n_symbols)`。
- `hp_hmm_train` 接受 `T >= 2`；收敛失败不返回错误，`rep.converged == 0`。
- `hp_hmm_filter` 输出因果滤波概率；`hp_hmm_posterior` 输出平滑后验，交易信号不要用后者。
- 状态编号是 0-based。

### 2.2 GMM-HMM

```c
hp_ghmm *m = hp_ghmm_new(3, 2, 1);    /* 3 states, 2 features, 1 mixture */
hp_ghmm_train_opts opt;
hp_ghmm_train_opts_default(&opt);
opt.seed = 42;

hp_ghmm_report rep;
st = hp_ghmm_train(x, T, 2, 3, 1, &opt, m, &rep);

double *filt = malloc((size_t)T * 3 * sizeof *filt);
hp_ghmm_filter(m, x, T, filt, NULL);  /* row t: P(q_t | x_1..t) */

hp_ghmm_free(m);
free(filt);
```

- `x` 是 `T*dim` 行主序，必须有限；特征缺失应在应用层处理。
- 方差下限默认 `1e-4`，在标准化特征上使用；如果特征尺度很大，调 `opt.var_floor`。
- `hp_ghmm_select` 做 BIC 选阶，`table` 由调用方分配。

### 2.3 渗流

```c
hp_lattice *lat = hp_lattice_new(128, HP_PERC_BOND, 0); /* L, kind, periodic */
hp_rng rng;
hp_rng_seed_stream(&rng, 1, 0);
hp_lattice_bernoulli(lat, 0.55, &rng);

hp_lattice_stats st;
hp_lattice_analyze(lat, &st);
printf("S=%.4f chi=%.4f spanning=%d\n", st.s_max, st.susceptibility, st.spanning);

size_t n = 0;
size_t *sizes = hp_lattice_cluster_sizes(lat, &n);
hp_free(sizes);
hp_lattice_free(lat);
```

有限尺寸标定：

```c
hp_perc_exp exp = {.L = 64, .kind = HP_PERC_BOND, .periodic = 0,
                   .reps = 200, .seed = 1};
double pc = 0, span = 0;
hp_perc_crossing(&exp, 0.35, 0.65, 16, &pc, &span);
```

相关图：

```c
double *corr = malloc((size_t)N * N * sizeof *corr);
hp_corr_matrix(returns, T, N, corr);

hp_corr_critical_result crit;
hp_corr_critical(corr, N, 0.10, 0.99, 45, &crit);
/* 应用层用 crit.theta_half 作为 theta_c */
```

`hp_perc_config_giant(pk, k_max, T, &S, &r0, &mean)` 计算配置模型的树近似；`pk` 必须归一化。

### 2.4 面板与回测

```c
hp_status st = HP_OK;
hp_panel *p = hp_panel_load_csv("data.csv", &st);
double *ew = malloc((size_t)p->T * sizeof *ew);
hp_panel_ew_returns(p, ew);

hp_bt_opts bo;
hp_bt_opts_default(&bo);
bo.cost_bps = 5.0;

hp_bt_result res;
hp_bt_metrics met;
hp_backtest(ew, p->T, signal, &bo, &res, &met);

hp_bt_result_free(&res);
hp_panel_free(p);
free(ew);
```

回测约定：`signal[t]` 在 t 日收盘确定，作用于 `ret[t+1]`。`hp_backtest` 自动检查 `ret[1..]` 有限，不检查 `ret[0]`。

## 3. CLI

### 3.1 `hp_synth`

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `--days` | 3000 | 工作日数 |
| `--assets` | 10 | 资产数 |
| `--regimes` | 3 | 2 或 3 |
| `--seed` | 20240101 | 随机种子 |
| `--out` | `data/synthetic_panel.csv` | 价格宽表 |
| `--truth` | `data/synthetic_truth.csv` | 状态/市场收益 |

### 3.2 `hp_regime`

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `--csv` | 必填 | 面板 |
| `--asset` | `EW` | 资产名或等权指数 |
| `--states` | 3 | 0 = BIC 在 2..5 中选 |
| `--mix` | 1 | GMM 分量数 |
| `--vol-window` | 20 | 滚动波动窗口 |
| `--train-frac` | 0.6 | 训练比例 |
| `--refit` | 250 | 走步重估间隔 |
| `--max-iter` | 150 | EM 上限 |
| `--restarts` | 4 | 多起点数 |
| `--seed` | 42 | 种子 |
| `--cost-bps` | 5 | 换手成本 |
| `--no-perc` | off | 关闭渗流闸门 |
| `--conn-window` | 60 | 相关窗口 |
| `--conn-theta-min/max/steps` | 0.10/0.99/45 | 阈值扫描 |
| `--conn-theta-ref` | 0.90 | 参考阈值 |
| `--risk-quantile` | 0.80 | 闸门分位 |
| `--risk-gate` | 0.5 | 触发时的仓位乘数 |
| `--outdir` | `out/regime` | 输出目录 |

输出：

- `states.csv`：日期、收益、波动、平滑后验、Viterbi、仓位、闸门、`theta_c`。
- `equity.csv`：样本外三条净值。
- `bic.csv`：状态数/BIC。
- `metrics.json`：配置、模型、状态画像、回测。
- `report.md`：中文摘要。
- `model_train.txt` / `model_final.txt`：文本模型。

### 3.3 `hp_connectedness`

| 参数 | 默认 |
| --- | --- |
| `--csv` | 必填 |
| `--window` | 60 |
| `--step` | 1 |
| `--theta-min/max/steps` | 0.10/0.99/45 |
| `--theta-ref` | 0.90 |
| `--outdir` | `out/connectedness` |

输出 `connectedness.csv`、`summary.json`、`corr_last.csv`。`theta_c` 是 0.5 巨分量交叉点，`theta_chi` 是 susceptibility 峰值。

### 3.4 `hp_perc_lab`

| 模式 | 主要参数 | 输出 |
| --- | --- | --- |
| `bond` / `site` | `--L --p --reps --seed` | 单点统计 JSON |
| `scan` | `--L --pmin --pmax --steps --reps` | p 扫描 CSV |
| `scaling` | `--sizes 32,64,128,256 --reps --iters` | 交叉点 CSV + 外推 JSON |
| `meanfield` | `--cmin --cmax --steps` | S(c) CSV |

共用的 `--out`、`--json`、`--periodic`、`--kind` 见上表；完整解析逻辑在 `apps/hp_perc_lab.cpp`，没有额外的 `--help`。

## 4. 模型文件

离散：

```
hp-hmm 1
states 3
symbols 3
pi
0.4 0.35 0.25
a
...
b
...
```

GMM：

```
hp-ghmm 1
states 3
dim 2
nmix 1
pi / a / mix / mean / var
```

数值用 `%.17g`，空白分隔。加载时校验行和与方差；格式升级改版本行。

## 5. Python

```bash
make shared
export PYTHONPATH=bindings/python
```

```python
import hpmarket as hp

h = hp.DiscreteHMM(2, 2)
logalpha, logL = h.forward([0, 1, 0])

m = hp.GaussianHMM(2, 1)
rep = m.train([[-2.0], [-1.8], [1.9], [2.1]] * 100, seed=1)
post, ll = m.posterior([[-2.0], [1.9]])

lat = hp.Lattice(64, "bond")
lat.bernoulli(0.55, hp.RNG(1))
print(lat.analyze())

crit = hp.corr_critical(hp.corr_matrix(returns_window))
print(crit["theta_half"])
```

对象需要 `close()`；也支持 `with`。`HPMARKET_LIB` 可覆盖库路径；`hpmarket.__version__` 从 C 的 `hp_version()` 读取。

## 6. Julia

```julia
include("bindings/julia/HP.jl")
using .HP

m = GaussianHMM(2, 1)
rep = train!(m, x; seed=42)
post, ll = posterior(m, x)        # 2 × T
path, lp = viterbi(m, x)          # 0-based states

lat = Lattice(64, :bond)
bernoulli!(lat, 0.55, RNG(1))
println(analyze(lat))
```

`HP_LIB` 可覆盖库路径；`HP.lib_version()` 从 C 读取版本。Julia 的 `params(m).mean` 是 `(states, mix, dim)`；`state_mean(m, j)` 的 `j` 从 0 开始。

## 7. 退出码

CLI 正常返回 0；参数/文件/数值错误返回 1，错误信息到 stderr。`hp_tests` 失败返回 1，`hp_tests_app` 同样。
