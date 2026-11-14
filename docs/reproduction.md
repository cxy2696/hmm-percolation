# 复现

这一页只写命令、种子和检查点。结果解释见 [results.md](results.md)。

## 1. 环境

本机验证环境：

| 项 | 版本 |
| --- | --- |
| macOS | 27.2 arm64 |
| Apple clang | 21.0.0 |
| Python | 3.13（绑定只用标准库） |
| Julia | 1.13.4（juliaup） |
| Make | GNU Make 3.81 |
| CMake | 3.31（仅用于验证 CMake 构建） |

C11 + libm，C++17 + 标准库。没有 BLAS、Eigen、numpy、网络依赖。不同 libm 在 `log/exp` 上可能有 1e-15 级差异；测试用相对容差，收益/回撤按两位小数核对。

## 2. 一键跑完

```bash
./scripts/reproduce.sh
```

脚本内容：

1. `make -j`、`make shared`、`make check`；
2. 生成 `data/synthetic_panel.csv` 与 truth；
3. 跑 connectedness；
4. 跑 3 / 2 / BIC / 无闸门 / 30 日窗口五组 regime；
5. 跑渗流 scan、键/点 scaling、meanfield；
6. Python 测试；
7. Julia 测试（`julia` 可执行时）；
8. 画图（如果 `build/venv` 存在）；
9. 打印摘要。

本机约 30 秒。

## 3. 分步

### 3.1 构建与测试

```bash
make clean
make -j
make check
```

期望：

```
30 passed, 0 failed      # C
5 passed, 0 failed       # C++ 应用层
```

开启 AddressSanitizer（可选）：

```bash
cc -std=c11 -g -O1 -fsanitize=address -Iinclude -Itests \
  src/core/*.c src/hmm/*.c src/perc/*.c src/market/*.c src/cross/*.c \
  tests/test_main.c tests/test_core.c tests/test_hmm.c tests/test_ghmm.c \
  tests/test_perc.c tests/test_market.c -lm -o build/hp_tests_asan
ASAN_OPTIONS=detect_leaks=0 ./build/hp_tests_asan
```

### 3.2 数据

```bash
./build/hp_synth --days 3000 --assets 10 --regimes 3 --seed 20240101 \
  --out data/synthetic_panel.csv --truth data/synthetic_truth.csv
```

面板 3000 行 × 10 资产，日期 2012-01-02 至 2023-06-30（工作日，不含节假日）。状态占比约为 61.9% / 22.3% / 15.8%。

### 3.3 横截面

```bash
./build/hp_connectedness --csv data/synthetic_panel.csv --window 60 \
  --outdir out/connectedness
```

期望（±0.002）：

```
theta_c mean ≈ 0.873, median ≈ 0.875, max ≈ 0.970 (2015-05-06)
2940 points
```

### 3.4 主模型

```bash
./build/hp_regime --csv data/synthetic_panel.csv --states 3 \
  --train-frac 0.6 --refit 250 --seed 42 --outdir out/regime
```

期望：

```
3 states, train 1788 days, OOS 1191 days
state 0 弱势  ann_ret -15.7%  ann_vol 34.2%
state 1 震荡  ann_ret   2.5%  ann_vol 13.0%
state 2 强势  ann_ret  14.5%  ann_vol 18.6%
buy&hold: 92.77%  Sharpe 0.85  MDD -33.83%
HMM     : 53.40%  Sharpe 0.77  MDD -22.40%
HMM+perc: 60.65%  Sharpe 0.88  MDD -16.90%
```

变体：

```bash
./build/hp_regime --states 2 --outdir out/regime_s2
./build/hp_regime --states 0 --outdir out/regime_s0    # BIC 2..5
./build/hp_regime --states 3 --no-perc --outdir out/regime_noperc
./build/hp_regime --states 3 --conn-window 30 --outdir out/regime_w30
```

### 3.5 渗流

```bash
./build/hp_perc_lab --mode scan --kind bond --L 128 \
  --pmin 0.35 --pmax 0.65 --steps 31 --reps 200 \
  --out out/perc_scan_bond.csv

./build/hp_perc_lab --mode scaling --kind bond --sizes 32,64,128,256 \
  --reps 200 --iters 16 --out out/perc_scaling_bond.csv \
  --json out/perc_scaling_bond.json

./build/hp_perc_lab --mode scaling --kind site --sizes 32,64,128,256 \
  --reps 200 --iters 16 --out out/perc_scaling_site.csv \
  --json out/perc_scaling_site.json

./build/hp_perc_lab --mode meanfield --cmin 0.5 --cmax 3 --steps 26 \
  --out out/perc_meanfield.csv
```

期望：

| 实验 | 结果 |
| --- | --- |
| 键渗流外推 | 0.4997（R²=0.9936） |
| 点渗流外推 | 0.5948（R²=0.9705） |
| ER c=2 巨分量 | 0.796812 |

### 3.6 Python / Julia

```bash
make shared
PYTHONPATH=bindings/python python3 -m unittest discover -s bindings/python/tests -v
PYTHONPATH=bindings/python python3 bindings/python/examples/regime_demo.py

julia --startup-file=no bindings/julia/test/runtests.jl
julia --startup-file=no bindings/julia/examples/regime_demo.jl
```

Python 12 项、Julia 22 项。Python 会复算 C++ 的 `model_train.txt` 对数似然（容差 2e-5）和最后一个窗口的 `theta_c`（容差 1e-10）。

如果 juliaup 启动器报 `Could not create lockfile`，给它一个可写 HOME：

```bash
mkdir -p /tmp/hp-julia-home
HOME=/tmp/hp-julia-home julia --startup-file=no bindings/julia/test/runtests.jl
```

### 3.7 图

```bash
python3 -m venv build/venv
build/venv/bin/pip install matplotlib
MPLCONFIGDIR=build/mplcache build/venv/bin/python scripts/plot_results.py
```

输出 `docs/figures/{regimes,equity,percolation,correlation}.png`。matplotlib 只用于展示。

## 4. 产物

| 路径 | 内容 |
| --- | --- |
| `out/regime/states.csv` | 日期、收益、波动、平滑后验、Viterbi、仓位、闸门、`theta_c` |
| `out/regime/equity.csv` | 样本外净值与仓位 |
| `out/regime/bic.csv` | 状态数/BIC |
| `out/regime/metrics.json` | 配置、模型、状态画像、回测、渗流摘要 |
| `out/regime/report.md` | 中文报告 |
| `out/regime/model_train.txt` | 训练段模型 |
| `out/regime/model_final.txt` | 最后一个走步窗口的模型 |
| `out/connectedness/` | 连通性序列、摘要、最后窗口相关矩阵 |
| `out/perc_*` | 渗流实验 CSV/JSON |

`model_train.txt` 和 `metrics.json.model` 是同一组参数，跨语言测试用它们做对照。

## 5. 检查点

| 文件 | 字段 | 期望 |
| --- | --- | --- |
| `out/perc_scaling_bond.json` | `pc_infinite` | 0.4997 ±0.001 |
| `out/perc_scaling_site.json` | `pc_infinite` | 0.5948 ±0.002 |
| `out/regime/metrics.json` | `model.loglik` | -3169.06 ±1e-4 |
| 同上 | `backtest.hmm_perc.max_drawdown` | -0.169 |
| 同上 | `data.oos_start` | 1809 |
| `out/connectedness/summary.json` | `theta_c_median` | 0.875 ±0.002 |

偏差超过容差时，先确认数据是否重建、种子是否一致、共享库是否与 CLI 同一次构建。

## 6. 清理

```bash
make clean                  # 删除 build/
rm -rf out                  # 删除参考产物（可选）
```

`data/` 可以保留，也可以用 `hp_synth` 重建。
