# 扩展点

这份文档只讲“在哪儿改”，不重复算法推导。数学见 [mathematics.md](mathematics.md)，整体结构见 [design.md](design.md)。

## 1. 核心库与应用层的边界

- `core/`、`hmm/`、`perc/` 不知道金融数据。新算法优先放在这三层。
- `market/` 是数据/回测适配层，可以做领域假设。
- `app/` 是编排层，可以为了一个策略写得很具体；不要把策略常数塞进核心库。

判断标准：如果把 `market/` 和 `app/` 整个删掉，`hmm/` 和 `perc/` 还能编译、测试、被别的语言调用，说明边界是干净的。

## 2. 新增一种发射分布

以 full-covariance Gaussian 为例，改动集中在 `hmm/ghmm.c`：

1. 在 `hp_ghmm` 里加协方差参数，或加一个指向扩展参数的指针；
2. `emission_row` / `compute_tables` 换成带 `logdet(Σ)` 和 Cholesky 求解的密度；
3. EM 的 M 步更新完整协方差，并加 ridge 或特征值下限；
4. `hp_ghmm_validate` 检查正定性；
5. 模型文件加版本行（例如 `hp-ghmm 2`），`hp_ghmm_load` 按版本分派。

推理函数签名（`forward/posterior/filter/viterbi`）不用变。测试至少要有：

- 非对角为 0 时与对角模型数值一致；
- 协方差正定性/病态输入；
- 小规模暴力路径似然；
- 保存/加载往返。

如果新发射族无法塞进 `hp_ghmm`，另开 `hmm/xxx.h`，不要改旧结构体的字段顺序。旧结构体只允许在末尾追加字段。

## 3. 新增渗流观测量

以“团簇尺寸直方图”为例：

1. `hp_lattice_stats` 末尾追加字段，或另写 `hp_lattice_histogram()`；
2. 在 `hp_lattice_analyze` 中计算；
3. 在 `hp_perc_lab` 的 CSV 中加列；
4. 用手工小格点写一个断言测试。

如果要做三角格点/三维格点，建议新开 `perc/lattice_tri.h`，复用 `hp_uf`，不要在一堆 `if (kind == ...)` 里改现有二维实现。

图上的点渗流（按顶点删除）目前没有实现；如果要做，基于 `hp_graph` 加一个“顶点存活”版本，而不是把格点逻辑硬套到 CSR 图。

## 4. 新增策略

策略分两种：

- **只改仓位映射**：`src/app/regime.cpp` 里的 `expo` 是最小改动点。它现在对 $n$ 个状态线性铺开 `0.1..1.0`；改成按状态波动/Sharpe 映射也可以。
- **新增一个策略模块**：在 `market/` 或 `app/` 加一个函数/类，输入滤波概率、特征、成本，输出目标仓位。不要动 HMM 推理。

无论哪种，必须满足：

- 只用 `filtered`，不用 `posterior`；
- 信号 $t$ 只作用于收益 $t+1$；
- 加一个“篡改未来收益，比较此前仓位”的测试（参考 `tests/test_app.cpp::test_regime_no_lookahead`）；
- 在 `metrics.json` 里记录映射参数和成本。

如果策略需要多资产，当前 `hp_regime` 的单资产/EW 路径不够用；应该新写一个组合层，而不是在 `run_regime` 里堆分支。

## 5. 新增 CLI 或绑定

### CLI

1. 复制 `apps/hp_connectedness.cpp` 的参数解析/输出骨架；
2. 用 `hp_csv_writer` 写表，`Json` 写摘要；
3. 在 `Makefile` 的 `APPS` 和 `CMakeLists.txt` 的 `foreach` 里登记（两处都要改，否则两套构建不一致）；
4. 在 `docs/usage.md` 补参数和输出。

### Python / Julia

1. 在 C 头文件里加 `extern "C"` 函数；
2. Python 在 `_capi.py` 声明 `argtypes/restype`，在 `core.py` 加 RAII 包装；
3. Julia 用 `ccall`，结构体布局要与 C 对齐；
4. 加一个“已知数值”测试，而不是只测“能跑”；
5. 更新 `docs/usage.md` 的 API 表。

当前绑定是手写 ABI 镜像。如果接口数量继续增长，可以考虑用 `cbindgen` 生成头文件，但不要自动生成包装层——包装层的生命周期/数组转换需要人来审。

## 6. 兼容性规则

- **C 公共结构体**：只允许末尾追加字段；改字段顺序/删除字段是 ABI break。
- **函数签名**：不要改已有参数顺序；新功能用新函数，或在末尾加可空参数（C 里意味着新符号）。
- **模型文件**：首行带版本；解析器按版本分派；旧版本必须还能读。
- **版本号**：`include/hp/core/version.h`、`CMakeLists.txt`、`README.md` 三处保持一致；Python/Julia 从 `hp_version()` 读，不重复。
- **测试**：改动核心数值路径后跑 `make check`、Python 测试、`scripts/reproduce.sh`；改内存布局后再跑 ASan。
- **数据/结果**：改默认参数会改变 `out/`，必须重新生成并更新 `docs/results.md` 的检查点。

## 7. 不要这样扩

- 不要在 `hp_hmm` 里塞市场特征或仓位字段。
- 不要把二维方格的特例（$p_c=1/2$）当成一般图的结论。
- 不要用 `posterior` 生成交易信号。
- 不要为了一个新策略修改 `hp_backtest` 的时序定义；先确认现有 `signal[t] -> ret[t+1]` 是否够用。
- 不要在没有测试的情况下替换 `hp_perc_config_giant` 的树近似；有限图结论需要模拟对照。

## 8. 新功能的最小验收

- [ ] 公共头文件注释说明输入/输出/所有权；
- [ ] 一个正常用例 + 一个边界用例；
- [ ] 一个已知数值/暴力枚举对照；
- [ ] `make check` 和绑定测试通过；
- [ ] 文档里写清楚“证明/引用/经验”；
- [ ] 没有新增第三方依赖。
