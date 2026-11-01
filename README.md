# Multitarget-Tracking · 双站协同机动多目标跟踪

双平台协同、方位量测（bearing-only）的机动多目标跟踪仿真平台，MATLAB R2024a 实现。

## 方法概要

- **IMM-UKF**：CV / CT / CA 三模型滤波，共用 7 维状态 `x = [px; py; vx; vy; ω; ax; ay]`，模型混合与融合都在这一状态空间内完成
- **JPDA**：多目标数据关联，特征谱相似度作为似然的乘性项（`L_total = L_kinematic × L_feature`），特征不进入滤波状态
- **双平台融合**：各平台独立关联后，逐模型做协方差交叉（CI）融合，模型概率按归一化几何平均融合
- **航迹管理**：三角定位初始化新航迹，两平台谱型一致性校验（Pearson ρ ≥ 0.4）抑制虚假航迹
- **ECM**：平台机动控制，最大化期望信息增益
- **GOSPA**：评价指标，可分解为定位误差 / 漏检 / 虚假航迹三部分

## 目录结构

```
matlab/     全部仿真代码（.m，扁平）
data/       舰船辐射噪声样本（wav 特征模式使用）
docs/       设计与实现说明
archive/    重构前的 .bak 回滚备份，代码不引用
CLAUDE.md   给 AI 编程助手的架构说明
```

## 运行

在仓库根目录执行：

```bash
# 完整仿真（机动场景）
matlab -batch "addpath('matlab'); main"

# CV 单模型 vs IMM 的蒙特卡洛对比
matlab -batch "addpath('matlab'); compare_imm"

# 四套预设场景总览
matlab -batch "addpath('matlab'); run_scenarios"

# 杂波密度扫描
matlab -batch "addpath('matlab'); sweep_clutter"

# 两艘舰船 LOFAR+DEMON 频谱对比
matlab -batch "addpath('matlab'); compare_ship_spectra"
```

也可以先 `cd matlab` 再 `matlab -batch "main"`。在 MATLAB GUI 里只需 `addpath('matlab')` 一次。

## 依赖

- MATLAB R2024a
- Statistics and Machine Learning Toolbox（`chi2inv`、`poissrnd`）
- Signal Processing Toolbox（`spectrogram`、`hilbert`、`butter`、`filtfilt`，仅 `use_features=true` 且从 wav 导入时需要）

## 常用开关

| 参数 | 位置 | 说明 |
|---|---|---|
| `params.scenario` | `matlab/main.m` 顶部 | `maneuver` / `parallel` / `crossing` / `speed` |
| `params.use_features` | `matlab/init_parameters.m` | `false` 退回纯方位模式，关闭全部特征路径 |
| `params.imm_models` | `matlab/init_parameters.m` | `{'CV'}` 为单模型基线，`{'CV','CT','CA'}` 为完整 IMM |
| `params.lambda_c` | `matlab/init_parameters.m` | 每帧每平台期望杂波数 |
| `params.sigma_bearing` | `matlab/init_parameters.m` | 方位量测噪声标准差（默认 1°） |

## 数据路径约定

`data/` 下的 wav 文件一律通过 `matlab/repo_path.m` 解析为绝对路径：

```matlab
wav = repo_path('data', 'Cargo_40s.wav');
```

该函数从自身所在位置反推仓库根目录，与当前工作目录无关。**不要写成相对路径**——`generate_ship_spectrum` 在文件缺失时会回退到仿真频谱，路径写错会让特征实验在不知情的情况下失真（现已改为告警而非静默）。
