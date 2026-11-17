# 数据关联算法对比

双站纯方位多目标跟踪系统（IMM-UKF + 双站 CI 融合）下的数据关联算法对比代码与结果。
实验结论见 [`实验结果.md`](实验结果.md)。

## 目录

```
matlab/数据关联算法对比/
├── association_likelihood.m      公共波门 + 运动学/特征似然
├── assoc_nn.m                    局部最近邻 (NN)
├── assoc_gnn.m                   全局最近邻 (GNN, 匈牙利)
├── assoc_pda.m                   逐航迹概率数据关联 (PDA)
├── associate_measurements.m      关联统一入口 (按 params.assoc_method 分发)
├── assoc_process_platform.m      非 JPDA 方法的单平台流程
├── eval_assoc.m                  关联层指标
├── compare_assoc.m               对比实验驱动
├── test_assoc.m                  一致性 / 回归测试
├── results/                      实验 .mat 与 .png/.svg
└── 实验结果.md                   实验结论
```

核心流程的入口在 `matlab/process_platform.m`：`assoc_method='jpda'` 时直接调用原有
`jpda_process_platform.m`（默认路径数值不变），其余方法走本目录的
`assoc_process_platform.m`。`process_platform.m` 会把本目录加入 MATLAB 搜索路径。

## 运行

从仓库根目录：

```matlab
addpath('matlab');
addpath('matlab/数据关联算法对比');

test_assoc                                  % 一致性/回归测试

% 主实验: 4 方法 x 特征开/关, 固定平台, 30 种子 (默认场景见 DEFAULT_SCENARIO = 'parallel')
r = compare_assoc(struct('mode','main','n_mc',30, ...
    'save_mat','matlab/数据关联算法对比/results/main_parallel.mat', ...
    'save_fig','matlab/数据关联算法对比/results/main_parallel'));

% 其他模式
r = compare_assoc(struct('mode','seeded','n_mc',30));     % 真值播种、禁止新生 (只比关联+滤波)
r = compare_assoc(struct('mode','gamma','n_mc',20));      % 特征强度 gamma 扫描 (默认 crossing)
r = compare_assoc(struct('mode','scenario','n_mc',10));   % 4 场景
r = compare_assoc(struct('mode','clutter','n_mc',20));    % 杂波扫描
r = compare_assoc(struct('mode','ecm','n_mc',10));        % 机动平台端到端
r = compare_assoc(struct('mode','complexity'));           % JPDA 事件枚举复杂度
```

也可以直接切换单次仿真：

```matlab
params = init_parameters(struct('assoc_method','gnn','use_features',false));
[pr, vr, gos] = run_simulation(params);
```

## 方法

| 方法 | `assoc_method` | β 形式 |
|---|---|---|
| 局部最近邻 | `'nn'` | one-hot, 航迹间可复用同一量测 |
| 全局最近邻 | `'gnn'` | one-hot, 匈牙利互斥分配 |
| 逐航迹 PDA | `'pda'` | 门内量测后验概率, 航迹间独立 |
| 联合 PDA | `'jpda'` | 联合事件边缘化 (默认) |

特征似然 `L_total = L_kin * L_feature` 由所有方法共用；`use_features` 控制其开关，
同时也会影响 `manage_tracks` 的频谱出生校验（系统原有行为）。

## 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `params.assoc_method` | `'jpda'` | 关联算法 |
| `params.assoc_unassoc_threshold` | 0.5 | 量测被任一航迹"占用"的判定阈值 |
| `params.assoc_jpda_max_events` | 50000 | JPDA 事件数上界; 超出退化为 PDA, `Inf` 关闭 |

### 场景切换

`compare_assoc.m` 顶部有一个集中参数：

```matlab
DEFAULT_SCENARIO = 'parallel';   % 'maneuver' | 'parallel' | 'crossing' | 'speed'
```

它控制 `main` / `seeded` / `clutter` / `ecm` 模式的默认场景（当前为 `parallel`）。
单次运行也可直接覆盖：

```matlab
r = compare_assoc(struct('scenario','crossing'));   % 覆盖默认场景
r = compare_assoc(struct('mode','gamma','scenario','parallel')); % gamma 默认 crossing
```

## 实验设计要点

- **配对比较**：主实验固定平台 (`mobile_platforms=false`)，同种子下量测序列与关联算法无关，
  差异只来自关联层。
- **seeded 模式**：第 1 帧用真值初始化航迹并设 `max_new_tracks=0`，排除出生/鬼影，
  只比较关联与滤波本身；与 main 模式对照可区分"关联质量"与"航迹管理后果"。
- **gamma 模式**：固定 `use_features=true`，只扫描 `feature_gamma`（默认 crossing 场景），
  展示特征从弱到强时四种方法的 GOSPA 如何从分离走向重合。
- **ECM 模式**：平台机动由航迹估计决定，不同方法的量测序列不同，只报均值±方差，不做配对检验。
- **指标**：GOSPA 及 loc/missed/false 分解、位置 RMSE、关联正确率 (剔除 β0 的归一化质量)、
  ID 切换、逐帧关联耗时。
- 结果 `.mat` 可用 `compare_assoc(struct('load_mat','...'))` 直接重绘。
