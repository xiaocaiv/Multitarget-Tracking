function res = compare_track_mgmt(cfg)
%COMPARE_TRACK_MGMT 航迹管理配置的配对蒙特卡洛对比 / 场景-特征扫描
%
%   res = COMPARE_TRACK_MGMT()                         % 配对模式 (speed, 特征关, 100 种子)
%   res = COMPARE_TRACK_MGMT(struct('mode','sweep'))   % 扫描模式 (4 场景 × 特征开/关)
%
%   cfg 字段 (全部可省略):
%     .mode          'paired' (默认) | 'sweep'
%     .scenario      配对模式场景名        默认 'speed'
%     .use_features  配对模式是否启用特征  默认 false
%     .n_mc          种子个数              默认 100 (sweep 默认 10)
%     .seeds         种子向量              默认 1:n_mc
%     .mgmt          [n_cfg x 3] 每行 [max_new_tracks, confirm_hits, max_misses]
%                    默认 [3 3 3; 2 3 2; 1 3 2]  (第一行为改动前的基线)
%     .labels        配置名 cell           默认由 mgmt 自动生成
%     .plot          是否绘图              默认 true
%     .save_fig      图片路径 (不含扩展名) 默认 '' (不保存)
%     .save_mat      结果 .mat 保存路径     默认 '' (不保存)
%     .load_mat      从已保存的 .mat 载入结果并直接重绘 (跳过仿真)
%
%   配对原理: 每个种子下依次运行各配置, 每次运行前重置随机数。量测生成只依赖真值与
%   检测/杂波实现, 与航迹数无关, 因此各配置面对完全相同的量测序列, 差异只来自航迹管理,
%   可直接做配对比较 (无需独立同分布假设)。
%
%   输出 res 结构体:
%     配对模式: .cfg / .met (各指标 n_mc x n_cfg) / .stat (结构体) / .sign_test
%     扫描模式: .cfg / .table (表格 cell) / .scen / .feat / .labels
%
%   示例:
%     compare_track_mgmt();                                  % 跑一遍看结果
%     compare_track_mgmt(struct('n_mc', 30, 'save_fig', 'docs/figs/track_mgmt_mc'));
%     compare_track_mgmt(struct('mode', 'sweep'));

if nargin < 1, cfg = struct(); end
cfg = fill_defaults(cfg);

switch lower(cfg.mode)
    case 'paired'
        res = run_paired(cfg);
    case 'sweep'
        res = run_sweep(cfg);
    otherwise
        error('compare_track_mgmt:badMode', '未知 mode: %s (应为 paired / sweep)', cfg.mode);
end
end

% ========================================================================
function cfg = fill_defaults(cfg)
if ~isfield(cfg, 'mode'),         cfg.mode = 'paired';        end
if ~isfield(cfg, 'scenario'),     cfg.scenario = 'speed';     end
if ~isfield(cfg, 'use_features'), cfg.use_features = false;   end
if strcmpi(cfg.mode, 'sweep')
    if ~isfield(cfg, 'n_mc'), cfg.n_mc = 10; end
else
    if ~isfield(cfg, 'n_mc'), cfg.n_mc = 100; end
end
if ~isfield(cfg, 'seeds'), cfg.seeds = 1:cfg.n_mc; end
if ~isfield(cfg, 'mgmt'),  cfg.mgmt = [3 3 3; 2 3 2; 1 3 2]; end
if ~isfield(cfg, 'plot'),  cfg.plot = true; end
if ~isfield(cfg, 'save_fig'), cfg.save_fig = ''; end
if ~isfield(cfg, 'labels')
    n = size(cfg.mgmt, 1);
    cfg.labels = cell(1, n);
    for i = 1:n
        cfg.labels{i} = sprintf('新增%d/确认%d/删除%d', ...
            cfg.mgmt(i, 1), cfg.mgmt(i, 2), cfg.mgmt(i, 3));
    end
end
end

% ========================================================================
function p = make_params(cfg, c)
p = init_parameters(struct( ...
    'scenario',       cfg.scenario, ...
    'use_features',   cfg.use_features, ...
    'max_new_tracks', cfg.mgmt(c, 1), ...
    'confirm_hits',   cfg.mgmt(c, 2), ...
    'max_misses',     cfg.mgmt(c, 3)));
end

% ========================================================================
function met = run_one(cfg, c, seed)
rng(seed);
p = make_params(cfg, c);
[pr, ~, gm, sim] = run_simulation(p);
g = sim.gos;
met = struct();
met.gos    = gm;
met.loc    = g.mean_loc;
met.miss   = g.mean_missed;
met.false  = g.mean_false;
met.nfalse = g.mean_n_false;
pr = pr(isfinite(pr));
if isempty(pr), met.prmse = NaN; else, met.prmse = mean(pr); end
met.ntrk = mean(cellfun(@(x) numel(x), sim.tracks_history));
end

% ========================================================================
function res = run_paired(cfg)
names = {'gos','loc','miss','false','nfalse','prmse','ntrk'};

if isfield(cfg, 'load_mat') && ~isempty(cfg.load_mat)
    % ---- 从缓存结果重绘: 跳过仿真, 只重算统计与出图 ----
    S = load(cfg.load_mat);
    met = S.res.met;
    if isfield(S.res.cfg, 'mgmt'),   cfg.mgmt   = S.res.cfg.mgmt;   end
    if isfield(S.res.cfg, 'labels'), cfg.labels = S.res.cfg.labels; end
    fprintf('\n从 %s 载入结果 (跳过仿真)\n', cfg.load_mat);
else
    n_mc  = numel(cfg.seeds);
    n_cfg = size(cfg.mgmt, 1);
    met = struct();
    for i = 1:numel(names), met.(names{i}) = nan(n_mc, n_cfg); end

    fprintf('\n===== 航迹管理配对对比 =====\n');
    fprintf('场景: %s | 特征: %s | 种子: %d 个 (%s)\n', cfg.scenario, ...
        onoff(cfg.use_features), n_mc, mat2str(cfg.seeds(1:min(3, end))));

    t0 = tic;
    for i = 1:n_mc
        for c = 1:n_cfg
            m = run_one(cfg, c, cfg.seeds(i));
            for k = 1:numel(names)
                met.(names{k})(i, c) = m.(names{k});
            end
        end
        if mod(i, max(1, round(n_mc / 10))) == 0
            fprintf('  %3d/%d 种子  已用 %.0f s\n', i, n_mc, toc(t0));
        end
    end
end
n_cfg = size(cfg.mgmt, 1);

%% ---- 统计 ----
stat = struct();
for c = 1:n_cfg
    g = met.gos(:, c);
    stat(c).label   = cfg.labels{c};
    stat(c).gos_mean = mean(g);
    stat(c).gos_std  = std(g, 1);
    stat(c).gos_med  = median(g);
    stat(c).gos_p90  = prctile_simple(g, 90);
    stat(c).loc_mean = mean(met.loc(:, c));
    stat(c).false_mean = mean(met.false(:, c));
    stat(c).nfalse_mean = mean(met.nfalse(:, c));
    stat(c).prmse_mean = mean(met.prmse(:, c), 'omitnan');
    stat(c).ntrk_mean  = mean(met.ntrk(:, c));
end

fprintf('\n%-22s %14s %9s %9s %8s %8s %8s\n', '配置', 'GOSPA 均值±std', ...
    '中位数', '90 分位', 'loc', 'false', '航迹/帧');
for c = 1:n_cfg
    fprintf('%-22s %7.2f±%-6.2f %9.2f %9.2f %8.2f %8.2f %8.2f\n', cfg.labels{c}, ...
        stat(c).gos_mean, stat(c).gos_std, stat(c).gos_med, stat(c).gos_p90, ...
        stat(c).loc_mean, stat(c).false_mean, stat(c).ntrk_mean);
end

%% ---- 与基线(第 1 行)的配对差异 + 符号检验 ----
base = met.gos(:, 1);
sign_test = struct();
fprintf('\n配对差异 (相对 %s):\n', cfg.labels{1});
for c = 2:n_cfg
    d = met.gos(:, c) - base;
    win = sum(d < 0);            % 该配置更优(更小)的种子数
    n   = sum(isfinite(d));
    % 符号检验 (正态近似 + 连续性校正)
    if n > 0
        z = (abs(win - n / 2) - 0.5) / sqrt(n / 4);
        pval = 2 * (1 - normcdf(max(z, 0)));
    else
        pval = NaN;
    end
    sign_test(c).label = cfg.labels{c};
    sign_test(c).mean_diff = mean(d, 'omitnan');
    sign_test(c).std_diff  = std(d(~isnan(d)), 1);
    sign_test(c).win_rate  = win / n;
    sign_test(c).p_value   = pval;
    fprintf('  %-22s 平均差 %+8.2f ± %6.2f | 更优 %3d/%d 个种子 (%.0f%%) | 符号检验 p = %.3g\n', ...
        cfg.labels{c}, sign_test(c).mean_diff, sign_test(c).std_diff, ...
        win, n, 100 * sign_test(c).win_rate, pval);
end

res = struct('cfg', cfg, 'met', met, 'stat', stat, 'sign_test', sign_test);

if isfield(cfg, 'save_mat') && ~isempty(cfg.save_mat)
    save(cfg.save_mat, 'res');
    fprintf('\n结果已保存: %s\n', cfg.save_mat);
end

if cfg.plot
    plot_paired(cfg, met, stat);
end
end

% ========================================================================
function res = run_sweep(cfg)
scens = {'speed', 'maneuver', 'parallel', 'crossing'};
feats = [false, true];
n_cfg = size(cfg.mgmt, 1);
n_mc  = numel(cfg.seeds);

G = nan(numel(scens), numel(feats), n_cfg);

fprintf('\n===== 场景 × 特征 扫描 (%d 种子/格) =====\n', n_mc);
t0 = tic;
for s = 1:numel(scens)
    for f = 1:numel(feats)
        for c = 1:n_cfg
            g = nan(n_mc, 1);
            for i = 1:n_mc
                rng(cfg.seeds(i));
                p = init_parameters(struct( ...
                    'scenario', scens{s}, 'use_features', feats(f), ...
                    'max_new_tracks', cfg.mgmt(c, 1), ...
                    'confirm_hits',   cfg.mgmt(c, 2), ...
                    'max_misses',     cfg.mgmt(c, 3)));
                [~, ~, gm] = run_simulation(p);
                g(i) = gm;
            end
            G(s, f, c) = mean(g);
            fprintf('  %-9s 特征=%-3s %-18s GOSPA=%6.2f  (%.0f s)\n', ...
                scens{s}, onoff(feats(f)), cfg.labels{c}, G(s, f, c), toc(t0));
        end
    end
end

fprintf('\n%-10s %-6s', '场景', '特征');
for c = 1:n_cfg, fprintf(' %20s', cfg.labels{c}); end
fprintf('\n');
for s = 1:numel(scens)
    for f = 1:numel(feats)
        fprintf('%-10s %-6s', scens{s}, onoff(feats(f)));
        for c = 1:n_cfg
            fprintf(' %20.2f', G(s, f, c));
        end
        fprintf('\n');
    end
end

res = struct('cfg', cfg, 'scen', {scens}, 'feat', feats, 'G', G, ...
    'labels', {cfg.labels});
end

% ========================================================================
function plot_paired(cfg, met, stat)
n_cfg = size(cfg.mgmt, 1);
col = [0.35 0.35 0.35; 0.85 0.55 0.20; 0.15 0.45 0.75;
       0.25 0.60 0.40; 0.55 0.35 0.70];
col = col(1:n_cfg, :);

fig = figure('Color', 'w', 'Position', [80 80 1180 380], 'Name', '航迹管理配对对比');

% ---- (a) 配对散点: 基线 vs 最优配置 ----
ax1 = subplot(1, 3, 1);
base = met.gos(:, 1);
c = n_cfg;
lo = min([base; met.gos(:, c)]) * 0.8;
hi = max([base; met.gos(:, c)]) * 1.2;
loglog(base, met.gos(:, c), 'o', 'MarkerSize', 5, 'MarkerFaceColor', col(c, :), ...
    'MarkerEdgeColor', 'none');
hold on;
plot([lo hi], [lo hi], 'k--', 'LineWidth', 0.9);
set(ax1, 'XScale', 'log', 'YScale', 'log', 'XLim', [lo hi], 'YLim', [lo hi]);
grid on; box on;
xlabel(sprintf('基线 GOSPA (m): %s', cfg.labels{1}), 'FontName', 'SimSun');
ylabel(sprintf('GOSPA (m): %s', cfg.labels{c}), 'FontName', 'SimSun');
title(sprintf('配对散点 (%d 种子, 对角线下方为更优)', size(base, 1)), 'FontName', 'SimSun');
n_better = sum(met.gos(:, c) < base);
text(lo * 1.15, hi / 1.6, sprintf('%d/%d 个种子更优', n_better, size(base, 1)), ...
    'FontName', 'SimSun', 'FontSize', 9, 'Color', col(c, :));

% ---- (b) GOSPA^2 分量分解 (平方和可加) ----
ax2 = subplot(1, 3, 2);
comp = zeros(n_cfg, 3);
for k = 1:n_cfg
    comp(k, 1) = mean(met.loc(:, k) .^ 2);
    comp(k, 2) = mean(met.miss(:, k) .^ 2);
    comp(k, 3) = mean(met.false(:, k) .^ 2);
end
hb = bar(comp, 'stacked');
hb(1).FaceColor = [0.30 0.55 0.75];
hb(2).FaceColor = [0.90 0.75 0.35];
hb(3).FaceColor = [0.80 0.35 0.35];
short_lbl = arrayfun(@(i) sprintf('%d/%d/%d', cfg.mgmt(i, 1), cfg.mgmt(i, 2), cfg.mgmt(i, 3)), ...
    1:n_cfg, 'UniformOutput', false);
set(ax2, 'XTick', 1:n_cfg, 'XTickLabel', short_lbl, 'FontName', 'SimSun');
xlabel(ax2, '新增 / 确认 / 删除', 'FontName', 'SimSun');
ylabel('GOSPA^2 分量均值 (m^2)', 'FontName', 'SimSun');
legend({'定位 loc^2', '漏检 miss^2', '虚假 false^2'}, 'Location', 'northwest', ...
    'FontName', 'SimSun');
title('误差平方分解', 'FontName', 'SimSun');
grid on; box on;

% ---- (c) 经验 CDF ----
ax3 = subplot(1, 3, 3);
hold on;
for k = 1:n_cfg
    g = sort(met.gos(:, k));
    stairs(g, (1:numel(g))' / numel(g), 'Color', col(k, :), 'LineWidth', 1.4);
end
set(ax3, 'XScale', 'log', 'FontName', 'SimSun');
grid on; box on;
xlabel('GOSPA (m)', 'FontName', 'SimSun');
ylabel('经验累积概率', 'FontName', 'SimSun');
legend(cfg.labels, 'Location', 'southeast', 'FontName', 'SimSun');
title('GOSPA 经验 CDF', 'FontName', 'SimSun');

if ~isempty(cfg.save_fig)
    out = cfg.save_fig;
    d = fileparts(out);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    print(fig, [out '.png'], '-dpng', '-r300');
    print(fig, [out '.svg'], '-dsvg');
    fprintf('\n图片已保存: %s.png / .svg\n', out);
end
end

% ========================================================================
function v = prctile_simple(x, p)
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
idx = max(1, min(numel(x), ceil(p / 100 * numel(x))));
v = x(idx);
end

function s = onoff(tf)
if tf, s = '开'; else, s = '关'; end
end
