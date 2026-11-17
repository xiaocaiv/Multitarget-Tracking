function res = compare_assoc(cfg)
%COMPARE_ASSOC 数据关联算法对比实验 (配对蒙特卡洛 / 场景 / 杂波 / ECM / 复杂度)
%
%   res = COMPARE_ASSOC()                       % 主实验: 4 方法 x 特征开/关
%   res = COMPARE_ASSOC(struct('mode','seeded'))    % 真值播种+禁止新生, 只比关联+滤波
%   res = COMPARE_ASSOC(struct('mode','gamma'))     % 特征强度 gamma 扫描 (默认 crossing)
%   res = COMPARE_ASSOC(struct('mode','scenario'))  % 4 场景 x 4 方法
%   res = COMPARE_ASSOC(struct('mode','clutter'))   % 杂波密度扫描
%   res = COMPARE_ASSOC(struct('mode','ecm'))       % 机动平台端到端
%   res = COMPARE_ASSOC(struct('mode','complexity'))% JPDA 事件数增长
%
%   cfg 字段 (均可省略):
%     .mode          'main' (默认) | 'seeded' | 'gamma' | 'scenario' | 'clutter' | 'ecm' | 'complexity'
%     .methods       {'nn','gnn','pda','jpda'} (默认)
%     .features      特征开/关向量; main 默认 [false true], 其余按模式
%     .scenario      main/seeded/clutter/ecm 场景 (默认 DEFAULT_SCENARIO='parallel');
%                    gamma 模式默认 'crossing' (可在 cfg 中覆盖)
%     .scenarios     scenario 模式的场景列表
%     .lambdas       clutter 模式的 lambda_c 列表
%     .gammas        gamma 模式的特征锐化因子列表 (默认 [0 1 2 3 4 6 8])
%     .feature_gamma 非 gamma 模式下固定的 gamma (默认 8)
%     .n_mc          种子数; main 30 / seeded 30 / gamma 20 / scenario 10 / clutter 20 / ecm 10
%     .seeds         种子向量 (默认 1:n_mc)
%     .mobile_platforms  是否启用 ECM; 默认仅 ecm 模式为 true
%     .seed_truth    是否用真值初始化航迹并禁止新生; 默认仅 seeded 模式为 true
%     .baseline      配对检验基线方法 (默认 'jpda')
%     .baseline_features 基线是否带特征 (默认 false, ecm 模式 true)
%     .plot          是否绘图 (默认 true)
%     .save_fig      图片路径前缀 (不含扩展名)
%     .save_mat      结果 .mat 路径
%     .load_mat      从 .mat 载入并重绘
%     .include_extreme  complexity 是否包含 8x20 极限点 (默认 false)
%
%   配对原理: 固定平台 (mobile_platforms=false) 时, 同种子下量测序列只依赖
%   真值/参数/随机数, 与关联算法无关, 因此各方法面对完全相同的量测; 差异
%   只来自关联层, 可做配对符号检验.
%
%   输出 res: .cfg / .confs / .met / .stat / .sign_test / .mode
%
%   示例:
%     compare_assoc(struct('n_mc',10,'save_mat','matlab/数据关联算法对比/results/main.mat'));
%     compare_assoc(struct('mode','clutter','n_mc',20));

% ============================================================
% 默认场景开关: 改这一行即可切换所有非 gamma 实验的场景
%   'maneuver' | 'parallel' | 'crossing' | 'speed'
DEFAULT_SCENARIO = 'parallel';
% ============================================================

if nargin < 1, cfg = struct(); end
if ~isfield(cfg, 'scenario') && ...
        ~(isfield(cfg, 'mode') && strcmpi(cfg.mode, 'gamma'))
    cfg.scenario = DEFAULT_SCENARIO;
end
cfg = fill_defaults(cfg);

if ~isempty(cfg.load_mat)
    S = load(cfg.load_mat);
    res = S.res;
    fprintf('\n从 %s 载入结果 (跳过仿真)\n', cfg.load_mat);
    if ~isempty(cfg.save_fig)
        res.cfg.save_fig = cfg.save_fig;
    end
    if cfg.plot, plot_results_assoc(res); end
    return;
end

switch lower(cfg.mode)
    case {'main', 'seeded', 'gamma', 'scenario', 'clutter', 'ecm'}
        res = run_grid(cfg);
    case 'complexity'
        res = run_complexity(cfg);
    otherwise
        error('compare_assoc:badMode', ...
            '未知 mode: %s (应为 main/scenario/clutter/ecm/complexity)', cfg.mode);
end

if ~isempty(cfg.save_mat)
    d = fileparts(cfg.save_mat);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    save(cfg.save_mat, 'res');
    fprintf('\n结果已保存: %s\n', cfg.save_mat);
end

if cfg.plot
    plot_results_assoc(res);
end
end

% ========================================================================
function cfg = fill_defaults(cfg)
if ~isfield(cfg, 'mode'), cfg.mode = 'main'; end
if ~isfield(cfg, 'scenario')
    if strcmpi(cfg.mode, 'gamma')
        cfg.scenario = 'crossing';   % 交叉场景最能体现关联歧义
    else
        cfg.scenario = 'maneuver';
    end
end
if ~isfield(cfg, 'methods'), cfg.methods = {'nn', 'gnn', 'pda', 'jpda'}; end
if ~isfield(cfg, 'lambdas'), cfg.lambdas = [0, 1, 2, 3, 5]; end
if ~isfield(cfg, 'gammas'), cfg.gammas = [0, 1, 2, 3, 4, 6, 8]; end
if ~isfield(cfg, 'feature_gamma'), cfg.feature_gamma = 8; end
if ~isfield(cfg, 'scenarios'), cfg.scenarios = {'maneuver', 'parallel', 'crossing', 'speed'}; end
if ~isfield(cfg, 'lambda_c'), cfg.lambda_c = 2; end

if ~isfield(cfg, 'features')
    switch lower(cfg.mode)
        case 'main',  cfg.features = [false, true];
        case {'ecm', 'gamma'}, cfg.features = true;
        otherwise,    cfg.features = false;
    end
end
if ~isfield(cfg, 'mobile_platforms')
    cfg.mobile_platforms = strcmpi(cfg.mode, 'ecm');
end
if ~isfield(cfg, 'seed_truth')
    cfg.seed_truth = strcmpi(cfg.mode, 'seeded');
end
if ~isfield(cfg, 'n_mc')
    switch lower(cfg.mode)
        case {'main', 'seeded'}, cfg.n_mc = 30;
        case 'gamma',    cfg.n_mc = 20;
        case 'scenario', cfg.n_mc = 10;
        case 'clutter',  cfg.n_mc = 20;
        case 'ecm',      cfg.n_mc = 10;
        otherwise,       cfg.n_mc = 3;
    end
end
if ~isfield(cfg, 'seeds'), cfg.seeds = 1:cfg.n_mc; end
if ~isfield(cfg, 'plot'), cfg.plot = true; end
if ~isfield(cfg, 'save_fig'), cfg.save_fig = ''; end
if ~isfield(cfg, 'save_mat'), cfg.save_mat = ''; end
if ~isfield(cfg, 'load_mat'), cfg.load_mat = ''; end
if ~isfield(cfg, 'baseline'), cfg.baseline = 'jpda'; end
if ~isfield(cfg, 'baseline_gamma'), cfg.baseline_gamma = max(cfg.gammas); end
if ~isfield(cfg, 'include_extreme'), cfg.include_extreme = false; end
if ~isfield(cfg, 'baseline_features')
    cfg.baseline_features = any(strcmpi(cfg.mode, {'ecm', 'gamma'}));
end
end

% ========================================================================
function confs = build_configs(cfg)
confs = struct('method', {}, 'use_features', {}, 'scenario', {}, ...
    'lambda_c', {}, 'feature_gamma', {}, 'label', {});

switch lower(cfg.mode)
    case {'main', 'seeded'}
        for mi = 1:numel(cfg.methods)
            for f = cfg.features
                confs(end + 1) = make_conf(cfg.methods{mi}, f, ...
                    cfg.scenario, cfg.lambda_c, cfg.feature_gamma); %#ok<AGROW>
            end
        end
    case 'gamma'
        for gi = 1:numel(cfg.gammas)
            for mi = 1:numel(cfg.methods)
                c = make_conf(cfg.methods{mi}, true, ...
                    cfg.scenario, cfg.lambda_c, cfg.gammas(gi));
                c.label = sprintf('%s/γ%g', upper(cfg.methods{mi}), cfg.gammas(gi));
                confs(end + 1) = c; %#ok<AGROW>
            end
        end
    case 'scenario'
        for si = 1:numel(cfg.scenarios)
            for mi = 1:numel(cfg.methods)
                c = make_conf(cfg.methods{mi}, ...
                    cfg.features(1), cfg.scenarios{si}, cfg.lambda_c, cfg.feature_gamma);
                c.label = sprintf('%s/%s', cfg.scenarios{si}, upper(cfg.methods{mi}));
                confs(end + 1) = c; %#ok<AGROW>
            end
        end
    case 'clutter'
        for li = 1:numel(cfg.lambdas)
            for mi = 1:numel(cfg.methods)
                c = make_conf(cfg.methods{mi}, ...
                    cfg.features(1), cfg.scenario, cfg.lambdas(li), cfg.feature_gamma);
                c.label = sprintf('%s/λ%d', upper(cfg.methods{mi}), cfg.lambdas(li));
                confs(end + 1) = c; %#ok<AGROW>
            end
        end
    case 'ecm'
        for mi = 1:numel(cfg.methods)
            confs(end + 1) = make_conf(cfg.methods{mi}, ...
                cfg.features(1), cfg.scenario, cfg.lambda_c, cfg.feature_gamma); %#ok<AGROW>
        end
end
end

function c = make_conf(method, use_features, scenario, lambda_c, feature_gamma)
c = struct();
c.method = method;
c.use_features = use_features;
c.scenario = scenario;
c.lambda_c = lambda_c;
c.feature_gamma = feature_gamma;
if use_features, ft = '特征'; else, ft = '无特征'; end
c.label = sprintf('%s (%s)', upper(method), ft);
end

% ========================================================================
function res = run_grid(cfg)
confs = build_configs(cfg);
n_cfg = numel(confs);
n_mc = numel(cfg.seeds);

names = {'gos', 'loc', 'miss', 'false', 'prmse', 'vrmse', ...
    'p_correct', 'p_correct_norm', 'p_wrong', 'p_clutter', 'p_miss', ...
    'hard', 'idsw', 'ntrk', 't_assoc', 'max_events', 'n_fallback'};
met = struct();
for i = 1:numel(names)
    met.(names{i}) = NaN(n_mc, n_cfg);
end
n_fail = zeros(1, n_cfg);

fprintf('\n===== 数据关联对比: %s =====\n', upper(cfg.mode));
fprintf('场景: %s | 特征: %s | 平台: %s | 种子: %d 个 | 配置: %d 个\n', ...
    cfg.scenario, mat2str(cfg.features), ...
    ternary(cfg.mobile_platforms, '机动(ECM)', '固定'), n_mc, n_cfg);

t0 = tic;
for i = 1:n_mc
    for c = 1:n_cfg
        try
            m = run_one(cfg, confs(c), cfg.seeds(i));
            for k = 1:numel(names)
                met.(names{k})(i, c) = m.(names{k});
            end
        catch ME
            n_fail(c) = n_fail(c) + 1;
            fprintf('\n  [失败] 种子 %d / %s: %s\n', cfg.seeds(i), ...
                confs(c).label, ME.message);
        end
    end
    if mod(i, max(1, round(n_mc / 10))) == 0
        fprintf('  %3d/%d 种子  已用 %.0f s\n', i, n_mc, toc(t0));
    end
end

stat = summarize(met, confs);
sign_test = paired_sign_test(met, confs, cfg);

print_table(stat, n_fail);
print_sign_test(sign_test);

res = struct('cfg', cfg, 'confs', confs, 'met', met, ...
    'stat', stat, 'sign_test', sign_test, 'mode', lower(cfg.mode), ...
    'n_fail', n_fail);
end

% ========================================================================
function m = run_one(cfg, conf, seed)
rng(seed);
over = struct( ...
    'scenario', conf.scenario, ...
    'use_features', conf.use_features, ...
    'assoc_method', conf.method, ...
    'mobile_platforms', cfg.mobile_platforms, ...
    'lambda_c', conf.lambda_c, ...
    'feature_gamma', conf.feature_gamma, ...
    'assoc_seed_truth', cfg.seed_truth);
if cfg.seed_truth
    over.max_new_tracks = 0;   % 禁止新生, 只比较关联+滤波
end
p = init_parameters(over);

[pr, vr, ~, sim] = run_simulation(p);
a = eval_assoc(sim.truth, sim.tracks_history, sim.meas_src, ...
    sim.assoc_history, sim.params, false);

m = struct();
m.gos = sim.gos.mean;
m.loc = sim.gos.mean_loc;
m.miss = sim.gos.mean_missed;
m.false = sim.gos.mean_false;
prf = pr(isfinite(pr));
if isempty(prf), m.prmse = NaN; else, m.prmse = mean(prf); end
vrf = vr(isfinite(vr));
if isempty(vrf), m.vrmse = NaN; else, m.vrmse = mean(vrf); end
m.p_correct = a.p_correct;
m.p_correct_norm = a.p_correct_norm;
m.p_wrong = a.p_wrong;
m.p_clutter = a.p_clutter;
m.p_miss = a.p_miss;
m.hard = a.hard_correct;
m.idsw = a.id_switches;
m.ntrk = mean(cellfun(@(x) numel(x), sim.tracks_history));

t_assoc = [];
events = [];
fb = 0;
for k = 1:sim.params.K
    for pp = 1:2
        info = sim.assoc_history{k, pp};
        t_assoc(end + 1) = info.t_assoc; %#ok<AGROW>
        events(end + 1) = info.n_events; %#ok<AGROW>
        fb = fb + info.fallback;
    end
end
m.t_assoc = mean(t_assoc);
m.max_events = max(events);
m.n_fallback = fb;
end

% ========================================================================
function stat = summarize(met, confs)
n_cfg = numel(confs);
stat = struct([]);
for c = 1:n_cfg
    g = met.gos(:, c);
    s = struct();
    s.label = confs(c).label;
    s.method = confs(c).method;
    s.use_features = confs(c).use_features;
    s.gos_mean = mean(g, 'omitnan');
    s.gos_std = std(g, 0, 1, 'omitnan');
    s.gos_med = median(g, 'omitnan');
    s.gos_p90 = prctile_simple(g, 90);
    s.loc_mean = mean(met.loc(:, c), 'omitnan');
    s.miss_mean = mean(met.miss(:, c), 'omitnan');
    s.false_mean = mean(met.false(:, c), 'omitnan');
    s.prmse_mean = mean(met.prmse(:, c), 'omitnan');
    s.p_correct_mean = mean(met.p_correct(:, c), 'omitnan');
    s.p_correct_norm_mean = mean(met.p_correct_norm(:, c), 'omitnan');
    s.hard_mean = mean(met.hard(:, c), 'omitnan');
    s.idsw_mean = mean(met.idsw(:, c), 'omitnan');
    s.t_assoc_mean = mean(met.t_assoc(:, c), 'omitnan');
    s.max_events = max(met.max_events(:, c));
    s.fallback_mean = mean(met.n_fallback(:, c), 'omitnan');
    stat = [stat, s]; %#ok<AGROW>
end
end

% ========================================================================
function st = paired_sign_test(met, confs, cfg)
n_cfg = numel(confs);
if strcmpi(cfg.mode, 'gamma')
    % gamma 模式的基线: 指定 gamma 处的 baseline 方法 (默认最强特征)
    base = find(strcmpi({confs.method}, cfg.baseline) & ...
        [confs.feature_gamma] == cfg.baseline_gamma, 1);
else
    base = find(strcmpi({confs.method}, cfg.baseline) & ...
        [confs.use_features] == cfg.baseline_features, 1);
end
if isempty(base), base = 1; end

st = struct([]);
for c = 1:n_cfg
    if c == base
        continue;
    end
    d = met.gos(:, c) - met.gos(:, base);
    win = sum(d < 0);              % 该配置更优(更小)的种子数
    n = sum(isfinite(d));
    if n > 0
        z = (abs(win - n / 2) - 0.5) / sqrt(n / 4);
        pval = 2 * (1 - normcdf(max(z, 0)));
    else
        pval = NaN;
    end
    s = struct();
    s.label = confs(c).label;
    s.baseline = confs(base).label;
    s.mean_diff = mean(d, 'omitnan');
    s.std_diff = std(d(isfinite(d)), 0, 1);
    s.win_rate = win / max(n, 1);
    s.p_value = pval;
    st = [st, s]; %#ok<AGROW>
end
end

% ========================================================================
function print_table(stat, n_fail)
fprintf('\n%-24s %15s %9s %8s %8s %9s %8s %8s\n', ...
    '配置', 'GOSPA 均值±std', 'loc', 'false', 'RMSE', '正确率', 'ID切换', '关联ms');
fprintf('%s\n', repmat('-', 1, 96));
for c = 1:numel(stat)
    s = stat(c);
    fprintf('%-24s %7.2f±%-6.2f %9.2f %8.2f %8.2f %9.3f %8.2f %8.2f', ...
        s.label, s.gos_mean, s.gos_std, s.loc_mean, s.false_mean, ...
        s.prmse_mean, s.p_correct_norm_mean, s.idsw_mean, 1000 * s.t_assoc_mean);
    if n_fail(c) > 0
        fprintf('  失败%d', n_fail(c));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('-', 1, 96));
end

% ========================================================================
function print_sign_test(st)
if isempty(st)
    return;
end
fprintf('\n配对差异 (相对基线, GOSPA 越小越好):\n');
for i = 1:numel(st)
    fprintf('  %-24s 平均差 %+8.2f ± %6.2f | 更优 %3.0f%% | 符号检验 p = %.3g\n', ...
        st(i).label, st(i).mean_diff, st(i).std_diff, ...
        100 * st(i).win_rate, st(i).p_value);
end
end

% ========================================================================
function plot_results_assoc(res)
if ~isfield(res, 'confs') || ~isfield(res, 'met')
    return;
end
cfg = res.cfg;
confs = res.confs;
met = res.met;
n_cfg = numel(confs);
labels = {confs.label};

if strcmpi(cfg.mode, 'gamma')
    plot_gamma(res);
    return;
end

fig = figure('Color', 'w', 'Position', [80, 80, 1200, 760], ...
    'Name', '数据关联算法对比');

% (a) GOSPA
subplot(2, 2, 1);
mu = mean(met.gos, 1, 'omitnan');
sd = std(met.gos, 0, 1, 'omitnan');
bar(mu); hold on;
errorbar(1:n_cfg, mu, sd, 'k.', 'LineWidth', 1);
set(gca, 'XTick', 1:n_cfg, 'XTickLabel', labels, ...
    'XTickLabelRotation', 25, 'FontName', 'SimSun');
ylabel('GOSPA (m)'); title('GOSPA'); grid on; box on;

% (b) 关联正确率
subplot(2, 2, 2);
bar(mean(met.p_correct_norm, 1, 'omitnan'));
set(gca, 'XTick', 1:n_cfg, 'XTickLabel', labels, ...
    'XTickLabelRotation', 25, 'FontName', 'SimSun');
ylabel('正确关联概率质量'); title('关联正确率'); grid on; box on;

% (c) 逐帧关联耗时
subplot(2, 2, 3);
bar(1000 * mean(met.t_assoc, 1, 'omitnan'));
set(gca, 'XTick', 1:n_cfg, 'XTickLabel', labels, ...
    'XTickLabelRotation', 25, 'FontName', 'SimSun');
ylabel('关联耗时 (ms/帧)'); title('关联层耗时'); grid on; box on;

% (d) 精度-耗时散点 (颜色=方法, 实心=特征开, 空心=特征关)
subplot(2, 2, 4);
methods = unique({confs.method}, 'stable');
cols = lines(max(3, numel(methods)));
hold on;
for c = 1:n_cfg
    mi = find(strcmp(methods, confs(c).method), 1);
    mk = 'o';
    if confs(c).use_features
        mfc = cols(mi, :);
    else
        mfc = 'none';
    end
    plot(1000 * mean(met.t_assoc(:, c), 'omitnan'), ...
        mean(met.gos(:, c), 'omitnan'), mk, 'MarkerSize', 9, ...
        'MarkerFaceColor', mfc, 'MarkerEdgeColor', cols(mi, :), ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end

% 图例: 每个方法一个实心标记 + 一个空心标记说明特征开关
h = gobjects(numel(methods), 1);
for mi = 1:numel(methods)
    h(mi) = plot(NaN, NaN, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', cols(mi, :), 'MarkerEdgeColor', cols(mi, :), ...
        'LineWidth', 1.2, 'DisplayName', upper(methods{mi}));
end
h_open = plot(NaN, NaN, 'o', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.35 0.35 0.35], ...
    'LineWidth', 1.2, 'DisplayName', '空心=特征关');

xlabel('关联耗时 (ms/帧)'); ylabel('GOSPA (m)');
title('精度-耗时 (实心=特征开, 空心=特征关)');
legend([h(:); h_open], 'Location', 'best', 'FontName', 'SimSun');
grid on; box on;

% 杂波模式: 额外画 GOSPA vs lambda
fig2 = [];
if strcmpi(cfg.mode, 'clutter')
    fig2 = figure('Color', 'w', 'Position', [120, 120, 720, 460], ...
        'Name', '杂波密度扫描');
    methods = unique({confs.method}, 'stable');
    hold on;
    for mi = 1:numel(methods)
        sel = find(strcmp({confs.method}, methods{mi}));
        x = [confs(sel).lambda_c];
        y = mean(met.gos(:, sel), 1, 'omitnan');
        plot(x, y, 'o-', 'LineWidth', 1.5, 'DisplayName', upper(methods{mi}));
    end
    xlabel('杂波密度 \lambda_c'); ylabel('GOSPA (m)');
    title('不同杂波密度下的 GOSPA'); legend('Location', 'best');
    grid on; box on;
end

if ~isempty(cfg.save_fig)
    d = fileparts(cfg.save_fig);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    print(fig, [cfg.save_fig '.png'], '-dpng', '-r300');
    print(fig, [cfg.save_fig '.svg'], '-dsvg');
    if ~isempty(fig2)
        print(fig2, [cfg.save_fig '_lambda.png'], '-dpng', '-r300');
        print(fig2, [cfg.save_fig '_lambda.svg'], '-dsvg');
    end
    fprintf('\n图片已保存: %s.png / .svg\n', cfg.save_fig);
end
end

% ========================================================================
function plot_gamma(res)
%PLOT_GAMMA 特征强度扫描专用图: GOSPA / loc / false / 正确率 vs gamma
cfg = res.cfg;
confs = res.confs;
met = res.met;

methods = unique({confs.method}, 'stable');
cols = lines(max(3, numel(methods)));
metrics = {'gos', 'loc', 'false', 'p_correct_norm'};
titles = {'GOSPA', '定位误差 loc', '虚假航迹 false', '关联正确率'};
ylabs = {'GOSPA (m)', 'loc (m)', 'false (m)', '正确率'};

fig = figure('Color', 'w', 'Position', [80, 80, 1000, 760], ...
    'Name', '特征强度扫描');
for p = 1:numel(metrics)
    subplot(2, 2, p); hold on;
    for mi = 1:numel(methods)
        sel = find(strcmp({confs.method}, methods{mi}));
        x = [confs(sel).feature_gamma];
        y = mean(met.(metrics{p})(:, sel), 1, 'omitnan');
        plot(x, y, 'o-', 'Color', cols(mi, :), 'LineWidth', 1.6, ...
            'MarkerSize', 7, 'DisplayName', upper(methods{mi}));
    end
    xlabel('feature\_gamma'); ylabel(ylabs{p}); title(titles{p});
    grid on; box on;
    if p == 1
        legend('Location', 'best');
    end
    if p == 3
        v = mean(met.false, 1, 'omitnan');
        v = v(isfinite(v));
        if isempty(v) || max(v) <= 0
            ylim([0, 0.1]);
        else
            ylim([0, 1.2 * max(v)]);
        end
    end
end

if ~isempty(cfg.save_fig)
    d = fileparts(cfg.save_fig);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    print(fig, [cfg.save_fig '.png'], '-dpng', '-r300');
    print(fig, [cfg.save_fig '.svg'], '-dsvg');
    fprintf('\n图片已保存: %s.png / .svg\n', cfg.save_fig);
end
end

% ========================================================================
function res = run_complexity(cfg)
pairs = [2, 5; 3, 8; 4, 10; 5, 12; 6, 15];
if cfg.include_extreme
    pairs = [pairs; 8, 20];
end
gate_p = 0.6;

n_ev = zeros(size(pairs, 1), 1);
t_ev = zeros(size(pairs, 1), 1);
fprintf('\n===== JPDA 事件枚举复杂度基准 (门概率 %.2f) =====\n', gate_p);
fprintf('%-8s %-8s %-14s %s\n', '航迹', '量测', '可行事件数', '耗时(s)');
for i = 1:size(pairs, 1)
    nt = pairs(i, 1);
    nm = pairs(i, 2);
    rng(1);
    V = rand(nt, nm) < gate_p;
    tic;
    ev = feasible_events(V);
    t_ev(i) = toc;
    n_ev(i) = numel(ev);
    fprintf('%-8d %-8d %-14d %.3f\n', nt, nm, n_ev(i), t_ev(i));
end

figure('Color', 'w', 'Position', [120, 120, 720, 460], ...
    'Name', 'JPDA 复杂度');
loglog(n_ev, t_ev, 'o-', 'LineWidth', 1.6, 'MarkerSize', 8);
xlabel('可行联合事件数'); ylabel('枚举耗时 (s)');
title('JPDA 事件枚举复杂度'); grid on; box on;

if ~isempty(cfg.save_fig)
    d = fileparts(cfg.save_fig);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    print(gcf, [cfg.save_fig '.png'], '-dpng', '-r300');
    print(gcf, [cfg.save_fig '.svg'], '-dsvg');
end

res = struct('cfg', cfg, 'mode', 'complexity', 'pairs', pairs, ...
    'n_events', n_ev, 'time', t_ev, 'gate_p', gate_p);
end

% ========================================================================
function v = prctile_simple(x, p)
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
idx = max(1, min(numel(x), ceil(p / 100 * numel(x))));
v = x(idx);
end

function s = ternary(tf, a, b)
if tf, s = a; else, s = b; end
end
