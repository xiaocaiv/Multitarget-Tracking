% run_scenarios.m — 多场景对比: 匀速平行 / 匀速交叉 / 一匀速一变加速
%
% 对每个场景运行 run_simulation, 统计位置 RMSE 与 GOSPA, 并用一张图
% 对比三组"真实轨迹 vs 估计轨迹"及 GOSPA 柱状图.
%
% 场景定义在 scenario_defs.m; 修改下方 scenario_list 可增删
% (例如把 'maneuver' 加进来).
%
% 运行: >> run_scenarios

clear; close all;

USE_FEATURES = true;                    % 是否启用频谱特征辅助
N_MC = 5;                               % 每个场景的蒙特卡洛次数
scenario_list = {'parallel', 'crossing', 'speed'};

n_sc = numel(scenario_list);
sc_label = cell(1, n_sc);
pos_mean = NaN(2, n_sc);
gos_mean = NaN(1, n_sc);
sims = cell(1, n_sc);

fprintf('===== 多场景对比 (N=%d, 特征辅助=%d) =====\n', N_MC, USE_FEATURES);

for s = 1:n_sc
    pr_all = NaN(2, N_MC);
    g_all = NaN(1, N_MC);

    for mc = 1:N_MC
        rng(mc);
        params = init_parameters(struct('scenario', scenario_list{s}, ...
            'use_features', USE_FEATURES));
        [pr, ~, g, sim] = run_simulation(params);
        pr_all(:, mc) = pr;
        g_all(mc) = g;
        if mc == 1
            sims{s} = sim;
            sc_label{s} = params.scenario_label;
        end
    end

    pos_mean(:, s) = mean(pr_all, 2, 'omitnan');
    gos_mean(s) = mean(g_all, 'omitnan');

    fprintf('%-16s (%-9s): 位置RMSE 目标1 %5.2f m, 目标2 %5.2f m, GOSPA %6.2f m\n', ...
        sc_label{s}, scenario_list{s}, pos_mean(1, s), pos_mean(2, s), gos_mean(s));
end

%% 轨迹与 GOSPA 对比图
figure('Name', '场景对比', 'Position', [80, 80, 1100, 820]);
colors = lines(2);

for s = 1:n_sc
    subplot(2, 2, s);
    plot_scene(sims{s}, colors);

    ttl = sprintf('%s  (GOSPA %.1f m)', sc_label{s}, gos_mean(s));
    if strcmp(scenario_list{s}, 'speed')
        v2 = sims{s}.truth.speed{2};
        ttl = {ttl, sprintf('目标2 速度 %.1f \\rightarrow %.1f \\rightarrow %.1f m/s', ...
            v2(1), max(v2), v2(end))};
    end
    title(ttl);

    axis equal; grid on;
    xlim(sims{s}.params.surveil_region(1:2));
    ylim(sims{s}.params.surveil_region(3:4));
    xlabel('X (m)'); ylabel('Y (m)');

    if s == 1
        text(0.03, 0.95, '实线: 真实   虚线: 估计   黑线: 平台', ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    end
end

subplot(2, 2, 4);
bar(gos_mean, 'FaceColor', [0.35 0.55 0.75]);
set(gca, 'XTickLabel', sc_label, 'XTickLabelRotation', 12);
ylabel('GOSPA (m)');
title(sprintf('各场景平均 GOSPA (N=%d, c=%g m)', N_MC, ...
    sims{1}.params.gospa_c));
grid on;

fprintf('完成。\n');

% ========================================================================
function plot_scene(sim, colors)
%PLOT_SCENE 在子图中画平台轨迹 + 真实轨迹 + 估计轨迹
params = sim.params;
truth = sim.truth;
K = params.K;
fh = sim.tracks_history;

hold on;

% 平台轨迹
if isfield(params, 'mobile_platforms') && params.mobile_platforms ...
        && isfield(params, 'plat1_pos_hist')
    plot(params.plat1_pos_hist(1, :), params.plat1_pos_hist(2, :), ...
        'k-', 'LineWidth', 1.2);
    plot(params.plat2_pos_hist(1, :), params.plat2_pos_hist(2, :), ...
        'k--', 'LineWidth', 1.2);
else
    plot(params.plat1_pos(1), params.plat1_pos(2), 'ks', ...
        'MarkerSize', 9, 'MarkerFaceColor', 'k');
    plot(params.plat2_pos(1), params.plat2_pos(2), 'kd', ...
        'MarkerSize', 9, 'MarkerFaceColor', 'k');
end

% 末帧按最近邻给每个目标分配一条确认航迹 (互不重复)
n_tgt = truth.num_targets;
ids = NaN(1, n_tgt);
used = false(1, numel(fh{K}));
for tgt = 1:n_tgt
    best = inf; best_idx = 0;
    for t = 1:numel(fh{K})
        tr = fh{K}{t};
        if used(t) || ~tr.confirmed, continue; end
        d = norm(tr.x(1:2) - truth.pos{tgt}(:, K));
        if d < best
            best = d; best_idx = t;
        end
    end
    if best_idx > 0
        used(best_idx) = true;
        ids(tgt) = fh{K}{best_idx}.id;
    end
end

% 真值 + 估计
for tgt = 1:n_tgt
    c = colors(tgt, :);
    plot(truth.pos{tgt}(1, :), truth.pos{tgt}(2, :), '-', ...
        'Color', c, 'LineWidth', 2);

    est = track_history(fh, ids(tgt), K);
    valid = ~isnan(est(1, :));
    if any(valid)
        plot(est(1, valid), est(2, valid), '--', 'Color', c, 'LineWidth', 1.2);
    end
end
end

% ========================================================================
function est = track_history(fh, id, K)
%TRACK_HISTORY 按航迹 ID 回溯位置历史
est = NaN(2, K);
if isnan(id), return; end
for k = 1:K
    for t = 1:numel(fh{k})
        if fh{k}{t}.id == id
            est(:, k) = fh{k}{t}.x(1:2);
            break;
        end
    end
end
end
