function plot_results(truth, tracks_history, measurements_history, params, plat_history)
%PLOT_RESULTS 生成全部可视化图表
%   PLOT_RESULTS(truth, tracks_history, measurements_history, params, plat_history)
%
%   图表:
%     Fig 1: 真实轨迹 vs 估计轨迹 (二维俯视图)
%     Fig 2: 方位角量测时间序列
%     Fig 3: 位置 RMSE 随时间变化
%     Fig 4: 速度 RMSE 随时间变化 (Vx / Vy 子图)
%     Fig 5: GOSPA 指标时间序列 (总量 + 分解分量 + 目标/航迹数)
%     Fig 6: IMM 模型概率与转弯率估计 (仅多模型时)
%     Fig 7: 频谱模板对比 (真值 vs 估计, 仅 use_features=true 时)

K = params.K;

%% 提取估计航迹 (航迹到目标的关联基于最近邻)
% 收集被确认为对应真实目标的航迹
confirmed_tracks = cell(K, 1);
for k = 1:K
    trk = tracks_history{k};
    confirmed_tracks{k} = {};
    for t = 1:length(trk)
        if trk{t}.confirmed
            confirmed_tracks{k}{end + 1} = trk{t};
        end
    end
end

% 将航迹与目标关联: 在末时刻选择距离最近的两个确认航迹
final_tracks = confirmed_tracks{K};
est_x = zeros(params.num_targets, K);
est_y = zeros(params.num_targets, K);
est_vx = zeros(params.num_targets, K);
est_vy = zeros(params.num_targets, K);
track_assigned = false(1, length(final_tracks));
tgt_track_id = zeros(1, params.num_targets);   % 每个目标对应的航迹 ID

for tgt = 1:params.num_targets
    best_dist = inf;
    best_trk = 0;
    for t = 1:length(final_tracks)
        if track_assigned(t)
            continue;
        end
        dist = norm(final_tracks{t}.x(1:2) - truth.pos{tgt}(:, K));
        if dist < best_dist
            best_dist = dist;
            best_trk = t;
        end
    end
    if best_trk > 0
        track_assigned(best_trk) = true;
        % 回溯该航迹在所有时刻的状态
        trk_id = final_tracks{best_trk}.id;
        tgt_track_id(tgt) = trk_id;
        for k = 1:K
            found = false;
            trk_list = tracks_history{k};
            for t = 1:length(trk_list)
                if trk_list{t}.id == trk_id
                    est_x(tgt, k)  = trk_list{t}.x(1);
                    est_y(tgt, k)  = trk_list{t}.x(2);
                    est_vx(tgt, k) = trk_list{t}.x(3);
                    est_vy(tgt, k) = trk_list{t}.x(4);
                    found = true;
                    break;
                end
            end
            if ~found
                est_x(tgt, k)  = NaN;
                est_y(tgt, k)  = NaN;
                est_vx(tgt, k) = NaN;
                est_vy(tgt, k) = NaN;
            end
        end
    end
end

%% 图1: 真实轨迹 vs 估计轨迹
figure('Name', '轨迹对比', 'Position', [100, 100, 900, 700]);
hold on;

colors = lines(params.num_targets);
h_legend = gobjects(0);  % 收集图例句柄
leg_labels = {};

% 平台 (移动轨迹 或 静态标记)
if nargin >= 5 && isfield(params, 'mobile_platforms') && params.mobile_platforms
    % 移动平台: 绘制完整轨迹
    h_legend(end+1) = plot(plat_history.pos1(1, :), plat_history.pos1(2, :), 'k-', 'LineWidth', 1.5);
    leg_labels{end+1} = '平台1';
    plot(plat_history.pos1(1, 1), plat_history.pos1(2, 1), 'ks', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(plat_history.pos1(1, end), plat_history.pos1(2, end), 'ko', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'g');

    h_legend(end+1) = plot(plat_history.pos2(1, :), plat_history.pos2(2, :), 'k--', 'LineWidth', 1.5);
    leg_labels{end+1} = '平台2';
    plot(plat_history.pos2(1, 1), plat_history.pos2(2, 1), 'kd', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(plat_history.pos2(1, end), plat_history.pos2(2, end), 'ko', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'g');
else
    % 静态平台: 固定标记
    h_legend(end+1) = plot(params.plat1_pos(1), params.plat1_pos(2), 'ks', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'k');
    leg_labels{end+1} = '平台1';
    h_legend(end+1) = plot(params.plat2_pos(1), params.plat2_pos(2), 'kd', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'k');
    leg_labels{end+1} = '平台2';
end

% 方位线 (视线): 间隔若干步从两平台指向目标, 展示双站三角定位几何
if ~isfield(params, 'plot_bearing_lines') || params.plot_bearing_lines
    if isfield(params, 'plot_bearing_every')
        step_line = max(1, round(params.plot_bearing_every));
    else
        step_line = 20;
    end
    frames = unique([step_line:step_line:K, K]);

    c_ray = [0.62, 0.62, 0.62];
    h_ray1 = gobjects(0);
    h_ray2 = gobjects(0);

    for k = frames
        if nargin >= 5 && isfield(params, 'mobile_platforms') && params.mobile_platforms
            p1 = plat_history.pos1(:, k);
            p2 = plat_history.pos2(:, k);
        else
            p1 = params.plat1_pos;
            p2 = params.plat2_pos;
        end

        for tgt = 1:params.num_targets
            q = truth.pos{tgt}(:, k);
            % 略微延长, 画成"方位线"而不是线段
            q1 = q + 0.05 * (q - p1);
            q2 = q + 0.05 * (q - p2);

            h = plot([p1(1), q1(1)], [p1(2), q1(2)], '-', ...
                'Color', c_ray, 'LineWidth', 0.6);
            if isempty(h_ray1), h_ray1 = h; end

            h = plot([p2(1), q2(1)], [p2(2), q2(2)], '-.', ...
                'Color', c_ray, 'LineWidth', 0.6);
            if isempty(h_ray2), h_ray2 = h; end
        end
    end

    h_legend(end+1) = h_ray1;
    leg_labels{end+1} = '方位线 (平台1)';
    h_legend(end+1) = h_ray2;
    leg_labels{end+1} = '方位线 (平台2)';
end

% 起点/终点标记 (仅第1个目标各画一次入图例, 其余目标只画不入图例)
h_start = gobjects(1, params.num_targets);
h_end   = gobjects(1, params.num_targets);

for tgt = 1:params.num_targets
    % 真实轨迹
    h_legend(end+1) = plot(truth.pos{tgt}(1, :), truth.pos{tgt}(2, :), '-', ...
        'Color', colors(tgt, :), 'LineWidth', 2);
    leg_labels{end+1} = sprintf('目标%d 真实', tgt);

    % 起点 (圆圈), 终点 (方块)
    h_start(tgt) = plot(truth.pos{tgt}(1, 1), truth.pos{tgt}(2, 1), 'o', ...
        'Color', colors(tgt, :), 'MarkerSize', 8, 'MarkerFaceColor', colors(tgt, :));
    h_end(tgt) = plot(truth.pos{tgt}(1, end), truth.pos{tgt}(2, end), 's', ...
        'Color', colors(tgt, :), 'MarkerSize', 10, 'MarkerFaceColor', colors(tgt, :));

    % 估计轨迹
    valid = ~isnan(est_x(tgt, :));
    if any(valid)
        h_legend(end+1) = plot(est_x(tgt, valid), est_y(tgt, valid), '--', ...
            'Color', colors(tgt, :), 'LineWidth', 1.5);
        leg_labels{end+1} = sprintf('目标%d 估计', tgt);
    end
end

% 图例加入起点/终点标记
h_legend(end+1) = h_start(1);
leg_labels{end+1} = '起点';
h_legend(end+1) = h_end(1);
leg_labels{end+1} = '终点';

legend(h_legend, leg_labels, 'Location', 'best'); grid on; axis equal;
xlim(params.surveil_region(1:2));
ylim(params.surveil_region(3:4));
xlabel('X轴 (m)'); ylabel('Y轴 (m)');
if isfield(params, 'scenario_label')
    title(sprintf('轨迹对比 — %s', params.scenario_label));
end

%% 图2: 方位角量测时间序列
figure('Name', '方位角量测', 'Position', [100, 100, 900, 400]);

time = params.dt * (1:K);

subplot(2, 1, 1);
hold on;
% 真实方位角 (无噪声)
for tgt = 1:params.num_targets
    true_bearing = zeros(1, K);
    for k = 1:K
        tx = truth.X{tgt}(1, k);
        ty = truth.X{tgt}(2, k);
        true_bearing(k) = atan2(ty - params.plat1_pos(2), tx - params.plat1_pos(1));
    end
    plot(time, true_bearing * 180 / pi, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
end
% 量测
for k = 1:K
    z1 = measurements_history{k, 1};
    if ~isempty(z1)
        plot(time(k) * ones(1, length(z1)), [z1.bearing] * 180 / pi, 'k.', 'MarkerSize', 4);
    end
end
xlabel('时间 (s)'); ylabel('方位角 (deg)');
title('平台1 方位角量测'); grid on;
legend_str1 = cell(1, params.num_targets);
for tgt = 1:params.num_targets
    legend_str1{tgt} = sprintf('目标%d 真实', tgt);
end
legend([legend_str1, {'量测'}], 'Location', 'best');

subplot(2, 1, 2);
hold on;
for tgt = 1:params.num_targets
    true_bearing = zeros(1, K);
    for k = 1:K
        tx = truth.X{tgt}(1, k);
        ty = truth.X{tgt}(2, k);
        true_bearing(k) = atan2(ty - params.plat2_pos(2), tx - params.plat2_pos(1));
    end
    plot(time, true_bearing * 180 / pi, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
end
for k = 1:K
    z2 = measurements_history{k, 2};
    if ~isempty(z2)
        plot(time(k) * ones(1, length(z2)), [z2.bearing] * 180 / pi, 'k.', 'MarkerSize', 4);
    end
end
xlabel('时间 (s)'); ylabel('方位角 (deg)');
title('平台2 方位角量测'); grid on;

%% 图3: 位置 RMSE
figure('Name', '位置 RMSE', 'Position', [100, 100, 600, 400]);

for tgt = 1:params.num_targets
    rmse = sqrt((est_x(tgt, :) - truth.pos{tgt}(1, :)).^2 + ...
                (est_y(tgt, :) - truth.pos{tgt}(2, :)).^2);
    plot(time, rmse, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
    hold on;
end
xlabel('时间 (s)'); ylabel('位置 RMSE (m)');
title('位置估计 RMSE');
legend_str = cell(1, params.num_targets);
for tgt = 1:params.num_targets
    legend_str{tgt} = sprintf('目标%d', tgt);
end
legend(legend_str, 'Location', 'best');
grid on;

%% 图4: 速度 RMSE
figure('Name', '速度 RMSE', 'Position', [100, 100, 700, 500]);

for tgt = 1:params.num_targets
    rmse_vx = abs(est_vx(tgt, :) - truth.X{tgt}(3, :));
    rmse_vy = abs(est_vy(tgt, :) - truth.X{tgt}(4, :));

    subplot(2, 1, 1);
    plot(time, rmse_vx, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
    hold on;

    subplot(2, 1, 2);
    plot(time, rmse_vy, '-', 'Color', colors(tgt, :), 'LineWidth', 1.5);
    hold on;
end

subplot(2, 1, 1);
xlabel('时间 (s)'); ylabel('V_x RMSE (m/s)');
title('速度估计 RMSE — X 方向');
legend_str = cell(1, params.num_targets);
for tgt = 1:params.num_targets
    legend_str{tgt} = sprintf('目标%d', tgt);
end
legend(legend_str, 'Location', 'best'); grid on;

subplot(2, 1, 2);
xlabel('时间 (s)'); ylabel('V_y RMSE (m/s)');
title('速度估计 RMSE — Y 方向');
legend(legend_str, 'Location', 'best'); grid on;

%% 图5: GOSPA 指标 (多目标跟踪综合评价)
% 与 RMSE 不同, GOSPA 同时惩罚定位误差、漏检目标和虚假航迹,
% 且无需事先指定"哪条航迹对应哪个目标" (最优关联由算法自动给出).
gos = eval_gospa(truth, tracks_history, params);
c_cut = gos.cfg.c;

figure('Name', 'GOSPA 指标', 'Position', [100, 100, 900, 750]);

% 子图1: 总 GOSPA
subplot(3, 1, 1);
plot(time, gos.gospa, 'k-', 'LineWidth', 1.8); hold on;
plot([time(1), time(end)], [c_cut, c_cut], 'r--', 'LineWidth', 1.2);
ylabel('GOSPA (m)');
title(sprintf('GOSPA (p=%g, \\alpha=%g, c=%g m) — 均值 %.2f m, 均方根 %.2f m', ...
    gos.cfg.p, gos.cfg.alpha, c_cut, gos.mean, gos.rms));
legend({'GOSPA', '截断参数 c'}, 'Location', 'best');
grid on; xlim([time(1), time(end)]);
ylim([0, max([gos.gospa(:); c_cut]) * 1.1 + 1]);

% 子图2: 分解分量 (定位误差 / 漏检 / 虚假航迹)
subplot(3, 1, 2);
plot(time, gos.loc, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4); hold on;
plot(time, gos.missed, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
plot(time, gos.false, ':', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8);
ylabel('分量 (m)');
title('GOSPA 分解: 定位误差 / 漏检目标 / 虚假航迹');
legend({'定位误差', '漏检目标', '虚假航迹'}, 'Location', 'best');
grid on; xlim([time(1), time(end)]);

% 子图3: 目标数与航迹数
subplot(3, 1, 3);
plot(time, gos.n_truth, 'k-', 'LineWidth', 1.4); hold on;
plot(time, gos.n_est, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4);
plot(time, gos.n_false, 'r--', 'LineWidth', 1.4);
xlabel('时间 (s)'); ylabel('个数');
title(sprintf('目标数 / 航迹数 (平均虚假航迹 %.2f 条)', gos.mean_n_false));
legend({'真实目标', '参与评价航迹', '虚假航迹'}, 'Location', 'best');
grid on; xlim([time(1), time(end)]);

%% 图6: IMM 模型概率与转弯率估计 (仅多模型时)
if isfield(params, 'mm') && numel(params.mm) > 1
    n_models = numel(params.mm);
    model_names = cell(1, n_models);
    for j = 1:n_models
        model_names{j} = params.mm(j).name;
    end

    % 真实机动段: 用于阴影标注
    turns_all = cell(1, params.num_targets);
    for tgt = 1:params.num_targets
        if isfield(params, 'target_turns') && ...
                (~isfield(params, 'maneuver_enable') || params.maneuver_enable)
            turns_all{tgt} = params.target_turns{tgt};
        else
            turns_all{tgt} = zeros(0, 3);
        end
    end

    figure('Name', 'IMM 模型概率', 'Position', [100, 100, 900, 780]);

    % ---- 上两幅: 各目标的模型概率 ----
    for tgt = 1:params.num_targets
        subplot(3, 1, tgt);
        hold on;

        % 机动段阴影
        for s = 1:size(turns_all{tgt}, 1)
            patch([turns_all{tgt}(s,1), turns_all{tgt}(s,2), ...
                   turns_all{tgt}(s,2), turns_all{tgt}(s,1)], ...
                  [0, 0, 1, 1], [0.85 0.85 0.85], 'EdgeColor', 'none');
        end

        % 该目标对应航迹的模型概率历史
        mu_hist = NaN(n_models, K);
        for k = 1:K
            trk_list = tracks_history{k};
            for t = 1:length(trk_list)
                if trk_list{t}.id == tgt_track_id(tgt) && isfield(trk_list{t}, 'imm')
                    mu_hist(:, k) = trk_list{t}.imm.mu;
                    break;
                end
            end
        end

        styles = {'-', '--', ':', '-.'};
        h = gobjects(1, n_models);
        for j = 1:n_models
            h(j) = plot(time, mu_hist(j, :), styles{min(j, numel(styles))}, ...
                'Color', colors(tgt, :), 'LineWidth', 1.8);
        end

        xlim([time(1), time(end)]); ylim([0, 1]);
        ylabel('模型概率');
        title(sprintf('目标%d 模型概率 (阴影 = 真实机动段)', tgt));
        grid on;
        if tgt == 1
            legend(h, model_names, 'Location', 'best');
        end
    end

    % ---- 第三幅: 转弯率估计 vs 真值 ----
    subplot(3, 1, 3);
    hold on;
    h_om = gobjects(0);
    leg_om = {};
    for tgt = 1:params.num_targets
        omega_est = NaN(1, K);
        for k = 1:K
            trk_list = tracks_history{k};
            for t = 1:length(trk_list)
                if trk_list{t}.id == tgt_track_id(tgt)
                    omega_est(k) = trk_list{t}.x(5);
                    break;
                end
            end
        end

        h_om(end+1) = plot(time, omega_est * 180 / pi, '-', ...
            'Color', colors(tgt, :), 'LineWidth', 1.8);
        leg_om{end+1} = sprintf('目标%d 估计', tgt);

        if isfield(truth, 'omega')
            h_om(end+1) = plot(time, truth.omega{tgt} * 180 / pi, '--', ...
                'Color', colors(tgt, :), 'LineWidth', 1.1);
            leg_om{end+1} = sprintf('目标%d 真实', tgt);
        end
    end
    xlabel('时间 (s)'); ylabel('\omega (deg/s)');
    title('转弯率估计 (实线) vs 真值 (虚线)');
    legend(h_om, leg_om, 'Location', 'best');
    grid on; xlim([time(1), time(end)]);
end

%% 图7: 频谱模板对比 (真值 vs 估计)
if isfield(params, 'use_features') && params.use_features && isfield(truth, 'spectrum')

    figure('Name', '频谱对比', 'Position', [100, 100, 1000, 250 * params.num_targets]);

    for tgt = 1:params.num_targets
        F_true = truth.spectrum{tgt};

        % 找该目标对应的估计频谱
        F_est = NaN(params.N_features, 1);
        final_trk = tracks_history{K};
        for t = 1:length(final_trk)
            if final_trk{t}.confirmed && isfield(final_trk{t}, 'feature_init') ...
                    && final_trk{t}.feature_init
                dist = norm(final_trk{t}.x(1:2) - truth.pos{tgt}(:, K));
                if dist < 300
                    F_est = final_trk{t}.spectrum;
                    break;
                end
            end
        end

        % LOFAR 子图
        subplot(params.num_targets, 2, 2*tgt - 1);
        fl = params.lofar_freqs;
        if any(~isnan(F_est))
            plot(fl, F_est(1:params.N_lofar), 'b-', 'LineWidth', 1.5); hold on;
        end
        plot(fl, F_true(1:params.N_lofar), 'r--', 'LineWidth', 1.5);
        set(gca, 'XScale', 'log');
        xlabel('频率 (Hz)'); ylabel('归一化幅度');
        title(sprintf('目标%d LOFAR (蓝=估计, 红=真值)', tgt));
        grid on;

        % DEMON 子图
        subplot(params.num_targets, 2, 2*tgt);
        fd = params.demon_freqs;
        demon_idx = params.N_lofar + (1:params.N_demon);
        if any(~isnan(F_est))
            plot(fd, F_est(demon_idx), 'b-', 'LineWidth', 1.5); hold on;
        end
        plot(fd, F_true(demon_idx), 'r--', 'LineWidth', 1.5);
        xlabel('频率 (Hz)'); ylabel('归一化幅度');
        title(sprintf('目标%d DEMON (蓝=估计, 红=真值)', tgt));
        grid on;
    end
end

%% 打印平均 RMSE 汇总
fprintf('\n========== 平均 RMSE 汇总 ==========\n');
fprintf('%-8s %-16s %-16s %-16s\n', '目标', '位置RMSE(m)', '速度Vx RMSE(m/s)', '速度Vy RMSE(m/s)');
fprintf('-----------------------------------------------\n');

for tgt = 1:params.num_targets
    pos_rmse = sqrt((est_x(tgt, :) - truth.X{tgt}(1, :)).^2 + ...
                    (est_y(tgt, :) - truth.X{tgt}(2, :)).^2);
    vx_rmse  = abs(est_vx(tgt, :) - truth.X{tgt}(3, :));
    vy_rmse  = abs(est_vy(tgt, :) - truth.X{tgt}(4, :));

    % 仅用有效估计值 (航迹存在的时间段)
    valid_pos = ~isnan(pos_rmse);
    valid_vx  = ~isnan(vx_rmse);
    valid_vy  = ~isnan(vy_rmse);

    avg_pos = mean(pos_rmse(valid_pos));
    avg_vx  = mean(vx_rmse(valid_vx));
    avg_vy  = mean(vy_rmse(valid_vy));

    fprintf('目标%-3d  %-16.2f %-16.2f %-16.2f\n', tgt, avg_pos, avg_vx, avg_vy);
end

% 总体平均值
all_pos = []; all_vx = []; all_vy = [];
for tgt = 1:params.num_targets
    pos_rmse = sqrt((est_x(tgt, :) - truth.X{tgt}(1, :)).^2 + ...
                    (est_y(tgt, :) - truth.X{tgt}(2, :)).^2);
    vx_rmse  = abs(est_vx(tgt, :) - truth.X{tgt}(3, :));
    vy_rmse  = abs(est_vy(tgt, :) - truth.X{tgt}(4, :));
    all_pos = [all_pos, pos_rmse(~isnan(pos_rmse))];  %#ok<AGROW>
    all_vx  = [all_vx,  vx_rmse(~isnan(vx_rmse))];    %#ok<AGROW>
    all_vy  = [all_vy,  vy_rmse(~isnan(vy_rmse))];    %#ok<AGROW>
end
fprintf('-----------------------------------------------\n');
fprintf('%-8s %-16.2f %-16.2f %-16.2f\n', '总体平均', mean(all_pos), mean(all_vx), mean(all_vy));
fprintf('===============================================\n');

%% 打印 GOSPA 汇总 (与 RMSE 互为补充: 含漏检/虚警代价)
fprintf('\n========== GOSPA 指标汇总 (p=%g, alpha=%g, c=%g m) ==========\n', ...
    gos.cfg.p, gos.cfg.alpha, gos.cfg.c);
fprintf('平均 GOSPA (m)                %10.2f\n', gos.mean);
fprintf('  |- 定位误差分量 (m)         %10.2f\n', gos.mean_loc);
fprintf('  |- 漏检目标分量 (m)         %10.2f\n', gos.mean_missed);
fprintf('  |- 虚假航迹分量 (m)         %10.2f\n', gos.mean_false);
fprintf('GOSPA 均方根 / 最大值 (m)     %10.2f / %.2f\n', gos.rms, gos.max);
fprintf('平均漏检目标数 / 虚假航迹数   %10.2f / %.2f\n', ...
    gos.mean_n_missed, gos.mean_n_false);
fprintf('评价航迹: %s\n', ternary(gos.cfg.use_confirmed, ...
    '仅已确认航迹', '全部航迹 (含暂定)'));
fprintf('============================================================\n');

fprintf('可视化完成。\n');
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

