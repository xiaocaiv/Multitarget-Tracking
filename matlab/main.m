% main.m — 双平台协同多目标跟踪系统主脚本
% 基于纯方位角量测, 使用 IMM-UKF + JPDA + 三角测量初始化
%
% 运行方式:
%   >> main
%
% 模块文件依赖:
%   init_parameters.m, generate_truth.m, generate_measurements.m
%   generate_sigma_points.m, motion_cv.m, motion_ct.m, build_motion_models.m
%   imm_init.m, imm_predict.m, imm_measurement_predict.m, imm_update.m,
%   imm_combine.m, ci_fuse.m
%   measurement_bearing.m, feasible_events.m, jpda_associate.m
%   jpda_process_platform.m, triangulate_init.m, manage_tracks.m
%   gospar.m, eval_gospa.m, plot_results.m
%
% 评价指标:
%   RMSE  — 单目标定位/测速精度 (需先按最近邻指定航迹-目标对应)
%   GOSPA — 多目标跟踪综合指标 (定位误差 + 漏检 + 虚假航迹, 无需指定对应)

clear; close all; rng(42);
fprintf('===== 双平台协同多目标跟踪系统 =====\n');

%% ---- 场景选择 ----
% 'maneuver' — 双目标协调转弯 (默认)
% 'parallel' — 两目标匀速平行运动
% 'crossing' — 两目标匀速交叉运动
% 'speed'    — 一目标匀速, 另一目标变速
SCENARIO = 'crossing';

%% ---- 特征辅助开关 ----
% true  = 启用 LOFAR+DEMON 频谱特征辅助 JPDA 关联
% false = 纯方位角跟踪
USE_FEATURES = true;

%% 初始化
params = init_parameters(struct('scenario', SCENARIO, ...
    'use_features', USE_FEATURES));
truth = generate_truth(params);
fprintf('场景: %s (%s)\n', params.scenario_label, params.scenario);

%% 移动平台初始化
if params.mobile_platforms
    params.plat1_heading_hist = zeros(1, params.K);
    params.plat2_heading_hist = zeros(1, params.K);
    params.plat1_pos_hist = zeros(2, params.K);
    params.plat2_pos_hist = zeros(2, params.K);
    params.plat1_heading_hist(1) = params.plat1_init_heading;
    params.plat2_heading_hist(1) = params.plat2_init_heading;
    params.plat1_pos_hist(:, 1) = params.plat1_pos;
    params.plat2_pos_hist(:, 1) = params.plat2_pos;
end

fprintf('仿真时长: %d s, 步长: %d s, 总步数: %d\n', ...
    params.T_sim, params.dt, params.K);
fprintf('运动模型: IMM [%s], 转移矩阵对角 %.2f\n', ...
    strjoin(params.imm_models, '+'), params.imm_transition(1, 1));
if params.mobile_platforms
    fprintf('目标数: %d, 平台初始位置: [%.0f, %.0f], [%.0f, %.0f]\n', ...
        params.num_targets, params.plat1_pos, params.plat2_pos);
else
    fprintf('目标数: %d, 平台位置: [%.0f, %.0f], [%.0f, %.0f]\n', ...
        params.num_targets, params.plat1_pos, params.plat2_pos);
end

% 航迹容器
tracks = {};
next_track_id = 1;

% 历史记录 (可视化用)
tracks_history = cell(params.K, 1);
measurements_history = cell(params.K, 2);
plat_history = struct();
plat_history.pos1 = zeros(2, params.K);
plat_history.pos2 = zeros(2, params.K);
plat_history.head1 = zeros(1, params.K);
plat_history.head2 = zeros(1, params.K);

%% 主仿真循环
fprintf('\n仿真进行中...\n');

for k = 1:params.K

    % ---------- 0. 移动平台: 覆写当前步位置 ----------
    if params.mobile_platforms
        params.plat1_pos = params.plat1_pos_hist(:, k);
        params.plat2_pos = params.plat2_pos_hist(:, k);
    end

    % ---------- 1. 生成量测 ----------
    [z1, z2] = generate_measurements(truth, params, k);
    measurements_history{k, 1} = z1;
    measurements_history{k, 2} = z2;

    % ---------- 2. IMM 预测所有航迹 (保存预测模型供融合使用) ----------
    num_trk = length(tracks);
    pred_imm = cell(num_trk, 1);
    for t = 1:num_trk
        tracks{t}.imm = imm_predict(tracks{t}.imm, params);
        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);
        pred_imm{t} = tracks{t}.imm;
    end

    % ---------- 3. 深拷贝航迹供平台2并行处理 ----------
    tracks2 = cell(num_trk, 1);
    for t = 1:num_trk
        tracks2{t} = tracks{t};
    end

    % ---------- 4. 双平台并行 JPDA + IMM-UKF 更新 ----------
    [tracks,  unassoc1] = jpda_process_platform(tracks,  z1, params.plat1_pos, params);
    [tracks2, unassoc2] = jpda_process_platform(tracks2, z2, params.plat2_pos, params);

    % ---------- 5. 协方差交叉融合 (逐模型) ----------
    n_models = numel(params.mm);
    for t = 1:num_trk
        % 模型概率: 两平台后验的几何平均 (相当于两平台似然各用一次)
        mu1 = tracks{t}.imm.mu;
        mu2 = tracks2{t}.imm.mu;
        mu = sqrt(max(mu1, 1e-12) .* max(mu2, 1e-12));
        tracks{t}.imm.mu = mu / sum(mu);

        % 逐模型融合: P^-1 = P1^-1 + P2^-1 - P_pred^-1
        for j = 1:n_models
            [x_f, P_f] = ci_fuse(...
                tracks{t}.imm.x(:, j), tracks{t}.imm.P(:, :, j), ...
                tracks2{t}.imm.x(:, j), tracks2{t}.imm.P(:, :, j), ...
                pred_imm{t}.x(:, j), pred_imm{t}.P(:, :, j));
            tracks{t}.imm.x(:, j) = x_f;
            tracks{t}.imm.P(:, :, j) = P_f;
        end

        % 组合估计
        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);

        % 合并 hit/miss (OR 逻辑: 任一平台检测到即算命中)
        det1 = tracks{t}.beta0_last < 0.5;
        det2 = tracks2{t}.beta0_last < 0.5;
        if det1 || det2
            tracks{t}.hit_count = tracks{t}.hit_count + 1;
            tracks{t}.miss_count = 0;
        else
            tracks{t}.miss_count = tracks{t}.miss_count + 1;
        end

        % 合并频谱 (仅在特征辅助时生效)
        if params.use_features && tracks{t}.feature_init && tracks2{t}.feature_init
            spec_avg = (tracks{t}.spectrum + tracks2{t}.spectrum) / 2;
            tracks{t}.spectrum = spec_avg / norm(spec_avg);
        end

        % 清除临时字段
        tracks{t} = rmfield(tracks{t}, 'beta0_last');
    end

    % ---------- 6. 航迹管理 ----------
    [tracks, next_track_id] = manage_tracks(tracks, z1, z2, ...
        unassoc1, unassoc2, params, next_track_id);

    % ---------- 7. ECM 优化: 决定下一步平台位置 ----------
    if params.mobile_platforms && k < params.K
        h1 = params.plat1_heading_hist(k);
        [h1_next, p1_next] = platform_ecm_optimize(tracks, ...
            params.plat1_pos, h1, params.plat1_speed, params);
        params.plat1_heading_hist(k+1) = h1_next;
        params.plat1_pos_hist(:, k+1) = p1_next;

        h2 = params.plat2_heading_hist(k);
        [h2_next, p2_next] = platform_ecm_optimize(tracks, ...
            params.plat2_pos, h2, params.plat2_speed, params);
        params.plat2_heading_hist(k+1) = h2_next;
        params.plat2_pos_hist(:, k+1) = p2_next;
    end

    % ---------- 8. 记录 ----------
    tracks_history{k} = tracks;
    plat_history.pos1(:, k) = params.plat1_pos;
    plat_history.pos2(:, k) = params.plat2_pos;
    if params.mobile_platforms
        plat_history.head1(k) = params.plat1_heading_hist(k);
        plat_history.head2(k) = params.plat2_heading_hist(k);
    end

    if mod(k, 10) == 0
        num_conf = sum(cellfun(@(t) t.confirmed, tracks));
        if params.use_features
            num_feat = sum(cellfun(@(t) isfield(t, 'feature_init') && t.feature_init, tracks));
        else
            num_feat = 0;
        end
        fprintf('  步 %3d/%d | 航迹: %d (已确认: %d, 特征: %d)\n', ...
            k, params.K, length(tracks), num_conf, num_feat);
    end
end

fprintf('仿真完成。\n');

%% 可视化
plot_results(truth, tracks_history, measurements_history, params, plat_history);
