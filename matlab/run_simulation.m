function [pos_rmse, vel_rmse, gos_mean, sim] = run_simulation(params)
%RUN_SIMULATION 运行一次仿真, 返回两目标的位置/速度 RMSE 与 GOSPA
%   [pos_rmse, vel_rmse] = RUN_SIMULATION(params)
%   [pos_rmse, vel_rmse, gos_mean] = RUN_SIMULATION(params)
%   [pos_rmse, vel_rmse, gos_mean, sim] = RUN_SIMULATION(params)
%   pos_rmse: [2x1] 每个目标的位置 RMSE (m)
%   vel_rmse: [2x1] 每个目标的速度 RMSE (m/s)
%   gos_mean: 标量, 逐帧 GOSPA 的时间平均 (m), 见 eval_gospa.m
%   sim:      仿真内部数据 (tracks_history / truth / params), 供调试与画图

truth = generate_truth(params);

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

tracks = {};
next_track_id = 1;
tracks_history = cell(params.K, 1);
assoc_history = cell(params.K, 2);   % 每帧每平台的关联记录 (beta / beta0 / 耗时)
meas_src = cell(params.K, 2);        % 每帧每平台量测的来源目标 (0=杂波)

% 关联层受控实验: 用真值初始化航迹 (配合 max_new_tracks=0 禁止新生)
if isfield(params, 'assoc_seed_truth') && params.assoc_seed_truth
    tracks = seed_tracks_from_truth(truth, params);
    next_track_id = numel(tracks) + 1;
end

for k = 1:params.K
    % 0. 移动平台: 覆写当前步位置
    if params.mobile_platforms
        params.plat1_pos = params.plat1_pos_hist(:, k);
        params.plat2_pos = params.plat2_pos_hist(:, k);
    end

    % 1. 生成量测
    [z1, z2] = generate_measurements(truth, params, k);

    % 2. IMM 预测 + 保存预测模型
    num_trk = length(tracks);
    pred_imm = cell(num_trk, 1);
    for t = 1:num_trk
        tracks{t}.imm = imm_predict(tracks{t}.imm, params);
        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);
        pred_imm{t} = tracks{t}.imm;
    end

    % 3. 深拷贝航迹
    tracks2 = cell(num_trk, 1);
    for t = 1:num_trk
        tracks2{t} = tracks{t};
    end

    % 4. 双平台并行处理 (关联算法由 params.assoc_method 选择)
    [tracks,  unassoc1, info1] = process_platform(tracks,  z1, params.plat1_pos, params);
    [tracks2, unassoc2, info2] = process_platform(tracks2, z2, params.plat2_pos, params);

    % 5. 逐模型协方差交叉融合
    n_models = numel(params.mm);
    for t = 1:num_trk
        mu1 = tracks{t}.imm.mu;
        mu2 = tracks2{t}.imm.mu;
        mu = sqrt(max(mu1, 1e-12) .* max(mu2, 1e-12));
        tracks{t}.imm.mu = mu / sum(mu);

        for j = 1:n_models
            [x_f, P_f] = ci_fuse(...
                tracks{t}.imm.x(:, j), tracks{t}.imm.P(:, :, j), ...
                tracks2{t}.imm.x(:, j), tracks2{t}.imm.P(:, :, j), ...
                pred_imm{t}.x(:, j), pred_imm{t}.P(:, :, j));
            tracks{t}.imm.x(:, j) = x_f;
            tracks{t}.imm.P(:, :, j) = P_f;
        end

        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);

        det1 = tracks{t}.beta0_last < 0.5;
        det2 = tracks2{t}.beta0_last < 0.5;
        if det1 || det2
            tracks{t}.hit_count = tracks{t}.hit_count + 1;
            tracks{t}.miss_count = 0;
        else
            tracks{t}.miss_count = tracks{t}.miss_count + 1;
        end

        if params.use_features && tracks{t}.feature_init && tracks2{t}.feature_init
            spec_avg = (tracks{t}.spectrum + tracks2{t}.spectrum) / 2;
            tracks{t}.spectrum = spec_avg / norm(spec_avg);
        end

        tracks{t} = rmfield(tracks{t}, 'beta0_last');
    end

    % 6. 航迹管理
    [tracks, next_track_id] = manage_tracks(tracks, z1, z2, ...
        unassoc1, unassoc2, params, next_track_id);

    tracks_history{k} = tracks;
    assoc_history{k, 1} = info1;
    assoc_history{k, 2} = info2;
    meas_src{k, 1} = reshape([z1.src], 1, []);
    meas_src{k, 2} = reshape([z2.src], 1, []);

    % ECM 优化: 决定下一步平台位置
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
end

%% 计算 RMSE: 末帧最近邻匹配航迹到目标
pos_rmse = NaN(params.num_targets, 1);
vel_rmse = NaN(params.num_targets, 1);

final_trk = tracks_history{params.K};
confirmed_trk = {};
for t = 1:length(final_trk)
    if final_trk{t}.confirmed
        confirmed_trk{end + 1} = final_trk{t};
    end
end

assigned = false(1, length(confirmed_trk));
for tgt = 1:params.num_targets
    best_dist = inf;
    best_idx = 0;
    for t = 1:length(confirmed_trk)
        if assigned(t), continue; end
        d = norm(confirmed_trk{t}.x(1:2) - truth.pos{tgt}(:, params.K));
        if d < best_dist
            best_dist = d;
            best_idx = t;
        end
    end
    if best_idx > 0
        assigned(best_idx) = true;
        trk_id = confirmed_trk{best_idx}.id;

        est_x = NaN(4, params.K);
        for k = 1:params.K
            trk_list = tracks_history{k};
            for t = 1:length(trk_list)
                if trk_list{t}.id == trk_id
                    est_x(:, k) = trk_list{t}.x(1:4);
                    break;
                end
            end
        end

        valid = ~any(isnan(est_x), 1);
        if any(valid)
            pos_err = sqrt((est_x(1, valid) - truth.X{tgt}(1, valid)).^2 ...
                         + (est_x(2, valid) - truth.X{tgt}(2, valid)).^2);
            vel_err = sqrt((est_x(3, valid) - truth.X{tgt}(3, valid)).^2 ...
                         + (est_x(4, valid) - truth.X{tgt}(4, valid)).^2);
            pos_rmse(tgt) = mean(pos_err);
            vel_rmse(tgt) = mean(vel_err);
        end
    end
end

%% 计算 GOSPA (逐帧, 时间平均; 同时计入漏检与虚假航迹)
gos = eval_gospa(truth, tracks_history, params);
gos_mean = gos.mean;

%% 内部数据 (可选输出)
if nargout >= 4
    sim = struct();
    sim.tracks_history = tracks_history;
    sim.truth = truth;
    sim.params = params;
    sim.gos = gos;
    sim.assoc_history = assoc_history;
    sim.meas_src = meas_src;
end
end

% ========================================================================
function tracks = seed_tracks_from_truth(truth, params)
%SEED_TRACKS_FROM_TRUTH 用第 1 帧真值初始化航迹 (关联层受控实验)
%   配合 max_new_tracks=0 使用: 排除出生/鬼影, 只比较关联+滤波.

if isfield(params, 'assoc_seed_cov')
    P = params.assoc_seed_cov;
else
    P = params.init_cov;
end

tracks = {};
for tgt = 1:truth.num_targets
    x = zeros(7, 1);
    x(1:4) = truth.X{tgt}(1:4, 1);
    x(5) = truth.X{tgt}(5, 1);

    tr = struct();
    tr.x = x;
    tr.P = P;
    tr.id = tgt;
    tr.imm = imm_init(x, P, params);
    tr.hit_count = params.confirm_hits;
    tr.miss_count = 0;
    tr.confirmed = true;

    if params.use_features
        s = truth.spectrum{tgt};
        s = s(:) / norm(s);
        tr.spectrum = s;
        tr.feature_init = true;
    else
        tr.spectrum = NaN;
        tr.feature_init = false;
    end

    tracks{end + 1} = tr; %#ok<AGROW>
end
end

