function res = run_case(params, truth, bank)
%RUN_CASE 执行一次闭环仿真（无绘图），返回评价指标
%
%   res = run_case(params, truth, bank)
%
%   与 main.m 的计算流程一致：IMM 预测 → 双站 process_platform → 逐模型 CI 融合
%   → hit/miss 合并 → 声纹模板合并 → manage_tracks。区别只是不出图、返回指标。
%
%   res 字段
%     .rmse        1×num_targets 位置 RMSE（收敛段，最近邻匹配）
%     .missed      平均漏检目标数/帧
%     .nfalse      平均虚假航迹数/帧
%     .gos         eval_gospa 的返回结构
%     .tracks_history

if nargin < 3, bank = []; end
if ~isempty(bank), params.audio_bank = bank; end
if ~isfield(params, 'k'), params.k = 1; end

K = params.K;

% 移动平台: 初始化航向/位置历史（与 main.m 相同）
if params.mobile_platforms && ~isfield(params, 'plat1_pos_hist')
    params.plat1_heading_hist = zeros(1, K);
    params.plat2_heading_hist = zeros(1, K);
    params.plat1_pos_hist = zeros(2, K);
    params.plat2_pos_hist = zeros(2, K);
    params.plat1_heading_hist(1) = params.plat1_init_heading;
    params.plat2_heading_hist(1) = params.plat2_init_heading;
    params.plat1_pos_hist(:, 1) = params.plat1_pos;
    params.plat2_pos_hist(:, 1) = params.plat2_pos;
end

tracks = {};
next_id = 1;
tracks_history = cell(K, 1);

for k = 1:K
    params.k = k;

    if params.mobile_platforms
        params.plat1_pos = params.plat1_pos_hist(:, k);
        params.plat2_pos = params.plat2_pos_hist(:, k);
    end

    % 1) 量测
    [z1, z2] = generate_measurements(truth, params, k);

    % 2) IMM 预测
    num_trk = numel(tracks);
    pred_imm = cell(num_trk, 1);
    for t = 1:num_trk
        tracks{t}.imm = imm_predict(tracks{t}.imm, params);
        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);
        pred_imm{t} = tracks{t}.imm;
    end

    % 3) 深拷贝
    tracks2 = cell(num_trk, 1);
    for t = 1:num_trk, tracks2{t} = tracks{t}; end

    % 4) 双站并行关联 + 更新
    [tracks,  unassoc1] = process_platform(tracks,  z1, params.plat1_pos, params);
    [tracks2, unassoc2] = process_platform(tracks2, z2, params.plat2_pos, params);

    % 5) 逐模型 CI 融合 + hit/miss + 模板合并
    n_models = numel(params.mm);
    for t = 1:num_trk
        mu1 = tracks{t}.imm.mu;
        mu2 = tracks2{t}.imm.mu;
        mu = sqrt(max(mu1, 1e-12) .* max(mu2, 1e-12));
        tracks{t}.imm.mu = mu / sum(mu);

        for j = 1:n_models
            [x_f, P_f] = ci_fuse( ...
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
            if isfield(tracks{t}, 'tpl_line') && ~isempty(tracks{t}.tpl_line) ...
                    && isfield(tracks2{t}, 'tpl_line') && ~isempty(tracks2{t}.tpl_line)
                w1 = max(tracks{t}.q, 1e-3);
                w2 = max(tracks2{t}.q, 1e-3);
                tl = w1 * tracks{t}.tpl_line + w2 * tracks2{t}.tpl_line;
                tb = w1 * tracks{t}.tpl_band + w2 * tracks2{t}.tpl_band;
                tracks{t}.tpl_line = tl / (norm(tl) + 1e-12);
                tracks{t}.tpl_band = tb / (norm(tb) + 1e-12);
            elseif ~isempty(tracks{t}.spectrum) && ~isempty(tracks2{t}.spectrum)
                spec_avg = (tracks{t}.spectrum + tracks2{t}.spectrum) / 2;
                tracks{t}.spectrum = spec_avg / norm(spec_avg);
            end
        end

        tracks{t} = rmfield(tracks{t}, 'beta0_last');
    end

    % 6) 航迹管理
    [tracks, next_id] = manage_tracks(tracks, z1, z2, unassoc1, unassoc2, params, next_id);

    % 7) 移动平台: 下一步位置
    if params.mobile_platforms && k < K
        h1 = params.plat1_heading_hist(k);
        [h1n, p1n] = platform_ecm_optimize(tracks, params.plat1_pos, h1, ...
            params.plat1_speed, params);
        params.plat1_heading_hist(k + 1) = h1n;
        params.plat1_pos_hist(:, k + 1) = p1n;

        h2 = params.plat2_heading_hist(k);
        [h2n, p2n] = platform_ecm_optimize(tracks, params.plat2_pos, h2, ...
            params.plat2_speed, params);
        params.plat2_heading_hist(k + 1) = h2n;
        params.plat2_pos_hist(:, k + 1) = p2n;
    end

    tracks_history{k} = tracks;
end

% ---- 评价 ----
[rmse, missed, nfalse] = eval_nn_rmse(tracks_history, truth, params);
res = struct();
res.rmse = rmse;
res.missed = missed;
res.nfalse = nfalse;
res.gos = eval_gospa(truth, tracks_history, params);
res.tracks_history = tracks_history;
end

% ========================================================================
function [rmse, missed, nfalse] = eval_nn_rmse(tracks_history, truth, params)
% 逐帧最近邻匹配：每个真值目标匹配 200 m 内最近的已确认航迹
K = numel(tracks_history);
conv = getf(params, 'rmse_conv_step', 20);
err2 = cell(1, truth.num_targets);
missed = 0; nfalse = 0; frames = 0;

for k = conv:K
    X = zeros(2, truth.num_targets);
    for t = 1:truth.num_targets, X(:, t) = truth.pos{t}(:, k); end
    Y = zeros(2, 0);
    for t = 1:numel(tracks_history{k})
        if tracks_history{k}{t}.confirmed
            Y = [Y, tracks_history{k}{t}.x(1:2)]; %#ok<AGROW>
        end
    end
    frames = frames + 1;
    used = false(size(Y, 2), 1);
    for i = 1:truth.num_targets
        if isempty(Y)
            missed = missed + 1;
            continue;
        end
        d = sqrt(sum((Y - X(:, i)) .^ 2, 1));
        [dmin, j] = min(d);
        if dmin < 200
            err2{i}(end + 1) = dmin ^ 2; %#ok<AGROW>
            used(j) = true;
        else
            missed = missed + 1;
        end
    end
    nfalse = nfalse + sum(~used);
end

rmse = zeros(1, truth.num_targets);
for i = 1:truth.num_targets
    if isempty(err2{i})
        rmse(i) = NaN;
    else
        rmse(i) = sqrt(mean(err2{i}));
    end
end
missed = missed / max(frames, 1);
nfalse = nfalse / max(frames, 1);
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
