function a = eval_assoc(truth, tracks_history, meas_src, assoc_history, params, verbose)
%EVAL_ASSOC 关联层评价: 正确/错误/杂波关联质量、硬判决正确率、ID 切换
%   a = EVAL_ASSOC(truth, tracks_history, meas_src, assoc_history, params)
%   a = EVAL_ASSOC(..., verbose)  verbose=false 时不打印汇总
%
%   输入:
%     truth           generate_truth 输出
%     tracks_history  K×1 cell, 每帧航迹 cell 数组
%     meas_src        K×2 cell, 每帧每平台量测来源 (0=杂波, k=目标 k)
%     assoc_history   K×2 cell, 每帧每平台关联记录 (info 结构体)
%     params          使用 .gospa_c 作为目标-航迹配对阈值 (默认 100 m)
%
%   输出 a 结构体:
%     .n_opp        参与统计的 track-frame 数 (目标在该平台被检测到)
%     .n_hard       硬判决次数 (beta0 < 0.5)
%     .p_correct    正确目标量测的期望关联概率质量
%     .p_correct_norm 正确质量 / (正确+错误+杂波) (剔除 beta0, 便于软硬方法对比)
%     .p_wrong      错误目标量测的期望关联概率质量
%     .p_clutter    杂波量测的期望关联概率质量
%     .p_miss       beta0 均值 (无关联质量)
%     .hard_correct 硬判决正确率
%     .id_switches  目标-航迹 ID 切换次数 (连续帧都匹配时才计)
%     .id_per_target 各目标 ID 切换次数
%     .match_dist   配对阈值
%
%   说明: 目标-航迹配对每帧用匈牙利算法 (距离阈值内) 完成; 交叉场景可能歧义,
%   因此 GOSPA 仍是主指标, 本函数用于解释关联层差异来源.
%
%   参见: compare_assoc.m, hungarian_square.m

if nargin < 6 || isempty(verbose)
    verbose = true;
end

K = numel(tracks_history);
n_tgt = truth.num_targets;

match_dist = 100;
if isfield(params, 'gospa_c') && ~isempty(params.gospa_c)
    match_dist = params.gospa_c;
end

n_opp = 0;
n_hard = 0;
n_hard_correct = 0;
sum_correct = 0;
sum_wrong = 0;
sum_clutter = 0;
sum_beta0 = 0;

id_hist = NaN(K, n_tgt);

for k = 1:K
    trk_all = tracks_history{k};
    if isempty(trk_all)
        continue;
    end
    keep = cellfun(@(t) t.confirmed, trk_all);
    trk = trk_all(keep);
    if isempty(trk)
        continue;
    end

    X = zeros(2, n_tgt);
    for tgt = 1:n_tgt
        X(:, tgt) = truth.pos{tgt}(:, k);
    end
    Y = zeros(2, numel(trk));
    for t = 1:numel(trk)
        Y(:, t) = trk{t}.x(1:2);
    end

    map = match_targets_to_tracks(X, Y, match_dist);
    for tgt = 1:n_tgt
        if map(tgt) > 0
            id_hist(k, tgt) = trk{map(tgt)}.id;
        end
    end

    for p = 1:2
        if isempty(assoc_history) || size(assoc_history, 2) < p
            continue;
        end
        info = assoc_history{k, p};
        if isempty(info) || ~isfield(info, 'beta')
            continue;
        end

        src = meas_src{k, p};
        if isempty(src)
            continue;
        end

        for tgt = 1:n_tgt
            ti = map(tgt);
            if ti == 0
                continue;
            end

            trk_id = trk{ti}.id;
            bi = find(info.track_ids == trk_id, 1);
            if isempty(bi)
                continue;
            end

            b = info.beta{bi};
            bs = b(2:end);
            if numel(bs) ~= numel(src)
                continue;
            end

            % 目标在该平台未被检测到: 不存在关联机会, 跳过
            if ~any(src == tgt)
                continue;
            end

            n_opp = n_opp + 1;
            sum_correct = sum_correct + sum(bs(src == tgt));
            sum_wrong = sum_wrong + sum(bs(src > 0 & src ~= tgt));
            sum_clutter = sum_clutter + sum(bs(src == 0));
            sum_beta0 = sum_beta0 + b(1);

            if b(1) < 0.5
                n_hard = n_hard + 1;
                [~, jmax] = max(bs);
                if src(jmax) == tgt
                    n_hard_correct = n_hard_correct + 1;
                end
            end
        end
    end
end

a = struct();
a.n_opp = n_opp;
a.n_hard = n_hard;
a.match_dist = match_dist;

if n_opp > 0
    a.p_correct = sum_correct / n_opp;
    a.p_wrong = sum_wrong / n_opp;
    a.p_clutter = sum_clutter / n_opp;
    a.p_miss = sum_beta0 / n_opp;
    assoc_mass = sum_correct + sum_wrong + sum_clutter;
    if assoc_mass > 0
        a.p_correct_norm = sum_correct / assoc_mass;
    else
        a.p_correct_norm = NaN;
    end
else
    a.p_correct = NaN;
    a.p_correct_norm = NaN;
    a.p_wrong = NaN;
    a.p_clutter = NaN;
    a.p_miss = NaN;
end

if n_hard > 0
    a.hard_correct = n_hard_correct / n_hard;
else
    a.hard_correct = NaN;
end

a.id_per_target = zeros(1, n_tgt);
for tgt = 1:n_tgt
    ids = id_hist(:, tgt);
    for k = 2:K
        if ~isnan(ids(k)) && ~isnan(ids(k - 1)) && ids(k) ~= ids(k - 1)
            a.id_per_target(tgt) = a.id_per_target(tgt) + 1;
        end
    end
end
a.id_switches = sum(a.id_per_target);

if verbose
    fprintf(['  关联层: 机会 %d | 正确 %.3f 错误 %.3f 杂波 %.3f beta0 %.3f | ' ...
        '硬判决 %.3f | ID 切换 %d\n'], ...
        a.n_opp, a.p_correct, a.p_wrong, a.p_clutter, a.p_miss, ...
        a.hard_correct, a.id_switches);
end
end

% ========================================================================
function map = match_targets_to_tracks(X, Y, dist)
% 目标->航迹的匈牙利配对 (距离阈值内); 未匹配返回 0
n = size(X, 2);
m = size(Y, 2);
map = zeros(1, n);
if n == 0 || m == 0
    return;
end

BIG = 1e6;
A = BIG * ones(n + m, n + m);
for i = 1:n
    for j = 1:m
        d = norm(X(:, i) - Y(:, j));
        if d < dist
            A(i, j) = d;
        end
    end
end
A(1:n, m+1:m+n) = dist / 2;
A(n+1:n+m, 1:m) = dist / 2;
A(n+1:n+m, m+1:m+n) = 0;

assign = hungarian_square(A);
for i = 1:n
    j = assign(i);
    if j <= m && A(i, j) < BIG / 2
        map(i) = j;
    end
end
end
