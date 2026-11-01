function gos = eval_gospa(truth, tracks_history, params)
%EVAL_GOSPA 由航迹历史逐帧计算 GOSPA 指标及其分解
%
%   gos = EVAL_GOSPA(truth, tracks_history, params)
%
%   评价集合:
%     真值集合 X_k  = 全部真实目标位置 (truth.pos)
%     估计集合 Y_k  = 已确认航迹的位置 (默认; 由 params.gospa_use_confirmed
%                     控制, 置 false 时使用全部航迹, 含未确认的暂定航迹)
%
%   输入:
%     truth            - generate_truth 的输出
%     tracks_history   - K×1 cell, 每帧的航迹 cell 数组
%     params           - 参数结构体, 使用以下字段:
%                          gospa_p            距离阶数 (默认 2)
%                          gospa_c            截断参数 (默认 100 m)
%                          gospa_alpha        基数惩罚因子 (默认 2)
%                          gospa_use_confirmed 仅统计确认航迹 (默认 true)
%
%   输出 gos 结构体 (各字段均为 1×K 行向量, 除汇总标量外):
%     .time        时间轴 (s)
%     .gospa       逐帧 GOSPA (m)
%     .loc         定位误差分量 (m)
%     .missed      漏检目标分量 (m)
%     .false       虚假航迹分量 (m)
%     .n_truth     真值目标数
%     .n_est       参与评价的航迹数
%     .n_matched   成功关联对数
%     .n_missed    漏检目标数
%     .n_false     虚假航迹数
%     .mean        平均 GOSPA (m)
%     .rms         GOSPA 均方根 (m), 对大误差更敏感
%     .max         最大 GOSPA (m)
%     .mean_loc / .mean_missed / .mean_false   各分量平均值 (m)
%     .mean_n_false / .mean_n_missed           平均虚假航迹数 / 漏检数
%     .cfg         本次评价使用的 (p, c, alpha, use_confirmed)

K = params.K;

% ---------- 评价参数 (缺省时使用推荐值) ----------
p     = get_param(params, 'gospa_p', 2);
c     = get_param(params, 'gospa_c', 100);
alpha = get_param(params, 'gospa_alpha', 2);
use_confirmed = get_param(params, 'gospa_use_confirmed', true);

n_tgt = truth.num_targets;

gos = struct();
gos.time      = (1:K) * params.dt;
gos.gospa     = NaN(1, K);
gos.loc       = NaN(1, K);
gos.missed    = NaN(1, K);
gos.false     = NaN(1, K);
gos.n_truth   = zeros(1, K);
gos.n_est     = zeros(1, K);
gos.n_matched = zeros(1, K);
gos.n_missed  = zeros(1, K);
gos.n_false   = zeros(1, K);

for k = 1:K
    % 真值集合
    X = zeros(2, n_tgt);
    for tgt = 1:n_tgt
        X(:, tgt) = truth.pos{tgt}(:, k);
    end

    % 估计集合: 确认航迹 (或全部航迹)
    trk = tracks_history{k};
    if isempty(trk)
        keep = [];
    elseif use_confirmed
        keep = cellfun(@(t) t.confirmed, trk);
    else
        keep = true(1, numel(trk));
    end
    trk = trk(keep);

    Y = zeros(2, numel(trk));
    for t = 1:numel(trk)
        Y(:, t) = trk{t}.x(1:2);
    end

    [d, comp] = gospar(X, Y, p, c, alpha);

    gos.gospa(k)     = d;
    gos.loc(k)       = comp.loc;
    gos.missed(k)    = comp.missed;
    gos.false(k)     = comp.false;
    gos.n_truth(k)   = n_tgt;
    gos.n_est(k)     = numel(trk);
    gos.n_matched(k) = comp.n_matched;
    gos.n_missed(k)  = comp.n_missed;
    gos.n_false(k)   = comp.n_false;
end

% ---------- 汇总统计 ----------
valid = ~isnan(gos.gospa);
gos.mean      = mean(gos.gospa(valid));
gos.rms       = sqrt(mean(gos.gospa(valid).^2));
gos.max       = max(gos.gospa(valid));
gos.mean_loc  = mean(gos.loc(valid));
gos.mean_missed = mean(gos.missed(valid));
gos.mean_false  = mean(gos.false(valid));
gos.mean_n_false  = mean(gos.n_false(valid));
gos.mean_n_missed = mean(gos.n_missed(valid));
gos.cfg = struct('p', p, 'c', c, 'alpha', alpha, ...
    'use_confirmed', use_confirmed);
end


function val = get_param(params, name, default)
%GET_PARAM 读取参数, 不存在或为空时返回默认值
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default;
end
end
