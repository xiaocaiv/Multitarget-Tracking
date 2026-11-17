function [v, S, L, gated] = association_likelihood(z_preds, Ss, measurements, tracks, params)
%ASSOCIATION_LIKELIHOOD 公共波门与关联似然 (方位角 + 可选频谱特征)
%   [v, S, L, gated] = ASSOCIATION_LIKELIHOOD(z_preds, Ss, measurements, tracks, params)
%
%   输入:
%     z_preds      各航迹预测量测 (元胞, 标量)
%     Ss           各航迹新息方差 (元胞, 标量)
%     measurements 当前平台量测 (结构体数组: .bearing, .spectrum, .src)
%     tracks       航迹元胞数组 (用于特征似然: .spectrum / .feature_init)
%     params       参数结构体 (.gate_gamma, .use_features, ...)
%
%   输出:
%     v      [n x m] 环绕后的新息 (rad)
%     S      [n x 1] 新息方差 (rad^2)
%     L      [n x m] 关联似然 (门外为 0; 含特征因子)
%     gated  [n x m] 波门逻辑矩阵
%
%   运动学部分与 jpda_associate.m 完全一致, 特征部分调用
%   spectrum_likelihood.m, 保证 NN/GNN/PDA/JPDA 的可比性.
%
%   参见: assoc_nn.m, assoc_gnn.m, assoc_pda.m, jpda_associate.m

n = numel(z_preds);
m = numel(measurements);

v = zeros(n, m);
S = zeros(n, 1);
L = zeros(n, m);
gated = false(n, m);

if n == 0 || m == 0
    return;
end

use_feat = isfield(params, 'use_features') && params.use_features;

for t = 1:n
    S(t) = max(Ss{t}, realmin);
    for j = 1:m
        d = mod(measurements(j).bearing - z_preds{t} + pi, 2 * pi) - pi;
        v(t, j) = d;

        if d^2 / S(t) < params.gate_gamma
            gated(t, j) = true;
            L(t, j) = exp(-0.5 * d^2 / S(t)) / sqrt(2 * pi * S(t));

            if use_feat && isfield(tracks{t}, 'feature_init') && tracks{t}.feature_init
                L(t, j) = L(t, j) * spectrum_likelihood(tracks{t}, measurements(j), params);
            end
        end
    end
end
end
