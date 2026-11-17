function [beta, beta0] = assoc_pda(z_preds, Ss, measurements, tracks, params)
%ASSOC_PDA 逐航迹概率数据关联 (Probabilistic Data Association)
%   [beta, beta0] = ASSOC_PDA(z_preds, Ss, measurements, tracks, params)
%
%   对每条航迹独立计算门内量测的后验关联概率:
%     b0 = (1 - Pd) * lambda_c / (2*pi)
%     e_j = Pd * L_ij
%     beta0 = b0 / (b0 + sum e_j),  beta_j = e_j / (b0 + sum e_j)
%
%   该式是现有 JPDA 事件模型在单航迹时的精确边缘化 (杂波因子 lambda^m 相消),
%   因此单航迹下 PDA 与 JPDA 的 beta 应完全一致.
%   多航迹时各航迹独立计算, 允许同一量测被多条航迹复用.
%
%   参见: jpda_associate.m, assoc_gnn.m, association_likelihood.m

n = numel(z_preds);
m = numel(measurements);

beta = cell(n, 1);
beta0 = ones(n, 1);

if n == 0
    return;
end

[~, ~, L, ~] = association_likelihood(z_preds, Ss, measurements, tracks, params);

lambda_spatial = params.lambda_c / (2 * pi);
b0 = (1 - params.Pd) * lambda_spatial;

for t = 1:n
    e = params.Pd * L(t, :);
    se = sum(e);

    if se > 0
        den = b0 + se;
        beta{t} = [b0 / den; (e / den)'];
        beta0(t) = b0 / den;
    else
        beta{t} = [1; zeros(m, 1)];
        beta0(t) = 1;
    end
end
end
