function [beta, beta0] = assoc_nn(z_preds, Ss, measurements, tracks, params)
%ASSOC_NN 局部最近邻 (Nearest Neighbor) 数据关联
%   [beta, beta0] = ASSOC_NN(z_preds, Ss, measurements, tracks, params)
%
%   每条航迹独立地在波门内取关联似然最大的量测, 输出 one-hot beta.
%   不同航迹允许选择同一量测 (这是 NN 与 GNN 的核心区别).
%   无门内量测时 beta0 = 1.
%
%   参见: assoc_gnn.m, association_likelihood.m, associate_measurements.m

n = numel(z_preds);
m = numel(measurements);

beta = cell(n, 1);
beta0 = ones(n, 1);

if n == 0
    return;
end

[~, ~, L, gated] = association_likelihood(z_preds, Ss, measurements, tracks, params);

for t = 1:n
    beta{t} = [1; zeros(m, 1)];
    js = find(gated(t, :));
    if isempty(js)
        continue;
    end

    % max 返回首个最大值索引 -> 平局时取最小量测索引, 结果确定
    [~, k] = max(L(t, js));
    j = js(k);

    beta{t} = zeros(m + 1, 1);
    beta{t}(j + 1) = 1;
    beta0(t) = 0;
end
end
