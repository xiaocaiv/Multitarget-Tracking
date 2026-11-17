function [beta, beta0] = assoc_gnn(z_preds, Ss, measurements, tracks, params)
%ASSOC_GNN 全局最近邻 (Global Nearest Neighbor) 数据关联
%   [beta, beta0] = ASSOC_GNN(z_preds, Ss, measurements, tracks, params)
%
%   以负对数似然为代价, 构造 (n+m)x(n+m) 的"航迹 + 虚拟行 / 量测 + 虚拟列"
%   方阵, 用匈牙利算法求全局最小代价的互斥分配:
%     c_ij      = -ln(Pd * L_ij)        实配对
%     c_miss    = -ln(1 - Pd)           航迹漏检
%     c_fa      = -ln(lambda_c / 2/pi)  量测为杂波
%   虚拟行列保证完美匹配总存在; 实配对代价 >= BIG/2 视为未关联.
%
%   输出 one-hot beta; 同一量测至多被一条航迹占用.
%
%   参见: assoc_nn.m, hungarian_square.m, association_likelihood.m

n = numel(z_preds);
m = numel(measurements);

beta = cell(n, 1);
beta0 = ones(n, 1);

if n == 0
    return;
end

for t = 1:n
    beta{t} = [1; zeros(m, 1)];
end
if m == 0
    return;
end

[~, ~, L, gated] = association_likelihood(z_preds, Ss, measurements, tracks, params);

BIG = 1e9;

% 实配对代价
c = BIG * ones(n, m);
c(gated) = -log(max(params.Pd * L(gated), 1e-300));

% 漏检 / 杂波代价 (Pd=1 或 lambda_c=0 时用有限大值截断)
c_miss = -log(max(1 - params.Pd, 1e-12));
c_fa = -log(max(params.lambda_c / (2 * pi), 1e-12));

A = zeros(n + m, n + m);
A(1:n, 1:m) = c;
A(1:n, m+1:m+n) = c_miss;
A(n+1:n+m, 1:m) = c_fa;
A(n+1:n+m, m+1:m+n) = 0;

assign = hungarian_square(A);

for t = 1:n
    j = assign(t);
    if j <= m && c(t, j) < BIG / 2
        beta{t} = zeros(m + 1, 1);
        beta{t}(j + 1) = 1;
        beta0(t) = 0;
    end
end
end
