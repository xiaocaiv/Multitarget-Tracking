function imm = imm_predict(imm, params)
%IMM_PREDICT IMM 一步预测: 模型混合 + 各模型时间更新
%   imm = IMM_PREDICT(imm, params)
%
%   1) 用马尔可夫转移矩阵 Pi 计算预测模型概率 c_j = sum_i Pi(i,j) mu_i
%   2) 计算混合概率 mu_{i|j} = Pi(i,j) mu_i / c_j
%   3) 混合初始条件: x0_j = sum_i mu_{i|j} x_i
%                     P0_j = sum_i mu_{i|j} [P_i + (x_i-x0_j)(x_i-x0_j)']
%   4) 各模型按自身状态转移函数做 UKF 时间更新 (CT 为非线性, 用 sigma 点)
%
%   返回的 imm 中 .mu 与 .mu_prior 均为预测模型概率 c.

n_models = numel(params.mm);
nx = imm.nx;

Pi = params.imm_transition;
if ~isequal(size(Pi), [n_models, n_models])
    Pi = eye(n_models);           % 转移矩阵与模型数不匹配时退化为不切换
end

mu = imm.mu(:);
mu = mu / sum(mu);

c = Pi' * mu;                     % 预测模型概率
c = max(c, 1e-12);
c = c / sum(c);

x_new = zeros(nx, n_models);
P_new = zeros(nx, nx, n_models);

for j = 1:n_models
    % ---- 模型混合 ----
    w = Pi(:, j) .* mu / c(j);

    x0 = zeros(nx, 1);
    for i = 1:n_models
        x0 = x0 + w(i) * imm.x(:, i);
    end

    P0 = zeros(nx);
    for i = 1:n_models
        dx = imm.x(:, i) - x0;
        P0 = P0 + w(i) * (imm.P(:, :, i) + dx * dx');
    end
    P0 = (P0 + P0') / 2;

    % ---- 模型条件预测 ----
    [x_new(:, j), P_new(:, :, j)] = predict_model(x0, P0, params.mm(j), params);
end

imm.x = x_new;
imm.P = P_new;
imm.mu = c;
imm.mu_prior = c;
end


% ========================================================================
function [x_pred, P_pred] = predict_model(x, P, model, params)
%PREDICT_MODEL 单个模型的 UKF 时间更新
%   线性模型直接用 F 传播 sigma 点; 非线性模型 (CT) 逐点调用 f.

if ~isreal(x) || ~isreal(P) || any(~isfinite(x)) || any(~isfinite(P(:)))
    error('imm_predict:invalidState', ...
        '模型 %s 的混合状态非实数/非有限 (x real=%d, P real=%d)', ...
        model.name, isreal(x), isreal(P));
end

[X, Wm, Wc] = generate_sigma_points(x, P, params.alpha, params.beta, params.kappa);
n = numel(x);
n_sigma = size(X, 2);

if ~isempty(model.F)
    X_pred = model.F * X;
else
    X_pred = zeros(n, n_sigma);
    for i = 1:n_sigma
        X_pred(:, i) = model.f(X(:, i), params.dt);
    end
end

x_pred = zeros(n, 1);
for i = 1:n_sigma
    x_pred = x_pred + Wm(i) * X_pred(:, i);
end

P_pred = zeros(n);
for i = 1:n_sigma
    d = X_pred(:, i) - x_pred;
    P_pred = P_pred + Wc(i) * (d * d');
end

P_pred = P_pred + model.Q;
P_pred = stabilize_cov(P_pred);
end
