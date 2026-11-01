function [x, P] = imm_combine(imm)
%IMM_COMBINE 由模型条件估计得到 IMM 组合估计 (全概率公式)
%   [x, P] = IMM_COMBINE(imm)
%
%   x = sum_j mu_j * x_j
%   P = sum_j mu_j * (P_j + (x_j-x)(x_j-x)')

mu = imm.mu(:);
mu = mu / sum(mu);
n_models = numel(mu);

x = zeros(imm.nx, 1);
for j = 1:n_models
    x = x + mu(j) * imm.x(:, j);
end

P = zeros(imm.nx);
for j = 1:n_models
    dx = imm.x(:, j) - x;
    P = P + mu(j) * (imm.P(:, :, j) + dx * dx');
end
P = (P + P') / 2;
end
