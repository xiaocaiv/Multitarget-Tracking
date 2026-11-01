function imm = imm_init(x, P, params)
%IMM_INIT 由单个高斯估计初始化 IMM 结构
%   imm = IMM_INIT(x, P, params)
%
%   所有模型取相同初值, 模型概率由 params.imm_mu0 给出 (缺省等概率).
%
%   imm 结构体:
%     .nx        状态维数
%     .mu        [M x 1] 模型概率 (后验)
%     .mu_prior  [M x 1] 预测模型概率 (混合后、量测更新前)
%     .x         [nx x M] 模型条件状态
%     .P         [nx x nx x M] 模型条件协方差

n_models = numel(params.mm);
x = x(:);

mu0 = params.imm_mu0(:);
if numel(mu0) ~= n_models
    % 模型数与先验概率不匹配时退化为等概率
    mu0 = ones(n_models, 1) / n_models;
end
mu0 = mu0 / sum(mu0);

imm = struct();
imm.nx = numel(x);
imm.mu = mu0;
imm.mu_prior = mu0;
imm.x = repmat(x, 1, n_models);
imm.P = repmat(P, [1, 1, n_models]);
end
