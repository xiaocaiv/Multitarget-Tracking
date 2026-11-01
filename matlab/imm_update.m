function imm = imm_update(imm, measurements, mp, beta, params)
%IMM_UPDATE IMM 量测更新 (模型条件 UKF 更新 + 模型概率更新)
%   imm = IMM_UPDATE(imm, measurements, mp, beta, params)
%
%   输入:
%     measurements - 当前平台量测 (结构体数组, 字段 .bearing)
%     mp           - imm_measurement_predict 的输出
%     beta         - JPDA 边缘关联概率 [beta0; beta_1; ...; beta_m]
%                    (由组合预测 z_comb/S_comb 计算, 各模型共用)
%
%   各模型用统一的 beta 加权新息做 UKF 更新, 并以该新息的似然更新模型概率
%   (IMM-JPDA 的常用近似, 避免对每个模型重复枚举联合事件).
%
%   模型概率更新:
%     mu_j = c_j * L_j / sum(c_j * L_j),  c = 预测模型概率 (.mu_prior)
%     L_j = N(nu_j; 0, S_j)
%   当 beta0 接近 1 (目标漏检/量测均落在波门外) 时, 不更新模型概率.

beta = beta(:);
beta0 = beta(1);
n_meas = numel(measurements);
n_models = numel(params.mm);

logL = zeros(n_models, 1);
use_likelihood = (beta0 < 0.9) && (n_meas > 0);

for j = 1:n_models
    % ---- 该模型的加权新息与散布 ----
    nu = 0;
    spread = 0;
    for m = 1:n_meas
        d = mod(measurements(m).bearing - mp.z(j) + pi, 2 * pi) - pi;
        nu = nu + beta(m + 1) * d;
        spread = spread + beta(m + 1) * d^2;
    end
    spread = max(spread - nu^2, 0);

    % ---- 模型条件 UKF 更新 ----
    K = mp.P_xz(:, j) / mp.S(j);

    imm.x(:, j) = imm.x(:, j) + K * nu;
    imm.P(:, :, j) = imm.P(:, :, j) ...
        - (1 - beta0) * (K * mp.S(j) * K') ...
        + K * spread * K';
    imm.P(:, :, j) = stabilize_cov(imm.P(:, :, j));

    % ---- 模型似然 ----
    if use_likelihood
        logL(j) = -0.5 * nu^2 / mp.S(j) - 0.5 * log(2 * pi * mp.S(j));
    end
end

% ---- 模型概率更新 ----
if use_likelihood
    logL = logL - max(logL);          % 防止下溢
    mu = imm.mu_prior .* exp(logL);
else
    mu = imm.mu_prior;
end
mu = real(max(mu, 1e-12));            % 数值保护
imm.mu = mu / sum(mu);
end
