function mp = imm_measurement_predict(imm, plat_pos, params)
%IMM_MEASUREMENT_PREDICT IMM 各模型的量测预测
%   mp = IMM_MEASUREMENT_PREDICT(imm, plat_pos, params)
%
%   输出 mp 结构体:
%     .z      [1 x M]  各模型预测量测 (方位角, rad)
%     .S      [1 x M]  各模型新息协方差 (rad^2)
%     .P_xz   [nx x M] 各模型状态-量测互协方差
%     .z_comb          模型概率加权的组合预测量测 (供 JPDA 波门/关联使用)
%     .S_comb          组合新息协方差 (含模型间离散项)
%
%   说明: JPDA 的波门与边缘关联概率用组合预测 (z_comb, S_comb) 计算一次,
%   各模型再用同一组 beta 做模型条件更新, 这是 IMM-JPDA 的常用近似.

n_models = numel(params.mm);
mu = imm.mu(:)';
mu = mu / sum(mu);

mp = struct();
mp.z = zeros(1, n_models);
mp.S = zeros(1, n_models);
mp.P_xz = zeros(imm.nx, n_models);

for j = 1:n_models
    [mp.z(j), mp.S(j), mp.P_xz(:, j)] = ukf_measurement_predict(...
        imm.x(:, j), imm.P(:, :, j), params.R, plat_pos, ...
        params.alpha, params.beta, params.kappa);

    % 数值保护: 理论上 S = H P H' + R >= R, 出现更小值说明 UKF 权重
    % (alpha=1e-3 时中心点权重很大) 引起数值误差, 这里以 R 为下界
    mp.S(j) = max(mp.S(j), params.R);
end

% ---- 组合预测 (角度需环绕处理) ----
z0 = sum(mu .* mp.z);
z0 = mod(z0 + pi, 2 * pi) - pi;

dz = mod(mp.z - z0 + pi, 2 * pi) - pi;
S0 = sum(mu .* (mp.S + dz.^2));

mp.z_comb = z0;
mp.S_comb = max(S0, params.R);
end
