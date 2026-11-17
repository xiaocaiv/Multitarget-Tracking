function spec_lik = spectrum_likelihood(track, meas, params)
%SPECTRUM_LIKELIHOOD 频谱特征似然 (LOFAR + DEMON 分段 Pearson)
%   spec_lik = SPECTRUM_LIKELIHOOD(track, meas, params)
%
%   原为 jpda_associate.m 的局部函数, 现独立成文件, 供 JPDA 与
%   NN/GNN/PDA 关联器共用, 保证各方法的特征因子完全一致.
%
%   分段计算 LOFAR/DEMON 的 Pearson 相关, 联合似然:
%     L_feature = exp(gamma * (max(0,rho_lofar)^2 + max(0,rho_demon)^2))
%   避免 LOFAR 的高相关性淹没 DEMON 的区分能力.
%
%   参见: jpda_associate.m, association_likelihood.m

if ~track.feature_init || isempty(meas.spectrum) || isempty(track.spectrum)
    spec_lik = 1.0;
    return;
end

n_lofar = params.N_lofar;
lofar_idx = 1:n_lofar;
demon_idx = n_lofar + 1 : params.N_features;

% LOFAR 段 Pearson
xl = track.spectrum(lofar_idx); xl = xl - mean(xl);
yl = meas.spectrum(lofar_idx);  yl = yl - mean(yl);
rho_lofar = (xl' * yl) / (norm(xl) * norm(yl) + 1e-10);

% DEMON 段 Pearson
xd = track.spectrum(demon_idx); xd = xd - mean(xd);
yd = meas.spectrum(demon_idx);  yd = yd - mean(yd);
rho_demon = (xd' * yd) / (norm(xd) * norm(yd) + 1e-10);

spec_lik = exp(params.feature_gamma * (max(0, rho_lofar)^2 + max(0, rho_demon)^2));
end
