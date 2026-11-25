function spec_lik = spectrum_likelihood(track, meas, params)
%SPECTRUM_LIKELIHOOD 频谱/声纹特征似然
%   spec_lik = SPECTRUM_LIKELIHOOD(track, meas, params)
%
%   两种实现，按数据结构自动选择：
%     (A) 新式（真实音频）：航迹存 tpl_line/tpl_band，量测带 line/band，
%         params 含 llr -> 用校准 LLR
%           L = exp( ell * q * scale * exp(-age/tau) )
%         ell 来自 feature_llr(...,'track')，q 为质量分，scale 为温度缩放，
%         age 为模板距上次更新的帧数（需要 params.k）。
%     (B) 旧式（96 维合成谱）：LOFAR/DEMON 分段 Pearson
%           L = exp(gamma * (max(0,rho_lofar)^2 + max(0,rho_demon)^2))
%
%   原为 jpda_associate.m 的局部函数, 现独立成文件, 供 JPDA 与
%   NN/GNN/PDA 关联器共用, 保证各方法的特征因子完全一致.
%
%   参见: jpda_associate.m, association_likelihood.m, feature_llr.m

%% ---------- (A) 新式：校准 LLR ----------
if isfield(track, 'tpl_line') && ~isempty(track.tpl_line) && ...
        isfield(meas, 'line') && ~isempty(meas.line) && isfield(params, 'llr')
    f1 = struct('line', track.tpl_line, 'band', track.tpl_band);
    f2 = struct('line', meas.line, 'band', meas.band);
    ell = feature_llr(f1, f2, params.llr, 'track');

    q = min(getf(track, 'q', 1), getf(meas, 'q', 1));
    age = 0;
    if isfield(params, 'k') && isfield(track, 'feat_step') && ~isempty(track.feat_step)
        age = max(params.k - track.feat_step, 0);
    end
    tau   = getf(params, 'feature_tau', 10);
    scale = getf(params, 'feature_scale', 0.7);

    spec_lik = exp(ell * q * scale * exp(-age / tau));
    return;
end

%% ---------- (B) 旧式：96 维谱 ----------
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

% ========================================================================
function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
