function track = update_track_spectrum(track, measurements, beta, params)
%UPDATE_TRACK_SPECTRUM 利用关联权重更新航迹声纹/频谱模板
%   track = UPDATE_TRACK_SPECTRUM(track, measurements, beta, params)
%
%   两种模板，按数据结构自动选择：
%     (A) 新式：tpl_line / tpl_band（真实音频链路），指数平滑 alpha=sig_alpha
%     (B) 旧式：spectrum（96 维合成谱），指数平滑 alpha=0.3
%
%   两者都会同步更新质量分 q 与特征时间戳 feat_step（供 spectrum_likelihood
%   做时间老化）。
%
%   原为 jpda_process_platform.m 的局部函数, 现独立成文件,
%   供 JPDA 与非 JPDA 关联流程共用.
%
%   参见: jpda_process_platform.m, assoc_process_platform.m, spectrum_likelihood.m

num_meas = length(measurements);
if num_meas == 0, return; end
alpha = getf(params, 'sig_alpha', 0.3);

%% ---------- (A) 新式：line / band ----------
if isfield(measurements(1), 'line') && ~isempty(measurements(1).line)
    L = zeros(size(measurements(1).line));
    B = zeros(size(measurements(1).band));
    tw = 0; qsum = 0;
    for j = 1:num_meas
        w = beta(j + 1);                     % beta(1)=beta0
        if w > 0.01 && ~isempty(measurements(j).line)
            L = L + w * measurements(j).line;
            B = B + w * measurements(j).band;
            qsum = qsum + w * getf(measurements(j), 'q', 1);
            tw = tw + w;
        end
    end
    if tw < 0.01, return; end
    L = L / tw; B = B / tw;
    L = L / (norm(L) + 1e-12);
    B = B / (norm(B) + 1e-12);

    if ~getf(track, 'feature_init', false) || isempty(getf(track, 'tpl_line', []))
        track.tpl_line = L;
        track.tpl_band = B;
        track.feature_init = true;
    else
        track.tpl_line = (1 - alpha) * track.tpl_line + alpha * L;
        track.tpl_line = track.tpl_line / (norm(track.tpl_line) + 1e-12);
        track.tpl_band = (1 - alpha) * track.tpl_band + alpha * B;
        track.tpl_band = track.tpl_band / (norm(track.tpl_band) + 1e-12);
    end
    qnew = qsum / tw;
    track.q = (1 - alpha) * getf(track, 'q', qnew) + alpha * qnew;
    track.feat_step = getf(params, 'k', NaN);
    return;
end

%% ---------- (B) 旧式：96 维 spectrum ----------
if ~isfield(measurements(1), 'spectrum') || isempty(measurements(1).spectrum)
    return;
end

spec_weighted = zeros(params.N_features, 1);
total_w = 0;
for j = 1:num_meas
    w = beta(j + 1);
    if w > 0.01
        spec_weighted = spec_weighted + w * measurements(j).spectrum;
        total_w = total_w + w;
    end
end
if total_w < 0.01, return; end

spec_new = spec_weighted / total_w;
spec_new = spec_new / norm(spec_new);

if ~getf(track, 'feature_init', false)
    track.spectrum = spec_new;
    track.feature_init = true;
else
    track.spectrum = (1 - 0.3) * track.spectrum + 0.3 * spec_new;
    track.spectrum = track.spectrum / norm(track.spectrum);
end
track.feat_step = getf(params, 'k', NaN);
end

% ========================================================================
function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
