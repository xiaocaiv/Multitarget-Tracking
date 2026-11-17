function track = update_track_spectrum(track, measurements, beta, params)
%UPDATE_TRACK_SPECTRUM 利用关联权重更新航迹频谱模板
%   track = UPDATE_TRACK_SPECTRUM(track, measurements, beta, params)
%
%   按 beta 对量测频谱做加权平均, 再与旧模板指数平滑 (alpha=0.3).
%   原为 jpda_process_platform.m 的局部函数, 现独立成文件,
%   供 JPDA 与非 JPDA 关联流程共用.
%
%   参见: jpda_process_platform.m, assoc_process_platform.m

num_meas = length(measurements);
if num_meas == 0, return; end

spec_weighted = zeros(params.N_features, 1);
total_w = 0;
for j = 1:num_meas
    w = beta(j + 1);  % beta(1)=beta0, beta(2:end) 对应各量测
    if w > 0.01
        spec_weighted = spec_weighted + w * measurements(j).spectrum;
        total_w = total_w + w;
    end
end
if total_w < 0.01, return; end

spec_new = spec_weighted / total_w;
spec_new = spec_new / norm(spec_new);

if ~track.feature_init
    track.spectrum = spec_new;
    track.feature_init = true;
else
    alpha = 0.3;
    track.spectrum = (1 - alpha) * track.spectrum + alpha * spec_new;
    track.spectrum = track.spectrum / norm(track.spectrum);
end
end
