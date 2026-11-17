function [tracks, unassociated, info] = assoc_process_platform(tracks, measurements, plat_pos, params)
%ASSOC_PROCESS_PLATFORM 非 JPDA 关联算法的单平台完整流程
%   [tracks, unassociated, info] = ASSOC_PROCESS_PLATFORM(tracks, measurements, plat_pos, params)
%
%   与 jpda_process_platform.m 同构:
%     1. IMM 各模型量测预测 (组合值用于关联, 模型值用于条件更新)
%     2. associate_measurements 产生 beta (nn / gnn / pda / jpda)
%     3. 各模型用同一 beta 做模型条件 UKF 更新
%     4. 由模型概率加权得到组合状态
%     5. 频谱模板更新 (仅特征辅助时)
%
%   参见: associate_measurements.m, process_platform.m, imm_update.m

num_tracks = length(tracks);
num_meas = length(measurements);
unassociated = true(1, max(num_meas, 0));

info = struct();
info.method = 'jpda';
if isfield(params, 'assoc_method') && ~isempty(params.assoc_method)
    info.method = lower(char(params.assoc_method));
end
info.beta = cell(num_tracks, 1);
info.beta0 = ones(num_tracks, 1);
info.track_ids = reshape(cellfun(@(t) t.id, tracks), 1, []);
info.t_assoc = 0;
info.n_events = 0;
info.fallback = false;
info.unassoc = unassociated;

if num_tracks == 0
    return;
end

% ---------- 第1步: 各航迹 IMM 量测预测 ----------
mps = cell(num_tracks, 1);
z_preds = cell(num_tracks, 1);
Ss = cell(num_tracks, 1);

for t = 1:num_tracks
    mps{t} = imm_measurement_predict(tracks{t}.imm, plat_pos, params);
    z_preds{t} = mps{t}.z_comb;
    Ss{t} = mps{t}.S_comb;
end

% ---------- 第2步: 数据关联 ----------
[beta, unassociated, info] = associate_measurements( ...
    z_preds, Ss, measurements, tracks, params);

% ---------- 第3步: 模型条件更新 + 组合 ----------
for t = 1:num_tracks
    if num_meas == 0
        tracks{t}.beta0_last = 1.0;
        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);
        continue;
    end

    tracks{t}.imm = imm_update(tracks{t}.imm, measurements, mps{t}, beta{t}, params);
    [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);

    tracks{t}.beta0_last = info.beta0(t);

    if params.use_features
        tracks{t} = update_track_spectrum(tracks{t}, measurements, beta{t}, params);
    end
end
end
