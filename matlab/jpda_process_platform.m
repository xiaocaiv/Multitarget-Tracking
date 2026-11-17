function [tracks, unassociated, info] = jpda_process_platform(tracks, measurements, plat_pos, params)
%JPDA_PROCESS_PLATFORM 对单个平台执行完整 JPDA + IMM-UKF 更新流程
%   [tracks, unassociated] = JPDA_PROCESS_PLATFORM(tracks, measurements, plat_pos, params)
%   [tracks, unassociated, info] = JPDA_PROCESS_PLATFORM(...)
%
%   unassociated: 未关联量测索引 (逻辑数组)
%   info:         关联层记录 (beta / beta0 / track_ids / 耗时 / 事件数),
%                 供 run_simulation 汇总与 eval_assoc 评价
%
%   对每个航迹:
%     1. IMM 各模型量测预测 (组合值用于波门/关联, 模型值用于条件更新)
%     2. 用组合预测运行 JPDA 获取边缘关联概率 beta
%     3. 各模型用同一 beta 做模型条件 UKF 更新, 并按新息似然更新模型概率
%     4. 由模型概率加权得到组合状态 (track.x / track.P)
%     5. 保存 beta0_last 供调用方做融合后的 hit/miss 合并

num_tracks = length(tracks);
num_meas = length(measurements);
unassociated = true(1, max(num_meas, 0));

info = struct();
info.method = 'jpda';
info.beta = cell(num_tracks, 1);
info.beta0 = ones(num_tracks, 1);
info.track_ids = reshape(cellfun(@(t) t.id, tracks), 1, []);
info.t_assoc = 0;
info.n_events = 0;
info.fallback = false;

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

% ---------- 第2步: JPDA 关联 (组合预测) ----------
tic;
[beta, ~, ~, beta0, unassociated, jinfo] = jpda_associate(...
    z_preds, Ss, cell(num_tracks, 1), measurements, tracks, plat_pos, params);
info.t_assoc = toc;
info.beta = beta;
info.beta0 = beta0;
info.n_events = jinfo.n_events;
info.fallback = jinfo.fallback;

% ---------- 第3步: 模型条件更新 + 组合 ----------
for t = 1:num_tracks
    if num_meas == 0
        tracks{t}.beta0_last = 1.0;
        [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);
        continue;
    end

    tracks{t}.imm = imm_update(tracks{t}.imm, measurements, mps{t}, beta{t}, params);
    [tracks{t}.x, tracks{t}.P] = imm_combine(tracks{t}.imm);

    % 保存检测状态 (hit/miss 由调用方融合后统一处理)
    tracks{t}.beta0_last = beta0(t);

    % 更新频谱模板 (仅在特征辅助时生效)
    if params.use_features
        tracks{t} = update_track_spectrum(tracks{t}, measurements, beta{t}, params);
    end
end
end
