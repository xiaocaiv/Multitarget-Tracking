function [beta, unassociated, info] = associate_measurements(z_preds, Ss, measurements, tracks, params)
%ASSOCIATE_MEASUREMENTS 数据关联统一入口 (按 params.assoc_method 分发)
%   [beta, unassociated, info] = ASSOCIATE_MEASUREMENTS(z_preds, Ss, measurements, tracks, params)
%
%   方法: 'nn' | 'gnn' | 'pda' | 'jpda' (缺省/未知字段时用 jpda)
%   输出:
%     beta          n×1 元胞, 每个 [beta0; beta_1; ...; beta_m]
%     unassociated  1×m 逻辑, 量测被任一航迹占用的总 beta < 阈值
%     info          关联层记录: .method / .beta / .beta0 / .track_ids /
%                   .t_assoc / .n_events / .fallback / .unassoc
%
%   各方法输出同格式 beta, 因此下游 imm_update / 融合 / 航迹管理无需区分.
%
%   参见: assoc_nn.m, assoc_gnn.m, assoc_pda.m, jpda_associate.m,
%         process_platform.m

n = numel(z_preds);
m = numel(measurements);

method = 'jpda';
if isfield(params, 'assoc_method') && ~isempty(params.assoc_method)
    method = lower(char(params.assoc_method));
end

info = struct();
info.method = method;
info.beta = cell(n, 1);
info.beta0 = ones(n, 1);
info.track_ids = reshape(cellfun(@(t) t.id, tracks), 1, []);
info.t_assoc = 0;
info.n_events = 0;
info.fallback = false;

% 无航迹或无量测: 所有量测视为杂波, 航迹保持预测
if n == 0 || m == 0
    beta = cell(n, 1);
    for t = 1:n
        beta{t} = [1; zeros(m, 1)];
    end
    unassociated = true(1, m);
    info.beta = beta;
    info.unassoc = unassociated;
    return;
end

tic;
switch method
    case 'nn'
        [beta, beta0] = assoc_nn(z_preds, Ss, measurements, tracks, params);
    case 'gnn'
        [beta, beta0] = assoc_gnn(z_preds, Ss, measurements, tracks, params);
    case 'pda'
        [beta, beta0] = assoc_pda(z_preds, Ss, measurements, tracks, params);
    case 'jpda'
        [beta, ~, ~, beta0, ~, jinfo] = jpda_associate( ...
            z_preds, Ss, cell(n, 1), measurements, tracks, [], params);
        info.n_events = jinfo.n_events;
        info.fallback = jinfo.fallback;
    otherwise
        error('associate_measurements:badMethod', ...
            '未知关联算法 "%s" (应为 nn / gnn / pda / jpda)', method);
end
info.t_assoc = toc;

% 未关联量测: 总关联概率低于阈值
unassoc_thr = 0.5;
if isfield(params, 'assoc_unassoc_threshold') && ~isempty(params.assoc_unassoc_threshold)
    unassoc_thr = params.assoc_unassoc_threshold;
end

unassociated = false(1, m);
for j = 1:m
    total_beta = 0;
    for t = 1:n
        total_beta = total_beta + beta{t}(j + 1);
    end
    unassociated(j) = (total_beta < unassoc_thr);
end

info.beta = beta;
info.beta0 = beta0;
info.unassoc = unassociated;
end
