function truth = generate_truth(params)
%GENERATE_TRUTH 生成目标真实轨迹 (分段匀速 / 协调转弯 / 切向变速, 无随机扰动)
%   truth = GENERATE_TRUTH(params)
%
%   目标状态为 5 维 [x; y; vx; vy; omega]:
%     - 按 params.target_turns 指定段做协调转弯 ([起始步, 结束步, deg/s])
%     - 按 params.target_accel 指定段做切向变速 ([起始步, 结束步, m/s^2])
%   其余时刻匀速直线. 场景由 params.scenario 选择 (见 scenario_defs.m).
%
%   truth 结构体:
%     truth.X:      元胞数组 {num_targets}, 每个 [5 x K] 完整状态历史
%     truth.pos:    元胞数组 {num_targets}, 每个 [2 x K] 位置历史
%     truth.omega:  元胞数组 {num_targets}, 每个 [1 x K] 转弯率历史 (rad/s)
%     truth.speed:  元胞数组 {num_targets}, 每个 [1 x K] 速率历史 (m/s)
%     truth.num_targets: 目标数

K = params.K;
num_targets = params.num_targets;

truth = struct();
truth.num_targets = num_targets;
truth.X = cell(num_targets, 1);
truth.pos = cell(num_targets, 1);
truth.omega = cell(num_targets, 1);
truth.speed = cell(num_targets, 1);

% 目标初始状态
targets_init = {params.target1_init, params.target2_init};

for tgt = 1:num_targets
    X_hist = zeros(5, K);
    X_hist(1:4, 1) = targets_init{tgt}(:);

    turns = params.target_turns{tgt};
    accel = params.target_accel{tgt};
    if ~params.maneuver_enable
        turns = zeros(0, 3);
        accel = zeros(0, 3);
    end

    omega_hist = zeros(1, K);

    for k = 2:K
        % 当前区间 (k-1, k] 的转弯率与切向加速度
        w = segment_value(turns, k, pi / 180);     % deg/s -> rad/s
        a = segment_value(accel, k, 1);

        v = X_hist(3:4, k - 1);
        if a ~= 0 && norm(v) > 0
            v_new = v + a * params.dt * (v / norm(v));
            if dot(v_new, v) <= 0
                v_new = [0; 0];                    % 防止减速到反向
            end
            v = v_new;
        end

        % 真值状态保持 5 维 [x;y;vx;vy;omega]; 用 7 维包装调用 kinematics
        % (加速度在区间内已直接作用于速度, 这里状态中的 a 置零)
        x_prev = [X_hist(1:2, k - 1); v; w; 0; 0];
        x_next = motion_ct(x_prev, params.dt);
        X_hist(:, k) = x_next(1:5);
        omega_hist(k) = w;
    end

    truth.X{tgt} = X_hist;
    truth.pos{tgt} = X_hist(1:2, :);
    truth.omega{tgt} = omega_hist;
    truth.speed{tgt} = sqrt(sum(X_hist(3:4, :).^2, 1));
end

%% 频谱特征真值 (调用 generate_ship_spectrum 生成一次, 后续帧从缓存读取)
truth.spectrum = cell(num_targets, 1);
truth.spec_info = cell(num_targets, 1);
for tgt = 1:num_targets
    [truth.spectrum{tgt}, truth.spec_info{tgt}] = ...
        generate_ship_spectrum(params, params.target_specs(tgt), ...
        params.target_specs(tgt).wav_file);
end
end

% ========================================================================
function val = segment_value(segments, k, scale)
%SEGMENT_VALUE 返回第 k 步所属区间的设定值 (无匹配段则为 0)
val = 0;
for s = 1:size(segments, 1)
    if k >= segments(s, 1) && k <= segments(s, 2)
        val = segments(s, 3) * scale;
        return;
    end
end
end
