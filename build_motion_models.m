function mm = build_motion_models(params)
%BUILD_MOTION_MODELS 构建 IMM 的模型集合
%   mm = BUILD_MOTION_MODELS(params)
%
%   返回结构体数组, 每个元素:
%     .name  模型名称 ('CV' / 'CT' / 'CA')
%     .F     线性状态转移矩阵 (非线性模型为空矩阵)
%     .Q     过程噪声协方差 [7x7]
%     .f     状态转移函数句柄 f(x, dt)
%
%   支持的模型:
%     CV — 匀速模型, 线性外推
%     CT — 协调转弯模型, 转弯率 omega 为状态量, 非线性 (由 UKF 处理)
%     CA — 匀加速模型, 加速度 (ax, ay) 为状态量
%
%   状态向量: x = [px; py; vx; vy; omega; ax; ay]

names = params.imm_models;
n_models = numel(names);

mm = struct('name', {}, 'F', {}, 'Q', {}, 'f', {});

for j = 1:n_models
    name = upper(names{j});
    switch name
        case 'CV'
            F = eye(7);
            F(1, 3) = params.dt;
            F(2, 4) = params.dt;
            Q = accel_noise(params.dt, params.sigma_a_cv);
            Q(5, 5) = params.sigma_omega_cv^2;
            Q(6, 6) = params.sigma_accel_rw^2;
            Q(7, 7) = params.sigma_accel_rw^2;
            f = @motion_cv;
        case 'CT'
            F = [];
            Q = accel_noise(params.dt, params.sigma_a_ct);
            Q(5, 5) = params.sigma_omega_ct^2;
            Q(6, 6) = params.sigma_accel_rw^2;
            Q(7, 7) = params.sigma_accel_rw^2;
            f = @motion_ct;
        case 'CA'
            F = [];
            Q = jerk_noise(params.dt, params.sigma_jerk);
            Q(5, 5) = params.sigma_omega_ca^2;
            f = @motion_ca;
        otherwise
            error('build_motion_models:badModel', '未知运动模型: %s', name);
    end

    mm(j).name = name;
    mm(j).F = F;
    mm(j).Q = Q;
    mm(j).f = f;
end
end


function Q = accel_noise(dt, sigma_a)
%ACCEL_NOISE 离散白噪声加速度模型 (x/y 两轴独立), 填充位置/速度 4 维
q11 = dt^3 / 3 * sigma_a^2;
q12 = dt^2 / 2 * sigma_a^2;
q22 = dt * sigma_a^2;

Q = zeros(7);
Q([1 2], [1 2]) = diag([q11, q11]);
Q([1 2], [3 4]) = diag([q12, q12]);
Q([3 4], [1 2]) = diag([q12, q12]);
Q([3 4], [3 4]) = diag([q22, q22]);
end


function Q = jerk_noise(dt, sigma_j)
%JERK_NOISE 离散白噪声 jerk 模型 (匀加速/CA), x/y 两轴独立
%   状态块为 [位置; 速度; 加速度], 噪声标准差 sigma_j (m/s^3)
q = sigma_j^2 * [dt^5/20, dt^4/8, dt^3/6;
                 dt^4/8,  dt^3/3, dt^2/2;
                 dt^3/6,  dt^2/2, dt];

Q = zeros(7);
Q([1 3 6], [1 3 6]) = q;      % x 轴: px, vx, ax
Q([2 4 7], [2 4 7]) = q;      % y 轴: py, vy, ay
end
