function x_next = motion_cv(x, dt)
%MOTION_CV 匀速 (CV) 模型状态转移
%   x_next = MOTION_CV(x, dt)
%
%   状态向量: x = [px; py; vx; vy; omega; ax; ay]
%   位置按速度线性外推, 速度/转弯率/加速度保持不变.
%   CV 模型不使用 omega 与加速度 (分别由 CT / CA 模型估计并在 IMM 中混合),
%   这里保持其原值, 使各模型共享同一状态维度.

x_next = [x(1) + dt * x(3);
          x(2) + dt * x(4);
          x(3);
          x(4);
          x(5);
          x(6);
          x(7)];
end
