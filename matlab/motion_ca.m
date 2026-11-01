function x_next = motion_ca(x, dt)
%MOTION_CA 匀加速 (CA) 模型状态转移
%   x_next = MOTION_CA(x, dt)
%
%   状态向量: x = [px; py; vx; vy; omega; ax; ay]
%   位置按当前速度与加速度外推, 速度按加速度更新, 加速度保持不变.
%   CA 模型不转弯, omega 保持原值 (不使用, 由 CT 模型负责).
%
%     px' = px + T*vx + T^2/2*ax
%     py' = py + T*vy + T^2/2*ay
%     vx' = vx + T*ax
%     vy' = vy + T*ay
%     ax' = ax,  ay' = ay,  omega' = omega

T = dt;

x_next = [x(1) + T * x(3) + 0.5 * T^2 * x(6);
          x(2) + T * x(4) + 0.5 * T^2 * x(7);
          x(3) + T * x(6);
          x(4) + T * x(7);
          x(5);
          x(6);
          x(7)];
end
