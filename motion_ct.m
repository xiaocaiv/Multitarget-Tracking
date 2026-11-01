function x_next = motion_ct(x, dt)
%MOTION_CT 协调转弯 (CT) 模型状态转移
%   x_next = MOTION_CT(x, dt)
%
%   状态向量: x = [px; py; vx; vy; omega; ax; ay]
%   omega > 0: 逆时针 (左转);  omega < 0: 顺时针 (右转).
%   |omega*dt| 很小时退化为匀速直线, 避免除以零.
%   CT 模型不使用加速度 (由 CA 模型负责), 保持其原值.
%
%   解析解 (匀速 + 恒定转弯率):
%     px' = px + sin(wT)/w * vx - (1-cos(wT))/w * vy
%     py' = py + (1-cos(wT))/w * vx + sin(wT)/w * vy
%     vx' = cos(wT) * vx - sin(wT) * vy
%     vy' = sin(wT) * vx + cos(wT) * vy
%     w'  = w

w = x(5);
wt = w * dt;

if abs(wt) < 1e-9
    % 极限情形: 匀速直线运动
    x_next = [x(1) + dt * x(3);
              x(2) + dt * x(4);
              x(3);
              x(4);
              w;
              x(6);
              x(7)];
else
    s = sin(wt) / w;            % = T * sinc(wT/pi)
    c = (1 - cos(wt)) / w;
    cw = cos(wt);
    sw = sin(wt);

    x_next = [x(1) + s * x(3) - c * x(4);
              x(2) + c * x(3) + s * x(4);
              cw * x(3) - sw * x(4);
              sw * x(3) + cw * x(4);
              w;
              x(6);
              x(7)];
end
end
