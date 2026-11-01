function [x, P] = ci_fuse(x1, P1, x2, P2, x_pred, P_pred)
%CI_FUSE 协方差交叉融合 (信息矩阵形式)
%   [x, P] = CI_FUSE(x1, P1, x2, P2, x_pred, P_pred)
%
%   利用同一先验 (预测) 去掉两个平台重复使用的信息:
%     P^-1 = P1^-1 + P2^-1 - P_pred^-1
%     x    = P * (P1^-1 x1 + P2^-1 x2 - P_pred^-1 x_pred)
%
%   若信息矩阵因重复计数而非正定, 退化为直接相加信息 (更保守),
%   此时不再减去先验项.

I1 = inv(P1);
I2 = inv(P2);
Ip = inv(P_pred);

I_new = I1 + I2 - Ip;
I_new = (I_new + I_new') / 2;

b = I1 * x1 + I2 * x2;

[~, flag] = chol(I_new);
if flag ~= 0
    % 非正定: 去掉先验项, 避免出现非物理的协方差
    I_new = I1 + I2;
else
    b = b - Ip * x_pred;
end

P = inv(I_new);
P = stabilize_cov(P);
x = I_new \ b;
end
