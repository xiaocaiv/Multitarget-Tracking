function P = stabilize_cov(P)
%STABILIZE_COV 对称化并把协方差矩阵拉回正半定 (数值保护)
%   P = STABILIZE_COV(P)
%
%   UKF 的加权协方差与 JPDA 的协方差更新在强非线性/病态情形下可能出现
%   微小负特征值; 这里做特征值截断, 避免后续 Cholesky 分解失败.
%   截断阈值取矩阵量级的 1e-10, 不影响正常协方差.

P = (P + P') / 2;
if any(~isfinite(P(:)))
    return;
end

[V, D] = eig(P);
d = diag(D);

tol = 1e-10 * max([1; abs(d)]);
if any(d < tol)
    d = max(d, tol);
    P = V * diag(d) * V';
    P = (P + P') / 2;
end
end
