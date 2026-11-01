function [d, comp] = gospar(X, Y, p, c, alpha, base_metric)
%GOSPAR 广义最优子模式分配 (GOSPA) 距离 —— 多目标跟踪性能评价指标
%
%   d = GOSPAR(X, Y)
%   [d, comp] = GOSPAR(X, Y, p, c, alpha, base_metric)
%
%   输入:
%     X, Y        - d×n / d×m 点集矩阵, 每列一个目标 (或航迹) 的位置
%     p           - 距离阶数, 默认 2 (有限正实数, 常用 1 或 2)
%     c           - 截断参数 (m), 默认 100
%     alpha       - 未关联项的基数惩罚因子, 默认 2 (推荐值)
%     base_metric - 基点距离函数句柄, 默认欧氏距离 @(a,b)norm(a-b)
%
%   输出:
%     d    - GOSPA 距离 (标量, 单位与基点距离相同, 即 m)
%     comp - 分解结构体:
%              .loc       定位误差分量 (m)
%              .missed    漏检目标分量 (m)
%              .false     虚假航迹分量 (m)
%              .n_matched 成功关联的目标-航迹对数
%              .n_missed  漏检目标数
%              .n_false   虚假航迹数
%              .cost_p    GOSPA^p (未开 p 次方, 便于跨帧累加/平均)
%
%   定义 (Rahmathullah, Garcia-Fernandez, Svensson,
%   "Generalized optimal sub-pattern assignment metric", FUSION 2017):
%
%     d_p(X,Y)^p = min_pi { sum_{i} d(x_i, y_pi(i))^p
%                           + (c^p/alpha) * (n + m - 2*|pi|) }
%
%   其中 pi 为目标-航迹的部分关联 (每个目标至多关联一条航迹, 反之亦然),
%   |pi| 为成功关联对数; 每个未关联的目标或航迹各罚 c^p/alpha.
%   因此当某对目标-航迹的距离 d 满足 d^p < 2*c^p/alpha 时才值得关联,
%   alpha=2 时的临界距离正好是 c (超过 c 的误差按 c 计, 即定位误差饱和).
%
%   分解关系 (GOSPA 的核心优点): GOSPA^p = loc^p + missed^p + false^p.
%   alpha=2 时 GOSPA 关于 X、Y 对称, 且无需事先指定目标-航迹的对应关系.
%
%   参考文献:
%     A. S. Rahmathullah, A. F. Garcia-Fernandez, L. Svensson,
%     "Generalized optimal sub-pattern assignment metric",
%     20th Int. Conf. on Information Fusion (FUSION), 2017, pp. 1-8.

if nargin < 3 || isempty(p),           p = 2;              end
if nargin < 4 || isempty(c),           c = 100;            end
if nargin < 5 || isempty(alpha),       alpha = 2;          end
if nargin < 6,                         base_metric = [];   end

if ~(isscalar(p) && p > 0 && isfinite(p))
    error('gospar:badP', 'p 必须为有限正实数 (常用 1 或 2).');
end
if ~(isscalar(c) && c > 0 && isfinite(c))
    error('gospar:badC', 'c 必须为有限正实数.');
end
if ~(isscalar(alpha) && alpha >= 1)
    error('gospar:badAlpha', 'alpha 必须 >= 1 (推荐 2).');
end

n = size(X, 2);   % 真值目标数
m = size(Y, 2);   % 估计航迹数

% 论文定义以较小集合为 "目标集"; alpha ~= 2 时度量不对称, 需交换后计算.
if n > m
    [d, comp] = gospar(Y, X, p, c, alpha, base_metric);
    tmp = comp.missed;   comp.missed   = comp.false;    comp.false   = tmp;
    tmp = comp.n_missed; comp.n_missed = comp.n_false;  comp.n_false = tmp;
    return;
end

cu = c^p / alpha;            % 单个未关联项 (漏检目标 或 虚假航迹) 的代价

% ---------- 1. 目标-航迹距离的 p 次方矩阵 ----------
D = zeros(n, m);
if isempty(base_metric)
    % 欧氏距离: 用内积展开, 避免逐点循环
    if n > 0 && m > 0
        D2 = sum(X.^2, 1)' + sum(Y.^2, 1) - 2 * (X' * Y);
        D2(D2 < 0) = 0;                 % 数值误差保护
        D = D2 .^ (p / 2);
    end
else
    for i = 1:n
        for j = 1:m
            D(i, j) = base_metric(X(:, i), Y(:, j))^p;
        end
    end
end

% ---------- 2. 部分关联的最优求解 ----------
% 构造 (n+m)×(n+m) 方阵:
%   行 = n 个目标 + m 个"空槽行"(吸收未关联的航迹)
%   列 = m 条航迹 + n 个"空槽列"(吸收未关联的目标)
% 目标-航迹配对的代价取 D(i,j)-2*cu, 其余全部为 0; 于是
%   总代价 = cu*(n+m) + 最优部分关联代价
% 其中"未关联"与"空槽-空槽"都以 0 计, 恰好等价于罚 cu.
% 该构造不需要大 M 罚值, 数值上比增广矩阵更稳健.
A = zeros(n + m, n + m);
A(1:n, 1:m) = D - 2 * cu;

[assign, ~] = hungarian_square(A);

n_matched = sum(assign(1:n) <= m);      % 目标行落到航迹列 = 成功关联
n_missed  = n - n_matched;
n_false   = m - n_matched;

% 关联对的定位误差平方和
loc_p = 0;
rows = find(assign(1:n) <= m);
for i = rows(:)'
    loc_p = loc_p + D(i, assign(i));
end

% ---------- 3. GOSPA 及其分解 ----------
cost_p = cu * (n + m) + sum(A(sub2ind(size(A), (1:n + m)', assign)));
cost_p = max(cost_p, 0);                % 数值误差保护

d = cost_p^(1 / p);

comp = struct();
comp.loc       = loc_p^(1 / p);
comp.missed    = c * (n_missed / alpha)^(1 / p);
comp.false     = c * (n_false  / alpha)^(1 / p);
comp.n_matched = n_matched;
comp.n_missed  = n_missed;
comp.n_false   = n_false;
comp.cost_p    = cost_p;
end


function [assign, cost] = hungarian_square(C)
%HUNGARIAN_SQUARE 方阵最小权完美匹配 (匈牙利算法, O(n^3))
%   [assign, cost] = HUNGARIAN_SQUARE(C)
%     assign - n×1, assign(i) 为第 i 行匹配的列号
%     cost   - 最小总代价 (列势 + 行势)
%
%   实现: 逐行增广 (最短增广路 + 对偶势更新), 允许负数/零代价,
%   不支持 Inf/NaN. 仅在本文件内部使用, 避免依赖 matchpairs (R2019a+).

n = size(C, 1);
u = zeros(n, 1);        % 行势
v = zeros(n, 1);        % 列势
pcol = zeros(n, 1);     % pcol(j): 第 j 列匹配的行号, 0 = 未匹配
way = zeros(n, 1);

for i = 1:n
    p0 = i;             % p[0]: 本轮起点行
    j0 = 0;             % 0 = 虚拟列
    minv = inf(n, 1);
    used = false(n, 1);

    while true
        if j0 > 0
            used(j0) = true;
            i0 = pcol(j0);
        else
            i0 = p0;
        end

        delta = inf;
        j1 = 0;
        for j = 1:n
            if ~used(j)
                cur = C(i0, j) - u(i0) - v(j);
                if cur < minv(j)
                    minv(j) = cur;
                    way(j) = j0;
                end
                if minv(j) < delta
                    delta = minv(j);
                    j1 = j;
                end
            end
        end

        % 对偶势更新 (虚拟列 0 始终处于"已用"状态)
        u(p0) = u(p0) + delta;
        for j = 1:n
            if used(j)
                u(pcol(j)) = u(pcol(j)) + delta;
                v(j) = v(j) - delta;
            else
                minv(j) = minv(j) - delta;
            end
        end

        j0 = j1;
        if pcol(j0) == 0
            break;
        end
    end

    % 沿增广路回溯
    while j0 ~= 0
        j1 = way(j0);
        if j1 == 0
            pcol(j0) = p0;
        else
            pcol(j0) = pcol(j1);
        end
        j0 = j1;
    end
end

assign = zeros(n, 1);
for j = 1:n
    if pcol(j) > 0
        assign(pcol(j)) = j;
    end
end
cost = sum(u) + sum(v);
end
