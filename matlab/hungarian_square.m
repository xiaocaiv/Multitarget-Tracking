function [assign, cost] = hungarian_square(C)
%HUNGARIAN_SQUARE 方阵最小权完美匹配 (匈牙利算法, O(n^3))
%   [assign, cost] = HUNGARIAN_SQUARE(C)
%     assign - n×1, assign(i) 为第 i 行匹配的列号
%     cost   - 最小总代价 (列势 + 行势)
%
%   实现: 逐行增广 (最短增广路 + 对偶势更新), 允许负数/零代价,
%   不支持 Inf/NaN. 原为 gospar.m 的局部函数, 现独立成文件以便
%   GNN 关联与关联层评价复用 (见 assoc_gnn.m / eval_assoc.m).
%
%   参见: gospar.m, assoc_gnn.m, eval_assoc.m

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
