function [tracks, next_id] = manage_tracks(tracks, z1, z2, unassoc1, unassoc2, ...
    params, next_id)
%MANAGE_TRACKS 航迹生命周期管理（删除 + 三角测量起批）
%
%   [tracks, next_id] = MANAGE_TRACKS(tracks, z1, z2, unassoc1, unassoc2, params, next_id)
%
%   删除条件
%     - 连续丢失次数 >= max_misses
%     - 航迹位置超出监视区域
%   起批条件
%     - 两平台未关联量测配对后三角测量，位置须在监视区内、两站前方、距离合理
%     - 不与现有航迹/本次候选重复
%     - 若量测带声纹（.line/.band）且 params.llr 存在：
%         用校准 LLR 打分，做**最大权匹配**择优，取代硬阈值门控
%       否则退回旧逻辑（96 维 spectrum 的 Pearson 硬阈值 / 纯几何）
%
%   与旧版的区别
%     旧版：rho_demon>=0.3 && rho_lofar>=0.6 的硬门控 + 双重循环先到先得，
%           2 站 2 目标时会随机锁死一种配对（可能是幽灵）。
%     新版：枚举所有几何可行的配对，按 LLR 求全局最优指派。

%% ---------- 删除丢失航迹 ----------
keep_mask = true(length(tracks), 1);
for t = 1:length(tracks)
    if tracks{t}.miss_count >= params.max_misses
        keep_mask(t) = false;
        continue;
    end
    x_region = params.surveil_region;
    if tracks{t}.x(1) < x_region(1) || tracks{t}.x(1) > x_region(2) ...
        || tracks{t}.x(2) < x_region(3) || tracks{t}.x(2) > x_region(4)
        keep_mask(t) = false;
    end
end
tracks = tracks(keep_mask);

%% ---------- 确认现有航迹 ----------
for t = 1:length(tracks)
    if tracks{t}.hit_count >= params.confirm_hits
        tracks{t}.confirmed = true;
    end
end

%% ---------- 起批 ----------
u1 = z1(unassoc1);
u2 = z2(unassoc2);
if isempty(u1) || isempty(u2)
    return;
end

has_llr = isfield(params, 'llr') && isfield(u1(1), 'line') && ~isempty(u1(1).line) ...
          && isfield(u2(1), 'line') && ~isempty(u2(1).line);

if has_llr
    [tracks, next_id] = init_by_llr(tracks, u1, u2, params, next_id);
elseif params.use_features
    [tracks, next_id] = init_legacy_spectrum(tracks, u1, u2, params, next_id);
else
    [tracks, next_id] = init_geometry(tracks, u1, u2, params, next_id);
end
end

%% ========================================================================
function [tracks, next_id] = init_by_llr(tracks, u1, u2, params, next_id)
% 声纹 LLR 打分 + 最大权匹配
ell_gate = getdef(params, 'ell_gate', 0);
dup_dist = getdef(params, 'dup_dist', 200);
r_min    = getdef(params, 'init_r_min', 50);
r_max    = getdef(params, 'init_r_max', 5000);
max_new  = getdef(params, 'max_new_tracks', 1);

n1 = numel(u1); n2 = numel(u2);
L    = -inf(n1, n2);
Xp   = nan(size(params.init_cov, 1), n1, n2);
feas = false(n1, n2);

for i = 1:n1
    for j = 1:n2
        [x_init, ~, ok] = tri_point(u1(i).bearing, u2(j).bearing, params, r_min, r_max);
        if ~ok, continue; end
        if is_dup(x_init, tracks, [], dup_dist), continue; end

        ell = feature_llr(u1(i), u2(j), params.llr, 'init');
        q   = min(getq(u1(i)), getq(u2(j)));
        L(i, j) = ell * q;
        Xp(:, i, j) = x_init;
        feas(i, j) = L(i, j) > ell_gate;
    end
end
if ~any(feas(:)), return; end

pairs = solve_assignment(L, feas);

created = zeros(size(params.init_cov, 1), 0);
for p = 1:size(pairs, 1)
    if size(created, 2) >= max_new, break; end
    i = pairs(p, 1); j = pairs(p, 2);
    x_init = Xp(:, i, j);
    if is_dup(x_init, tracks, created, dup_dist), continue; end

    nt = make_track(x_init, params, next_id);
    next_id = next_id + 1;

    qa = getq(u1(i)); qb = getq(u2(j));
    tl = qa * u1(i).line + qb * u2(j).line;
    tb = qa * u1(i).band + qb * u2(j).band;
    nt.tpl_line = tl / (norm(tl) + 1e-12);
    nt.tpl_band = tb / (norm(tb) + 1e-12);
    nt.q = min(qa, qb);
    nt.feature_init = true;
    nt.feat_step = getdef(params, 'k', NaN);
    nt.ell0 = L(i, j);
    nt.eval_src = getdef(u1(i), 'src', 0);     % 仅供评估

    tracks{end + 1} = nt; %#ok<AGROW>
    created(:, end + 1) = x_init; %#ok<AGROW>
end
end

%% ========================================================================
function [tracks, next_id] = init_legacy_spectrum(tracks, u1, u2, params, next_id)
% 旧逻辑：96 维 spectrum 的 Pearson 硬阈值（保留以便回归对比）
n_lofar = params.N_lofar;
new_count = 0;
candidates = [];
for i = 1:numel(u1)
    for j = 1:numel(u2)
        if new_count >= params.max_new_tracks, break; end
        [x_init, P_init] = triangulate_init(u1(i).bearing, u2(j).bearing, ...
            params.plat1_pos, params.plat2_pos, params.init_cov);
        xr = params.surveil_region;
        if x_init(1) < xr(1) || x_init(1) > xr(2) || ...
           x_init(2) < xr(3) || x_init(2) > xr(4)
            continue;
        end
        if is_dup(x_init, tracks, candidates, 200), continue; end

        xl = u1(i).spectrum(1:n_lofar); xl = xl - mean(xl);
        yl = u2(j).spectrum(1:n_lofar); yl = yl - mean(yl);
        rho_lofar = (xl' * yl) / (norm(xl) * norm(yl) + 1e-10);
        xd = u1(i).spectrum(n_lofar+1:end); xd = xd - mean(xd);
        yd = u2(j).spectrum(n_lofar+1:end); yd = yd - mean(yd);
        rho_demon = (xd' * yd) / (norm(xd) * norm(yd) + 1e-10);
        if rho_demon < 0.3 || rho_lofar < 0.6, continue; end

        nt = make_track(x_init, params, next_id);
        next_id = next_id + 1;
        spec_avg = (u1(i).spectrum + u2(j).spectrum) / 2;
        nt.spectrum = spec_avg / norm(spec_avg);
        nt.feature_init = true;
        nt.eval_src = getdef(u1(i), 'src', 0);
        tracks{end + 1} = nt; %#ok<AGROW>
        candidates(:, end + 1) = x_init(1:2); %#ok<AGROW>
        new_count = new_count + 1;
    end
    if new_count >= params.max_new_tracks, break; end
end
end

%% ========================================================================
function [tracks, next_id] = init_geometry(tracks, u1, u2, params, next_id)
% 纯几何：只能用"合法指派"，但在 2 种可行方案之间无法区分（抛硬币）。
% 为保证基线不固定偏好某一种配对，这里先随机打乱量测顺序再做匹配。
dup_dist = getdef(params, 'dup_dist', 200);
r_min = getdef(params, 'init_r_min', 50);
r_max = getdef(params, 'init_r_max', 5000);
max_new = getdef(params, 'max_new_tracks', 1);

n1 = numel(u1); n2 = numel(u2);
feas = false(n1, n2);
Xp = nan(size(params.init_cov, 1), n1, n2);
for i = 1:n1
    for j = 1:n2
        [x_init, ~, ok] = tri_point(u1(i).bearing, u2(j).bearing, params, r_min, r_max);
        if ~ok, continue; end
        if is_dup(x_init, tracks, [], dup_dist), continue; end
        feas(i, j) = true;
        Xp(:, i, j) = x_init;
    end
end
if ~any(feas(:)), return; end

pr = randperm(n1); pc = randperm(n2);
Fr = feas(pr, pc);
pairs_r = solve_assignment(double(Fr), Fr);
pairs = [pr(pairs_r(:, 1))', pc(pairs_r(:, 2))'];

created = zeros(size(params.init_cov, 1), 0);
for p = 1:size(pairs, 1)
    if size(created, 2) >= max_new, break; end
    i = pairs(p, 1); j = pairs(p, 2);
    x_init = Xp(:, i, j);
    if is_dup(x_init, tracks, created, dup_dist), continue; end
    nt = make_track(x_init, params, next_id);
    next_id = next_id + 1;
    nt.feature_init = false;
    nt.eval_src = getdef(u1(i), 'src', 0);
    tracks{end + 1} = nt; %#ok<AGROW>
    created(:, end + 1) = x_init; %#ok<AGROW>
end
end

%% ========================================================================
function nt = make_track(x_init, params, id)
nt = struct();
nt.x = x_init;
nt.P = params.init_cov;
nt.id = id;
nt.hit_count = 0;
nt.miss_count = 0;
nt.confirmed = false;
nt.imm = imm_init(x_init, params.init_cov, params);   % IMM 初始化
nt.tpl_line = [];
nt.tpl_band = [];
nt.q = 0;
nt.feature_init = false;
nt.feat_step = NaN;
nt.ell0 = NaN;
nt.eval_src = 0;
nt.spectrum = [];
end

function [x, P, ok] = tri_point(b1, b2, params, r_min, r_max)
% 三角定位 + 几何合理性检查
ok = false;
x = zeros(size(params.init_cov, 1), 1);
P = params.init_cov;
if abs(wrap_pi(b1 - b2)) < 1e-3 || abs(abs(wrap_pi(b1 - b2)) - pi) < 1e-3
    return;                                   % 近平行，无法定位
end
[x_init, P_init] = triangulate_init(b1, b2, params.plat1_pos, params.plat2_pos, ...
    params.init_cov);
p  = x_init(1:2);
p1 = params.plat1_pos(:); p2 = params.plat2_pos(:);
r1 = norm(p - p1); r2 = norm(p - p2);
d1 = [cos(b1); sin(b1)];
d2 = [cos(b2); sin(b2)];
in_front = dot(p - p1, d1) > 0 && dot(p - p2, d2) > 0;
xr = params.surveil_region;
in_region = p(1) >= xr(1) && p(1) <= xr(2) && p(2) >= xr(3) && p(2) <= xr(4);
if in_front && in_region && r1 >= r_min && r1 <= r_max && r2 >= r_min && r2 <= r_max
    x = x_init;
    P = P_init;
    ok = true;
end
end

function tf = is_dup(x, tracks, extra, dmin)
tf = false;
p = x(1:2);
for t = 1:numel(tracks)
    if norm(p - tracks{t}.x(1:2)) < dmin, tf = true; return; end
end
if ~isempty(extra)
    for c = 1:size(extra, 2)
        if norm(p - extra(1:2, c)) < dmin, tf = true; return; end
    end
end
end

function pairs = solve_assignment(L, feas)
% 最大权匹配；有 matchpairs 就用它，否则退化为贪心
[n1, n2] = size(L);
pairs = zeros(0, 2);
if exist('matchpairs', 'file') == 2
    Cmax = max(L(feas)) + 1;
    BIG  = 1e6;
    cost = BIG * ones(n1, n2);
    cost(feas) = Cmax - L(feas);
    costUnmatched = Cmax;
    try
        m = matchpairs(cost, costUnmatched);
        keep = false(size(m, 1), 1);
        for r = 1:size(m, 1)
            keep(r) = feas(m(r, 1), m(r, 2));
        end
        m = m(keep, :);
        [~, ord] = sort(L(sub2ind(size(L), m(:, 1), m(:, 2))), 'descend');
        pairs = m(ord, :);
        return;
    catch
        % 落到贪心
    end
end
[~, idx] = sort(L(feas), 'descend');
[ii, jj] = find(feas);
used1 = false(n1, 1); used2 = false(n2, 1);
for q = 1:numel(idx)
    i = ii(idx(q)); j = jj(idx(q));
    if used1(i) || used2(j), continue; end
    used1(i) = true; used2(j) = true;
    pairs(end + 1, :) = [i, j]; %#ok<AGROW>
end
end

function v = getdef(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function q = getq(m)
q = getdef(m, 'q', 1);
if isnan(q) || q <= 0, q = 0.1; end
q = min(max(q, 0), 1);
end

function th = wrap_pi(th)
th = mod(th + pi, 2 * pi) - pi;
end
