% test_assoc.m — 数据关联算法一致性 / 回归测试
%
% 运行 (仓库根目录):
%   matlab -batch "addpath('matlab'); addpath('matlab/数据关联算法对比'); test_assoc"
%
% 覆盖:
%   1. 各方法 beta 合法性 (非负、归一、长度)
%   2. 空输入 / 全部门外
%   3. 单航迹 PDA 与 JPDA 一致
%   4. GNN 互斥性 / NN 复用性
%   5. 特征因子判别
%   6. 确定性
%   7. JPDA 事件上限回退
%   8. GOSPA 回归 (匈牙利抽取)
%   9. 默认路径 run_simulation 回归 (bit-for-bit)
%  10. 路径自举 + GNN 冒烟

clear; close all;

here = fileparts(mfilename('fullpath'));
parent = fileparts(here);
addpath(parent);
addpath(here);

fprintf('===== 数据关联算法测试 =====\n');
n_pass = 0;
k = 0;

params = init_parameters(struct('use_features', false));

%% 1. beta 合法性
methods = {'nn', 'gnn', 'pda', 'jpda'};
z = {0.10, 0.12, -0.40};
S = {0.01, 0.01, 0.01};
ms = make_meas([0.11, 0.13, 0.50], [1, 2, 0], {[], [], []});
tracks = {make_track(1, params), make_track(2, params), make_track(3, params)};

for mi = 1:numel(methods)
    p = params;
    p.assoc_method = methods{mi};
    [beta, u, info] = associate_measurements(z, S, ms, tracks, p);
    check_beta(beta, numel(z), numel(ms), methods{mi});
    assert(numel(u) == numel(ms));
    assert(strcmp(info.method, methods{mi}));
    n_pass = n_pass + 1;
end
k = k + 1;
fprintf('  [%d] beta 合法性\n', k);

%% 2. 空输入
p = params; p.assoc_method = 'gnn';
[beta, u] = associate_measurements({}, {}, ms, {}, p);
assert(isempty(beta) && all(u));
[beta, u] = associate_measurements(z, S, make_meas([], [], {}), tracks, p);
assert(numel(beta) == numel(tracks));
for t = 1:numel(beta)
    assert(isequal(beta{t}, 1));
end
assert(isempty(u));
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 空输入\n', k);

%% 3. 全部量测在门外
p.assoc_method = 'gnn';
ms_far = make_meas([-2.0, 2.0], [1, 0], {[], []});
[beta, u] = associate_measurements({0.0}, {0.01}, ms_far, {make_track(1, params)}, p);
assert(abs(beta{1}(1) - 1) < 1e-12);
assert(all(u));
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 全部门外\n', k);

%% 4. 单航迹: PDA 与 JPDA 一致
p = params; p.assoc_method = 'pda';
ms2 = make_meas([0.1, 0.2], [1, 0], {[], []});
[b_pda, ~] = associate_measurements({0.1}, {0.01}, ms2, {make_track(1, params)}, p);
p.assoc_method = 'jpda';
[b_jpda, ~, info_j] = associate_measurements({0.1}, {0.01}, ms2, {make_track(1, params)}, p);
assert(max(abs(b_pda{1} - b_jpda{1})) < 1e-12);
assert(info_j.n_events > 0 && ~info_j.fallback);
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 单航迹 PDA == JPDA\n', k);

%% 5. GNN 互斥 / NN 复用
ms_shared = make_meas([0.10], [1], {[]});
tracks2 = {make_track(1, params), make_track(2, params)};
p.assoc_method = 'gnn';
[b_gnn, ~] = associate_measurements({0.10, 0.10}, {0.01, 0.01}, ms_shared, tracks2, p);
assert(abs(sum(cellfun(@(b) b(1), b_gnn)) - 1) < 1e-12);
assert(abs(sum(cellfun(@(b) sum(b(2:end)), b_gnn)) - 1) < 1e-12);

p.assoc_method = 'nn';
[b_nn, ~] = associate_measurements({0.10, 0.10}, {0.01, 0.01}, ms_shared, tracks2, p);
assert(abs(sum(cellfun(@(b) sum(b(2:end)), b_nn)) - 2) < 1e-12);
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] GNN 互斥 / NN 复用\n', k);

%% 6. 特征因子判别
pf = init_parameters(struct('use_features', true));
rng(3);
spec_t = randn(pf.N_features, 1); spec_t = spec_t / norm(spec_t);
spec_a = spec_t + 0.01 * randn(pf.N_features, 1);
spec_b = randn(pf.N_features, 1);
tr_f = make_track(1, pf);
tr_f.spectrum = spec_t;
tr_f.feature_init = true;
ms_f = make_meas([0.10, 0.10], [1, 0], {spec_a, spec_b});
for mi = 1:numel(methods)
    p = pf; p.assoc_method = methods{mi};
    [b, ~] = associate_measurements({0.10}, {0.0025}, ms_f, {tr_f}, p);
    assert(b{1}(2) > b{1}(3));
end
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 特征因子判别\n', k);

%% 7. 确定性
p.assoc_method = 'gnn';
[b1, u1] = associate_measurements(z, S, ms, tracks, p);
[b2, u2] = associate_measurements(z, S, ms, tracks, p);
assert(isequal(b1, b2) && isequal(u1, u2));
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 确定性\n', k);

%% 8. JPDA 事件上限回退
p = params; p.assoc_method = 'jpda';
p.assoc_jpda_max_events = 50;
nt = 8; nm = 20;
z_big = repmat({0}, 1, nt);
S_big = repmat({1}, 1, nt);
ms_big = make_meas(linspace(-2.5, 2.5, nm), zeros(1, nm), repmat({[]}, 1, nm));
tracks_big = cell(1, nt);
for t = 1:nt
    tracks_big{t} = make_track(t, params);
end
[b_fb, ~, info_fb] = associate_measurements(z_big, S_big, ms_big, tracks_big, p);
assert(info_fb.fallback);
p.assoc_method = 'pda';
[b_pda_big, ~] = associate_measurements(z_big, S_big, ms_big, tracks_big, p);
for t = 1:nt
    assert(max(abs(b_fb{t} - b_pda_big{t})) < 1e-12);
end
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] JPDA 事件上限回退\n', k);

%% 9. GOSPA 回归 (匈牙利抽取前后)
rng(7);
X = rand(2, 3) * 100;
Y = rand(2, 4) * 100;
[d, c] = gospar(X, Y, 2, 100, 2);
assert(abs(d - 106.05412966196911) < 1e-9);
assert(abs(c.loc - 79.040992012738286) < 1e-9);
assert(abs(c.false - 70.710678118654755) < 1e-9);
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] GOSPA 回归\n', k);

%% 10. 默认路径 run_simulation 回归
rng(42);
p = init_parameters(struct('scenario', 'maneuver', 'use_features', false));
[pr, vr, g] = run_simulation(p);
assert(abs(g - 37.976556817217748) < 1e-9);
assert(abs(pr(1) - 7.9212781526990055) < 1e-9);
assert(abs(pr(2) - 9.2124480046336092) < 1e-9);
assert(abs(vr(1) - 1.5670254458956132) < 1e-9);
assert(abs(vr(2) - 2.2117032573867759) < 1e-9);
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 默认 JPDA 路径回归\n', k);

%% 11. 路径自举 + GNN 冒烟
rmpath(here);
rng(1);
p = init_parameters(struct('scenario', 'maneuver', 'use_features', false, ...
    'assoc_method', 'gnn'));
[~, ~, g_gnn] = run_simulation(p);
assert(isfinite(g_gnn));
assert(contains(path, here));
addpath(here);
n_pass = n_pass + 1; k = k + 1;
fprintf('  [%d] 路径自举 + GNN 冒烟\n', k);

fprintf('\n全部 %d 项测试通过。\n', k);

% ========================================================================
function tr = make_track(id, params)
tr = struct();
tr.x = zeros(7, 1);
tr.P = eye(7);
tr.id = id;
tr.imm = imm_init(tr.x, tr.P, params);
tr.spectrum = [];
tr.feature_init = false;
tr.confirmed = true;
tr.hit_count = 3;
tr.miss_count = 0;
end

function ms = make_meas(bearings, srcs, spectra)
ms = struct('bearing', {}, 'spectrum', {}, 'src', {});
for i = 1:numel(bearings)
    m.bearing = bearings(i);
    m.spectrum = spectra{i};
    m.src = srcs(i);
    ms(end + 1) = m; %#ok<AGROW>
end
end

function check_beta(beta, n, m, label)
assert(numel(beta) == n, '%s: beta 数量错误', label);
for t = 1:n
    b = beta{t};
    assert(numel(b) == m + 1, '%s: beta 长度错误', label);
    assert(all(b >= -1e-12), '%s: beta 出现负值', label);
    assert(abs(sum(b) - 1) < 1e-12, '%s: beta 未归一', label);
end
end
