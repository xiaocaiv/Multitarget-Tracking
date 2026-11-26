function S = sweep_audio_03(nSeed, T_sim, scen)
%SWEEP_AUDIO_03  多随机种子对比：纯几何 vs 声纹辅助（真实音频量测）
%
%   S = sweep_audio_03()              6 个种子，60 s，scenario='speed'
%   S = sweep_audio_03(20, 60, 'crossing')
%
%   两种配置除 use_features 外完全相同，因此差异可归因于声纹。
%   成功判据：两个目标的位置 RMSE 都 < 50 m（即没有锁到幽灵上）。
%
%   输出 S 为结构体数组，字段：setting / rmse / missed / nfalse / gos /
%   ok / succ_rate / ci

if nargin < 1, nSeed = 50; end
if nargin < 2, T_sim = 60; end
if nargin < 3, scen = 'speed'; end

fprintf('===== 真实音频闭环对比：%d 个种子，T=%d s，scenario=%s =====\n', ...
    nSeed, T_sim, scen);

% 声纹库只加载一次
fprintf('加载声纹库...\n');
bank = audio_bank('source', 'origin', 'verbose', false);

settings = {'geometry', 'feature'};
S = struct('setting', {}, 'seed', {}, 'rmse', {}, 'missed', {}, ...
           'nfalse', {}, 'gos', {}, 'ok', {});

for si = 1:numel(settings)
    isFeat = strcmp(settings{si}, 'feature');
    fprintf('\n--- %s（use_features=%d）---\n', settings{si}, isFeat);
    fprintf('%5s %9s %9s %7s %7s %10s  %s\n', ...
        'seed', 'RMSE1', 'RMSE2', 'miss', 'false', 'GOSPA', '判定');

    for sd = 1:nSeed
        rng(sd);

        over = struct('scenario', scen, 'mobile_platforms', false, ...
            'T_sim', T_sim, 'K', T_sim, 'use_audio', true, ...
            'use_features', isFeat, 'assoc_method', 'jpda');
        params = init_parameters(over);

        params.use_audio   = true;
        params.use_features = isFeat;
        params.audio_bank  = bank;
        params.llr         = bank.llr;
        params.k           = 1;
        params.max_new_tracks = 2;      % 2 目标，一次建满便于看配对
        for t = 1:numel(params.target_specs)
            params.target_specs(t).wav_file = '';   % 不需要合成频谱
        end

        truth = generate_truth(params);
        res = run_case(params, truth, bank);

        ok = all(isfinite(res.rmse)) && all(res.rmse < 50);
        fprintf('%5d %9.1f %9.1f %7.2f %7.2f %10.1f  %s\n', ...
            sd, res.rmse(1), res.rmse(2), res.missed, res.nfalse, ...
            res.gos.mean, ternary(ok, 'OK', 'FAIL'));

        S(end + 1) = struct('setting', settings{si}, 'seed', sd, ...
            'rmse', res.rmse, 'missed', res.missed, 'nfalse', res.nfalse, ...
            'gos', res.gos.mean, 'ok', ok); %#ok<AGROW>
    end
end

%% ---------- 汇总 ----------
fprintf('\n===== 汇总（成功判据：两目标 RMSE<50 m）=====\n');
fprintf('%-10s %8s %10s %12s %12s %12s\n', ...
    '配置', '成功', '成功率', '95% CI', '成功RMSE均值', 'GOSPA均值');
for si = 1:numel(settings)
    idx = strcmp({S.setting}, settings{si});
    okv = [S(idx).ok];
    k = sum(okv); n = numel(okv);
    [lo, hi] = wilson(k, n, 1.96);
    rm = reshape([S(idx).rmse], 2, [])';
    if any(okv)
        rmean  = mean(rm(okv, 1), 'omitnan');
        rmean2 = mean(rm(okv, 2), 'omitnan');
    else
        rmean = NaN; rmean2 = NaN;
    end
    gmean = mean([S(idx).gos]);
    fprintf('%-10s %5d/%-2d %9.0f%% [%4.0f%%,%4.0f%%]  %.1f/%.1f m  %10.1f m\n', ...
        settings{si}, k, n, 100 * k / n, 100 * lo, 100 * hi, ...
        rmean, rmean2, gmean);
end

% 配对比较：同一 seed 下 feature 是否优于 geometry
fprintf('\n===== 逐种子配对对比（RMSE1+RMSE2）=====\n');
fprintf('%5s %18s %18s %8s\n', 'seed', 'geometry', 'feature', '改善');
for sd = 1:nSeed
    ig = find(strcmp({S.setting}, 'geometry') & [S.seed] == sd, 1);
    if_ = find(strcmp({S.setting}, 'feature') & [S.seed] == sd, 1);
    if isempty(ig) || isempty(if_), continue; end
    g = sum(S(ig).rmse); f = sum(S(if_).rmse);
    fprintf('%5d %18.1f %18.1f %8.1f\n', sd, g, f, g - f);
end
end

% ========================================================================
function [lo, hi] = wilson(k, n, z)
% Wilson 置信区间（比例）
if n == 0, lo = NaN; hi = NaN; return; end
p = k / n;
den = 1 + z ^ 2 / n;
c = p + z ^ 2 / (2 * n);
r = z * sqrt(p * (1 - p) / n + z ^ 2 / (4 * n ^ 2));
lo = max(0, (c - r) / den);
hi = min(1, (c + r) / den);
end

function y = ternary(c, a, b)
if c, y = a; else, y = b; end
end
