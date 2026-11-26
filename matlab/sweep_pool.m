function R = sweep_pool(nSeed, T_sim, pools, scens, opts)
%SWEEP_POOL  P2 闭环：按目标池 × 场景 × 配置做配对蒙特卡洛
%
%   R = sweep_pool()                                    % 10 种子 / 60 s / A+B / speed
%   R = sweep_pool(20, 60, {'B_same'}, {'speed','crossing'})
%   R = sweep_pool(20, 60, {'A_hetero','B_same'}, {'speed'}, struct('save','x.mat'))
%
%   每个 (池, 场景, 种子) 下跑两种配置（geometry / feature）。同一种子下
%   量测序列与关联算法完全相同，因此差异可归因于声纹。
%
%   目标池取自 pools.json（见 scripts/pools.py）：
%     A_hetero  异型：Mar de Cangas + Zodiac
%     B_same    同型：Mar de Cangas + Mar de Onza
%     C_same3 / D_hetero3 为三目标，需要先把 scenario_defs 扩到 3 目标
%
%   成功判据：所有目标的位置 RMSE < 50 m（即没有锁到鬼影上）
%
%   输出 R 为结构体数组，字段：pool / scen / setting / seed / rmse / missed /
%   nfalse / gos / ok

if nargin < 1, nSeed  = 10;   end
if nargin < 2, T_sim  = 60;   end
if nargin < 3, pools  = {'A_hetero', 'B_same'}; end
if nargin < 4, scens  = {'speed'}; end
if nargin < 5, opts   = struct(); end
if ~isfield(opts, 'assoc'), opts.assoc = 'jpda'; end
if ~isfield(opts, 'snr_db'), opts.snr_db = [12, 6]; end
if ~isfield(opts, 'save'), opts.save = sprintf('sweep_pool_%dseeds.mat', nSeed); end

settings = {'geometry', 'feature'};
R = struct('pool', {}, 'scen', {}, 'setting', {}, 'seed', {}, ...
           'rmse', {}, 'missed', {}, 'nfalse', {}, 'gos', {}, 'ok', {});

fprintf('===== P2 闭环配对扫描：%d 种子，T=%d s，关联=%s =====\n', ...
    nSeed, T_sim, opts.assoc);

for pi_ = 1:numel(pools)
    poolName = pools{pi_};
    bank = audio_bank(struct('pool', poolName, 'verbose', false));
    nT = numel(bank.names);
    fprintf('\n########## 池 %s：%d 个目标（%s）##########\n', ...
        poolName, nT, strjoin(bank.names, ', '));
    if nT > 2
        error(['池 %s 有 %d 个目标。scenario_defs 目前只定义 2 目标，' ...
               '请先把场景扩到 3 目标（计划 §6 E10）再用此池。'], poolName, nT);
    end

    for si = 1:numel(scens)
        scen = scens{si};
        fprintf('\n--- 池 %s / 场景 %s ---\n', poolName, scen);
        fprintf('%5s %9s %9s %7s %7s %10s  %-8s %-8s\n', ...
            'seed', 'RMSE1', 'RMSE2', 'miss', 'false', 'GOSPA', '几何', '声纹');

        for sd = 1:nSeed
            pair_ = struct();
            for ti = 1:numel(settings)
                isFeat = strcmp(settings{ti}, 'feature');
                rng(sd);                                  % 两种配置共用同一随机流

                over = struct('scenario', scen, 'mobile_platforms', false, ...
                    'T_sim', T_sim, 'K', T_sim, 'use_audio', true, ...
                    'use_features', isFeat, 'assoc_method', opts.assoc);
                params = init_parameters(over);
                params.use_audio    = true;
                params.use_features = isFeat;
                params.audio_bank   = bank;
                params.llr          = bank.llr;
                params.k            = 1;
                params.max_new_tracks = nT;
                params.audio.vessel = bank.names;         % 目标顺序 = 池顺序
                params.audio.snr_db = opts.snr_db;
                for t = 1:numel(params.target_specs)
                    params.target_specs(t).wav_file = '';
                end

                truth = generate_truth(params);
                res = run_case(params, truth, bank);

                ok = all(isfinite(res.rmse)) && all(res.rmse < 50);
                pair_.(settings{ti}) = struct('rmse', res.rmse, ...
                    'missed', res.missed, 'nfalse', res.nfalse, ...
                    'gos', res.gos.mean, 'ok', ok);
                R(end + 1) = struct('pool', poolName, 'scen', scen, ...
                    'setting', settings{ti}, 'seed', sd, ...
                    'rmse', res.rmse, 'missed', res.missed, ...
                    'nfalse', res.nfalse, 'gos', res.gos.mean, ...
                    'ok', ok); %#ok<AGROW>
            end
            fprintf('%5d %9.1f %9.1f %7.2f %7.2f %10.1f  %-8s %-8s\n', sd, ...
                pair_.geometry.rmse(1), pair_.geometry.rmse(2), ...
                pair_.geometry.missed, pair_.geometry.nfalse, ...
                pair_.geometry.gos, tf(pair_.geometry.ok), ...
                tf(pair_.feature.ok));
        end
    end
end

%% ---------- 汇总 ----------
fprintf('\n\n===== 汇总（成功判据：所有目标 RMSE<50 m）=====\n');
fprintf('%-10s %-10s %8s %12s %14s %14s %10s\n', ...
    '池', '场景', '配置', '成功率', '95% CI', '成功RMSE均值', 'GOSPA均值');

keys = unique(strcat({R.pool}, '|', {R.scen}));
for ki = 1:numel(keys)
    kk = split(keys{ki}, '|');
    for ti = 1:numel(settings)
        idx = strcmp({R.pool}, kk{1}) & strcmp({R.scen}, kk{2}) & ...
              strcmp({R.setting}, settings{ti});
        okv = [R(idx).ok];
        k = sum(okv); n = numel(okv);
        [lo, hi] = wilson(k, n, 1.96);
        rm = reshape([R(idx).rmse], 2, [])';
        if any(okv)
            rmean = mean(rm(okv, :), 1, 'omitnan');   % 保持 1×2
        else
            rmean = [NaN, NaN];
        end
        if numel(rmean) < 2, rmean = [rmean(1), NaN]; end
        gmean = mean([R(idx).gos]);
        fprintf('%-10s %-10s %8s %5d/%-2d %5.0f%% [%3.0f%%,%3.0f%%]  %.1f/%.1f m  %10.1f m\n', ...
            kk{1}, kk{2}, settings{ti}, k, n, 100 * k / n, 100 * lo, 100 * hi, ...
            rmean(1), rmean(2), gmean);
    end
    % 配对检验
    b = 0; c = 0;
    for sd = 1:nSeed
        ig = find(strcmp({R.pool}, kk{1}) & strcmp({R.scen}, kk{2}) & ...
                  strcmp({R.setting}, 'geometry') & [R.seed] == sd, 1);
        iff = find(strcmp({R.pool}, kk{1}) & strcmp({R.scen}, kk{2}) & ...
                   strcmp({R.setting}, 'feature') & [R.seed] == sd, 1);
        if isempty(ig) || isempty(iff), continue; end
        if ~R(ig).ok && R(iff).ok, b = b + 1; end
        if R(ig).ok && ~R(iff).ok, c = c + 1; end
    end
    p = mcnemar_exact(b, c);
    fprintf('           配对：几何败/声纹成 %d 次，反向 %d 次，McNemar p=%.3g\n', b, c, p);
end

save(opts.save, 'R', 'pools', 'scens', 'settings', 'nSeed', 'T_sim', 'opts');
fprintf('\n结果已存 %s\n', fullfile(pwd, opts.save));
end

% ========================================================================
function [lo, hi] = wilson(k, n, z)
if n == 0, lo = NaN; hi = NaN; return; end
p = k / n;
den = 1 + z ^ 2 / n;
c = p + z ^ 2 / (2 * n);
r = z * sqrt(p * (1 - p) / n + z ^ 2 / (4 * n ^ 2));
lo = max(0, (c - r) / den);
hi = min(1, (c + r) / den);
end

function p = mcnemar_exact(b, c)
% 精确 McNemar：不一致对数 n=b+c，检验 b 与 c 是否对称
n = b + c;
if n == 0, p = 1; return; end
k = min(b, c);
p = 2 * sum(arrayfun(@(i) nchoosek(n, i), 0:k)) / 2 ^ n;
p = min(1, p);
end

function s = tf(v)
if v, s = 'OK'; else, s = 'FAIL'; end
end
