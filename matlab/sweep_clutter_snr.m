function R = sweep_clutter_snr(nSeed, T_sim)
%SWEEP_CLUTTER_SNR  杂波密度 × 接收 SNR 的性能扫描（真实音频量测）
%
%   R = sweep_clutter_snr()          4 个种子，T=45 s
%   R = sweep_clutter_snr(10, 60)
%
%   网格：lambda_c ∈ {0.5, 1, 2, 4}（每帧每站）
%         SNR(dB)  ∈ {[18 12], [12 6], [6 0]}（两站宽带 SNR）
%   两种配置：纯几何 / 声纹辅助
%   成功判据：两目标位置 RMSE < 50 m
%
%   输出 R 为结构体数组；同时保存 fig5_clutter_snr.png 与 sweep_clutter_snr.mat

if nargin < 1, nSeed = 4; end
if nargin < 2, T_sim = 45; end

lams = [0.5, 1, 2, 4];
snrs = {[18 12], [12 6], [6 0]};
sets = {'geometry', 'feature'};

fprintf('===== 杂波 × SNR 扫描：%d 种子，T=%d s =====\n', nSeed, T_sim);
bank = audio_bank('source', 'origin', 'verbose', false);

R = struct('lam', {}, 'si', {}, 'setting', {}, 'seed', {}, 'ok', {}, ...
           'rmse1', {}, 'rmse2', {}, 'gos', {}, 'nfalse', {});

for li = 1:numel(lams)
    for si = 1:numel(snrs)
        for ti = 1:numel(sets)
            isFeat = strcmp(sets{ti}, 'feature');
            nOK = 0;
            for sd = 1:nSeed
                rng(1000 * sd + 10 * li + si);

                over = struct('scenario', 'speed', 'mobile_platforms', false, ...
                    'T_sim', T_sim, 'K', T_sim, 'use_audio', true, ...
                    'use_features', isFeat, 'lambda_c', lams(li), ...
                    'assoc_method', 'jpda');
                params = init_parameters(over);
                params.audio.snr_db = snrs{si};
                params.audio_bank   = bank;
                params.llr          = bank.llr;
                params.k            = 1;
                params.max_new_tracks = 2;
                for t = 1:numel(params.target_specs)
                    params.target_specs(t).wav_file = '';
                end

                truth = generate_truth(params);
                res = run_case(params, truth, bank);

                ok = all(isfinite(res.rmse)) && all(res.rmse < 50);
                nOK = nOK + ok;
                R(end + 1) = struct('lam', lams(li), 'si', si, ...
                    'setting', sets{ti}, 'seed', sd, 'ok', ok, ...
                    'rmse1', res.rmse(1), 'rmse2', res.rmse(2), ...
                    'gos', res.gos.mean, 'nfalse', res.nfalse); %#ok<AGROW>
            end
            fprintf('  lam=%.1f snr=[%2d %2d] %-8s 成功 %d/%d\n', ...
                lams(li), snrs{si}(1), snrs{si}(2), sets{ti}, nOK, nSeed);
        end
    end
end

%% ---------- 汇总表 ----------
fprintf('\n%-6s %-12s %-10s %8s %10s\n', 'lam', 'SNR', '配置', '成功率', 'GOSPA');
for li = 1:numel(lams)
    for si = 1:numel(snrs)
        for ti = 1:numel(sets)
            idx = [R.lam] == lams(li) & [R.si] == si & strcmp({R.setting}, sets{ti});
            k = sum([R(idx).ok]); n = sum(idx);
            gm = mean([R(idx).gos]);
            fprintf('%-6.1f [%2d %2d]      %-10s %5.0f%% %10.1f\n', ...
                lams(li), snrs{si}(1), snrs{si}(2), sets{ti}, 100 * k / max(n, 1), gm);
        end
    end
end

%% ---------- 画图 ----------
fig = figure('Visible', 'off', 'Position', [100 100 1200 480]);

% (a) 成功率 vs 杂波密度
subplot(1, 2, 1); hold on; grid on;
cols = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19];
for si = 1:numel(snrs)
    for ti = 1:numel(sets)
        y = zeros(1, numel(lams));
        for li = 1:numel(lams)
            idx = [R.lam] == lams(li) & [R.si] == si & strcmp({R.setting}, sets{ti});
            y(li) = 100 * sum([R(idx).ok]) / max(sum(idx), 1);
        end
        ls = '-'; if ti == 1, ls = '--'; end
        plot(lams, y, [ls 'o'], 'Color', cols(si, :), 'LineWidth', 1.4, ...
            'MarkerSize', 5, 'MarkerFaceColor', cols(si, :));
    end
end
xlabel('杂波密度 \lambda_c（每帧每站）'); ylabel('成功率 (%)');
title('成功率 vs 杂波密度（实线=声纹，虚线=纯几何）');
legend(arrayfun(@(s) sprintf('SNR [%d %d]', s{1}(1), s{1}(2)), snrs, ...
    'UniformOutput', false), 'Location', 'southwest', 'FontSize', 8);
ylim([0 105]);

% (b) GOSPA vs 杂波密度
subplot(1, 2, 2); hold on; grid on;
for si = 1:numel(snrs)
    for ti = 1:numel(sets)
        y = zeros(1, numel(lams));
        for li = 1:numel(lams)
            idx = [R.lam] == lams(li) & [R.si] == si & strcmp({R.setting}, sets{ti});
            y(li) = mean([R(idx).gos]);
        end
        ls = '-'; if ti == 1, ls = '--'; end
        plot(lams, y, [ls 's'], 'Color', cols(si, :), 'LineWidth', 1.4, ...
            'MarkerSize', 5, 'MarkerFaceColor', cols(si, :));
    end
end
xlabel('杂波密度 \lambda_c（每帧每站）'); ylabel('GOSPA 均值 (m)');
title('GOSPA vs 杂波密度');
legend(arrayfun(@(s) sprintf('SNR [%d %d]', s{1}(1), s{1}(2)), snrs, ...
    'UniformOutput', false), 'Location', 'northwest', 'FontSize', 8);

outdir = 'D:\AI_Agent\05_目标噪声特征';
if ~isfolder(outdir), outdir = pwd; end
png = fullfile(outdir, 'fig5_clutter_snr.png');
saveas(fig, png);
fprintf('\n已保存 %s\n', png);

save(fullfile(outdir, 'sweep_clutter_snr.mat'), 'R', 'lams', 'snrs', 'sets');
fprintf('已保存 sweep_clutter_snr.mat\n');
end
