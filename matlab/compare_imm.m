% compare_imm.m — 单模型 CV 与 IMM (CV+CT) 的蒙特卡洛对比
%
% 在相同的量测序列下 (同一随机种子) 分别用单模型 CV 与 IMM 滤波,
% 比较位置/速度 RMSE 与 GOSPA (含分解分量), 检验机动条件下的跟踪性能.
%
% 运行: >> compare_imm
% 参数: N_MC 蒙特卡洛次数; USE_FEATURES 是否启用频谱特征辅助
%       机动段由 params.target_turns 控制 (可在 init_parameters.m 修改)

clear; close all;

USE_FEATURES = true;                    % 与 init_parameters 默认一致 (论文主场景)
SCENARIO = 'maneuver';                  % 'maneuver' | 'parallel' | 'crossing' | 'speed'
N_MC = 20;                              % 蒙特卡洛次数

names = {'CV', 'IMM (CV+CT)', 'IMM (CV+CT+CA)'};
model_sets = {{'CV'}, {'CV', 'CT'}, {'CV', 'CT', 'CA'}};
n_modes = numel(names);

%% 蒙特卡洛
fprintf('模型集对比 (场景=%s, N=%d, 特征辅助=%d)...\n', ...
    SCENARIO, N_MC, USE_FEATURES);

pos = NaN(N_MC, 2, n_modes);            % 位置 RMSE [N_MC x 目标 x 模式]
vel = NaN(N_MC, 2, n_modes);
gos_all = NaN(N_MC, n_modes);
gos_loc = NaN(N_MC, n_modes);
gos_miss = NaN(N_MC, n_modes);
gos_false = NaN(N_MC, n_modes);
n_fail = zeros(1, n_modes);

for m = 1:n_modes
    fprintf('  模式 %d/%d: %s [', m, n_modes, names{m});
    for mc = 1:N_MC
        rng(mc);
        params = init_parameters(struct('imm_models', {model_sets{m}}, ...
            'use_features', USE_FEATURES, 'scenario', SCENARIO));
        try
            [pr, vr, ~, sim] = run_simulation(params);
            pos(mc, :, m) = pr';
            vel(mc, :, m) = vr';
            gos_all(mc, m) = sim.gos.mean;
            gos_loc(mc, m) = sim.gos.mean_loc;
            gos_miss(mc, m) = sim.gos.mean_missed;
            gos_false(mc, m) = sim.gos.mean_false;
        catch
            n_fail(m) = n_fail(m) + 1;
        end
        if mod(mc, 5) == 0, fprintf('.'); end
    end
    fprintf('] 失败 %d 次\n', n_fail(m));
end

%% 汇总表
fprintf('\n================= 指标汇总 (均值 ± 标准差, N=%d) =================\n', N_MC);
fprintf('%-24s', '指标');
for m = 1:n_modes
    fprintf(' %-20s', names{m});
end
fprintf('\n%s\n', repmat('-', 1, 24 + 21 * n_modes));

rows = {
    '位置 RMSE 目标1 (m)',      squeeze(pos(:, 1, :));
    '位置 RMSE 目标2 (m)',      squeeze(pos(:, 2, :));
    '速度 RMSE 目标1 (m/s)',    squeeze(vel(:, 1, :));
    '速度 RMSE 目标2 (m/s)',    squeeze(vel(:, 2, :));
    'GOSPA (m)',                gos_all;
    '  |- 定位误差分量 (m)',    gos_loc;
    '  |- 漏检目标分量 (m)',    gos_miss;
    '  |- 虚假航迹分量 (m)',    gos_false
};

for r = 1:size(rows, 1)
    v = rows{r, 2};
    fprintf('%-24s', rows{r, 1});
    for m = 1:n_modes
        fprintf(' %-20s', sprintf('%.2f ± %.2f', ...
            mean(v(:, m), 'omitnan'), std(v(:, m), 'omitnan')));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('=', 1, 24 + 21 * n_modes));

%% 柱状图对比
figure('Name', '模型集对比', 'Position', [100, 100, 1000, 400]);

subplot(1, 2, 1);
data = [mean(gos_all, 'omitnan')', ...
        mean(gos_loc, 'omitnan')', ...
        mean(gos_miss, 'omitnan')', ...
        mean(gos_false, 'omitnan')'];
bar(data);
set(gca, 'XTickLabel', names, 'XTickLabelRotation', 12);
ylabel('GOSPA (m)');
legend({'GOSPA', '定位误差', '漏检', '虚假航迹'}, 'Location', 'best');
title('GOSPA 及分解分量'); grid on;

subplot(1, 2, 2);
d1 = squeeze(mean(pos(:, 1, :), 1, 'omitnan'));
d2 = squeeze(mean(pos(:, 2, :), 1, 'omitnan'));
bar([d1(:), d2(:)]);
set(gca, 'XTickLabel', names, 'XTickLabelRotation', 12);
ylabel('位置 RMSE (m)');
legend({'目标1', '目标2'}, 'Location', 'best');
title('位置 RMSE'); grid on;

fprintf('对比完成。\n');
