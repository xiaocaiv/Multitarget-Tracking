function test_audio_03()
%TEST_AUDIO_03  验证 03 工程里的音频量测 + LLR 起批接线
clc; rng(7);

fprintf('===== 载入声纹库 =====\n');
bank = audio_bank('source', 'origin');

fprintf('\n===== 参数与真值 =====\n');
params = init_parameters(struct('use_audio', true, 'mobile_platforms', false, ...
    'T_sim', 60, 'use_features', true, 'assoc_method', 'jpda'));
params.use_audio = true;
params.mobile_platforms = false;
params.use_features = true;
params.audio_bank = bank;
params.llr = bank.llr;
params.max_new_tracks = 2;          % 本测试一次建满两条，便于看配对
params.k = 1;

truth = generate_truth(params);
fprintf('场景 %s，目标数 %d，K=%d\n', params.scenario, params.num_targets, params.K);
for t = 1:params.num_targets
    fprintf('  目标%d 起点 (%.0f, %.0f)\n', t, truth.pos{t}(1,1), truth.pos{t}(2,1));
end

fprintf('\n===== 逐帧量测 + 起批 =====\n');
tracks = {}; next_id = 1;
for k = 1:3
    params.k = k;
    [z1, z2] = gen_meas_audio(bank, truth, params, k);
    fprintf('\n--- k=%d 量测 %d/%d ---\n', k, numel(z1), numel(z2));
    for j = 1:numel(z1)
        fprintf('  站1[%d] src=%d %7.2f° 线谱=%d q=%.2f\n', ...
            j, z1(j).src, rad2deg(z1(j).bearing), numel(z1(j).line), z1(j).q);
    end
    for j = 1:numel(z2)
        fprintf('  站2[%d] src=%d %7.2f° 线谱=%d q=%.2f\n', ...
            j, z2(j).src, rad2deg(z2(j).bearing), numel(z2(j).line), z2(j).q);
    end

    % 2x2 LLR（取真实目标量测）
    i1 = find([z1.src] > 0); i2 = find([z2.src] > 0);
    if numel(i1) >= 2 && numel(i2) >= 2
        M = zeros(numel(i1), numel(i2));
        for a = 1:numel(i1)
            for b = 1:numel(i2)
                M(a, b) = feature_llr(z1(i1(a)), z2(i2(b)), params.llr, 'init');
            end
        end
        fprintf('  LLR 矩阵（行=站1 src %s，列=站2 src %s）\n', ...
            mat2str([z1(i1).src]), mat2str([z2(i2).src]));
        disp(M);
    end

    u1 = true(numel(z1), 1); u2 = true(numel(z2), 1);
    [tracks, next_id] = manage_tracks(tracks, z1, z2, u1, u2, params, next_id);
    fprintf('  → 航迹 %d 条\n', numel(tracks));
end

fprintf('\n===== 航迹 =====\n');
for t = 1:numel(tracks)
    tr = tracks{t};
    fprintf('  航迹%d: (%.1f, %.1f)  src=%d  ell0=%+.2f  q=%.2f  imm模型=%d\n', ...
        tr.id, tr.x(1), tr.x(2), tr.eval_src, tr.ell0, tr.q, numel(tr.imm.mu));
end

% 频谱似然接口自检（新式）
if numel(tracks) >= 2
    meas = struct('line', tracks{1}.tpl_line, 'band', tracks{1}.tpl_band, 'q', 1);
    fprintf('\nspectrum_likelihood 自检: 与自身模板 = %.3f\n', ...
        spectrum_likelihood(tracks{1}, meas, params));
    meas2 = struct('line', tracks{2}.tpl_line, 'band', tracks{2}.tpl_band, 'q', 1);
    fprintf('                          与他人模板 = %.3f\n', ...
        spectrum_likelihood(tracks{1}, meas2, params));
end
end
