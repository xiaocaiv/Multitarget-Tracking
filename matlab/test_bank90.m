function test_bank90()
%TEST_BANK90 audio_bank 完整库模式的接线自检（不跑跟踪，只验取数与特征）
%
%   覆盖：默认池 / 四个池 / 旧 origin 模式 / 取窗越界与循环 / 特征维度 /
%         跨目标相似度是否明显小于 1 / LLR 是否为新标定参数

addpath(fileparts(mfilename('fullpath')));

fprintf('\n================ 1. 默认（library + B_same）================\n');
b = audio_bank();
assert(strcmp(b.source, 'library'), '默认应为 library');
assert(numel(b.names) == 2, 'B_same 应有 2 个目标，实得 %d', numel(b.names));
fprintf('llr 文件      : %s\n', b.llr_name);
fprintf('fs             : %d\n', b.fs);
fprintf('背景池         : %.1f s\n', b.noiseDur);
assert(isfield(b.llr.methods, 'track') && isfield(b.llr.methods, 'init'), ...
    'LLR 缺少 track/init');
fprintf('track.w        : [%.3f %.3f]  b=%.3f\n', ...
    b.llr.methods.track.w(1), b.llr.methods.track.w(2), b.llr.methods.track.b);
fprintf('init.w         : [%.3f %.3f]  b=%.3f\n', ...
    b.llr.methods.init.w(1), b.llr.methods.init.w(2), b.llr.methods.init.b);
% 标定必须与运行时特征同尺度：mu_line 应落在实测 H0~H1 的区间内
assert(b.llr.methods.track.mu(1) > 0.05 && b.llr.methods.track.mu(1) < 0.75, ...
    'track.mu(1)=%.3f 与 MATLAB 特征尺度不符（标定可能用了别的特征定义）', ...
    b.llr.methods.track.mu(1));

fprintf('\n---- 取窗与特征 ----\n');
for i = 1:numel(b.names)
    x = b.get_window(b.names{i}, 0);
    assert(numel(x) == round(b.win * b.fs), '窗长不对');
    [line, band, info] = b.get_feature(b.names{i}, 0);
    assert(numel(line) == 491 && numel(band) == 16, '特征维度不对');
    fprintf('  %-22s 窗 %d 点  线谱 %d 维(范数 %.3f)  带级 %d 维  线数 %d\n', ...
        b.names{i}, numel(x), numel(line), norm(line), numel(band), info.n_lines);
end

fprintf('\n---- 跨目标相似度（应为低值）----\n');
[L1, B1] = b.get_feature(b.names{1}, 10);
[L2, B2] = b.get_feature(b.names{2}, 10);
fprintf('  line cos = %+.3f   band cos = %+.3f\n', L1' * L2, B1' * B2);
assert(L1' * L2 < 0.6, '两个同型客渡的线谱相似度过高，检查是否取到同一段音频');

fprintf('\n---- 同目标不同时刻（应高）----\n');
[L1b, ~] = b.get_feature(b.names{1}, 20);
fprintf('  line cos(t=10 vs t=20) = %+.3f\n', L1' * L1b);

fprintf('\n---- 同目标 vs 跨目标 线谱相似度统计（决定特征是否可用）----\n');
ts = 5:5:60;
sin_ = [];
for i = 1:numel(b.names)
    for k = 1:numel(ts) - 1
        [La, ~] = b.get_feature(b.names{i}, ts(k));
        [Lb, ~] = b.get_feature(b.names{i}, ts(k + 1));
        sin_(end + 1) = La' * Lb; %#ok<AGROW>
    end
end
sout = [];
for k = 1:numel(ts) - 1
    [La, ~] = b.get_feature(b.names{1}, ts(k));
    [Lb, ~] = b.get_feature(b.names{2}, ts(k + 1));
    sout(end + 1) = La' * Lb; %#ok<AGROW>
end
fprintf('  同目标(间隔5s) mean=%.3f sd=%.3f n=%d\n', mean(sin_), std(sin_), numel(sin_));
fprintf('  跨目标(同刻)   mean=%.3f sd=%.3f n=%d\n', mean(sout), std(sout), numel(sout));
dp = (mean(sin_) - mean(sout)) / sqrt(0.5 * (var(sin_) + var(sout)));
fprintf('  对比度 = %.2fx   d'' = %.2f\n', mean(sin_) / mean(sout), dp);

fprintf('\n---- LLR 端到端（feature_llr + 当前 LLR 文件）----\n');
e1 = []; e0 = [];
for i = 1:numel(b.names)
    for k = 1:numel(ts) - 1
        [La, Ba] = b.get_feature(b.names{i}, ts(k));
        [Lb, Bb] = b.get_feature(b.names{i}, ts(k + 1));
        e1(end + 1) = feature_llr(struct('line', La, 'band', Ba), ...
            struct('line', Lb, 'band', Bb), b.llr, 'track'); %#ok<AGROW>
    end
end
for k = 1:numel(ts) - 1
    [La, Ba] = b.get_feature(b.names{1}, ts(k));
    [Lb, Bb] = b.get_feature(b.names{2}, ts(k + 1));
    e0(end + 1) = feature_llr(struct('line', La, 'band', Ba), ...
        struct('line', Lb, 'band', Bb), b.llr, 'track'); %#ok<AGROW>
end
fprintf('  E[ell|H1]=%+.2f  E[ell|H0]=%+.2f  分离=%.2f nat\n', ...
    mean(e1), mean(e0), mean(e1) - mean(e0));
assert(mean(e1) > 0 && mean(e0) < 0, ...
    'LLR 符号不对：H1 应为正、H0 应为负（检查特征与标定是否同源）');
fprintf('  符号正确（H1>0, H0<0）\n');

fprintf('\n---- 噪声取窗 ----\n');
xn = b.get_noise(0, b.win);
fprintf('  noise %d 点（期望 %d）  RMS=%.1f dBFS\n', numel(xn), ...
    round(b.win * b.fs), 20 * log10(sqrt(mean(xn .^ 2)) + 1e-20));

fprintf('\n================ 2. 四个池 ================\n');
for p = {'A_hetero', 'B_same', 'C_same3', 'D_hetero3'}
    bp = audio_bank(struct('pool', p{1}, 'verbose', false));
    fprintf('  %-10s 目标 %d 个：%s\n', p{1}, numel(bp.names), ...
        strjoin(bp.names, ', '));
    for i = 1:numel(bp.names)
        x = bp.get_window(bp.names{i}, 0);
        assert(numel(x) == round(bp.win * bp.fs), '%s 取窗失败', bp.names{i});
    end
end

fprintf('\n================ 3. 旧 origin 模式 ================\n');
bo = audio_bank('source', 'origin', 'verbose', false);
fprintf('  目标 %d 个：%s\n', numel(bo.names), strjoin(bo.names, ', '));
fprintf('  llr 文件: %s\n', bo.llr_name);
assert(numel(bo.tpl.line) > 0, 'origin 模式应带离线模板');
xo = bo.get_window(bo.names{1}, 0);
fprintf('  取窗 %d 点\n', numel(xo));

fprintf('\n================ 4. 越界与循环 ================\n');
[~, t0a] = b.get_window(b.names{1}, 0);
[~, t0b] = b.get_window(b.names{1}, 1e5);
fprintf('  t=0 -> %.1f s,  t=1e5 -> %.1f s （循环生效，均未越界）\n', t0a, t0b);

fprintf('\n全部断言通过。\n');
end
