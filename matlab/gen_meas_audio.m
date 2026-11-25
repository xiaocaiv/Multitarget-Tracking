function [z1, z2] = gen_meas_audio(bank, truth, params, k)
%GEN_MEAS_AUDIO  生成带真实声纹的双站量测（含杂波）
%
%   [z1, z2] = gen_meas_audio(bank, truth, params, k)
%
%   与 generate_measurements.m 的区别
%     特征不再由"模板 + 高斯噪声"合成，而是从真实音频经差算子后重算：
%       取测试段 5 s 窗 → 多普勒频率缩放 → 站间平滑频响 → 混入降雨噪声到指定 SNR
%       → LOFAR 线谱残差 + 1/3 倍频程带级
%
%   输入
%     bank    audio_bank 返回的结构体
%     truth   需含 .num_targets 与 .X{tgt}（前 4 维为 [x;y;vx;vy]，可为 5/7 维）
%     params  需含 .dt .Pd .sigma_bearing .lambda_c .clutter_range
%             .plat1_pos .plat2_pos .use_features
%             .audio.vessel   1×num_targets 元胞，指定每个目标用哪段录音
%             .audio.snr_db   1×2，两站的接收 SNR (dB)
%             .audio.eq_node_db  2×3，站间平滑频响节点 (dB)
%             .audio.eq_node_f   1×3，节点频率 (Hz)
%
%   输出（结构体数组，可能为空 0x0）
%     .bearing  方位 (rad)
%     .spectrum []（与旧量测字段对齐，音频链路不用）
%     .line     491x1 线谱残差（L2 归一化）
%     .band     16x1  1/3 倍频程（去均值 + 归一化）
%     .q        质量分 [0,1]
%     .snr      线谱峰值高出局部噪底 (dB)
%     .src      目标编号，杂波为 0（与 generate_measurements 的 .src 一致）

if nargin < 4, k = 1; end
a = audio_cfg(params, truth.num_targets);

z1 = empty_meas();
z2 = empty_meas();
tsec = (k - 1) * params.dt;
c = 1500;                                  % 声速 (m/s)

%% ---------- 目标量测 ----------
for tgt = 1:truth.num_targets
    X = truth.X{tgt};
    p = X(1:2, k);
    v = X(3:4, k);
    vessel = a.vessel{tgt};

    for s = 1:2
        if s == 1, ps = params.plat1_pos(:); else, ps = params.plat2_pos(:); end
        d  = p - ps;
        r  = norm(d) + 1e-9;
        th = atan2(d(2), d(1));
        rdot  = dot(v, d) / r;
        alpha = 1 - rdot / c;

        if rand() >= params.Pd, continue; end

        m = struct();
        m.bearing  = wrap_pi(th + params.sigma_bearing * randn());
        m.spectrum = [];
        m.src      = tgt;

        if params.use_features
            xw = bank.get_window(vessel, tsec);
            y  = apply_diff_operator(xw, bank.fs, alpha, s, a, bank, tsec);
            [ln, bd, inf1] = bank.features(y);
            m.line = ln;
            m.band = bd;
            m.snr  = inf1.snr;
            m.q    = quality(inf1);
        else
            m.line = []; m.band = []; m.snr = NaN; m.q = 0;
        end

        if s == 1, z1(end + 1) = m; else, z2(end + 1) = m; end %#ok<AGROW>
    end
end

%% ---------- 杂波 ----------
for s = 1:2
    nc = poissrnd(params.lambda_c);
    for j = 1:nc
        m = struct();
        m.bearing  = params.clutter_range(1) + ...
            (params.clutter_range(2) - params.clutter_range(1)) * rand();
        m.spectrum = [];
        m.src      = 0;
        if params.use_features
            xn = bank.get_noise(rand() * max(bank.noiseDur - bank.win, 0.1), bank.win);
            [ln, bd, inf1] = bank.features(xn);
            m.line = ln;
            m.band = bd;
            m.snr  = inf1.snr;
            m.q    = quality(inf1) * 0.5;      % 杂波质量分人为压低
        else
            m.line = []; m.band = []; m.snr = NaN; m.q = 0;
        end
        if s == 1, z1(end + 1) = m; else, z2(end + 1) = m; end %#ok<AGROW>
    end
end
end

%% ================= 局部函数 =================

function a = audio_cfg(params, nT)
a = struct();
a.vessel = {'MarDeCangas', 'Zodiac'};
a.snr_db = [12, 6];
a.eq_node_db = zeros(2, 3);
a.eq_node_f  = [10, 100, 500];
if isfield(params, 'audio')
    fn = fieldnames(params.audio);
    for i = 1:numel(fn), a.(fn{i}) = params.audio.(fn{i}); end
end
if numel(a.vessel) < nT
    error('params.audio.vessel 只有 %d 项，但目标有 %d 个', numel(a.vessel), nT);
end
end

function z = empty_meas()
z = struct('bearing', {}, 'spectrum', {}, 'line', {}, 'band', {}, ...
           'q', {}, 'snr', {}, 'src', {});
end

function th = wrap_pi(th)
th = mod(th + pi, 2 * pi) - pi;
end

function y = apply_diff_operator(x, fs, alpha, s, a, bank, tsec)
% 差算子：频率缩放（多普勒）+ 站间平滑频响 + 混噪
x = x(:);

% 1) 多普勒：y[n] = x(alpha*n)，输出长度与输入一致（保证频率栅格不变）
n   = numel(x);
pos = 1 + (0:n-1)' * alpha;
pos = min(max(pos, 1), n);
y   = interp1((1:n)', x, pos, 'linear');

% 2) 站间平滑频响
eq = a.eq_node_db(min(s, size(a.eq_node_db, 1)), :);
y  = apply_eq(y, fs, a.eq_node_f, eq);

% 3) 混入降雨噪声至指定 SNR
snr_db = a.snr_db(min(s, numel(a.snr_db)));
ps = mean(y .^ 2) + 1e-20;
nv = bank.get_noise(tsec, n / fs);
nv = nv(:) - mean(nv);
nv = nv / (sqrt(mean(nv .^ 2)) + 1e-20);
nv = nv * sqrt(ps / (10 ^ (snr_db / 10)));
y  = y + nv;
end

function y = apply_eq(x, fs, fnode, dnode)
n  = numel(x);
X  = fft(x);
f  = (0:n-1)' * fs / n;
f  = min(f, fs - f);                       % 折叠到 [0, fs/2]
lf = log10(max(f, 1));
lo = log10(fnode(1)); hi = log10(fnode(end));
lf = min(max(lf, lo), hi);                 % 两端平延
eq = interp1(log10(fnode(:)), dnode(:), lf, 'linear');
y  = real(ifft(X .* 10 .^ (eq / 20)));
end

function q = quality(info)
q = 1 / (1 + exp(-((info.n_lines - 3) / 2 + (info.snr - 6) / 4)));
q = min(max(q, 0.05), 1);
end
