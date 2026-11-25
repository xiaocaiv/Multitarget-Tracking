function bank = audio_bank(varargin)
%AUDIO_BANK  加载官方样例音频、声纹模板与校准参数，提供取窗与特征提取
%
%   bank = audio_bank()                    使用默认配置
%   bank = audio_bank(cfg)                 用结构体覆盖默认配置
%   bank = audio_bank('win', 5, ...)       用名值对覆盖默认配置
%
%   cfg 字段（均可省略）
%     .base     数据目录（manifest.csv / tpl_*.csv / clr_params.json）
%     .root     音频目录（Shipsear_origin）
%     .win      特征窗长 (s)，默认 5
%     .loop     测试段用尽后是否循环，默认 true
%     .verbose  是否打印加载信息，默认 true
%
%   bank 接口
%     [x, t0]            = bank.get_window(vessel, t)   取场景时刻 t 对应的 5 s 音频
%     [x, t0]            = bank.get_noise(t, dur)       取降雨噪声段
%     [line, band, info] = bank.features(x)             LOFAR 线谱残差 + 1/3 倍频程
%     [line, band, info] = bank.get_feature(vessel, t)  取窗并提特征
%
%   说明
%     只加载各目标"测试段"的音频（时间上与模板段不相交），因此不可能发生泄漏。
%     特征定义与 Python 标定脚本完全一致（5 s 窗、0.2 Hz 栅格、中值 11 bin、
%     插值到 1 Hz 的 491 维线谱残差；10–500 Hz 的 17 个 1/3 倍频程带级）。

if nargin == 1 && isstruct(varargin{1})
    cfg = varargin{1};
elseif nargin == 0
    cfg = struct();
elseif mod(nargin, 2) == 0
    cfg = struct();
    for i = 1:2:nargin
        cfg.(char(varargin{i})) = varargin{i + 1};
    end
else
    error('audio_bank：参数应为结构体，或成对的"名称, 值"');
end
base_default = local_data_dir();
def = struct( ...
    'base',    base_default, ...
    'root',    fullfile(base_default, 'Shipsear_origin'), ...
    'win',     5.0, ...
    'loop',    true, ...
    'verbose', true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
end

%% ---------- 读清单 ----------
mpath = fullfile(cfg.base, 'manifest.csv');
if ~isfile(mpath), error('找不到 %s', mpath); end
T = readtable(mpath, 'Delimiter', ',', 'TextType', 'string');

%% ---------- 读模板与校准参数 ----------
names = string(strtrim(splitlines(string(fileread(fullfile(cfg.base, 'tpl_names.txt'))))));
names = names(names ~= "");
nV = numel(names);
tplLine = readmatrix(fullfile(cfg.base, 'tpl_line.csv'));   % nV x 491
tplBand = readmatrix(fullfile(cfg.base, 'tpl_band.csv'));   % nV x 16
if size(tplLine, 1) ~= nV || size(tplBand, 1) ~= nV
    error('模板行数与 tpl_names.txt 不一致');
end
llr = jsondecode(fileread(fullfile(cfg.base, 'clr_params.json')));

%% ---------- 只加载"测试段"音频 ----------
sig  = cell(nV, 1);
segB = zeros(nV, 2);
fs   = [];
for i = 1:nV
    rows = T(T.vessel == names(i) & T.segment == "test", :);
    if isempty(rows)
        error('目标 %s 没有 test 段记录', names(i));
    end
    afile = fullfile(cfg.root, char(rows.file(1)));
    ai = audioinfo(afile);
    if isempty(fs), fs = ai.SampleRate; end
    if ai.SampleRate ~= fs
        error('采样率不一致：%s 为 %d Hz，前面为 %d Hz', char(names(i)), ai.SampleRate, fs);
    end
    t0 = min(rows.t_start);
    t1 = max(rows.t_end);
    s1 = round(t0 * fs) + 1;
    s2 = min(round(t1 * fs), ai.TotalSamples);
    x  = audioread(afile, [s1 s2]);
    if size(x, 2) > 1, x = mean(x, 2); end
    sig{i}  = x(:);
    segB(i, :) = [t0, t1];
end

%% ---------- 背景噪声池（降雨） ----------
nrows = T(T.vessel == "Rain", :);
if isempty(nrows)
    noise = zeros(round(cfg.win * fs), 1);
    noiseDur = cfg.win;
else
    nfile = fullfile(cfg.root, char(nrows.file(1)));
    ni = audioinfo(nfile);
    if ni.SampleRate ~= fs
        error('降雨文件采样率 %d 与目标 %d 不一致', ni.SampleRate, fs);
    end
    noise = audioread(nfile);
    if size(noise, 2) > 1, noise = mean(noise, 2); end
    noise = noise(:);
    noiseDur = numel(noise) / fs;
end

%% ---------- 组装 ----------
bank = struct();
bank.cfg     = cfg;
bank.fs      = fs;
bank.win     = cfg.win;
bank.names   = names;
bank.tpl     = struct('line', tplLine, 'band', tplBand);
bank.llr     = llr;
bank.seg     = segB;
bank.noise   = noise;
bank.noiseDur = noiseDur;
bank.idx     = containers.Map(cellstr(names), num2cell(1:nV));

bank.get_window   = @get_window_impl;
bank.get_noise    = @get_noise_impl;
bank.features     = @features_impl;
bank.get_feature  = @get_feature_impl;

if cfg.verbose
    fprintf('[audio_bank] fs=%d Hz, win=%.1f s, 目标 %d 个: %s\n', ...
        fs, cfg.win, nV, strjoin(cellstr(names), ', '));
    for i = 1:nV
        fprintf('             %-12s test 段 %6.1f–%6.1f s (%5.1f s)\n', ...
            names(i), segB(i, 1), segB(i, 2), segB(i, 2) - segB(i, 1));
    end
    fprintf('             背景池 %.1f s\n', noiseDur);
end

%% ================= 嵌套函数 =================

    function [x, t0] = get_window_impl(vessel, t)
        % 取场景时刻 t（秒，从 0 开始）对应的 5 s 音频
        key = char(string(vessel));
        if ~isKey(bank.idx, key)
            error('未知目标 "%s"；可选：%s', key, strjoin(cellstr(names), ', '));
        end
        i = bank.idx(key);
        a0 = segB(i, 1);
        a1 = segB(i, 2);
        span = a1 - a0 - cfg.win;
        if span <= 0
            t0 = a0;
        else
            tt = mod(t, span);          % 循环（cfg.loop=true 时）
            t0 = a0 + tt;
            if ~cfg.loop
                t0 = a0 + min(max(t, 0), span);
            end
        end
        i1 = round((t0 - a0) * fs) + 1;          % 相对段首的索引
        i2 = i1 + round(cfg.win * fs) - 1;
        if i1 < 1 || i2 > numel(sig{i})
            error('取窗越界: %s, t=%.2f s (段 %0.1f–%0.1f)', key, t, a0, a1);
        end
        x  = sig{i}(i1:i2);
        if numel(x) < round(cfg.win * fs)
            x = [x; zeros(round(cfg.win * fs) - numel(x), 1)];
        end
    end

    function [x, t0] = get_noise_impl(t, dur)
        if nargin < 2, dur = cfg.win; end
        n = round(dur * fs);
        if numel(noise) <= n
            x = [noise; zeros(n - numel(noise), 1)];
            t0 = 0; return;
        end
        t0 = mod(t, max(noiseDur - dur, 0));
        i1 = round(t0 * fs) + 1;
        x  = noise(i1:min(i1 + n - 1, numel(noise)));
        if numel(x) < n, x = [x; zeros(n - numel(x), 1)]; end
    end

    function [line, band, info] = get_feature_impl(vessel, t)
        x = get_window_impl(vessel, t);
        [line, band, info] = features_impl(x);
    end

    function [line, band, info] = features_impl(x)        % 与 Python 标定脚本一致的实现
        x = x(:);
        n = numel(x);
        nn = round(cfg.win * fs);
        if n < nn, x = [x; zeros(nn - n, 1)]; n = nn; end
        idx = (0:n-1)';
        w = 0.5 * (1 - cos(2 * pi * idx / max(n - 1, 1)));   % 对称 Hann，免工具箱
        X = fft(x .* w);
        nh = floor(n / 2) + 1;
        P = abs(X(1:nh)).^2 / (fs * sum(w.^2));
        fr = (0:nh-1)' * fs / n;
        m  = (fr >= 10) & (fr <= 500);
        fm = fr(m); Pm = P(m);
        PdB = 10 * log10(max(Pm, 1e-30));
        % 中值滤波（窗宽 11 bin ≈ 2.2 Hz，两端复制填充）
        pad = 5;
        ext = [repmat(PdB(1), pad, 1); PdB; repmat(PdB(end), pad, 1)];
        med = zeros(size(PdB));
        for q = 1:numel(PdB)
            med(q) = median(ext(q:q + 2 * pad));
        end
        R = max(PdB - med, 0);
        g = (10:500)';
        line = interp1(fm, R, g, 'linear', 0);   % 端点外插取 0，避免 NaN
        line = line / (norm(line) + 1e-12);
        % 1/3 倍频程
        edges = 10 * 2 .^ ((0:17) / 3);
        edges = edges(edges <= 500);
        nb = numel(edges) - 1;
        B  = zeros(nb, 1);
        for q = 1:nb
            sel = (fm >= edges(q)) & (fm < edges(q + 1));
            if any(sel)
                B(q) = 10 * log10(mean(Pm(sel)) + 1e-30);
            else
                B(q) = -120;
            end
        end
        band = B - mean(B);
        band = band / (norm(band) + 1e-12);
        % 质量指标
        loc = find(R(2:end-1) > R(1:end-2) & R(2:end-1) >= R(3:end) & R(2:end-1) > 6) + 1;
        info = struct('n_lines', numel(loc), 'snr', max(R), 't', []);
    end
end

% ========================================================================
function base = local_data_dir()
%LOCAL_DATA_DIR 定位 05_目标噪声特征 数据目录（相对优先，绝对兜底）
%   工程布局：<root>\03_课题代码\双站机动跟踪多目标\matlab\audio_bank.m
%   数据目录：<root>\05_目标噪声特征
base = '';
p = fileparts(mfilename('fullpath'));      % ...\matlab
for i = 1:3
    p = fileparts(p);                       % ...\双站机动跟踪多目标\03_课题代码\AI_Agent
end
cand = fullfile(p, '05_目标噪声特征');
if isfolder(cand)
    base = cand;
elseif isfolder('D:\AI_Agent\05_目标噪声特征')
    base = 'D:\AI_Agent\05_目标噪声特征';
else
    error('找不到数据目录 05_目标噪声特征（请用 audio_bank(struct(''base'', ...)) 指定）');
end
end
