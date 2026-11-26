function bank = audio_bank(varargin)
%AUDIO_BANK  加载 ShipEar 音频、目标池与校准 LLR，提供取窗与特征提取
%
%   两种数据源
%     'library'（默认）  完整库 shipsEar_AUDIOS + manifest_90.csv + pools.json
%     'origin'          旧的 4 段官方样例（Shipsear_origin + manifest.csv），用于复现旧结果
%
%   bank = audio_bank()                        完整库，默认目标池
%   bank = audio_bank(struct('pool','B_same')) 指定目标池（A_hetero/B_same/C_same3/D_hetero3）
%   bank = audio_bank('source','origin')       回到旧的 4 文件模式
%
%   cfg 字段（均可省略）
%     .base       数据目录（manifest_90.csv / pools.json / p1_llr.json），默认自动定位
%     .root       音频目录，默认 base/shipsEar_AUDIOS（library）或 base/Shipsear_origin（origin）
%     .source     'library' | 'origin'，默认 'library'
%     .pool       目标池名，默认 'B_same'（同型主线组）；空串表示用库内全部可用目标
%     .win        特征窗长 (s)，默认 5
%     .loop       测试段用尽后是否循环，默认 true
%     .verbose    是否打印加载信息，默认 true
%     .llr_file   LLR 参数文件名，默认先找 p1_llr.json，退回 clr_params.json
%     .max_noise_s 背景池最多加载多少秒，默认 120
%     .max_targets 不指定池时最多加载几个目标，默认 6
%
%   bank 接口
%     [x, t0]            = bank.get_window(vessel, t)   场景时刻 t 对应的 5 s 音频
%     [x, t0]            = bank.get_noise(t, dur)       环境噪声段
%     [line, band, info] = bank.features(x)             LOFAR 线谱残差 + 1/3 倍频程
%     [line, band, info] = bank.get_feature(vessel, t)  取窗并提特征
%     S                  = bank.pool_spec               目标装配表（vessel/id/file/test 段）
%
%   说明
%     只加载各目标的"测试段"音频（manifest_90.csv 的 test_start..test_end，与模板段
%     时间上不相交），因此不会发生泄漏。特征定义与 Python 标定脚本的 MATLAB 参考
%     实现一致：0.2 Hz 周期图、11 bin(≈2.2 Hz) 中值、整流后插值到 1 Hz 的 491 维
%     线谱残差；10–500 Hz 的 16 个 1/3 倍频程带级。
%
%     采样率：同一批目标必须同采样率（库内 pool A–D 均为 52.734 kHz）。混入 44.1 kHz
%     的录音会直接报错并列出冲突项，而不是偷偷重采样。

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
def = struct('base', base_default, 'root', '', 'source', 'library', ...
             'pool', 'B_same', 'win', 5.0, 'loop', true, 'verbose', true, ...
             'llr_file', '', 'max_noise_s', 120, 'max_targets', 6);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), cfg.(f{i}) = def.(f{i}); end
end
if isempty(cfg.root)
    if strcmpi(cfg.source, 'origin')
        cfg.root = rp(cfg.base, 'Shipsear_origin');
    else
        cfg.root = rp(cfg.base, 'shipsEar_AUDIOS');
    end
end

%% ---------- 取数 ----------
if strcmpi(cfg.source, 'origin')
    S = load_origin(cfg);
else
    S = load_library(cfg);
end
names = S.names;  segB = S.seg;  sig = S.sig;  fs = S.fs;
tpl   = S.tpl;    noise = S.noise;  noiseDur = S.noiseDur;
llr   = S.llr;    poolSpec = S.pool_spec;

%% ---------- 组装 ----------
bank = struct();
bank.cfg       = cfg;
bank.source    = cfg.source;
bank.fs        = fs;
bank.win       = cfg.win;
bank.names     = names;
bank.tpl       = tpl;
bank.llr       = llr;
bank.llr_name  = S.llr_name;
bank.seg       = segB;
bank.noise     = noise;
bank.noiseDur  = noiseDur;
bank.pool_spec = poolSpec;
bank.idx       = containers.Map(cellstr(names), num2cell(1:numel(names)));

bank.get_window  = @get_window_impl;
bank.get_noise   = @get_noise_impl;
bank.features    = @features_impl;
bank.get_feature = @get_feature_impl;

if cfg.verbose
    fprintf('[audio_bank] source=%s  fs=%d Hz  win=%.1f s  目标 %d 个\n', ...
        cfg.source, fs, cfg.win, numel(names));
    for i = 1:numel(names)
        fprintf('             %-22s test 段 %6.1f–%6.1f s (%5.1f s)\n', ...
            names{i}, segB(i, 1), segB(i, 2), segB(i, 2) - segB(i, 1));
    end
    fprintf('             背景池 %.1f s，LLR=%s\n', noiseDur, S.llr_name);
end

%% ================= 嵌套函数 =================

    function [x, t0] = get_window_impl(vessel, t)
        key = char(string(vessel));
        if ~isKey(bank.idx, key)
            error('未知目标 "%s"；可选：%s', key, strjoin(cellstr(names), ', '));
        end
        i  = bank.idx(key);
        a0 = segB(i, 1);
        a1 = segB(i, 2);
        span = a1 - a0 - cfg.win;
        if span <= 0
            t0 = a0;
        elseif cfg.loop
            t0 = a0 + mod(t, span);
        else
            t0 = a0 + min(max(t, 0), span);
        end
        i1 = round((t0 - a0) * fs) + 1;
        i2 = i1 + round(cfg.win * fs) - 1;
        if i1 < 1 || i2 > numel(sig{i})
            error('取窗越界: %s, t=%.2f s (段 %.1f–%.1f)', key, t, a0, a1);
        end
        x = sig{i}(i1:i2);
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

    function [line, band, info] = features_impl(x)     % 与 Python 标定的参考实现一致
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
        R = max(PdB - med, 0);                    % 整流：只保留高出本底的线
        g = (10:500)';
        line = interp1(fm, R, g, 'linear', 0);
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

%% ========================================================================
function S = load_library(cfg)
%LOAD_LIBRARY 从完整库加载目标与背景池
mpath = rp(cfg.base, 'manifest_90.csv');
if ~isfile(mpath), error('找不到 %s（先运行 scripts/make_manifest90.py）', mpath); end
T = readtable(mpath, 'Delimiter', ',', 'TextType', 'string');

sel = T.role == "target" & T.usable == 1;

% ---- 目标池筛选 ----
poolName = '';
if ~isempty(cfg.pool)
    ppath = rp(cfg.base, 'pools.json');
    if ~isfile(ppath), error('找不到 %s（先运行 scripts/pools.py）', ppath); end
    P = jsondecode(fileread(ppath));
    if ~isfield(P.pools, cfg.pool)
        error('pools.json 中没有池 "%s"；可选：%s', cfg.pool, ...
            strjoin(fieldnames(P.pools)', ', '));
    end
    ids = cellfun(@(s) s.id, num2cell(P.pools.(cfg.pool).targets));
    sel = sel & ismember(T.id, ids);
    poolName = cfg.pool;
end
rows = T(sel, :);
if isempty(rows), error('没有可用的目标行（pool="%s"）', poolName); end

% ---- 无池时限制数量，避免一次读进 GB 级音频 ----
if isempty(poolName) && height(rows) > cfg.max_targets
    error(['不指定 pool 时最多加载 %d 个目标（当前可用 %d 个）。' ...
           '请用 audio_bank(struct(''pool'',''B_same'')) 选池，' ...
           '或调大 cfg.max_targets。'], cfg.max_targets, height(rows));
end

% ---- 采样率一致性 ----
ufs = unique(rows.fs);
if numel(ufs) > 1
    if isempty(poolName)
        % 无池：只保留占多数的采样率，并明确告知排除了谁
        maj = mode(rows.fs);
        bad = rows(rows.fs ~= maj, :);
        warning('audio_bank:mixedFs', ...
            '目标采样率不一致，已只保留 %d Hz，排除：%s', maj, ...
            strjoin(cellstr(bad.vessel_key + "(" + string(bad.id) + ")"), ', '));
        rows = rows(rows.fs == maj, :);
        ufs = maj;
        if isempty(rows), error('排除后没有剩余目标'); end
    else
        bad = rows(rows.fs ~= ufs(1), :);
        error(['目标池 "%s" 内采样率不一致：%s。请改池定义或排除这些录音。'], ...
            poolName, strjoin(cellstr( ...
            bad.vessel_key + "(" + string(bad.id) + ":" + string(bad.fs) + ")"), ', '));
    end
end
rows = sortrows(rows, 'id');

nV  = height(rows);
fs  = double(ufs(1));
names   = cell(nV, 1);
segB    = zeros(nV, 2);
sig     = cell(nV, 1);
poolSpec = struct('vessel', {}, 'id', {}, 'file', {}, 'tpl', {}, ...
                  'guard', {}, 'test', {}, 'test_len', {});

for i = 1:nV
    afile = fullfile(cfg.root, char(rows.file(i)));
    if ~isfile(afile), error('找不到音频 %s', afile); end
    t0 = rows.test_start(i);
    t1 = rows.test_end(i);
    s1 = round(t0 * fs) + 1;
    s2 = round(t1 * fs);
    x  = audioread(afile, [s1 s2]);
    if size(x, 2) > 1, x = mean(x, 2); end
    names{i} = char(rows.vessel_key(i));
    sig{i}   = x(:);
    segB(i, :) = [t0, t1];
    poolSpec(end+1) = struct('vessel', names{i}, 'id', double(rows.id(i)), ...
        'file', char(rows.file(i)), ...
        'tpl', [rows.tpl_start(i), rows.tpl_end(i)], ...
        'guard', [rows.guard_start(i), rows.guard_end(i)], ...
        'test', [t0, t1], 'test_len', double(rows.test_len(i)));
end

% ---- 背景池 ----
noise = zeros(round(cfg.win * fs), 1);
noiseDur = cfg.win;
bg = T(T.role == "background" & T.fs == fs, :);
if ~isempty(bg)
    nfile = fullfile(cfg.root, char(bg.file(1)));
    ni = audioinfo(nfile);
    dur = min(cfg.max_noise_s, ni.Duration);
    nz = audioread(nfile, [1 round(dur * fs)]);
    if size(nz, 2) > 1, nz = mean(nz, 2); end
    noise = nz(:);
    noiseDur = numel(noise) / fs;
else
    warning('audio_bank:noNoise', '没有与目标同采样率的背景噪声，噪声池用零填充');
end

% ---- 模板（库模式下不使用离线模板，置空） ----
tpl = struct('line', [], 'band', []);

[llr, llr_name] = load_llr(cfg);

S = struct('names', {names}, 'seg', segB, 'sig', {sig}, 'fs', fs, ...
           'tpl', tpl, 'noise', noise, 'noiseDur', noiseDur, ...
           'llr', llr, 'llr_name', llr_name, 'pool_spec', poolSpec);
end

%% ========================================================================
function S = load_origin(cfg)
%LOAD_ORIGIN 旧路径：4 段官方样例（保持逐位兼容）
mpath = rp(cfg.base, 'manifest.csv');
if ~isfile(mpath), error('找不到 %s', mpath); end
T = readtable(mpath, 'Delimiter', ',', 'TextType', 'string');

names = string(strtrim(splitlines(string(fileread(rp(cfg.base, 'tpl_names.txt'))))));
names = names(names ~= "");
nV = numel(names);
tplLine = readmatrix(rp(cfg.base, 'tpl_line.csv'));
tplBand = readmatrix(rp(cfg.base, 'tpl_band.csv'));
if size(tplLine, 1) ~= nV || size(tplBand, 1) ~= nV
    error('模板行数与 tpl_names.txt 不一致');
end

sig = cell(nV, 1);
segB = zeros(nV, 2);
fs = [];
for i = 1:nV
    rows = T(T.vessel == names(i) & T.segment == "test", :);
    if isempty(rows), error('目标 %s 没有 test 段记录', names(i)); end
    afile = fullfile(cfg.root, char(rows.file(1)));
    ai = audioinfo(afile);
    if isempty(fs), fs = ai.SampleRate; end
    if ai.SampleRate ~= fs
        error('采样率不一致：%s 为 %d Hz，前面为 %d Hz', char(names(i)), ai.SampleRate, fs);
    end
    t0 = min(rows.t_start); t1 = max(rows.t_end);
    s1 = round(t0 * fs) + 1;
    s2 = min(round(t1 * fs), ai.TotalSamples);
    x  = audioread(afile, [s1 s2]);
    if size(x, 2) > 1, x = mean(x, 2); end
    sig{i} = x(:);
    segB(i, :) = [t0, t1];
end

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

[llr, llr_name] = load_llr(cfg);

S = struct('names', {cellstr(names)}, 'seg', segB, 'sig', {sig}, 'fs', fs, ...
           'tpl', struct('line', tplLine, 'band', tplBand), ...
           'noise', noise, 'noiseDur', noiseDur, 'llr', llr, ...
           'llr_name', llr_name, 'pool_spec', []);
end

%% ========================================================================
function p = rp(base, name)
%RP 解析数据文件路径：先找 base/<sub>/name，再回退 base/name
%   新版布局把接口文件放在 config/、音频库放在 data/；
%   旧布局平铺在 base 下。两种都能用，因此迁移不必一次性完成。
subs = {'config', 'data'};
for i = 1:numel(subs)
    cand = fullfile(base, subs{i}, name);
    if isfile(cand) || isfolder(cand)
        p = cand;
        return;
    end
end
p = fullfile(base, name);
end

%% ========================================================================
function [llr, nm] = load_llr(cfg)
%LOAD_LLR 按数据源选择标定文件
%   library → p1_llr.json（同型负样本重标定），回退 clr_params.json
%   origin  → clr_params.json（旧 4 文件标定），回退 p1_llr.json
%   这样旧结果（98%/24%）仍然可逐位复现，不会被新标定悄悄改掉。
if ~isempty(cfg.llr_file)
    cand = {cfg.llr_file};
elseif strcmpi(cfg.source, 'origin')
    cand = {'clr_params.json', 'p1_llr.json'};
else
    cand = {'p1_llr.json', 'clr_params.json'};
end
for i = 1:numel(cand)
    p = rp(cfg.base, cand{i});
    if isfile(p)
        llr = jsondecode(fileread(p));
        nm = cand{i};
        return;
    end
end
error('找不到 LLR 参数文件，尝试过：%s', strjoin(cand, ', '));
end

%% ========================================================================
function base = local_data_dir()
%LOCAL_DATA_DIR 定位 05_目标噪声特征 数据目录（相对优先，绝对兜底）
base = '';
p = fileparts(mfilename('fullpath'));      % ...\matlab
for i = 1:3
    p = fileparts(p);                       % ...\AI_Agent
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
