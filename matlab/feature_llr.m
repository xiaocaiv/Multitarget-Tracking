function [ell, s] = feature_llr(f1, f2, llr, mode)
%FEATURE_LLR  由两个特征算余弦相似度与对数似然比（calibrated LLR）
%
%   [ell, s] = feature_llr(f1, f2, llr, mode)
%
%   输入
%     f1, f2  特征结构体，字段 .line (491x1)、.band (16x1)；
%             也可以直接传向量：feature_llr(line1, line2, llr, mode) 仅用线谱
%     llr     audio_bank 返回的 bank.llr（含 methods.track / methods.init）
%     mode    'track'（量测 vs 航迹模板，默认）或 'init'（量测 vs 量测）
%
%   输出
%     ell     对数似然比 (nat)，已按参数截断
%     s       结构体，.line / .band 两路余弦相似度
%
%   说明
%     两套参数必须分开用：起批阶段两边都是单次观测，相似度分布与
%     "观测 vs 已平均的模板"不同，混用会让 LLR 系统性偏高。
%     融合形式：ell = w' * ((s - mu) ./ sd) + b，其中 s = [s_line, s_band]。

if nargin < 4, mode = 'track'; end
if ~isfield(llr.methods, mode)
    error('llr 中没有 mode = "%s" 的参数（可选：%s）', mode, ...
        strjoin(fieldnames(llr.methods), ', '));
end
p = llr.methods.(mode);

% ---- 取相似度 ----
if isstruct(f1) && isstruct(f2)
    s.line = cos_sim(f1.line, f2.line);
    s.band = cos_sim(f1.band, f2.band);
else
    s.line = cos_sim(f1, f2);
    s.band = 0;                      % 仅线谱时，带级项按均值处理（贡献为 0）
end

% ---- 标准化 + 融合 ----
w  = p.w(:);
mu = p.mu(:);
sd = p.sd(:);
z  = ([s.line; s.band] - mu) ./ sd;
ell = z' * w + p.b;
ell = min(max(ell, -p.clip), p.clip);
end

% ------------------------------------------------------------------------
function c = cos_sim(a, b)
a = a(:); b = b(:);
na = norm(a); nb = norm(b);
if na < 1e-12 || nb < 1e-12
    c = 0;
else
    c = dot(a, b) / (na * nb);
end
end
