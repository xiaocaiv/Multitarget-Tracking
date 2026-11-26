function dump_feature(outfile)
%DUMP_FEATURE 用 audio_bank 导出特征，供 Python 复刻版对拍
%   M 每 491+16 列为一组，顺序与 b.names 一致
b = audio_bank(struct('pool', 'B_same', 'verbose', false));
M = [];
for i = 1:numel(b.names)
    [line, band] = b.get_feature(b.names{i}, 10);
    M = [M, [line(:); band(:)]]; %#ok<AGROW>
end
writematrix(M, outfile);
fid = fopen([outfile '.info'], 'w');
fprintf(fid, 'fs=%d\nwin=%.6f\n', b.fs, b.win);
for i = 1:numel(b.names)
    fprintf(fid, 'name=%s\nt0=%.6f\n', b.names{i}, b.seg(i, 1) + 10);
end
fclose(fid);
fprintf('saved %s  %dx%d\n', outfile, size(M, 1), size(M, 2));
end
