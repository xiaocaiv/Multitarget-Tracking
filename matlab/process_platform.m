function [tracks, unassociated, info] = process_platform(tracks, measurements, plat_pos, params)
%PROCESS_PLATFORM 单平台处理总入口 (按 assoc_method 选择关联算法)
%   [tracks, unassociated, info] = PROCESS_PLATFORM(tracks, measurements, plat_pos, params)
%
%   jpda  -> 调用现有 jpda_process_platform (默认路径, 数值与改动前一致)
%   其余  -> 调用 数据关联算法对比/assoc_process_platform
%
%   本文件位于 matlab/ 根, 负责把 "数据关联算法对比/" 子文件夹加入
%   MATLAB 搜索路径; 这样 run_simulation / main 只需 addpath('matlab').
%
%   参见: init_parameters.m, run_simulation.m, main.m,
%         数据关联算法对比/assoc_process_platform.m

assoc_dir = fullfile(fileparts(mfilename('fullpath')), '数据关联算法对比');
if exist(assoc_dir, 'dir') && ~contains(path, assoc_dir)
    addpath(assoc_dir);
end

method = 'jpda';
if isfield(params, 'assoc_method') && ~isempty(params.assoc_method)
    method = lower(char(params.assoc_method));
end

if strcmp(method, 'jpda')
    [tracks, unassociated, info] = jpda_process_platform( ...
        tracks, measurements, plat_pos, params);
else
    [tracks, unassociated, info] = assoc_process_platform( ...
        tracks, measurements, plat_pos, params);
end
end
