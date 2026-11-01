function p = repo_path(varargin)
%REPO_PATH  拼接相对于仓库根目录的路径，与当前工作目录 (CWD) 无关。
%
%   p = repo_path('data', 'Cargo_40s.wav')
%
%   仓库根目录由本文件自身位置推出（本文件位于 <root>/matlab/repo_path.m），
%   因此无论从仓库根目录、matlab/ 子目录还是在 MATLAB GUI 里任意位置启动，
%   解析结果都一致。
%
%   存在的理由：generate_ship_spectrum 在 wav 文件不存在时会静默回退到仿真
%   频谱。若路径写成相对 CWD 的形式，代码搬家或换个目录启动就会在不知情的
%   情况下产出假的特征数据，故所有数据文件路径都统一由此函数解析。
%
%   参见: init_parameters.m, compare_ship_spectra.m

here = fileparts(mfilename('fullpath'));   % <root>/matlab
root = fileparts(here);                    % <root>
p = fullfile(root, varargin{:});

end
