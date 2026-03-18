function abaqus_odb_open(path, T, Phe)
% 打开指定 job 的 ODB（GUI）
% path : Abaqus 工作目录
% T    : 代数
% Phe  : 个体编号

matlabpath = pwd();
cd(path);

% 纯英文运行环境，避免 locale / unicode 问题
runtime_dir = fullfile(path, '_abq_runtime');
tmp_dir     = fullfile(runtime_dir, 'tmp');
home_dir    = fullfile(runtime_dir, 'home');

if ~exist(runtime_dir, 'dir'); mkdir(runtime_dir); end
if ~exist(tmp_dir, 'dir'); mkdir(tmp_dir); end
if ~exist(home_dir, 'dir'); mkdir(home_dir); end

setenv('TMP',  tmp_dir);
setenv('TEMP', tmp_dir);
setenv('HOME', home_dir);
setenv('PYTHONIOENCODING', 'utf-8');

save('mat-odb.mat','Phe','T');

ABAQUS_CMD = '"C:\SIMULIA\Commands\abq2024.bat"';

% 用 start 打开 GUI，不阻塞 MATLAB
cmd = sprintf('start "" cmd /c ""%s" cae script=python_ODB_open.py"', ABAQUS_CMD);
ret = system(cmd);

if ret ~= 0
    warning('abaqus_odb_open:LaunchFailed', ...
        '启动 Abaqus/CAE 打开 ODB 失败，返回码 = %d', ret);
end

cd(matlabpath);
end