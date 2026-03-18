function abaqus_run(path)
% 向 Abaqus 提交 inp 文件
% path : inp 计算文件的绝对路径

global GA FEA T

matlabpath = pwd();  % 记录当前 matlab 目录
fprintf('--------第%d代计算开始：--------\n', T);
time_abaqus = tic;

cd(path);

% 纯英文运行环境目录，避免 locale / unicode 问题
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

save('mat-cae.mat', 'FEA');

ABAQUS_CMD = '"C:\SIMULIA\Commands\abq2024.bat"';

% 第一步：只生成 inp，不在 CAE 里 submit
system([ABAQUS_CMD ' cae noGUI=python_CAE_Press_M.py']);

% 第二步：命令行逐个提交 inp
for n = 1:GA.N
    jobname = ['Job-', num2str(T), '-', num2str(n)];
    inpname = [jobname, '.inp'];

    if exist(inpname, 'file') == 2
        cmd = sprintf('%s job=%s input=%s cpus=%d interactive ask_delete=OFF', ...
            ABAQUS_CMD, jobname, inpname, FEA.cpus);
        fprintf('提交 %s\n', jobname);
        system(cmd);
    else
        fprintf('------%s 的 inp 未生成------\n', jobname);
    end
end

fprintf('-------第%d代计算完成-----------------\n', T);
toc(time_abaqus)
fprintf('运算情况：\n');

pause(3);

for n = 1:GA.N
    jobname = ['Job-', num2str(T), '-', num2str(n)];
    joblog  = [jobname, '.log'];

    if exist(joblog, 'file') == 2
        S = readlines(joblog);
        A = contains(S, 'COMPLETED');
        if any(A)
            fprintf('-%s 计算成功-\n', jobname);
        else
            fprintf('-----%s计算中断-----\n', jobname);
        end
    else
        fprintf('------%s创建失败------\n', jobname);
    end
end

cd(matlabpath);
end