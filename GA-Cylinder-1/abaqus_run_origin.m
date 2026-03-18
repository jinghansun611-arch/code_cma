function abaqus_run(path)
% 向 Abaqus 提交 inp 文件（代内并行版，Abaqus 2024）
% 逻辑：
%   1) 用 python_CAE_Press_M.py 仅生成本代所有 .inp
%   2) 命令行后台提交 job，最多同时运行 FEA.multitask 个
%   3) 轮询 .log/.sta/.lck 文件，直到本代全部结束
%
% 前提：python_CAE_Press_M.py 内部已经把 submit() 改成了 writeInput()

global GA FEA T

matlabpath = pwd();
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

% 建议改成你机器上的真实路径
ABAQUS_CMD = '"C:\SIMULIA\Commands\abq2024.bat"';

%% 第一步：只生成 inp（不在 CAE 里 submit）
ret = system([ABAQUS_CMD ' cae noGUI=python_CAE_Press_M.py']);
if ret ~= 0
    warning('第%d代：生成 inp 阶段返回码非 0（ret=%d）。', T, ret);
end

%% 第二步：后台并行提交 job
max_parallel = max(1, FEA.multitask);
job_queue = 1:GA.N;
running = struct('name', {}, 'bat', {}, 'idx', {}, 'launch_time', {});
next_ptr = 1;

fprintf('本代并行设置：最多同时运行 %d 个 job；每个 job 使用 %d 个 CPU。\n', max_parallel, FEA.cpus);

while (next_ptr <= numel(job_queue)) || ~isempty(running)

    % 补满并发槽位
    while (numel(running) < max_parallel) && (next_ptr <= numel(job_queue))
        n = job_queue(next_ptr);
        jobname = ['Job-', num2str(T), '-', num2str(n)];
        inpname = [jobname, '.inp'];

        if exist(inpname, 'file') ~= 2
            fprintf('------%s 的 inp 未生成，跳过------\n', jobname);
            next_ptr = next_ptr + 1;
            continue;
        end

        cleanup_old_jobfiles(jobname);
        batfile = write_job_bat(path, jobname, inpname, FEA.cpus, tmp_dir, home_dir, ABAQUS_CMD);

        launch_cmd = sprintf('start "" /B cmd /c call "%s"', batfile);
        fprintf('启动 %s\n', jobname);
        system(launch_cmd);

        s.name = jobname;
        s.bat = batfile;
        s.idx = n;
        s.launch_time = tic;
        running(end+1) = s; %#ok<AGROW>

        next_ptr = next_ptr + 1;
        pause(1.0); % 稍微错开一点启动时间，减少同时抢文件
    end

    if isempty(running)
        pause(2);
        continue;
    end

    pause(5); % 轮询间隔

    keep = true(1, numel(running));
    for i = 1:numel(running)
        jobname = running(i).name;
        state = get_job_state(jobname);
        switch state
            case 'done'
                fprintf('-%s 计算成功-\n', jobname);
                keep(i) = false;
            case 'fail'
                fprintf('-----%s 计算中断-----\n', jobname);
                keep(i) = false;
            otherwise
                % 仍在运行/排队，不输出，避免刷屏
        end
    end
    running = running(keep);
end

fprintf('-------第%d代计算完成-----------------\n', T);
toc(time_abaqus)

fprintf('运算情况：\n');
for n = 1:GA.N
    jobname = ['Job-', num2str(T), '-', num2str(n)];
    st = get_job_state(jobname);
    switch st
        case 'done'
            fprintf('-%s 计算成功-\n', jobname);
        case 'fail'
            fprintf('-----%s计算中断-----\n', jobname);
        otherwise
            fprintf('------%s状态未知------\n', jobname);
    end
end

cd(matlabpath);
end


function batfile = write_job_bat(path, jobname, inpname, cpus, tmp_dir, home_dir, ABAQUS_CMD)
% 为每个 job 生成一个独立 bat，便于后台并发启动
batfile = fullfile(path, ['run_', jobname, '.bat']);
fid = fopen(batfile, 'w');
if fid < 0
    error('无法创建 bat 文件：%s', batfile);
end
fprintf(fid, '@echo off\r\n');
fprintf(fid, 'set "TMP=%s"\r\n', tmp_dir);
fprintf(fid, 'set "TEMP=%s"\r\n', tmp_dir);
fprintf(fid, 'set "HOME=%s"\r\n', home_dir);
fprintf(fid, 'set "PYTHONIOENCODING=utf-8"\r\n');
fprintf(fid, 'cd /d "%s"\r\n', path);
fprintf(fid, 'call %s job=%s input=%s cpus=%d interactive ask_delete=OFF\r\n', ...
    ABAQUS_CMD, jobname, inpname, cpus);
fclose(fid);
end


function st = get_job_state(jobname)
% 通过 log/sta/lck 等文件判断 job 状态
% 返回：done / fail / running / pending
st = 'pending';

logf = [jobname, '.log'];
staf = [jobname, '.sta'];
lckf = [jobname, '.lck'];
msgf = [jobname, '.msg'];
datf = [jobname, '.dat'];
prtf = [jobname, '.prt'];
odbf = [jobname, '.odb'];

% 1) 先看 sta（最直接）
if exist(staf, 'file') == 2
    txt = safe_read_text(staf);
    utxt = upper(txt);
    if contains(utxt, 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY')
        st = 'done';
        return;
    end
    if contains(utxt, 'THE ANALYSIS HAS NOT BEEN COMPLETED') || ...
       contains(utxt, 'ERROR') || contains(utxt, 'A FATAL ERROR HAS OCCURRED')
        st = 'fail';
        return;
    end
end

% 2) 再看 log
if exist(logf, 'file') == 2
    txt = safe_read_text(logf);
    utxt = upper(txt);
    if contains(utxt, 'COMPLETED')
        st = 'done';
        return;
    end
    if contains(utxt, 'EXITED WITH ERRORS') || ...
       contains(utxt, 'A FATAL ERROR HAS OCCURRED') || ...
       contains(utxt, 'ABORTED')
        st = 'fail';
        return;
    end
end

% 3) dat 里也可能先出现错误
if exist(datf, 'file') == 2
    txt = safe_read_text(datf);
    utxt = upper(txt);
    if contains(utxt, 'A FATAL ERROR HAS OCCURRED') || ...
       contains(utxt, 'ERROR') || ...
       contains(utxt, 'THE ANALYSIS HAS NOT BEEN COMPLETED')
        st = 'fail';
        return;
    end
end

% 4) 还在跑的迹象
if exist(lckf, 'file') == 2 || exist(msgf, 'file') == 2 || ...
   exist(prtf, 'file') == 2 || exist(odbf, 'file') == 2
    st = 'running';
    return;
end
end


function txt = safe_read_text(fname)
try
    txt = fileread(fname);
catch
    try
        txt = strjoin(cellstr(readlines(fname)), newline);
    catch
        txt = '';
    end
end
end


function cleanup_old_jobfiles(jobname)
% 删除同名旧文件，避免读到上一轮残留状态
exts = {'.com','.dat','.ipm','.log','.msg','.odb','.prt','.sim','.sta','.lck'};
for i = 1:numel(exts)
    f = [jobname, exts{i}];
    if exist(f, 'file') == 2
        try
            delete(f);
        catch
        end
    end
end
end
