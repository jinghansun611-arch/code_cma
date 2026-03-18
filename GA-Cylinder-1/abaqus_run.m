function run_info = abaqus_run(path)
% 向 Abaqus 提交 inp 文件（代内并行版，严格成功判据）
% 成功判据：仅当 .sta 显示分析已成功完成时，才认为 DONE。
% 失败判据：.sta/.dat/launcher 明确显示失败，或 inp 未生成。
% 额外增加：若 .sta/.msg/.dat 长时间完全无更新，则判定为“卡死失败”。

global GA FEA T

matlabpath = pwd();
fprintf('--------第%d代计算开始：--------\n', T);
time_abaqus = tic;
cd(path);

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

abaqus_bat = 'D:\SIMULIA\Commands\abq2024.bat';
if exist(abaqus_bat, 'file') ~= 2
    error('找不到 Abaqus 启动文件：%s', abaqus_bat);
end
ABAQUS_CMD = ['"', abaqus_bat, '"'];

%% 第一步：只生成 inp
fprintf('[Abaqus-1/2] 生成本代 inp 文件...\n');
t_inp = tic;
ret = system([ABAQUS_CMD ' cae noGUI=python_CAE_Press_M.py']);
fprintf('Exit code: %d\n', ret);

inp_exists = false(1, GA.N);
missing_list = [];
for n = 1:GA.N
    inpname = sprintf('Job-%d-%d.inp', T, n);
    inpfull = fullfile(path, inpname);
    inp_exists(n) = isfile(inpfull);
    if ~inp_exists(n)
        missing_list(end+1) = n; %#ok<AGROW>
    end
end

fprintf('[Abaqus-1/2] inp 统计：已生成 %d / %d，耗时 %.2f s\n', nnz(inp_exists), GA.N, toc(t_inp));
if isempty(missing_list)
    fprintf('[Abaqus-1/2] 所有 inp 已生成。\n');
else
    fprintf('[Abaqus-1/2] 未生成 inp 的个体：%s \n', num2str(missing_list));
end

run_info = struct();
run_info.inp_generated = nnz(inp_exists);
run_info.total_jobs = GA.N;
run_info.status_map = strings(1, GA.N);
run_info.status_map(:) = "PENDING";
run_info.status_map(~inp_exists) = "NOINP";
run_info.solve_done = 0;
run_info.solve_fail = 0;
run_info.solve_other = 0;

if run_info.inp_generated == 0
    fprintf('[Abaqus-2/2] 跳过求解：本代未生成任何 inp。\n');
    run_info.solve_fail = nnz(run_info.status_map == "NOINP");
    cd(matlabpath);
    return;
end

running_file = fullfile(path, sprintf('running_jobs_T%d.txt', T));
status_file  = fullfile(path, sprintf('running_jobs_status_T%d.txt', T));
if exist(running_file, 'file') == 2, delete(running_file); end
if exist(status_file,  'file') == 2, delete(status_file);  end

max_parallel = max(1, FEA.multitask);
job_queue = 1:GA.N;

% ===== 轮询与卡死判据 =====
poll_pause_sec = 20;                 % 每轮轮询间隔（秒）
stall_minutes  = 10;                 % 连续多少分钟无文件进展，判为卡死
max_stall_polls = ceil(stall_minutes * 60 / poll_pause_sec);
min_runtime_before_stall_min = 10;   % 至少运行多久后才允许判卡死

running = struct('name', {}, 'bat', {}, 'idx', {}, ...
                 'launch_time', {}, 'last_sig', {}, 'stall_polls', {});
next_ptr = 1;
poll_id = 0;

fprintf('[Abaqus-2/2] 并行求解开始：最多同时运行 %d 个 job；每个 job 使用 %d 个 CPU。\n', max_parallel, FEA.cpus);
fprintf('[Abaqus-2/2] 卡死判据：连续 %.1f 分钟 .sta/.msg/.dat 无变化，则按失败处理。\n', stall_minutes);

while (next_ptr <= numel(job_queue)) || ~isempty(running)

    % ===== 往队列里补充新 job =====
    while (numel(running) < max_parallel) && (next_ptr <= numel(job_queue))
        n = job_queue(next_ptr);
        jobname = sprintf('Job-%d-%d', T, n);
        inpname = sprintf('%s.inp', jobname);
        inpfull = fullfile(path, inpname);

        if ~isfile(inpfull)
            fprintf('------%s 的 inp 未生成，跳过------\n', jobname);
            run_info.status_map(n) = "NOINP";
            next_ptr = next_ptr + 1;
            continue;
        end

        cleanup_old_jobfiles(path, jobname);
        batfile = write_job_bat(path, jobname, inpname, FEA.cpus, tmp_dir, home_dir, ABAQUS_CMD);
        launch_cmd = sprintf('start "" /B cmd /c call "%s"', batfile);

        fprintf('启动 %-12s (队列推进 %d/%d)\n', jobname, n, GA.N);
        system(launch_cmd);

        s.name = jobname;
        s.bat = batfile;
        s.idx = n;
        s.launch_time = tic;
        s.last_sig = get_job_progress_signature(path, jobname);
        s.stall_polls = 0;

        running(end+1) = s; %#ok<AGROW>
        run_info.status_map(n) = "RUNNING";

        fid = fopen(running_file, 'a');
        if fid > 0
            fprintf(fid, '%s\n', jobname);
            fclose(fid);
        end

        fid = fopen(status_file, 'a');
        if fid > 0
            fprintf(fid, '%s\tRUNNING\n', jobname);
            fclose(fid);
        end

        next_ptr = next_ptr + 1;
        pause(1.0);
    end

    if isempty(running)
        pause(2);
        continue;
    end

    pause(poll_pause_sec);
    poll_id = poll_id + 1;
    keep = true(1, numel(running));

    % ===== 检查当前运行中的 job =====
    for i = 1:numel(running)
        jobname = running(i).name;
        n = running(i).idx;

        [st, reason] = check_abaqus_job_status(path, jobname);

        switch st
            case "DONE"
                fprintf('-%s 真正完成（状态文件确认成功）-\n', jobname);
                run_info.status_map(n) = "DONE";
                keep(i) = false;
                append_status(status_file, jobname, "DONE", reason);

            case "FAILED"
                fprintf('-----%s 计算失败/中断（%s）-----\n', jobname, reason);
                run_info.status_map(n) = "FAILED";
                keep(i) = false;
                append_status(status_file, jobname, "FAILED", reason);

            otherwise
                % RUNNING / PENDING：检查是否“长时间无进展”
                new_sig = get_job_progress_signature(path, jobname);

                if strcmp(new_sig, running(i).last_sig)
                    running(i).stall_polls = running(i).stall_polls + 1;
                else
                    running(i).stall_polls = 0;
                    running(i).last_sig = new_sig;
                end

                wall_minutes = toc(running(i).launch_time) / 60;

                if running(i).stall_polls >= max_stall_polls && ...
                   wall_minutes >= min_runtime_before_stall_min

                    fprintf('-----%s 疑似卡死：%.1f 分钟无任何 .sta/.msg/.dat 进展，按失败处理-----\n', ...
                        jobname, wall_minutes);

                    run_info.status_map(n) = "FAILED";
                    keep(i) = false;
                    append_status(status_file, jobname, "FAILED", "STALLED: no sta/msg/dat progress");
                end
        end
    end

    running = running(keep);

    n_done   = nnz(run_info.status_map == "DONE");
    n_failed = nnz(run_info.status_map == "FAILED") + nnz(run_info.status_map == "NOINP");
    n_run    = nnz(run_info.status_map == "RUNNING") + nnz(run_info.status_map == "PENDING");
    n_left   = max(0, GA.N - next_ptr + 1);

    fprintf('[轮询 #%d] 队列剩余=%d | 运行中=%d | 成功=%d | 失败=%d\n', ...
        poll_id, n_left, n_run, n_done, n_failed);
end

% ===== 汇总 =====
run_info.solve_done = nnz(run_info.status_map == "DONE");
run_info.solve_fail = nnz(run_info.status_map == "FAILED") + nnz(run_info.status_map == "NOINP");
run_info.solve_other = GA.N - run_info.solve_done - run_info.solve_fail;

fprintf('-------第%d代计算完成-----------------\n', T);
toc(time_abaqus)
fprintf('运算情况（严格成功判据）：\n');

for n = 1:GA.N
    jobname = sprintf('Job-%d-%d', T, n);
    switch string(run_info.status_map(n))
        case "DONE"
            fprintf('-%s 真正完成-\n', jobname);
        case "FAILED"
            fprintf('-----%s 计算失败/中断-----\n', jobname);
        case "NOINP"
            fprintf('------%s 未生成 inp------\n', jobname);
        otherwise
            fprintf('------%s 状态未确认：%s------\n', jobname, string(run_info.status_map(n)));
    end
end

cd(matlabpath);
end


function batfile = write_job_bat(path, jobname, inpname, cpus, tmp_dir, home_dir, ABAQUS_CMD)
% 生成后台启动 bat，并把命令行输出写入 launcher 日志
batfile = fullfile(path, ['run_', jobname, '.bat']);
launcher_log = fullfile(path, ['launcher_', jobname, '.out']);

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

fprintf(fid, 'echo [%%date%% %%time%%] START %s > "%s"\r\n', jobname, launcher_log);
fprintf(fid, 'echo CMD: call %s job=%s input="%s" cpus=%d interactive ask_delete=OFF >> "%s"\r\n', ...
    ABAQUS_CMD, jobname, inpname, cpus, launcher_log);
fprintf(fid, 'call %s job=%s input="%s" cpus=%d interactive ask_delete=OFF >> "%s" 2>&1\r\n', ...
    ABAQUS_CMD, jobname, inpname, cpus, launcher_log);
fprintf(fid, 'echo EXITCODE=%%ERRORLEVEL%% >> "%s"\r\n', launcher_log);
fprintf(fid, 'echo [%%date%% %%time%%] END %s >> "%s"\r\n', jobname, launcher_log);

fclose(fid);
end


function append_status(status_file, jobname, st, reason)
fid = fopen(status_file, 'a');
if fid > 0
    fprintf(fid, '%s\t%s\t%s\n', jobname, st, reason);
    fclose(fid);
end
end


function cleanup_old_jobfiles(path, jobname)
exts = {'.com','.dat','.ipm','.log','.msg','.odb','.prt','.sim','.sta','.lck','.023','.mdl','.res'};
for i = 1:numel(exts)
    f = fullfile(path, [jobname, exts{i}]);
    if exist(f, 'file') == 2
        try, delete(f); catch, end
    end
end

launcher = fullfile(path, ['launcher_', jobname, '.out']);
if exist(launcher,'file') == 2
    try, delete(launcher); catch, end
end

bat = fullfile(path, ['run_', jobname, '.bat']);
if exist(bat,'file') == 2
    try, delete(bat); catch, end
end
end


function sig = get_job_progress_signature(path, jobname)
% 用 .sta/.msg/.dat 的“文件大小 + 修改时间”作为进展签名
% 只要这些文件持续变化，就说明 job 还在往前走

files = { ...
    fullfile(path, [jobname, '.sta']), ...
    fullfile(path, [jobname, '.msg']), ...
    fullfile(path, [jobname, '.dat']) ...
    };

sig = '';

for i = 1:numel(files)
    f = files{i};
    if exist(f, 'file') == 2
        d = dir(f);
        sig = [sig, sprintf('|%s|%d|%.12f', d.name, d.bytes, d.datenum)]; %#ok<AGROW>
    else
        sig = [sig, '|MISSING']; %#ok<AGROW>
    end
end
end