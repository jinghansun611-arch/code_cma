clear
clear global
clc

curdir = pwd;
[~, folder_name] = fileparts(curdir);

if ~strcmp(folder_name, 'GA-Cylinder-1')
    error('当前工作目录不是 GA-Cylinder-1，请先 cd 到 GA-Cylinder-1 后再运行主函数。\n当前目录：%s', curdir);
end

%记得如果换电脑要改abq2024.bat路径 要改三个地方 abaqus_odb_m.m abaqus_run.m abaqus_odb_open.m
% 所有单位均为 mm，与实际保持一致（PI 厚度 50 微米 = 0.05 mm）
%
% ================================
%   CMA-ES 主入口（替代 Fun_GA）
% ================================

%% ------------------------------ 参数设定区 ------------------------------
global GA FEA tar T Phe

GA.T_MAX = 40;
GA.N     = 24;
GA.TYPE  = 5;
GA.D     = [21,21,1,1,1];
GA.EXPECT = 0.005;
GA.BIT    = 10000;

GA.CMA.sigma0 = 0.28;
GA.CMA.x0     = [];
GA.CMA.elite_inject = true;
GA.CMA.elite_warmup = 0;

% if evalin('base','exist(''GA_restart_seed'',''var'')') %这是之前没有checkpoint版本 也想要继续往后算的情况下 使用的
%        GA = evalin('base','GA_restart_seed');
%    end

FEA.multitask = 4;
FEA.cpus      = 2;
FEA.load_disturb   = 5e-7;
FEA.load_pressall  = 5e-8;
FEA.load_pressbond = 5e-8;
FEA.mesh =[0.5;0.5];
FEA.t = 0.05;
FEA.step_stabilize = [1e-6;1e-6;1e-6;1e-6;1e-6];
FEA.step_time = [5;1;1;1;5];
FEA.N = GA.N;
FEA.penalty_1 = 1;

%% ------------------------------ 目标构型设置 ------------------------------
A_set_tar_equation;
Fun_cal_para_tar_eqation;

%% ------------------------------ CMA-ES 初始化 ------------------------------
map = Fun_CMAES_build_map();

if isfield(GA,'CMA') && isfield(GA.CMA,'x0') && ~isempty(GA.CMA.x0)
    x0 = GA.CMA.x0;
else
    x0 = [];
end
if isfield(GA,'CMA') && isfield(GA.CMA,'sigma0') && ~isempty(GA.CMA.sigma0)
    sigma0 = GA.CMA.sigma0;
else
    sigma0 = 0.25;
end

cma = Fun_CMAES_init(map.nvars, GA.N, x0, sigma0);
cma.elite_inject = false;   % 第1代采样不注入

%% ------------------------------ 个体库/历史量 ------------------------------
P = cell(GA.N, GA.TYPE, GA.T_MAX);
BEST_NUM = zeros(1,2);
BEST_P_G = P(1,:,1);
BEST_D_G = inf;
BEST_X_G = [];
HIS_D_P  = zeros(GA.T_MAX, GA.N);
HIS_D_G  = zeros(GA.T_MAX, 1);
HIS_D_TP = zeros(GA.T_MAX, 2);

[X, cma] = Fun_CMAES_ask(cma);
for pn = 1:GA.N
    P(pn,:,1) = Fun_CMAES_vec2PP(X(pn,:), map);
end

T = 1; Phe = 1;

%% ------------------------------ 存储设定 ------------------------------
ORDER = GA.ORDER;
Structure = GA.Structure;
name_data = ['mat-CMAES-', ORDER, Structure, '.mat'];
name_oldfile = 'A-GA-';
name_newfile = ['A-CMAES-', ORDER];

path_ga = pwd;         % 当前应为 GA-Cylinder-1
path_root = fileparts(path_ga);

% 续跑状态快照文件（保存在 GA-Cylinder-1 目录）
state_mat_file = fullfile(path_ga, ['state-CMAES-', ORDER, Structure, '.mat']);

% ===== fresh start 时，旧 checkpoint 也删掉，避免误续跑 =====
if exist(state_mat_file, 'file') == 2
    delete(state_mat_file);
    fprintf('旧 checkpoint 已删除：%s\n', state_mat_file);
end

if exist(fullfile(path_root, name_newfile), 'dir')
    fprintf('检测到旧结果目录已存在，正在删除：%s\n', fullfile(path_root, name_newfile));
    try
        rmdir(fullfile(path_root, name_newfile), 's');
        fprintf('旧结果目录删除完成。\n');
    catch ME
        error('无法删除旧结果目录 %s\n可能原因：有 Abaqus/Matlab 进程仍在占用其中的文件。\n原始错误：%s', ...
            fullfile(path_root, name_newfile), ME.message);
    end
end

copyfile(fullfile(path_root, name_oldfile), fullfile(path_root, name_newfile));
path_abaqus = fullfile(path_root, name_newfile);

fprintf('\n==================== CMA-ES 优化启动 ====================\n');
fprintf('启动时间      : %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf('ORDER / 结构  : %s  |  %s\n', GA.ORDER, GA.Structure);
fprintf('控制点 D      : %s\n', mat2str(GA.D));
fprintf('自由变量数    : %d\n', map.nvars);
fprintf('目标离散      : tar.n = %d, 条带数 = %d\n', tar.n, tar.line);
fprintf('CMA-ES 参数   : 种群 = %d, T_MAX = %d, sigma0 = %.4g\n', GA.N, GA.T_MAX, GA.CMA.sigma0);
fprintf('精英注入      : %d（前 %d 代关闭）\n', GA.CMA.elite_inject, GA.CMA.elite_warmup);
fprintf('Abaqus 并行   : multitask = %d, cpus/job = %d\n', FEA.multitask, FEA.cpus);
fprintf('工作目录      : %s\n', path_abaqus);
fprintf('=========================================================\n');

time_all = tic;
fprintf('初始化完成：');
toc(time_all)

progress_log = fullfile(path_abaqus, 'progress_log.txt');
if exist(progress_log,'file') == 2
    delete(progress_log);
end

%% =============================== 主循环 ===============================
for T = 1:GA.T_MAX
    fprintf('\n==================== 第 %d / %d 代开始 ====================\n', T, GA.T_MAX);
    fprintf('时间            : %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    if isfinite(BEST_D_G)
        fprintf('当前全局最优    : %.7g  (T=%d, Phe=%d)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));
    else
        fprintf('当前全局最优    : 暂无\n');
    end
    fprintf('当前 sigma      : %.6g\n', cma.sigma);
    fprintf('当前精英注入状态: %d | elite_count = %d | fbest = %.6g\n', cma.elite_inject, cma.elite_count, cma.fbest);
    fprintf('本代个体数      : %d\n', GA.N);
    fprintf('----------------------------------------------------------\n');

    time_T = tic;
    FEA.NUM_T = T;

    P(:,:,T) = Fun_set_special(P(:,:,T));
    X_feasible = zeros(GA.N, map.nvars);
    for pn = 1:GA.N
        X_feasible(pn,:) = Fun_CMAES_PP2vec(P(pn,:,T), map)';
    end
    X_T = X_feasible;

    %% 阶段1：生成本代 mat 文件
    fprintf('[阶段 1/3] 生成本代草图 mat 文件...\n');
    t_mat = tic;
    for Phe = 1:GA.N
        name_cae_mat = ['mat-', num2str(T), '-', num2str(Phe), '.mat'];
        sketch = Fun_cal_para_2D(P(Phe,:,T));
    
        % Python 脚本会读取 Phe / FEA / tar / sketch
        save(fullfile(path_abaqus, name_cae_mat), 'Phe', 'FEA', 'tar', 'sketch');
    
        fprintf('  [MAT %2d/%d] %s 已生成\n', Phe, GA.N, name_cae_mat);
    end
    time_mat = toc(t_mat);
    fprintf('[阶段 1/3] 完成：共生成 %d 个 mat，耗时 %.3f s\n', GA.N, time_mat);

    %% 阶段2：Abaqus 求解
    fprintf('[阶段 2/3] Abaqus 并行求解开始...\n');
    t_run = tic;
    run_info = abaqus_run(path_abaqus);
    time_run = toc(t_run);
    fprintf('[阶段 2/3] Abaqus 求解结束，耗时 %.3f s\n', time_run);

    %% 阶段3：严格按 DONE 的 job 读取 ODB
    fprintf('[阶段 3/3] 读取 ODB 并计算 DEL...\n');
    t_odb = tic;

    solve_done_cnt   = nnz(run_info.status_map == "DONE");
    solve_failed_cnt = nnz(run_info.status_map == "FAILED") + nnz(run_info.status_map == "NOINP");
    solve_other_cnt  = GA.N - solve_done_cnt - solve_failed_cnt;

    odb_ok_cnt = 0;
    odb_fail_cnt = 0;
    odb_skip_cnt = 0;

    for Phe = 1:GA.N
        jobname = sprintf('Job-%d-%d', T, Phe);
        job_state = string(run_info.status_map(Phe));

        if job_state ~= "DONE"
            del = 1;
            odb_skip_cnt = odb_skip_cnt + 1;
            fprintf('  [ODB %2d/%d] Phe=%d | job=%s | 求解状态=%s | 跳过 ODB | DEL=1\n', ...
                Phe, GA.N, Phe, jobname, job_state);
        else
            [result, odb_ok, odb_msg] = abaqus_odb_m(path_abaqus);
            if odb_ok && ~isempty(result.rx{1})
                result = Fun_cal_para_odb_m(result);
                %del = cal_delta_C_m(result);
                del = cal_delta_C_m_new(result);
                odb_ok_cnt = odb_ok_cnt + 1;
                fprintf('  [ODB %2d/%d] Phe=%d | job=%s | ODB读取成功 | DEL=%.6g\n', ...
                    Phe, GA.N, Phe, jobname, del);
            else
                del = 1;
                odb_fail_cnt = odb_fail_cnt + 1;
                fprintf('  [ODB %2d/%d] Phe=%d | job=%s | ODB读取失败：%s | DEL=1\n', ...
                    Phe, GA.N, Phe, jobname, odb_msg);
            end
        end

        HIS_D_P(T, Phe) = del;
        if del < BEST_D_G
            BEST_D_G = del;
            BEST_P_G = P(Phe,:,T);
            BEST_X_G = X_T(Phe,:).';
            BEST_NUM = [T, Phe];
        end
    end

    HIS_D_G(T) = BEST_D_G;
    HIS_D_TP(T,:) = BEST_NUM;
    time_odb = toc(t_odb);
    fprintf('[阶段 3/3] ODB 处理结束，耗时 %.3f s\n', time_odb);

    %% 本代统计
    del_vec = HIS_D_P(T,:);
    [best_gen, idx_best_gen] = min(del_vec);
    mean_gen = mean(del_vec);
    med_gen  = median(del_vec);

    fprintf('\n当前最佳误差 ：%g    最佳个体信息 ： JOB-%d-%d \n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));
    fprintf('本代(T=%d) 各个体 DEL：\n', T);
    for p = 1:GA.N
        fprintf('  %2d:% .6g', p, del_vec(p));
        if mod(p,6) == 0 || p == GA.N
            fprintf('\n');
        end
    end
    fprintf('本代最优 DEL = %.6g (Phe=%d)\n', best_gen, idx_best_gen);
    if T > 1
        prev_best = min(HIS_D_P(T-1,:));
        improve = prev_best - best_gen;
        if isfinite(prev_best) && prev_best > 0
            improve_pct = 100 * improve / prev_best;
            fprintf('较上一代最优(%.6g) 变化：%+.6g  (%+.2f%%)\n', prev_best, improve, improve_pct);
        else
            fprintf('较上一代最优(%.6g) 变化：%+.6g\n', prev_best, improve);
        end
    else
        fprintf('第一代：无上一代可对比。\n');
    end

    fprintf('-------------------- 本代汇总 --------------------\n');
    fprintf('求解成功 / 失败 / 其它 : %d / %d / %d\n', solve_done_cnt, solve_failed_cnt, solve_other_cnt);
    fprintf('ODB成功 / 失败 / 跳过 : %d / %d / %d\n', odb_ok_cnt, odb_fail_cnt, odb_skip_cnt);
    fprintf('本代均值 / 中位数     : %.6g / %.6g\n', mean_gen, med_gen);
    fprintf('当前全局最优          : %.6g  (T=%d, Phe=%d)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));
    fprintf('阶段耗时(s)           : mat=%.3f | run=%.3f | odb=%.3f | total=%.3f\n', ...
        time_mat, time_run, time_odb, toc(time_T));
    fprintf('--------------------------------------------------\n');

    fid = fopen(progress_log, 'a');
    if fid > 0
        fprintf(fid, 'T=%d best=%.8g mean=%.8g median=%.8g solve_done=%d solve_fail=%d solve_other=%d odb_ok=%d odb_fail=%d odb_skip=%d sigma=%.8g gen_time=%.3f\n', ...
            T, best_gen, mean_gen, med_gen, solve_done_cnt, solve_failed_cnt, solve_other_cnt, odb_ok_cnt, odb_fail_cnt, odb_skip_cnt, cma.sigma, toc(time_T));
        fclose(fid);
    end

    if BEST_D_G < GA.EXPECT
        fprintf('CMA-ES complete! T=%d. Best DEL=%g\n', T, BEST_D_G);
        break;
    end

    if T < GA.T_MAX
        cma = Fun_CMAES_tell(cma, X_feasible, HIS_D_P(T,:)');
        if (T+1) <= GA.CMA.elite_warmup
            cma.elite_inject = false;
        else
            cma.elite_inject = true;
        end

        [X, cma] = Fun_CMAES_ask(cma);
        if cma.elite_inject && isfinite(cma.fbest) && all(isfinite(cma.xbest))
            elite_diff = norm(X(1,:) - cma.xbest(:)');
            fprintf('精英注入检查（T=%d -> T=%d）: ||X(1,:)-xbest|| = %.3e\n', T, T+1, elite_diff);
        else
            fprintf('精英注入检查（T=%d -> T=%d）: 当前未启用或 xbest/fbest 无效\n', T, T+1);
        end
        for pn = 1:GA.N
            P(pn,:,T+1) = Fun_CMAES_vec2PP(X(pn,:), map);
        end
    end

    fprintf('本代历时 '); toc(time_T)
    fprintf('总运行时间：'); toc(time_all)
    fprintf('-----------------------------第%d代结束-----------------------------\n', T);
        % ===== 每代结束后保存一次完整 checkpoint，便于续跑 =====
    T_NOW = T;
    save(state_mat_file, ...
        'GA','FEA','tar','map','cma','P', ...
        'BEST_NUM','BEST_P_G','BEST_X_G','BEST_D_G', ...
        'HIS_D_P','HIS_D_G','HIS_D_TP', ...
        'T_NOW','path_abaqus','progress_log', ...
        '-v7.3');
    
    fprintf('Checkpoint 已保存：%s\n', state_mat_file);
end



%% 保存/绘图（仅对 best job 真正成功时读取 ODB）
save_sketch = Fun_cal_para_2D(P(BEST_NUM(2),:,BEST_NUM(1)));
T = BEST_NUM(1);
Phe = BEST_NUM(2);

best_jobname = sprintf('Job-%d-%d', BEST_NUM(1), BEST_NUM(2));
[best_state, ~] = check_abaqus_job_status(path_abaqus, best_jobname);

if best_state ~= "DONE"
    warning('最佳个体 %s 并未真正成功完成，跳过最终 ODB 读取与绘图。', best_jobname);
    best_odb_ok = false;
    save_result = [];
else
    [save_result, best_odb_ok, best_odb_msg] = abaqus_odb_m(path_abaqus);
    if ~best_odb_ok
        warning('最终最佳个体 ODB 读取失败：%s', best_odb_msg);
    else
        save_result = Fun_cal_para_odb_m(save_result);
    end
end

save_his_del = HIS_D_G(1:T_NOW);
save(fullfile(path_abaqus, name_data));

bestfile = fullfile(path_abaqus, 'best_stage1.mat');
save(bestfile, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','-v7.3');
fprintf('Saved best to: %s (DEL=%g at T=%d, Phe=%d)\n', bestfile, BEST_D_G, BEST_NUM(1), BEST_NUM(2));

fprintf('-----------------------------CMA-ES 结束！共计完成%d代！--------------------------\n', T_NOW);
toc(time_all)



if best_odb_ok && ~isempty(save_result.rx{1})
    Result_plot(save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM);
else
    fprintf('跳过 Result_plot：最佳个体 ODB 未成功读取。\n');
end

if best_odb_ok
    abaqus_odb_open(path_abaqus, BEST_NUM(1), BEST_NUM(2));
end