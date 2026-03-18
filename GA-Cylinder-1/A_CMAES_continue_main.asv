function A_CMAES_continue_main(state_mat_file, add_generations)
% 从已保存的 checkpoint 继续往下算
%
% 用法：
%   A_CMAES_continue_main( ...
%       'E:\code_cma_newnew\GA-Cylinder-1\state-CMAES-N3R-chiralDoublePeakSoft.mat', 20)
%
% 含义：
%   - 从 checkpoint 记录的最后一代继续
%   - 再额外跑 add_generations 代
%
% 前提：
%   1) 这个 checkpoint 是由 A_CMAES_muity_main.m 每代自动保存出来的
%   2) 结果目录 path_abaqus 还存在
%   3) 不会删除旧结果目录，而是直接在原目录续算

clear global
clc

global GA FEA tar T Phe

if nargin < 1 || isempty(state_mat_file)
    error('请提供 checkpoint 文件路径。');
end
if nargin < 2 || isempty(add_generations)
    add_generations = 10;
end

if exist(state_mat_file, 'file') ~= 2
    error('找不到 checkpoint 文件：%s', state_mat_file);
end

S = load(state_mat_file);

% ===== 检查关键变量 =====
need_vars = {'GA','FEA','tar','map','cma','P','BEST_NUM','BEST_P_G','BEST_X_G', ...
             'BEST_D_G','HIS_D_P','HIS_D_G','HIS_D_TP','T_NOW','path_abaqus'};
for i = 1:numel(need_vars)
    if ~isfield(S, need_vars{i})
        error('checkpoint 缺少关键变量：%s', need_vars{i});
    end
end

GA       = S.GA;
FEA      = S.FEA;
tar      = S.tar;
map      = S.map;
cma      = S.cma;
P        = S.P;
BEST_NUM = S.BEST_NUM;
BEST_P_G = S.BEST_P_G;
BEST_X_G = S.BEST_X_G;
BEST_D_G = S.BEST_D_G;
HIS_D_P  = S.HIS_D_P;
HIS_D_G  = S.HIS_D_G;
HIS_D_TP = S.HIS_D_TP;
T_NOW    = S.T_NOW;
path_abaqus = S.path_abaqus;

if exist(path_abaqus, 'dir') ~= 7
    error('结果目录不存在：%s', path_abaqus);
end

T_start = T_NOW + 1;
T_end   = T_NOW + add_generations;
old_TMAX = GA.T_MAX;
GA.T_MAX = T_end;

% ===== 如果 checkpoint 里还没有下一代种群，则用最后一代结果补做 tell + ask =====
if T_start <= size(P,3)
    need_rebuild_next_gen = isempty(P{1,1,T_start});
else
    need_rebuild_next_gen = true;
end

if need_rebuild_next_gen
    fprintf('\n检测到 checkpoint 中第 %d 代种群为空，开始根据第 %d 代结果补建下一代...\n', T_start, T_NOW);

    % 1) 把第 T_NOW 代的 PP 映射回 CMA-ES 向量
    X_last = zeros(GA.N, map.nvars);
    for pn = 1:GA.N
        X_last(pn,:) = Fun_CMAES_PP2vec(P(pn,:,T_NOW), map)';
    end

    % 2) 用第 T_NOW 代的 DEL 做 tell
    cma = Fun_CMAES_tell(cma, X_last, HIS_D_P(T_NOW,:)');

    % 3) 设置下一代的精英注入状态（和主函数逻辑保持一致）
    if (T_NOW + 1) <= GA.CMA.elite_warmup
        cma.elite_inject = false;
    else
        cma.elite_inject = true;
    end

    % 4) ask 出下一代
    [X_new, cma] = Fun_CMAES_ask(cma);

    % 5) 写入 P(:,:,T_NOW+1)
    for pn = 1:GA.N
        P(pn,:,T_NOW+1) = Fun_CMAES_vec2PP(X_new(pn,:), map);
    end

    fprintf('第 %d 代种群已补建完成。\n\n', T_start);
end

ORDER = GA.ORDER;
Structure = GA.Structure;
name_data = ['mat-CMAES-', ORDER, Structure, '.mat'];
progress_log = fullfile(path_abaqus, 'progress_log.txt');

fprintf('\n==================== CMA-ES 续跑启动 ====================\n');
fprintf('checkpoint      : %s\n', state_mat_file);
fprintf('结果目录        : %s\n', path_abaqus);
fprintf('已完成代数      : %d\n', T_NOW);
fprintf('续跑范围        : 第 %d 代 -> 第 %d 代\n', T_start, T_end);
fprintf('旧 T_MAX / 新 T_MAX : %d / %d\n', old_TMAX, GA.T_MAX);
fprintf('当前全局最优    : %.7g  (T=%d, Phe=%d)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));
fprintf('当前 sigma      : %.6g\n', cma.sigma);
fprintf('=========================================================\n');

% ===== 扩容 =====
if size(P,3) < GA.T_MAX
    P_new = cell(GA.N, GA.TYPE, GA.T_MAX);
    P_new(:,:,1:size(P,3)) = P;
    P = P_new;
end

if size(HIS_D_P,1) < GA.T_MAX
    tmp = zeros(GA.T_MAX, GA.N);
    tmp(1:size(HIS_D_P,1),:) = HIS_D_P;
    HIS_D_P = tmp;
end

if size(HIS_D_G,1) < GA.T_MAX
    tmp = zeros(GA.T_MAX, 1);
    tmp(1:size(HIS_D_G,1),:) = HIS_D_G;
    HIS_D_G = tmp;
end

if size(HIS_D_TP,1) < GA.T_MAX
    tmp = zeros(GA.T_MAX, 2);
    tmp(1:size(HIS_D_TP,1),:) = HIS_D_TP;
    HIS_D_TP = tmp;
end

time_all = tic;

%% =============================== 主循环（续跑） ===============================
for T = T_start:T_end
    fprintf('\n==================== 第 %d / %d 代开始（续跑） ====================\n', T, GA.T_MAX);
    fprintf('时间            : %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    fprintf('当前全局最优    : %.7g  (T=%d, Phe=%d)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));
    fprintf('当前 sigma      : %.6g\n', cma.sigma);
    fprintf('当前精英注入状态: %d | elite_count = %d | fbest = %.6g\n', ...
        cma.elite_inject, cma.elite_count, cma.fbest);
    fprintf('本代个体数      : %d\n', GA.N);
    fprintf('----------------------------------------------------------\n');

    time_T = tic;
    FEA.NUM_T = T;

    if isempty(P{1,1,T})
        error('P(:,:,T) 在第 %d 代为空，且补建下一代后仍为空，请检查 checkpoint 是否损坏。', T);
    end

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

    %% 阶段3：只对 DONE 的 job 读取 ODB
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
                del = cal_delta_C_m(result);
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

    prev_best = min(HIS_D_P(T-1,:));
    improve = prev_best - best_gen;
    if isfinite(prev_best) && prev_best > 0
        improve_pct = 100 * improve / prev_best;
        fprintf('较上一代最优(%.6g) 变化：%+.6g  (%+.2f%%)\n', prev_best, improve, improve_pct);
    else
        fprintf('较上一代最优(%.6g) 变化：%+.6g\n', prev_best, improve);
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

    % 每代结束都保存 checkpoint
    T_NOW = T;
    save(state_mat_file, ...
        'GA','FEA','tar','map','cma','P', ...
        'BEST_NUM','BEST_P_G','BEST_X_G','BEST_D_G', ...
        'HIS_D_P','HIS_D_G','HIS_D_TP', ...
        'T_NOW','path_abaqus','progress_log', ...
        '-v7.3');
    fprintf('Checkpoint 已更新：%s\n', state_mat_file);
end

T_NOW = T;
bestfile = fullfile(path_abaqus, 'best_stage1.mat');
save(bestfile, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','-v7.3');

save(state_mat_file, ...
    'GA','FEA','tar','map','cma','P', ...
    'BEST_NUM','BEST_P_G','BEST_X_G','BEST_D_G', ...
    'HIS_D_P','HIS_D_G','HIS_D_TP', ...
    'T_NOW','path_abaqus','progress_log', ...
    '-v7.3');

save(fullfile(path_abaqus, name_data));

fprintf('\n==================== CMA-ES 续跑结束 ====================\n');
fprintf('最终已完成到第 %d 代\n', T_NOW);
fprintf('当前全局最优 DEL = %.8g  (T=%d, Phe=%d)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));
fprintf('checkpoint 已保存：%s\n', state_mat_file);
fprintf('=========================================================\n');
end