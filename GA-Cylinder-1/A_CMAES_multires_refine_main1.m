clear
clear global
clc
% ============================================================
%   CMA-ES 多分辨率精修主入口（粗 -> 细）  v2
%
%   现在流程
%     - 用粗分辨率（例：GA.D=[11 11 1 1 1], tar.n=51）先找到大概构型
%     - 然后把控制点和插值点加密（例：GA.D=[21 21 1 1 1], tar.n=101）
%       从“粗搜最优”出发继续精修贴合。
%
%   本脚本的改动重点：
%     1) x0 = 由 stage1 best(PP) 重采样到新 D_new 后得到，并在新 tar.n 下重新离散目标。
%        （并额外做一次“可行化” set_special，确保注入的精英=实际送入仿真的精英）
%     2) 可选：先单独评估一次 seed（T=0, Phe=1）得到“新评价下的初始 del0”
%        用它初始化 cma.fbest，从第1代就能精英注入。
%     3) 输出信息增强：每代打印 del 列表、均值/中位数/最优、与上一代对比、全局最优等。
%
%   使用方式：
%     - 把本文件放到 GA-Cylinder-1 目录
%     - 修改下面【用户只需要改这里】的 3 个变量
%     - 直接运行本脚本
% ============================================================

%% ------------------------------ 全局变量 ------------------------------
global GA FEA tar T Phe

%% ========== 0) 【用户只需要改这里】=========
% (1) 细化后的控制点数量（每个 TYPE 一项，长度必须等于 GA.TYPE=5）
D_new = [21,21,1,1,1];

% (2) 细化后的插值点数量（建议 2~4 倍于粗搜的 tar.n，比如 51->101/151；最好用奇数）
n_new = 101;

% (3) Stage-1 的 best_stage1.mat 路径（你也可以写为空，让程序自动找 root_dir 下最新的一个）
% 例：'D:\\CMAES_inverse_design_pack\\code_cma_withtwist\\A-CMAES-7S-1-80-25-refine1\\best_stage1.mat'
bestfile_stage1 = 'D:\CMAES_inverse_design_pack\code_cma_withtwist\A-CMAES-7S-1-80-25-refine1\best_stage1.mat';
% ==========================================

%% ========== 1) 精修建议参数（按预算改） ==========
GA.T_MAX  = 60;            % 精修迭代数（建议 20~60）
GA.N      = 25;            % 每代个体数（并行多就调大；一般 10~25）
GA.TYPE   = 5;             % 不改
GA.D      = [11,11,1,1,1]; % 先占位；稍后会被 D_new 覆盖
GA.EXPECT = 0.01;          % 期望误差
GA.BIT    = 10000;

% 精英注入：你说你想“每一代都开”，那就：elite_warmup=0
GA.CMA.elite_inject = true;
GA.CMA.elite_warmup = 0;

% sigma0：先给默认，稍后用 multires_seed 的建议值覆盖
GA.CMA.sigma0 = 0.10;
GA.CMA.x0     = [];

% 是否先评估一次 seed（得到新评价下的 del0，并用于初始化 cma.fbest）
CFG.eval_seed_once = true;

% 每代是否打印所有个体 del（很直观，但比较长）
CFG.print_each_del = true;

%% --- 仿真计算参数（照搬你的默认） ---
FEA.multitask = 5;              % 同时运行 Job 个数
FEA.cpus      = 2;              % 每个 Job 使用的 CPU 核数
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
% 这一步定义 tar.funx/funy/funz、tar.range、tar.judge、tar.symmetry、GA.ORDER 等
A_set_tar_equation;

%% ------------------------------ 定位 stage1 bestfile ------------------------------
path_GA   = pwd;                 % 当前运行目录（通常是 ...\code_cma\GA-Cylinder-1）
root_dir  = fileparts(path_GA);  % 上一级目录（...\code_cma）

if isempty(bestfile_stage1)
    cand = dir(fullfile(root_dir,'A-CMAES-*','best_stage1.mat'));
    if isempty(cand)
        error(['未找到任何 Stage-1 best_stage1.mat。\n' ...
               '请先跑粗搜 A_CMAES_muity_main.m，或手动在本脚本里指定 bestfile_stage1。']);
    end
    [~,idx] = max([cand.datenum]);
    bestfile_stage1 = fullfile(cand(idx).folder, cand(idx).name);
end

if ~exist(bestfile_stage1,'file')
    error('找不到 Stage-1 best_stage1.mat：%s', bestfile_stage1);
end

S_stage1 = load(bestfile_stage1, 'GA','tar','BEST_D_G','BEST_NUM');

%% ------------------------------ 多分辨率 seed：粗 -> 细 ------------------------------
% 关键：用 stage1 的 tar 定义来保证“完全顺着之前的目标继续做”（避免你中途改了 set_tar_equation）
seed = Fun_CMAES_multires_seed(bestfile_stage1, D_new, n_new, 'UseStage1Tar', true);

GA.D          = seed.D_new;
tar.n         = seed.n_new;
map           = seed.map;
GA.CMA.x0     = seed.x0;
GA.CMA.sigma0 = seed.sigma0_suggest;

% 再做一次“可行化”：保证注入的精英 = 实际送进 Abaqus 的精英（避免 set_special 改写导致退化）
PP0 = Fun_CMAES_vec2PP(GA.CMA.x0(:)', map);
PP0 = Fun_set_special(PP0);
x0_feasible = Fun_CMAES_PP2vec(PP0, map);
GA.CMA.x0 = x0_feasible(:);

FEA.N = GA.N;

%% ------------------------------ 输出目录（refineMR） ------------------------------
ORDER = GA.ORDER;
Structure = GA.Structure;
name_data = ['mat-CMAES-', ORDER, Structure, '-multires-refine.mat'];

name_oldfile = 'A-GA-';
Dtag = regexprep(mat2str(GA.D),'[\[\] ,]','-');
Dtag = regexprep(Dtag,'-+','-');
Dtag = regexprep(Dtag,'^-|-$','');
name_newfile = sprintf('A-CMAES-%s-refineMR-D%s-n%d', ORDER, Dtag, tar.n);

path_matlab = cd('..');
if exist(name_newfile,'dir')
    rmdir(name_newfile,'s');
end
copyfile(name_oldfile, name_newfile);
cd(['.\', name_newfile]);
path_abaqus = cd(path_matlab);   % refineMR 的绝对路径

%% ------------------------------ 开场说明 ------------------------------
fprintf('\n================ 多分辨率精修（粗 → 细） ================\n');
fprintf('Stage1 文件: %s\n', bestfile_stage1);
if isfield(S_stage1,'GA') && isfield(S_stage1.GA,'D')
    fprintf('粗分辨率: GA.D=%s\n', mat2str(S_stage1.GA.D));
end
if isfield(S_stage1,'tar') && isfield(S_stage1.tar,'n')
    fprintf('粗分辨率: tar.n=%d\n', S_stage1.tar.n);
end
fprintf('细分辨率: GA.D=%s, tar.n=%d, nvars=%d\n', mat2str(GA.D), tar.n, map.nvars);
fprintf('CMA-ES: N=%d/代, T_MAX=%d, sigma0=%.4g, 精英注入=%d(每代)\n', GA.N, GA.T_MAX, GA.CMA.sigma0, GA.CMA.elite_inject);
fprintf('提示：插值点 tar.n 加密后，即使构型不变，del 也可能变大——这是“评价更细”导致的，属正常现象。\n');

fprintf('RefineMR 工作目录: %s\n', path_abaqus);
time_all = tic;

%% ------------------------------ 可选：先评估 seed 的 del0（新评价下） ------------------------------
SEED_DEL0 = NaN;
if CFG.eval_seed_once
    fprintf('\n[Seed评估] 先用新分辨率对 seed(x0) 单独评估一次，得到 del0（用于基线 & 精英注入）。\n');
    % 临时把 GA.N 改成 1，只算 1 个 Job
    GA_N_bak  = GA.N;
    FEA_N_bak = FEA.N;
    GA.N  = 1;
    FEA.N = 1;

    % 设定 T=0, Phe=1 对应的 mat / Job
    T = 0; FEA.NUM_T = T; Phe = 1;

    name_cae_mat = ['mat-', num2str(T), '-', num2str(Phe), '.mat'];
    sketch = Fun_cal_para_2D(PP0);
    cd(path_abaqus);
    save(name_cae_mat);
    cd(path_matlab);

    abaqus_run(path_abaqus);

    result = abaqus_odb_m(path_abaqus);
    if isempty(result.rx{1}) == 0
        result = Fun_cal_para_odb_m(result);
        SEED_DEL0 = cal_delta_C_m(result);
    else
        SEED_DEL0 = 1;
    end

    fprintf('[Seed评估] del0(新评价) = %.6g\n', SEED_DEL0);

    % 恢复 GA.N
    GA.N  = GA_N_bak;
    FEA.N = FEA_N_bak;
end

%% ------------------------------ CMA-ES 初始化 ------------------------------
x0 = GA.CMA.x0;
sigma0 = GA.CMA.sigma0;

cma = Fun_CMAES_init(map.nvars, GA.N, x0, sigma0);

% 让精英注入从第1代就生效：需要 cma.fbest 是有限数
cma.elite_inject = GA.CMA.elite_inject;
cma.xbest = x0;
if isfinite(SEED_DEL0)
    cma.fbest = SEED_DEL0;
else
    % 没评估 seed，也要给一个“有限的大数”来触发注入
    cma.fbest = 1e99;
end

%% ------------------------------ 个体库/历史量 ------------------------------
P = cell(GA.N, GA.TYPE, GA.T_MAX);

BEST_NUM = zeros(1,2);   % [T, Phe]
BEST_P_G = P(1,:,1);
BEST_D_G = inf;
BEST_X_G = [];

HIS_D_P  = zeros(GA.T_MAX, GA.N);
HIS_D_G  = zeros(GA.T_MAX, 1);
HIS_D_TP = zeros(GA.T_MAX, 2);

% ask 第1代
T = 1; Phe = 1;
[X_raw, cma] = Fun_CMAES_ask(cma);
for pn = 1:GA.N
    P(pn,:,1) = Fun_CMAES_vec2PP(X_raw(pn,:), map);
end

fprintf('\n初始化完成，总计耗时：');
toc(time_all)

%% =============================== 主循环 ===============================
for T = 1:GA.T_MAX

    fprintf('\n================ 第 %d / %d 代 ================\n', T, GA.T_MAX);
    time_T = tic;
    FEA.NUM_T = T;

    % 1) 可行化（对称/全局约束）
    P(:,:,T) = Fun_set_special(P(:,:,T));

    % 2) 抽取可行域内的 X（让 CMA-ES 在可行域里学习）
    X_feasible = zeros(GA.N, map.nvars);
    for pn = 1:GA.N
        X_feasible(pn,:) = Fun_CMAES_PP2vec(P(pn,:,T), map)';
    end

    % Debug：精英是否被 set_special 改写（改写越大，越可能出现你说的“精英注入但 del 变大”）
    elite_changed_norm = norm(X_raw(1,:) - X_feasible(1,:));
    fprintf('参数概览：sigma=%.4g | elite_inject=%d | elite改变幅度(norm)=%.3e\n', cma.sigma, cma.elite_inject, elite_changed_norm);

    %% -------- 3) 生成每个个体的 mat --------
    time_part = tic;
    for Phe = 1:GA.N
        name_cae_mat = ['mat-', num2str(T), '-', num2str(Phe), '.mat'];
        sketch = Fun_cal_para_2D(P(Phe,:,T));
        cd(path_abaqus);
        save(name_cae_mat);
        cd(path_matlab);
    end
    fprintf('生成 mat 完成（%d 个体），耗时：', GA.N);
    toc(time_part)

    %% -------- 4) Abaqus 批量计算 --------
    abaqus_run(path_abaqus);

    %% -------- 5) 读取 odb + 计算 del --------
    fprintf('开始读取 ODB 并计算 del...\n');
    time_part = tic;

    fail_cnt = 0;
    for Phe = 1:GA.N
        result = abaqus_odb_m(path_abaqus);
        if isempty(result.rx{1}) == 0
            result = Fun_cal_para_odb_m(result);
            del = cal_delta_C_m(result);
        else
            del = 1;
            fail_cnt = fail_cnt + 1;
        end

        HIS_D_P(T, Phe) = del;

        if del < BEST_D_G
            BEST_D_G = del;
            BEST_P_G = P(Phe,:,T);
            BEST_X_G = X_feasible(Phe,:).';
            BEST_NUM = [T, Phe];
        end
    end

    HIS_D_G(T) = BEST_D_G;
    HIS_D_TP(T,:) = BEST_NUM;

    fprintf('本代 ODB 处理耗时：');
    toc(time_part)

    %% -------- 6) 每代输出（更直观的中文信息） --------
    del_vec = HIS_D_P(T,:);
    [best_gen, idx_best_gen] = min(del_vec);
    mean_gen = mean(del_vec);
    med_gen  = median(del_vec);

    if CFG.print_each_del
        fprintf('\n本代各个体 del（越小越好）：\n');
        for p = 1:GA.N
            fprintf('  %2d:% .6g', p, del_vec(p));
            if mod(p,6) == 0 || p == GA.N
                fprintf('\n');
            end
        end
    end

    fprintf('本代统计：最优=%.6g (Phe=%d) | 均值=%.6g | 中位数=%.6g | 失败=%d/%d\n', ...
        best_gen, idx_best_gen, mean_gen, med_gen, fail_cnt, GA.N);

    % 与上一代对比（用“本代最优”对比“上一代最优”）
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

    % 全局最优
    fprintf('全局最优：BEST_D_G=%.6g  (出现于 第%d代 / 第%d个体)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));

    % 若评估了 seed，再给一个“相对 seed 的提升”
    if isfinite(SEED_DEL0)
        fprintf('相对 seed 基线(del0=%.6g) 改善：%+.6g\n', SEED_DEL0, SEED_DEL0 - BEST_D_G);
    end

    %% -------- 7) 提前停止 --------
    if BEST_D_G < GA.EXPECT
        fprintf('\n达到期望误差 GA.EXPECT=%.6g，提前停止。\n', GA.EXPECT);
        break;
    end

    %% -------- 8) CMA-ES 更新并 ask 下一代 --------
    if T < GA.T_MAX
        cma = Fun_CMAES_tell(cma, X_feasible, del_vec');

        % 你要“每代都精英注入”：这里就不做 warmup 逻辑
        cma.elite_inject = GA.CMA.elite_inject;

        fprintf('CMA-ES 更新完成：下一代采样前 sigma=%.4g | 当前历史最优 fbest=%.6g\n', cma.sigma, cma.fbest);

        [X_raw, cma] = Fun_CMAES_ask(cma);
        for pn = 1:GA.N
            P(pn,:,T+1) = Fun_CMAES_vec2PP(X_raw(pn,:), map);
        end
    end

    fprintf('本代总耗时：');
    toc(time_T)
    fprintf('累计耗时：');
    toc(time_all)
end

T_NOW = T;

%% ------------------------------ 保存最优解 ------------------------------
bestfile = fullfile(path_abaqus, 'best_stage1.mat');
save(bestfile, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','seed','SEED_DEL0','-v7.3');%这是为了兼容原来的流程

bestfile2 = fullfile(path_abaqus, 'best_multires_refine.mat');
save(bestfile2, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','seed','SEED_DEL0','-v7.3');%这是为了名字更明确 和上面的一模一样

fprintf('\n已保存最优解：\n  %s\n  %s\n', bestfile, bestfile2);

%% ------------------------------ 保存/绘图（保持原工程输出） ------------------------------
save_sketch = Fun_cal_para_2D(P(BEST_NUM(2),:,BEST_NUM(1)));
T = BEST_NUM(1);
Phe = BEST_NUM(2);
save_result = abaqus_odb_m(path_abaqus);
save_result = Fun_cal_para_odb_m(save_result);
save_his_del = HIS_D_G(1:T_NOW);
% ===== 追加保存绘图变量到 best 文件（以后只 load best 就能画图）=====
try
    save(bestfile,  'save_sketch','save_result','save_his_del','-append');
    save(bestfile2, 'save_sketch','save_result','save_his_del','-append');
    fprintf('已将 save_sketch/save_result/save_his_del 追加保存到 best 文件。\n');
catch ME
    warning('multires:AppendSaveFailed', '追加保存 save_* 到 best 文件失败：%s', ME.message);
end

path_matlab = pwd();
save(name_data);

fprintf('\n================ 多分辨率精修结束：共完成 %d 代 ================\n', T_NOW);
fprintf('最终 BEST_D_G=%.6g  (第%d代/第%d个体)\n', BEST_D_G, BEST_NUM(1), BEST_NUM(2));

Result_plot(save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM);
%如果想要单独画图 运行：
% load('D:\...\best_multires_refine.mat', ...
%      'save_sketch','save_result','save_his_del','BEST_D_G','BEST_NUM');
% Result_plot(save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM

if 1==1
    abaqus_odb_open(path_abaqus, BEST_NUM(1), BEST_NUM(2));
end
