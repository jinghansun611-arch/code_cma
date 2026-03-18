clear
clear global
clc
% ============================================================
%   CMA-ES 多分辨率精修主入口（粗 -> 细）
%   目标：
%     - 先用 Stage-1（少控制点/少插值点）得到大概构型
%     - 再把【控制点 GA.D】和【插值点 tar.n】加密，用 CMA-ES 局部精修贴合目标
%
%   你只需要：   
%     1) 把 Fun_CMAES_multires_seed.m 放在本目录
%     2) 复制本脚本（或直接用本脚本）
%     3) 修改下面"用户只需要改这里"的 3 个变量：D_new / n_new / bestfile_stage1(可选)
% ============================================================

%% ------------------------------ 参数设定区 ------------------------------
global GA FEA tar T Phe

% ========== 0) 【用户只需要改这里】=========
% (1) 细化后的控制点数量（每个 TYPE 一项，长度必须等于 GA.TYPE=5）
D_new = [21,21,1,1,1];

% (2) 细化后的插值点数量（建议 2~4 倍于粗搜的 tar.n=51，比如 101/151；最好是奇数）
n_new = 101;

% (3) Stage-1 的 best_stage1.mat 路径：
%     - 如果留空，本脚本会在根目录下自动寻找最新的 A-CMAES-*/best_stage1.mat
bestfile_stage1 = 'D:\CMAES_inverse_design_pack\code_cma_withtwist\A-CMAES-7S-1-80-25-refine1\best_stage1.mat';
% ==========================================

% ========== 1) 精修建议参数（可按预算改） ==========
GA.T_MAX  = 50;          % 精修迭代数（建议 20~60）
GA.N      = 20;          % 每代个体数（并行多就调大；一般 10~25）
GA.TYPE   = 5;           % 不改
GA.D      = [11,11,1,1,1];% 这里先给个默认，后面会被 D_new 覆盖
GA.EXPECT = 0.01;        % 期望误差
GA.BIT    = 10000;

GA.CMA.sigma0 = 0.10;    % 先给个默认，后面会被 multires 的建议覆盖
GA.CMA.x0     = [];
GA.CMA.elite_inject = false;
GA.CMA.elite_warmup = 0; % 这个数值是n意思就是 前n带不开启精英注入 如果设置为0就是一直开启

% --- 仿真计算参数（你也可以不动，A_set_tar_equation 里会覆盖一些） ---
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

%% ------------------------------ 自动定位 Stage-1 best file ------------------------------
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

%% ------------------------------ 多分辨率种子：粗 -> 细 ------------------------------
% 这一步会：
%   - 把 GA.D 改成 D_new
%   - 把 tar.n 改成 n_new 并重算 tar.x/y/z/s/K1/location/cross/center
%   - 把 stage1 的 BEST_P_G 重采样到新控制点数
%   - build 新 map，并输出 x0/sigma0 建议
seed = Fun_CMAES_multires_seed(bestfile_stage1, D_new, n_new);

GA.D          = seed.D_new;
tar.n         = seed.n_new;
map           = seed.map;
GA.CMA.x0     = seed.x0;
GA.CMA.sigma0 = seed.sigma0_suggest;

% 同步 FEA.N（因为你可能改了 GA.N）
FEA.N = GA.N;

fprintf('\n【MultiRes-Refine】Stage1: %s\n', bestfile_stage1);
fprintf('【MultiRes-Refine】New GA.D = %s, New tar.n = %d\n', mat2str(GA.D), tar.n);
fprintf('【MultiRes-Refine】New nvars = %d, x0 loaded, sigma0 = %.4g\n\n', map.nvars, GA.CMA.sigma0);

%% ------------------------------ CMA-ES 初始化 ------------------------------
x0 = GA.CMA.x0;
sigma0 = GA.CMA.sigma0;

cma = Fun_CMAES_init(map.nvars, GA.N, x0, sigma0);
cma.elite_inject = GA.CMA.elite_inject;

%% ------------------------------ 个体库/历史量 ------------------------------
P = cell(GA.N, GA.TYPE, GA.T_MAX);

BEST_NUM = zeros(1,2);   % [T, Phe]
BEST_P_G = P(1,:,1);
BEST_D_G = inf;
BEST_X_G = [];

HIS_D_P  = zeros(GA.T_MAX, GA.N);
HIS_D_G  = zeros(GA.T_MAX, 1);
HIS_D_TP = zeros(GA.T_MAX, 2);

% 先 ask 一代，填充 P(:,:,1)
[X, cma] = Fun_CMAES_ask(cma);
for pn = 1:GA.N
    P(pn,:,1) = Fun_CMAES_vec2PP(X(pn,:), map);
end

T = 1; Phe = 1;

%% ------------------------------ 存储设定 ------------------------------
ORDER = GA.ORDER;
Structure = GA.Structure;
name_data = ['mat-CMAES-', ORDER, Structure, '-multires-refine.mat'];

name_oldfile = 'A-GA-';

% 输出目录名：A-CMAES-<ORDER>-refineMR-Dxx-...-nxxx
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
path_abaqus = cd(path_matlab);

fprintf('------------CMA-ES 多分辨率精修开始-----------\n');
fprintf('Workdir : %s\n', path_abaqus);
time_all = tic;
fprintf('初始化完成:');
toc(time_all)

%% =============================== 主循环 ===============================
for T = 1:GA.T_MAX

    fprintf('------------第%d代迭代开始\n', T);
    time_T = tic;
    FEA.NUM_T = T;

    % 对称/全局变量规则
    P(:,:,T) = Fun_set_special(P(:,:,T));

    % 可行域 x
    X_feasible = zeros(GA.N, map.nvars);
    for pn = 1:GA.N
        X_feasible(pn,:) = Fun_CMAES_PP2vec(P(pn,:,T), map)';
    end
    X_T = X_feasible;

    %% -------- 1) 生成每个个体的 mat --------
    time_part = tic;
    for Phe = 1:GA.N
        name_cae_mat = ['mat-', num2str(T), '-', num2str(Phe), '.mat'];
        sketch = Fun_cal_para_2D(P(Phe,:,T));
        cd(path_abaqus);
        save(name_cae_mat);
        cd(path_matlab);
    end
    fprintf('[%d]草图参数计算完成 ：', Phe);
    toc(time_part)

    %% -------- 2) Abaqus 批量计算 --------
    abaqus_run(path_abaqus);

    %% -------- 3) 读取 odb + 计算误差 del --------
    fprintf('----------第%d代ODB处理开始----------\n', T);
    time_part = tic;

    for Phe = 1:GA.N
        result = abaqus_odb_m(path_abaqus);
        if isempty(result.rx{1}) == 0
            result = Fun_cal_para_odb_m(result);
            del = cal_delta_C_m(result);
        else
            del = 1;
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

    fprintf('----------第%d代ODB处理结束：', T);
    toc(time_part)

    fprintf('当前最佳误差 ：%g    ', BEST_D_G);
    fprintf('最佳个体信息 ： JOB-%d-%d \n', BEST_NUM(1), BEST_NUM(2));

    %% -------- 4) 提前停止 --------
    if BEST_D_G < GA.EXPECT
        fprintf('CMA-ES complete! T=%d. Best DEL=%g\n', T, BEST_D_G);
        break;
    end

    %% -------- 5) CMA-ES 更新并 ask 下一代 --------
    if T < GA.T_MAX

        cma = Fun_CMAES_tell(cma, X_feasible, HIS_D_P(T,:)');

        % 精英注入调度
        if (T+1) <= GA.CMA.elite_warmup
            cma.elite_inject = false;
        else
            cma.elite_inject = true;
        end

        fprintf('T=%d -> ask T+1, elite_inject=%d, sigma=%g\n', T, cma.elite_inject, cma.sigma);

        [X, cma] = Fun_CMAES_ask(cma);
        for pn = 1:GA.N
            P(pn,:,T+1) = Fun_CMAES_vec2PP(X(pn,:), map);
        end
    end

    fprintf('本代');
    toc(time_T)
    fprintf('总运行时间：');
    toc(time_all)
    fprintf('-----------------------------第%d代结束-----------------------------\n', T);
end

T_NOW = T;

%% ------------------------------ 保存最优解（保存到 refineMR 目录） ------------------------------
% 为了兼容你原来的 refine/续跑逻辑，这里仍保存为 best_stage1.mat（但在新的 refineMR 文件夹里）
bestfile = fullfile(path_abaqus, 'best_stage1.mat');
save(bestfile, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','-v7.3');

% 另外再保存一份更明确的名字
bestfile2 = fullfile(path_abaqus, 'best_multires_refine.mat');
save(bestfile2, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','seed','-v7.3');

fprintf('Saved best to: %s (DEL=%g at T=%d, Phe=%d)\n', bestfile, BEST_D_G, BEST_NUM(1), BEST_NUM(2));
fprintf('Saved best to: %s\n', bestfile2);

%% ------------------------------ 保存/绘图（保持你原工程输出） ------------------------------
save_sketch = Fun_cal_para_2D(P(BEST_NUM(2),:,BEST_NUM(1)));
T = BEST_NUM(1);
Phe = BEST_NUM(2);
save_result = abaqus_odb_m(path_abaqus);
save_result = Fun_cal_para_odb_m(save_result);
save_his_del = HIS_D_G(1:T_NOW);

path_matlab = pwd();
save(name_data);

fprintf('-----------------------------CMA-ES 多分辨率精修结束！共计完成%d代！--------------------------\n', T_NOW);
toc(time_all)

Result_plot(save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM);

if 1==1
    abaqus_odb_open(path_abaqus, BEST_NUM(1), BEST_NUM(2));
end