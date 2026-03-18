clear
clear global
clc
% 所有单位均为 mm，与实际保持一致（PI 厚度 50 微米 = 0.05 mm）
%
% ================================
%   CMA-ES 主入口（替代 Fun_GA）
% ================================
%
% 你原始流程：
%   生成 P → Fun_cal_para_2D → 写 mat → Abaqus 批量计算 → 读 odb → cal_delta_C_m
%   然后用 Fun_GA(轮盘赌+交叉+变异) 产生下一代。
%
% 本脚本做的事：
%   - 保留原来的“建模/提交/读取/误差”管线不变
%   - 把“产生下一代”从 GA 替换成 CMA-ES（连续变量自适应协方差进化策略）
%
% 为什么 CMA-ES 更适合你这个问题：
%   - 目标函数来自 Abaqus 黑盒、不可导、可能有噪声/失败个体
%   - 参数是连续的 [0,1]，且维度中等（几十维）
%   - CMA-ES 能用协方差学习变量耦合关系，通常比 GA/PSO 更省评估次数

%% ------------------------------ 参数设定区 ------------------------------
global GA FEA tar T Phe

% --- 优化器参数（CMA-ES 用 GA.N 作为种群规模 λ）
GA.T_MAX = 1;         % 最大迭代数 30to60
GA.N     = 25;         % 每代个体数（= CMA-ES 的 λ），不再要求偶数
GA.TYPE  = 4;          % 变量类型\数（与你工程一致）
%GA.D     = [11,11,1,1,1];% 每类变量的控制点个数
GA.D     = [21,21,1,1,1];
GA.EXPECT = 0.005;      % 期望误差，用于提前停止 0.05to0.08
GA.BIT    = 10000;     % 数字精度（保留你原工程习惯）

% --- CMA-ES 初始超参数（在 [0,1] 变量空间）
GA.CMA.sigma0 = 0.28;  % 初始步长（越大探索越广，越小收敛越稳；建议 0.15~0.35 试）0.25to0.16
GA.CMA.x0     = [];    % 初始均值（留空=0.5；如果你有“经验初值”，可以填 nvars×1 向量）
GA.CMA.elite_inject = true; % 是否把历史最优注入到每代种群（更稳）原来是true
GA.CMA.elite_warmup = 10;   % 前10代不注入精英（第11代开始注入）


% --- 仿真计算参数（如需针对某构型重置，在 A_set_tar_equation 里写覆盖）
FEA.multitask = 4;              % 同时运行 Job 个数
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
% 参数设置 + 目标质心线离散
A_set_tar_equation;
Fun_cal_para_tar_eqation;

%% ------------------------------ CMA-ES 初始化 ------------------------------
% 建立 x ↔ PP 的映射（自动识别 judge/symmetry/lambda 全局变量）
map = Fun_CMAES_build_map();

% 初始化 CMA-ES 状态
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
cma.elite_inject = GA.CMA.elite_inject;
% ===== 精英注入调度：第1~elite_warmup代不注入 =====
cma.elite_inject = false;   % 保证第1代采样不注入
% ===============================================

%% ------------------------------ 个体库/历史量 ------------------------------
P = cell(GA.N, GA.TYPE, GA.T_MAX);

BEST_NUM = zeros(1,2);   % [T, Phe]
BEST_P_G = P(1,:,1);
BEST_D_G = inf;
BEST_X_G = [];           % 对应的自由变量 x（0~1, nvars×1）

HIS_D_P  = zeros(GA.T_MAX, GA.N);
HIS_D_G  = zeros(GA.T_MAX, 1);
HIS_D_TP = zeros(GA.T_MAX, 2);

% 先 ask 一代，填充 P(:,:,1)
[X, cma] = Fun_CMAES_ask(cma);
for pn = 1:GA.N
    P(pn,:,1) = Fun_CMAES_vec2PP(X(pn,:), map);
end

T = 1; Phe = 1;

%% ------------------------------ 存储设定（保持你原工程逻辑） ------------------------------
ORDER = GA.ORDER;
Structure = GA.Structure;
name_data = ['mat-CMAES-', ORDER, Structure, '.mat'];
name_oldfile = 'A-GA-';
name_newfile = ['A-CMAES-', ORDER];

path_matlab = cd('..');
copyfile(name_oldfile, name_newfile);
cd(['.\', name_newfile]);
path_abaqus = cd(path_matlab);

fprintf('------------CMA-ES开始-----------\n');
time_all = tic;
fprintf('初始化完成:');
toc(time_all)

%% =============================== 主循环 ===============================
for T = 1:GA.T_MAX
    fprintf('------------第%d代迭代开始\n', T);
    time_T = tic;
    FEA.NUM_T = T;

    % 先把对称/全局规则应用到本代 P（保证送进 Abaqus 的都是可行解）
    P(:,:,T) = Fun_set_special(P(:,:,T));

    % 用于 CMA-ES 更新的“可行域 x”（在 Fun_set_special 之后再抽取）
    X_feasible = zeros(GA.N, map.nvars);
    for pn = 1:GA.N
        X_feasible(pn,:) = Fun_CMAES_PP2vec(P(pn,:,T), map)';
    end
    X_T = X_feasible;   % 备份本代用于评估的 X（用于记录BEST_X_G）

    %% -------- 1) 生成每个个体的 mat（sketch/bc 等） --------
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
            del = 1; % 失败个体直接给一个较差的目标值
        end

        HIS_D_P(T, Phe) = del;

        if del < BEST_D_G
            BEST_D_G = del;
            BEST_P_G = P(Phe,:,T);
            BEST_X_G = X_T(Phe,:).';   % ★新增：保存自由变量最优（nvars×1）
            BEST_NUM = [T, Phe];
        end
    end

    HIS_D_G(T) = BEST_D_G;
    HIS_D_TP(T,:) = BEST_NUM;

    fprintf('----------第%d代ODB处理结束：', T);
    toc(time_part)

    fprintf('当前最佳误差 ：%g    ', BEST_D_G);
    fprintf('最佳个体信息 ： JOB-%d-%d \n', BEST_NUM(1), BEST_NUM(2));

    %% ====== 新增：本代误差全列表 + 与上一代对比 ======
    del_vec = HIS_D_P(T,:);   % 本代每个体的 del（1×N）

    %打印本代所有个体 del（每行6个）
    fprintf('本代(T=%d) 各个体 DEL：\n', T);
    for p = 1:GA.N
        fprintf('  %2d:% .6g', p, del_vec(p));
        if mod(p,6) == 0 || p == GA.N
            fprintf('\n');
        end
    end
    
    %本代最优（只看本代，不是全局）
    [best_gen, idx_best_gen] = min(del_vec);
    fprintf('本代最优 DEL = %.6g (Phe=%d)\n', best_gen, idx_best_gen);
    
    %与上一代最优对比（提升=上一代最优 - 本代最优）
    if T > 1
        prev_best = min(HIS_D_P(T-1,:));
        improve = prev_best - best_gen;   % >0 表示变好
        if isfinite(prev_best) && prev_best > 0
            improve_pct = 100 * improve / prev_best;
            fprintf('较上一代最优(%.6g) 变化：%+.6g  (%+.2f%%)\n', prev_best, improve, improve_pct);
        else
            fprintf('较上一代最优(%.6g) 变化：%+.6g\n', prev_best, improve);
        end
    else
        fprintf('第一代：无上一代可对比。\n');
    end
    
    fprintf('----------------------------------------------------\n');
    %% =====================================================

    % -------- 4) 提前停止 --------
    if BEST_D_G < GA.EXPECT
        fprintf('CMA-ES complete! T=%d. Best DEL=%g\n', T, BEST_D_G);
        break;
    end

    % -------- 5) CMA-ES 更新并 ask 下一代 --------
    if T < GA.T_MAX
        
        cma = Fun_CMAES_tell(cma, X_feasible, HIS_D_P(T,:)');

        % ===== 精英注入调度：前 elite_warmup 代不注入，从下一代开始注入 =====
        % 当前 T 代评估完以后，会 ask 生成 T+1 代
        if (T+1) <= GA.CMA.elite_warmup
            cma.elite_inject = false;
        else
            cma.elite_inject = true;
        end
        % =====================================================================

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
end %一次大循环到这里结束

T_NOW = T;
% -------- 保存最优解，供精修阶段读取 --------
bestfile = fullfile(path_abaqus, 'best_stage1.mat');  % 存到 A-CMAES-XX
save(bestfile, 'BEST_P_G','BEST_X_G','BEST_D_G','BEST_NUM','GA','tar','map','-v7.3');
fprintf('Saved best to: %s (DEL=%g at T=%d, Phe=%d)\n', bestfile, BEST_D_G, BEST_NUM(1), BEST_NUM(2));


%% ------------------------------ 保存/绘图（保持你原工程输出） ------------------------------
save_sketch = Fun_cal_para_2D(P(BEST_NUM(2),:,BEST_NUM(1)));
T = BEST_NUM(1);
Phe = BEST_NUM(2);
save_result = abaqus_odb_m(path_abaqus);
save_result = Fun_cal_para_odb_m(save_result);
save_his_del = HIS_D_G(1:T_NOW);

path_matlab = pwd();
save(name_data);
fprintf('-----------------------------CMA-ES 结束！共计完成%d代！--------------------------\n', T_NOW);
toc(time_all)

Result_plot(save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM);

if 1==1
    abaqus_odb_open(path_abaqus, BEST_NUM(1), BEST_NUM(2));
end
