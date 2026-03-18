clear
clear global
clc
%所有单位均为mm，与实际保持一致（PI厚度 50微米 = 0.05mm ）
%%
%-------------------------------------------参数设定区
global GA FEA tar T Phe    
% -- lp 
% sketch被标记因为其仅在python中使用
%设置遗传算法参数 
    GA.T_MAX = 2;            %最大迭代数
    GA.N = 4;            %种群规模，必须为偶数，不再保留最优个体
    GA.TYPE = 4;          %待优化的变量类型个数，宽度与角度各算一类  【宽度、扭转角、随机旋转角度】
    GA.D = [ 11, 11, 1, 1];    %每种变量类型的粒子维度，即优化变量个数 这里添加了圆柱半径这个变量
    GA.Pc  = 0.6;     %交叉率
    GA.Pm  = 0.01;     %变异率
    GA.EXPECT =0.05;     %期望误差，用于判断是否终止迭代
    GA.BIT = 10000;  %数字的精度
%设置仿真计算参数    ： 给出的是较为通用的单条带、小周向跨度情况下的参数，如需修改，应在A_set_tar中设置修改                                                        
    FEA.multitask = 4;            %同时运行Job个数
    FEA.cpus      = 2 ;           %多核计算
    FEA.load_disturb   = 5e-7   ; %扰动载荷大小
    FEA.load_pressall  = 5e-8   ; %贴合载荷大小
    FEA.load_pressbond = 5e-8   ; %贴合载荷大小
    FEA.mesh =[ 0.5; 0.5];     %曲面和支撑条带的网格大小
    FEA.t = 0.05 ;   %条带厚度
    FEA.step_stabilize = [1e-6;1e-6;1e-6;1e-6;1e-6] ;   %稳定因子大小
    FEA.step_time = [ 5; 1; 1; 1; 5];  %存储各分析步的步长
    FEA.N = GA.N;
    FEA.penalty_1 = 1 ;
%%
%---------------------------------------------------------------目标构型设置
%参数设置，并计算目标构型的质心线插值点，误差仅计算坐标
%质心线形式：parameter equation
A_set_tar_equation; %定义：圆柱半径 r、条带数量、每条带的目标中心线参数方程、交叉/定位/键合信息、优化变量范围等
Fun_cal_para_tar_eqation;   %把目标参数方程离散成 tar.x/tar.y/tar.z，并计算弧长 tar.s、曲率 tar.K1 等
%离散点形式：discrete point
%A_set_tar_point_set;
%Fun_cal_para_tar_point_set;

%%
P = cell( GA.N, GA.TYPE, GA.T_MAX);   %建立个体库，存储每代、每个粒子的优化参数（位置）
%P构成：（个体数，变量类型数，总代数）：成员cell
%每个成员cell构成：（条带数，优化变量数）：成员float
for pn = 1:GA.N
    for ptype = 1:GA.TYPE
        P{ pn, ptype, 1} = round( rand( tar.line, GA.D(ptype) ) * GA.BIT) / GA.BIT;  %初始化P0
    end
end
T = 1; Phe = 1;%便于调试，初始化
%%
    fprintf('------------PSO开始-----------\n');  
    time_all = tic;

%初始化个体最优位置和最优值
BEST_P_P  = P(:,:,1);            %存储每个粒子迄今为止搜索到的最优位置  (个体极值)
BEST_D_P  = ones( GA.N, 1) *100;    %存储每个个体的最小误差(DEL)
%初始化全局最优位置和最优值
BEST_NUM  = zeros(1,2);       %存储最优个体信息 [代数T，个体号Phe]
BEST_P_G  = P(1,:,1)  ;       %存储整个种群迄今为止搜索到的最优位置  (全局极值)
BEST_D_G  = inf;              %存储全局最优的误差

HIS_D_P  = zeros( GA.T_MAX, GA.N); %存储所有个体的历史误差(DEL与DEL_theta)
HIS_D_G  = zeros( GA.T_MAX,1);   %存储适应度(DEL_G)变化趋势
HIS_D_TP = zeros( GA.T_MAX,2);   %存储各代对应的最佳个体TP值

%%
%-------------------------------------------存储设定
    ORDER = GA.ORDER ;      %存储路径编号
    Structure = GA.Structure;
    name_data =  ['mat-GA-', ORDER, Structure, '.mat'] ;
    name_oldfile = 'A-GA-' ;
    name_newfile = ['A-GA-',ORDER];
    
    path_matlab = cd('..') ; %进入上一级，同时将matlab目录存储至path_m
    copyfile(name_oldfile,name_newfile)     %复制存放python脚本的文件夹并用于计算文件存储
    cd([ '.\', name_newfile ] )
    path_abaqus = cd(path_matlab) ;  %返回matlab目录，同时将计算文件目录存储至path_abaqus
    
fprintf('初始化完成:');
toc(time_all)
%%
for T = 1 : GA.T_MAX     %                                                      T --- START
    %T=1; Phe = 1;   %单代调试
fprintf('------------第%d代迭代开始\n', T);  
time_T = tic;
    FEA.NUM_T = T;
    P(:,:,T) = Fun_set_special( P(:,:,T));
    %%
time_part = tic;  
    for Phe = 1:GA.N  %                                                       Phe_CAE --- START
    %% 
    name_cae_mat = ['mat-',num2str(T),'-',num2str(Phe),'.mat'];...                   mat-T-Phe
    %生成质心线2D信息和位移边界信息
    %生成Python所需的坐标信息（草图、路径）
    sketch = Fun_cal_para_2D( P(Phe,:,T) ) ; %sketch 这边获取abaqus建模需要的所有文件 存在sketch结构体中
    %%
    cd(path_abaqus);
    save(name_cae_mat);... 全存，仅用部分                                  mat-T-Phe 这边sketch的怎么给到mat呢
    cd(path_matlab);

    end   %                                                                Phe_CAE --- END
fprintf('[%d]草图参数计算完成 ：',Phe);  
toc(time_part)   
    %%
    %创建模型提交计算，并检查计算结果
    abaqus_run( path_abaqus );...                                   python_CAE(_Fold).py
    %%
fprintf('----------第%d代ODB处理开始----------\n',T);
time_part = tic;

    for Phe = 1:GA.N   %                                                      Phe_ODB --- START
        %%     
time_odb = tic;
        %提取odb数据，写入resule_i并输出读取状态
        result = abaqus_odb_m( path_abaqus );...python_ODB.py
toc(time_odb);
%%
        if isempty(result.rx{1}) == 0
            %计算相关参数
            %[ tar_xr, tar_yr, tar_zr, res_x, res_y, res_z] = ...
             result = Fun_cal_para_odb_m( result);...    Fun_cal_cir_odb_para            
            %计算误差函数δ
            %仅计算坐标
            del  = cal_delta_C_m( result);
            %del  = cal_delta_C( tar_xr, tar_yr, tar_zr, res_x, res_y, res_z, tar_center);
        else
            del = 1; %如果发生没有成功提取的情况发生 直接让del=1排除这种情况
        end

        %%
        HIS_D_P(T,Phe) = del; 
        %更新个体的极值与最优位置
        if del < BEST_D_P(Phe)
            BEST_D_P(Phe)    = del; %存下来最小的误差函数值
            BEST_P_P(Phe,:) = P(Phe,:,T); %存下来最好的误差函数值对应的基因
        end
        %更新全局的极值与最优位置，并记录最优个体的信息
        if del < BEST_D_G
            BEST_D_G   = del;
            BEST_P_G  = BEST_P_P(Phe,:);
            BEST_NUM = [T,Phe];
        end       
    end    %                                                               Phe_ODB --- END
    %%
    %记录适应度迭代历程
    HIS_D_G(T)    = BEST_D_G;
    HIS_D_TP(T,:) = BEST_NUM ;
    %判断是否达到要求。
    if BEST_D_G < GA.EXPECT && T>=30          
        fprintf('PSO complete! T=%d.\n',T);
        fprintf('Best DEL=:%g   \n ',BEST_D_G);
        fprintf('T & Phe : JOB-%d-%d \n',BEST_NUM(1),BEST_NUM(1));
        T_NOW = T;
        break;
    end
fprintf('----------第%d代ODB处理结束：',T);toc(time_part)
fprintf('进行第%d代更新\n',T);
    %HIS_DEL_P(Phe,:,T) = [del,0];
    %进行遗传操作--参数更新
    if T < GA.T_MAX                                   %未到最后一代则计算下一代参数
        P(:,:,T+1) = Fun_GA( HIS_D_P(T,:), P(:,:,T));
        %同时进行边界条件判断：所有的P位于 [0,1]
    end      
    %%
fprintf('-----------------------------第%d代结束-----------------------------\n',T);
fprintf('当前最佳误差 ：:%g    ',BEST_D_G);
fprintf('最佳个体信息 ： JOB-%d-%d \n', BEST_NUM(1), BEST_NUM(2));
fprintf('本代');
toc(time_T)
fprintf('总运行时间：');
toc(time_all)
end  ...                                                                   T --- END

%%
T_NOW = T;
%%
%提取最佳个体信息并存储
%tar.judge_type_line = ones( GA.TYPE, tar.line);     %补充
save_sketch = Fun_cal_para_2D( P(BEST_NUM(2),:,BEST_NUM(1)) ) ;
%读取对应odb
T = BEST_NUM(1);
Phe = BEST_NUM(2);
save_result = abaqus_odb_m( path_abaqus );...python_ODB.py                        
%计算相关参数
save_result = Fun_cal_para_odb_m( save_result);...    Fun_cal_cir_odb_para
save_his_del = HIS_D_G( 1 : T_NOW);
%--------------------------------------------------------------------------保存本次计算的数据
%%
%保存数据文件
path_matlab = pwd();
save(name_data);
fprintf('-----------------------------PSO 结束！共计完成%d代！--------------------------\n',T_NOW);
toc(time_all)
%%
Result_plot( save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM);
%%
if 1==1
    abaqus_odb_open( path_abaqus, BEST_NUM(1), BEST_NUM(2) );
end

%%
if 1 == 0
    fprintf('-----------------------------PSO 结束！共计完成%d代！--------------------------\n',GA.T_NOW);
    fprintf('Best DEL=:%g    ',BEST_D_G);
    fprintf('最佳个体信息 ： JOB-%d-%d \n', BEST_NUM(1), BEST_NUM(2));
    %绘制适应度曲线
    figure('color','w')
    plot   (save_his_del,'r')
    xlabel ('迭代次数')
    ylabel ('误差值')
    title  ('适应度进化曲线')
    set(gca,'yticklabel',get(gca,'ytick'));
end

