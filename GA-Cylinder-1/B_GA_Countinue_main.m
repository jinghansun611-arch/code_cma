clear
clc
%%
%填写需要拓展的算例mat文件名  和  拓展要求
Before_mat = 'mat-GA-1-1.mat';
T_NEW = 40;

%%
time_all = tic;
%进行数组扩充

path_m = pwd();%记录matlab目录
%%
T_OLD = T_MAX;
T_MAX = T_NEW;
clear('T_NEW');

P_OLD = P;
P = zeros( N, D, T_MAX);
P( :, :, 1:T_OLD) = P_OLD;

HIS_DEL_P_OLD = HIS_DEL_P;
HIS_DEL_P = zeros( T_MAX, N);
HIS_DEL_P(1:T_OLD ,:) = HIS_DEL_P_OLD;

HIS_DEL_G_OLD = HIS_DEL_G; 
HIS_DEL_G = zeros(T_MAX,1); 
HIS_DEL_G(1:T_OLD) = HIS_DEL_G_OLD; 
%%
%更新得到第 T_OLD + 1 的参数
    P(:,:,T_OLD+1) = Fun_GA( N, D, HIS_DEL_P(T_OLD,:), P(:,:,T_OLD), AG_Pc, AG_Pm );
    %同时进行边界条件判断：所有的P位于 [0,1] 
%%
fprintf('重载完成:');
toc(time_all)
%%
for T = T_OLD+1:T_MAX     %                                                      T --- START
    %T=31;   %单代调试
time_T = tic;
    charT        = num2str(T);
    %%
    for Phe = 1:N  %                                                       Phe_CAE --- START
    %%
fprintf('----Phe-%d ：',Phe);  time_part = tic;
    charPhe      = num2str(Phe);    
    name_cae_mat = ['mat-',charT,'-',charPhe,'.mat'];...                   mat-T-Phe
    %生成质心线2D信息和位移边界信息
    %生成Python所需的坐标信息（草图、路径）
    [ sketch_1, sketch_2, sketch_3, sketch_4, sketch_path, sketch_div, sketch_point, sketch_center, sketch_U, sketch_tarline, sketch_bonding] = ...
        Fun_cal_para_2D( N, D, tar_n, tar_x, tar_y, tar_z, tar_s, tar_r, tar_K, P(Phe,:,T), Range);
    cd(path);
    if Phe == 1
        save('mat-cae.mat','charT','N','multy','cpus');%减少重复存储步骤    mat-cae
    end
    %%
    save(name_cae_mat,'sketch_1','sketch_2','sketch_3','sketch_4','sketch_U','sketch_T',...
        'sketch_path','sketch_div','sketch_point','sketch_center','sketch_tarline','sketch_bonding',...
        'tar_r','charPhe','Load_disturb','Load_press','Mesh');...                         mat-T-Phe
    cd(path_m);
toc(time_part)  
    end   %                                                                Phe_CAE --- END
    %%
    %创建模型提交计算，并检查计算结果
    abaqus_run( path, T, N);...                                   python_CAE(_Fold).py
    %%
fprintf('----------第%d代ODB处理开始----------\n',T);
time_part = tic;
    for Phe = 1:N   %                                                      Phe_ODB --- START
        %%     
time_odb = tic;
        %提取odb数据，写入resule_i并输出读取状态
        [ result_x, result_y, result_z] = abaqus_odb( path, Phe, charT);...python_ODB.py
toc(time_odb);
        %%
        %计算相关参数
        [ tar_xr, tar_yr, tar_zr, res_x, res_y, res_z] = ...
            Fun_cal_para_odb...
            ( result_x, result_y, result_z, tar_x, tar_y, tar_z, tar_s);...    Fun_cal_cir_odb_para            
        %%
        %计算误差函数δ
        %尝试仅计算坐标
        del  = cal_delta_C( tar_xr, tar_yr, tar_zr, res_x, res_y, res_z, tar_center);
        %%
        HIS_DEL_P(T,Phe) = del; 
        HIS_DEL_TP(T,:) = A_BEST_TP ;
        %更新个体的极值与最优位置
        if del < P_DEL(Phe)
            P_DEL(Phe)    = del;
            P_BEST(Phe,:) = P(Phe,:,T);
        end
        %更新全局的极值与最优位置，并记录最优个体的信息
        if del < P_DEL_G
            P_DEL_G   = del;
            P_BEST_G  = P_BEST(Phe,:);
            A_BEST_TP = [T,Phe];
        end       
    end    %                                                               Phe_ODB --- END
    %%
    %记录适应度迭代历程
    HIS_DEL_G(T) = P_DEL_G;
    %判断是否达到要求。
    if P_DEL_G < EXPECT          
        fprintf('PSO complete! T=%d.\n',T);
        fprintf('Best DEL=:%g   \n ',P_DEL_G);
        fprintf('T & Phe : JOB-%d-%d \n',A_BEST_TP(1),A_BEST_TP(1));
        %break;
    end
fprintf('----------第%d代ODB处理结束：',T);toc(time_part)
fprintf('进行第%d代更新\n',T);
    %HIS_DEL_P(Phe,:,T) = [del,0];
    %进行遗传操作--参数更新
    if T < T_MAX                                   %未到最后一代则计算下一代参数
        P(:,:,T+1) = Fun_GA( N, D, HIS_DEL_P(T,:), P(:,:,T), AG_Pc, AG_Pm );
        %同时进行边界条件判断：所有的P位于 [0,1]
    end      
    %%
fprintf('-----------------------------第%d代结束-----------------------------\n',T);
fprintf('当前最佳误差 ：:%g    ',P_DEL_G);
fprintf('最佳个体信息 ： JOB-%d-%d \n', A_BEST_TP(1), A_BEST_TP(2));
fprintf('本代');
toc(time_T)
fprintf('总运行时间：');
toc(time_all)
end  ...                                                                   T --- END

%%
%提取最佳个体信息并存储
[ Z_2D_1, Z_2D_2, Z_2D_3, Z_2D_4, ~, ~, ~, Z_2D_C, ~] = ...
        Fun_cal_para_2D( N, D, tar_n, tar_x, tar_y, tar_z, tar_s, tar_r, tar_K, P(A_BEST_TP(2),:,A_BEST_TP(1)), Range);
%读取对应odb
[ result_x, result_y, result_z] = abaqus_odb( path, A_BEST_TP(2), num2str(A_BEST_TP(1)));...python_ODB.py
%计算相关参数
[ Z_3D_TX, Z_3D_TY, Z_3D_TZ, Z_3D_RX, Z_3D_RY, Z_3D_RZ] = ...
    Fun_cal_para_odb...
    ( result_x, result_y, result_z, tar_x, tar_y, tar_z, tar_s);...    Fun_cal_cir_odb_para 
%Z_result_del = HIS_DEL_G( 1:T-1);...未跑完
Z_result_del = HIS_DEL_G( 1:T_MAX);...全部完成
%--------------------------------------------------------------------------保存本次计算的数据
save(name_result_mat);

%%
fprintf('-----------------------------PSO 结束！共计完成%d代！--------------------------\n',T);
toc(time_all)

Result_plot( Z_2D_1, Z_2D_2, Z_2D_3, Z_2D_4 , Z_2D_C,Z_3D_TX, Z_3D_TY, Z_3D_TZ, Z_3D_RX, Z_3D_RY, Z_3D_RZ, Z_result_del, P_DEL_G, A_BEST_TP, T_MAX);
%%
if 1==0
    abaqus_odb_open( path, A_BEST_TP(2), num2str(A_BEST_TP(1)));
end

if 1==0
    fprintf('-----------------------------PSO 结束！共计完成%d代！--------------------------\n',T_MAX);
    toc(time_all)
    fprintf('Best DEL=:%g    ',P_DEL_G);
    fprintf('最佳个体信息 ： JOB-%d-%d \n', A_BEST_TP(1), A_BEST_TP(2));
    %绘制迭代过程曲线
    figure ('color','w')

    plot   (Z_result_del,'r')
    %ylim   ([0,inf])
    xlabel ('迭代次数')
    ylabel ('误差值')
    title  ('适应度进化曲线')
end


