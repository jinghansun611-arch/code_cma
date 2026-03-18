
%%
for T = 1 : GA.T_MAX     %                                                      T --- START
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
            del = 1;
        end

        %%
        HIS_D_P(T,Phe) = del; 
        %更新个体的极值与最优位置
        if del < BEST_D_P(Phe)
            BEST_D_P(Phe)    = del;
            BEST_P_P(Phe,:) = P(Phe,:,T);
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
    if BEST_D_G < GA.EXPECT          
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
%提取最佳个体信息并存储
%tar.judge_type_line = ones( GA.TYPE, tar.line);     %补充
save_sketch = Fun_cal_para_2D_multy_w( P(BEST_NUM(2),:,BEST_NUM(1)) ) ;
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
save(name_data);
fprintf('-----------------------------PSO 结束！共计完成%d代！--------------------------\n',T_NOW);
toc(time_all)
%%
Result_plot( save_sketch, save_result, save_his_del, BEST_D_G, BEST_NUM);
%%
if 1==1
    abaqus_odb_open( path_abaqus, BEST_NUM(1), BEST_NUM(2) );
end

if 1==0
    fprintf('-----------------------------PSO 结束！共计完成%d代！--------------------------\n',GA.T_NOW);
    toc(time_all)
    fprintf('Best DEL=:%g    ',BEST_D_G);
    fprintf('最佳个体信息 ： JOB-%d-%d \n', BEST_NUM(1), BEST_NUM(2));
    %绘制迭代过程曲线
    figure ('color','w')

    plot   (save_his_del,'r')
    %ylim   ([0,inf])
    xlabel ('迭代次数')
    ylabel ('误差值')
    title  ('适应度进化曲线')
end


