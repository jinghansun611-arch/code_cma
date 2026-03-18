function result = Fun_cal_para_odb_m( result)
global tar

tar_x = tar.x;
tar_y = tar.y;
tar_z = tar.z;
tar_s = tar.s;
tar_r = tar.r;
num_line = tar.line;

result_x = result.rx;
result_y = result.ry;
result_z = result.rz;


tar_data = cell(num_line,3);
res_data = cell(num_line,3);
for k = 1:num_line
    %以路径提取点的弧长坐标为准
    res_s = result_x{k}(:,1);
    
    %计算目标曲面在各点的数值
    tar_data{k,1} = spline( tar_s(:,k), tar_x(:,k), res_s);
    tar_data{k,2} = spline( tar_s(:,k), tar_y(:,k), res_s);
    tar_data{k,3} = spline( tar_s(:,k), tar_z(:,k), res_s);
    %%
    %-----------------------------------------------------------------------直接提取坐标
    res_data{k,1} = result_x{k}(:,2);
    res_data{k,2} = result_y{k}(:,2);
    res_data{k,3} = result_z{k}(:,2) - tar_r;
end

result.resdata = res_data;
result.tardata = tar_data;




end