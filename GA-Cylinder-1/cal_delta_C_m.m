function del = cal_delta_C_m( result)
%计算误差函数δC
%数据提取
global tar
center = tar.center;
num_line = tar.line;
%质心为center
cx = center(1);
cy = center(2);
cz = center(3);
res_data = result.resdata;
tar_data = result.tardata;


up   = 0;  %分子
down = 0;  %分母
for k = 1 : num_line
    for i = 1 : length(res_data{k,1})
        xr = res_data{k,1};
        yr = res_data{k,2};
        zr = res_data{k,3};
        xt = tar_data{k,1};
        yt = tar_data{k,2};
        zt = tar_data{k,3};
        
        up   = up   + ( ( xr(i) - xt(i) )^2 + ( yr(i) - yt(i) )^2 + ( zr(i) - zt(i) )^2 );
        down = down + ( ( xt(i) - cx)^2     + ( yt(i) - cy)^2     + ( zt(i) - cz)^2  );
    end
end
del   = up^0.5 / down^0.5;

end

