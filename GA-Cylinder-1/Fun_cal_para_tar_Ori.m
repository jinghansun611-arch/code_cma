function [tar_x, tar_y, tar_z, tar_s, center, tar_K] = ...
    Fun_cal_para_tar( tar_n , funx, funy, funz, t, t_lim)
%若不是直接给出质心线插值坐标，由该函数给出
%测试1：由t参数方程直接给定质心线方程，仅给出质心线起止点。

%计算插值参数
t1 = double(t_lim(1) );
t2 = double(t_lim(2) );
tt = linspace( t1, t2, tar_n)';
%由参数方程计算插值点
xx = zeros(tar_n, 1);
yy = zeros(tar_n, 1);
zz = zeros(tar_n, 1);
for i = 1 : tar_n
    xx(i) = double ( funx( tt(i) ) );
    yy(i) = double ( funy( tt(i) ) );
    zz(i) = double ( funz( tt(i) ) );
end
ave_t = linspace(1,tar_n,tar_n)';
xt = spline(ave_t, xx);
yt = spline(ave_t, yy);
zt = spline(ave_t, zz);
ss = cal_s_integral( xt, yt, zt);%数值积分，近似,对于该问题足够精确且迅速
%均匀化，便于划分网格
tar_s = linspace( 0, ss( tar_n), tar_n)';%s_tar即2D下的 x_2D
tar_x = spline( ss, xx, tar_s);
tar_y = spline( ss, yy, tar_s);
tar_z = spline( ss, zz, tar_s);
%令左Bonding点位于坐标原点
tar_x = tar_x - tar_x(1);
tar_y = tar_y - tar_y(1);

center = [ sum(tar_x,'all'); sum(tar_z,'all'); sum(tar_y,'all')] / tar_n;

%%
%---------------------------------------计算2D参数(部分)
%使用S对XYZ进行三次样条拟合
xs = spline( tar_s, tar_x);
ys = spline( tar_s, tar_y);
zs = spline( tar_s, tar_z);
%计算质心线曲线曲率，K即EML中的K1*
kx_2 = cal_k2( xs, tar_n);
ky_2 = cal_k2( ys, tar_n);
kz_2 = cal_k2( zs, tar_n);%计算分段三次拟合曲线在各个插值点处XYZ对S的二阶导数
tar_K = ( kx_2.^2 +ky_2.^2 +kz_2.^2).^0.5;

end
%%
%局部函数
    function s_tar = cal_s_integral( xt, yt, zt)
        %计算各插值点的弧长S
        %数值积分求各段弧长
        %[breaks,coefs,L,order,dim] = unmkpp(pp)
        [ ~, cx, L, ~, ~] = unmkpp(xt);
        [ ~, cy, ~, ~, ~] = unmkpp(yt);
        [ ~, cz, ~, ~, ~] = unmkpp(zt);
        syms tt ;
        ds(tt) = tt;%只是初始化预分配内存，等号后没有含义
        s_tar  = zeros(1,L+1);
        for i = 1:L
            breaks = tt - i;
            dx = 3*cx(i,1)*breaks^2 + 2*cx(i,2)*breaks + cx(i,3);
            dy = 3*cy(i,1)*breaks^2 + 2*cy(i,2)*breaks + cy(i,3);
            dz = 3*cz(i,1)*breaks^2 + 2*cz(i,2)*breaks + cz(i,3);
            ds(tt) = ( dx^2 + dy^2 + dz^2 )^0.5;
            %列出ds表达式，并由matlabFunction转为函数句柄
            DS     = matlabFunction(ds);
            fun_ds = integral(DS,i,i+1);
            s_tar(i+1) = s_tar(i) + fun_ds;    
        end
    end

    function k = cal_k2(pp,N_inter)
    %求第i个插值点处的r对s的二阶导
    %coefs(N_inter,4)的四列分别为3210阶的系数，k=6a(t-t1)+2b
        k=zeros(N_inter,1);
        bi=pp.coefs(:,2);
        for i=1:N_inter-1
            k(i)=2*bi(i);
        end
        d=pp.breaks(N_inter)-pp.breaks(N_inter-1);
        a=pp.coefs(N_inter-1,1);
        b=pp.coefs(N_inter-1,2);
        k(N_inter)=6*a*d+2*b;
    end