function Fun_cal_para_tar_point_set
%由t参数方程直接给定质心线方程，仅给出质心线起止点。
global tar
syms t
%数据导入：每条质心线为一个元胞，每个元胞为 ponit*3 矩阵，可以考虑Load 或者 readmatrix
%数据存储：每 列 为一条质心线  【方便查看】
%数据提取
num_line = tar.line;
linelist = tar.linelist;
location = tar.location;
tar_n = tar.n;
%%
xx = zeros( tar_n, num_line);
yy = zeros( tar_n, num_line);
zz = zeros( tar_n, num_line);
ss = zeros( tar_n, num_line);
s_old = cell(num_line,1);
K1 = zeros( tar_n, num_line);
for k = 1:num_line
    num_point = size(linelist{k},1);
    s = zeros(num_point, 1);
    cr = linelist{k};
    %计算初始弧长，直接计算直线距离
    for i = 1:num_point - 1
        ds = 0;
        for j = 1:3
            ds = ds + ( cr(i+1,j) - cr(i,j) )^2;%XYZ累加
        end
        s(i+1) = s(i) + (ds)^0.5;
    end
    xs = spline( s, cr(:,1));
    ys = spline( s, cr(:,2));
    zs = spline( s, cr(:,3));
    s = cal_s_integral( xs, ys, zs, s);
    s_old{k} = s;
    ss(:,k) = linspace( 0, s(end), tar_n);
    xx(:,k) = spline( s, cr(:,1), ss(:,k) );
    yy(:,k) = spline( s, cr(:,2), ss(:,k) );
    zz(:,k) = spline( s, cr(:,3), ss(:,k) );
    %使用S对XYZ进行三次样条拟合
    xs = spline( ss( :, k), xx( :, k) );
    ys = spline( ss( :, k), yy( :, k) );
    zs = spline( ss( :, k), zz( :, k) );
    %计算质心线曲线曲率，K即EML中的K1*
    kx_2 = cal_k2( xs, tar_n);
    ky_2 = cal_k2( ys, tar_n);
    kz_2 = cal_k2( zs, tar_n);%计算分段三次拟合曲线在各个插值点处XYZ对S的二阶导数
    K1(:,k) = ( kx_2.^2 +ky_2.^2 +kz_2.^2).^0.5;
end
tar.x  = xx;
tar.y  = yy;
tar.z  = zz;
tar.s  = ss;
tar.K1 = K1;
%计算结构整体的中心位置，用于误差计算
tar.center = [ sum(tar.x,'all'); sum(tar.z,'all'); sum(tar.y,'all')] / tar.n / num_line;

%计算定位点的参数,该点Y坐标为0
cc1 = location(1) ;
cc2 = location(2) ;
cc3 = location(3) ;
cc4 = location(4) ;
[~,location(5)] = min( abs( ss(:,k) - s_old{k}(cc2) ) );
[ ~, cc3] =  min( abs( ss(:,k) - s_old{k}(cc3) ) );
[ ~, cc4] =  min( abs( ss(:,k) - s_old{k}(cc4) ) );
if cc3 == cc4
    if cc3 == 1
        tt1 = 2;
        tt2 = 1;
    elseif cc3 == tar_n
        tt1 = tar_n  ;
        tt2 = tar_n - 1 ;
    else
        tt1 = cc4+1;
        tt2 = cc4-1;
    end
else
    tt1 = cc4; %后
    tt2 = cc3; %前
end
location(6) = atan2( tar.y( tt1,cc1) - tar.y( tt2,cc1) , tar.x( tt1,cc1) - tar.x( tt2,cc1) );
location(3) = tt2;
location(4) = tt1;
tar.location = location;
%%
if num_line == 1
    return
end
%计算交叉点的参数情况
num_cross = size(tar.cross,1);

for c = 1:num_cross
    c1 = tar.cross(c,1);
    c2 = tar.cross(c,2);
    t1 = tar.cross(c,3);
    t2 = tar.cross(c,4);
    [~,n1] = min( abs( ss(:,c1) - s_old{c1}(t1) ) );
    [~,n2] = min( abs( ss(:,c2) - s_old{c2}(t2) ) );
    tar.cross(c,5) = n1;
    tar.cross(c,6) = n2;
    if n2 == 1
        tt1 = 2;
        tt2 = 1;
    elseif n2 == tar_n
        tt1 = tar_n  ;
        tt2 = tar_n - 1 ;
    else
        tt1 = n2+1;
        tt2 = n2-1;
    end
    tar.cross(c,7) = atan2( tar.y( tt1,c2) - tar.y(tt2,c2) , tar.x( tt1,c2) - tar.x( tt2,c2) );    %   转角以XOY平面投影的转角为准。这种方式仅适用于交叉点位于 Y = 0的平面内！
end


end %主函数

%%
%局部函数
function s_tar = cal_s_integral( xt, yt, zt, ss)
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
            breaks = tt - ss(i);
            dx = 3*cx(i,1)*breaks^2 + 2*cx(i,2)*breaks + cx(i,3);
            dy = 3*cy(i,1)*breaks^2 + 2*cy(i,2)*breaks + cy(i,3);
            dz = 3*cz(i,1)*breaks^2 + 2*cz(i,2)*breaks + cz(i,3);
            ds(tt) = ( dx^2 + dy^2 + dz^2 )^0.5;
            %列出ds表达式，并由matlabFunction转为函数句柄
            DS     = matlabFunction(ds);
            fun_ds = integral( DS, ss(i), ss(i+1) );
            s_tar(i+1) = s_tar(i) + fun_ds;    
        end
end
%局部函数

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
%局部函数
