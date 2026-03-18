function Fun_cal_para_tar_eqation
%由t参数方程直接给定质心线方程，仅给出质心线起止点。
global tar
syms t
%数据导入：每 行 为一条质心线  【方便书写】
%数据存储：每 列 为一条质心线  【方便查看】
%数据提取
num_line = tar.line;
location = tar.location;
%针对不同曲线分别插值
tt = zeros( tar.n, num_line);
for i = 1:num_line
    tt( :, i) = linspace( tar.t_lim( i, 1), tar.t_lim( i, end), tar.n)';
end

%%
%确定各插值点所选用的 分段函数
shape = ones( tar.n, num_line);
for a = 1:num_line
    judge = 1;
    for b = 1:tar.n
        if tt( b, a) > tar.t_lim( a,judge + 1)
           judge = judge + 1;
        end
        shape( b, a) = judge;
    end
end

xx = zeros( tar.n, num_line);
yy = zeros( tar.n, num_line);
zz = zeros( tar.n, num_line);
ss = zeros( tar.n, num_line);

for k = 1:num_line
    for i = 1:tar.n
        fx(t) = tar.funx( k, shape(i, k) );
        fy(t) = tar.funy( k, shape(i, k) );
        fz(t) = tar.funz( k, shape(i, k) );
        %求坐标
        xx( i, k) = double ( fx( tt( i ,k) ) );
        yy( i, k) = double ( fy( tt( i ,k) ) );
        zz( i, k) = double ( fz( tt( i ,k) ) );   
        %求弧长
        fx1(t) = diff( tar.funx( k, shape(i, k) ) );
        fy1(t) = diff( tar.funy( k, shape(i, k) ) );
        fz1(t) = diff( tar.funz( k, shape(i, k) ) );
        ds = (fx1^2+fy1^2+fz1^2)^0.5;
        DS = matlabFunction(ds);
        if i>1
            dss = integral( DS, tt(i-1), tt(i) ,'ArrayValued',true);
            ss(i,k) = ss( i-1, k) + dss;
        end    
    end
end
%%

tar.x  = zeros( tar.n, num_line);
tar.y  = zeros( tar.n, num_line);
tar.z  = zeros( tar.n, num_line);
tar.s  = zeros( tar.n, num_line);
tar.K1 = zeros( tar.n, num_line);
%均匀化，便于划分网格
for k = 1:num_line
    tar.s( :, k) = linspace( 0, ss( end, k), tar.n)';%s_tar即2D下的 x_2D
    tar.x( :, k) = spline( ss( :, k), xx( :, k), tar.s( :, k));
    tar.y( :, k) = spline( ss( :, k), yy( :, k), tar.s( :, k));
    tar.z( :, k) = spline( ss( :, k), zz( :, k), tar.s( :, k));
    %使用S对XYZ进行三次样条拟合
    xs = spline( tar.s( :, k), tar.x( :, k) );
    ys = spline( tar.s( :, k), tar.y( :, k) );
    zs = spline( tar.s( :, k), tar.z( :, k) );
    %计算质心线曲线曲率，K即EML中的K1*
    kx_2 = cal_k2( xs, tar.n);
    ky_2 = cal_k2( ys, tar.n);
    kz_2 = cal_k2( zs, tar.n);%计算分段三次拟合曲线在各个插值点处XYZ对S的二阶导数
    tar.K1(:,k) = ( kx_2.^2 +ky_2.^2 +kz_2.^2).^0.5;
end

%计算结构整体的中心位置，用于误差计算
tar.center = [ sum(tar.x,'all'); sum(tar.z,'all'); sum(tar.y,'all')] / tar.n / num_line;

%计算定位点的参数,该点Y坐标为0
cc1 = location(1) ;
cc2 = location(2) ;
cc3 = location(3) ;
cc4 = location(4) ;
[~,location(5)] = min( abs( tt( :, cc1) - cc2) );
[ ~, cc3] =  min( abs( tt( :, cc1) - cc3) );
[ ~, cc4] =  min( abs( tt( :, cc1) - cc4) );
%  tt1为靠后点， tt2为靠前点
if cc3 == cc4
    if cc3 == 1
        tt1 = 2;
        tt2 = 1;
    elseif cc3 == tar.n
        tt1 = tar.n  ;
        tt2 = tar.n - 1 ;
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
%  对外轮廓线的特殊调整                                                         ##############
if isfield(tar,'inflection') == 1 %对称判断
    inf = tar.inflection;
    num_inf = size(inf,1);
    for kl = 1:num_inf
        line = inf(kl,1);
        kt = inf(kl,2);
        [~,tt] = min( abs( tt( :, line) - kt) );  % 将参数坐标替换为序号
        tar.inflection(kl,2) = tt;
    end
end


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
    [~,n1] = min( abs( tt( :, c1) - t1) );
    [~,n2] = min( abs( tt( :, c2) - t2) );
    tar.cross(c,5) = n1;  %角度所在插值点序号
    tar.cross(c,6) = n2;  %角度所在插值点序号
    %  ***1为靠后点， ***2为靠前点
    [ ct11, ct12] = judge_order(n1);
    [ ct21, ct22] = judge_order(n2);
    %首先旋转至顶端，正投影
    theta1 = pi/2 - atan2( tar.z(n1,c1), tar.y( n1,c1) );
    theta2 = pi/2 - atan2( tar.z(n2,c2), tar.y( n2,c2) );
    rotation1 = [ cos(theta1), -sin(theta1); sin(theta1), cos(theta1) ];
    rotation2 = [ cos(theta2), -sin(theta2); sin(theta2), cos(theta2) ];
    vector1 = [ tar.y( ct11,c1) - tar.y(ct12,c1) , tar.x( ct11,c1) - tar.x( ct12,c1)] * rotation1;
    vector2 = [ tar.y( ct21,c2) - tar.y(ct22,c2) , tar.x( ct21,c2) - tar.x( ct22,c2)] * rotation2;
    tar.cross(c,7) = atan2( vector1(1), vector1(2) );
    tar.cross(c,8) = atan2( vector2(1), vector2(2) );
end



end

%%
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

function [ order1, order2] = judge_order(num)
%判断角度计算点的位置，得到该点处相邻两点的序列
% order1 --> 靠后点
% order2 --> 靠前点
global tar
    if num == 1
        order1 = 2;
        order2 = 1;
    elseif num == tar.n
        order1 = tar.n  ;
        order2 = tar.n - 1 ;
    else
        order1 = num + 1;
        order2 = num - 1;
    end
    
end
