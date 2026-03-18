function sketch = Z_cal_scr( PP )
global tar GA FEA
%数据提取
num_line = tar.line;
bond = tar.bond;
K1 = tar.K1;
tar_n = tar.n;
tar_s = tar.s;
D = GA.D;
TYPE = GA.TYPE;
Cross = tar.cross;
Loc = tar.location;
JUDGE = tar.judge_type_line;  %格式为[ TYPE, LINE]， 分别判断每条带是否启用该优化参数（==1）。
%设置装配调整参数 （为使计算更容易收敛）
P1_AZ = FEA.P1_AZ;
P2_AZ = FEA.P2_AZ ;

BL = 2.5;  %键合区长度

%存储曲面条带的参数
sketch_tarline = cell(num_line,1); 
for k = 1:num_line
    sketch_tarline{k} = [ tar.x(:,k), tar.y(:,k), tar.z(:,k) + tar.r];  %建立目标曲线用
end

%计算2D坐标下的质心线坐标
range = tar.range ;
w_ave = zeros( tar_n, num_line);
t_ave = zeros( tar_n, num_line);

for k = 1:num_line
    w1 = range( 1, 1);
    w2 = range( 1, 2);
    w3 = PP{1}(k,:);
    rw = w1 + w2 * w3;
    tt = linspace( 0, tar_s( end, k), D(1));
    w_ave(:,k) = spline( tt, rw, tar_s(:,k) );
    %修正拟合值
    for i = 1:tar_n
        if w_ave(i,k) < range( 1, 1)
            w_ave(i,k) = range( 1,1);
        end
        if w_ave(i,k) > range( 1, 1) + range( 1, 2)
            w_ave(i,k) = range( 1, 1) + range( 1, 2) + 0.1 ;
        end
    end
    if TYPE > 1 && JUDGE( 2, k) == 1
        t1 = range( 2, 1);
        t2 = range( 2, 2);
        t3 = PP{2}(k,:);
        rt = t1 + t2 * t3;
        tt = linspace( 0, tar_s( end, k), D(2));
        t_ave(:,k) = spline( tt, rt, tar_s(:,k) );
        for i = 1:tar_n
            if t_ave(i,k) < range( 2, 1)
                t_ave(i,k) = range( 2, 1);
            end
            if t_ave(i,k) > range( 2, 1) + range( 2, 2)
                t_ave(i,k) = range( 2, 1) + range( 2, 2);
            end
        end
    end
end

%%
%生成纯直线的二维质心线坐标 , Y坐标均为0
%旋转角度，令其相对位置与三维下保持一致
%  [1] :上下轮廓线 
sketch_center = cell(num_line,1);    %Part时使用，用于创建质心线切分   zeros( tar_n, 2, num_line + 1)
sketch_1 = cell(num_line,1);  %上
sketch_2 = cell(num_line,1);  %下

%计算各条带所需要的  旋转角度
theta_rotation = zeros( num_line ,1); %旋转角度
theta_loc      = zeros( num_line ,4); %旋转角度计算的插值点位置，仅对TYPE>1启用

theta_rotation(Loc(1)) = Loc(6);  %第一条带的角度为定位角度
theta_loc(Loc(1),:) = [ Loc(3), Loc(4),  Loc(3), Loc(4)] ;   
line_list = [ Loc(1), Loc(1); Cross(:,1:2) ];
%如果为多条带，执行Cross计算(旋转)
if num_line > 1
    for c = 1:size(Cross,1)
        c2 = Cross(c,2); %Line2号
        t1 = Cross(c,7); %线1角度
        t2 = Cross(c,8); %线2角度
        theta_rotation(c2) = t2 - t1;    %存角度(夹角)
        theta_loc(c2,1:2) = Fun_Point_judge( Cross(c,5), tar_n); %存插值点位置：定位线
        theta_loc(c2,3:4) = Fun_Point_judge( Cross(c,6), tar_n); %存插值点位置：待定位线
    end
end

%如果为单条带且带有随机角度，进行随机角度修正【随机角度仅在单条带下启用】
if num_line ==1 && TYPE == 3
    theta_rand = range(3,1) + PP{3}(1) * range(3,2) ;
    theta_rotation(1) = theta_rotation(1) + theta_rand ;
end
%计算上下轮廓线  （TYPE<1时进行旋转）
if TYPE == 1
%############################################
    theta_right = zeros(num_line,1);
    theta_left  = zeros(num_line,1);
    for k = 1:num_line
        [ sketch_center{k}, sketch_1{k},  sketch_2{k}, theta_right(k), theta_left(k)] = Fun_TYPE_1( w_ave(:,k), tar_s(:,k), tar_n, theta_rotation(k)) ;
    end
%############################################
else
    Fun_TYPE_2;  %带有截面扭转角度参数的情况
end
%############################################

%定位（平移） - Location  
Loc_tar = sketch_center{Loc(1)}(Loc(5), :);
Loc_Ori = [ tar.x(Loc(5),Loc(1)), tar.y(Loc(5),Loc(1))];
LocXY = Loc_tar - Loc_Ori ;
sketch_center{Loc(1)} = sketch_center{Loc(1)} - LocXY ;
sketch_1{Loc(1)} = sketch_1{Loc(1)} - LocXY ;
sketch_2{Loc(1)} = sketch_2{Loc(1)} - LocXY ;

%定位（平移） - Cross
if num_line > 1 
    for c = 1:size(Cross,1)
        cc1 = Cross(c,1); %线号
        cc2 = Cross(c,2); %线号
        tt1 = Cross(c,5); %位置
        tt2 = Cross(c,6); %位置
        tarLoc = sketch_center{cc1}( tt1, :);
        oriLoc = sketch_center{cc2}( tt2, :);
        del = tarLoc - oriLoc;
        sketch_center{cc2} = sketch_center{cc2} + del;
        sketch_1{cc2}( :, :) = sketch_1{cc2}( :, :) + del;
        sketch_2{cc2}( :, :) = sketch_2{cc2}( :, :) + del;
        %根据cross调整center点
        sketch_center{cc2}( tt2, :) = sketch_center{cc1}( tt1, :);
    end  

end

%%
sketch_3 = cell( num_line, 1);
sketch_4 = cell( num_line, 1);
sketch_div  = cell(num_line,1);
sketch_ring = cell(num_line,1);
sketch_find_b = cell(num_line,1);
sketch_U    = cell(num_line,2);   %cell 的 行为条带，列为分析步 ； 内部两行分别为左右Bonding点的数据
sketch_path = cell(num_line,1);         %ODB处理时使用

sketch_find_f = zeros( 4, 2, num_line + 1);

for k = 1:num_line
%  [2] :   Bonding点信息
%方式：以原点为中心建立相对关系，先旋转至与边缘处相切（对应），再平移至对应位置
    Bonding_area  = [  0,  BL/2;
                     -BL,  BL/2;
                     -BL, -BL/2; 
                       0, -BL/2]; 
    left_location  = [sketch_center{k}(1,1)  , sketch_center{k}(1,2)  ];
    right_location = [sketch_center{k}(end,1), sketch_center{k}(end,2)];
    LEFT  = theta_left(k);
    RIGHT = theta_right(k);
    matrix_left  = [ cos( LEFT ), sin( LEFT ); -sin( LEFT ), cos( LEFT )];
    matrix_right = -[ cos( RIGHT ), sin( RIGHT ); -sin( RIGHT ), cos( RIGHT )];
    sketch_3{k} = Bonding_area * matrix_left  + left_location;
    sketch_4{k} = Bonding_area * matrix_right + right_location;
   
    %FindAt的点数据，取Bonding面的两部分
    sketch_find_b{k}( 1, :) = [ sketch_3{k}(1,1) + sketch_3{k}(3,1)  ,sketch_3{k}(1,2) + sketch_3{k}(3,2) ]/2 ;
    sketch_find_b{k}( 2, :) = [ sketch_4{k}(1,1) + sketch_4{k}(3,1)  ,sketch_4{k}(1,2) + sketch_4{k}(3,2) ]/2 ;
    
    %键合点数量判定,若不存在键合点则改为直接封口
    if bond(k) == -1
        sketch_4{k} = zeros(2,2);
        sketch_find_b{k}(2,:) = [];
    elseif bond(k) == 1
        sketch_3{k} = zeros(2,2);
        sketch_find_b{k}(1,:) = [];
    elseif bond(k) == 0
        sketch_3{k} = zeros(2,2);
        sketch_4{k} = zeros(2,2);
        sketch_find_b{k} = [];
    end
    

    
%  [4] : 用于定位Ring的数据
    sketch_ring{k} = [ sketch_center{k}( 1, 1);  sketch_center{k}(end, 1) ];
                   
%  [5] : 位移边界条件信息（左Bonding点始终为固定）
    theta_t1 = pi/2 - atan2( tar.z(   1, k) + tar.r , tar.y(   1, k) );
    theta_t2 = pi/2 - atan2( tar.z( end, k) + tar.r , tar.y( end, k) );
    theta_b1 = sketch_center{k}(   1, 2) / tar.r;
    theta_b2 = sketch_center{k}( end, 2) / tar.r;
    UX_L  = tar.x(   1, k) - sketch_center{k}(   1, 1);  %左位移X
    UX_R  = tar.x( end, k) - sketch_center{k}( end, 1);  %右位移X
    UXR_L = -( theta_t1 - theta_b1 );   %左角度XR
    UXR_R = -( theta_t2 - theta_b2 );   %右角度XR
    sketch_U{k,1} = [ UX_L, 0, 0, UXR_L/2, 0, 0; 
                      UX_R, 0, 0, UXR_R/2, 0, 0];
    sketch_U{k,2} = [ UX_L, 0, 0, UXR_L, 0, 0; 
                      UX_R, 0, 0, UXR_R, 0, 0];
    
%  [6] : 用于findAt 的点数据
    %取FREE的两部分
    sketch_find_f(1,:,k) = [ sketch_1{k}(1,1) + sketch_center{k}(2,1) , sketch_1{k}(1,2) + sketch_center{k}(2,2)]/2 ; 
    sketch_find_f(2,:,k) = [ sketch_2{k}(1,1) + sketch_center{k}(2,1) , sketch_2{k}(1,2) + sketch_center{k}(2,2)]/2 ;
    sketch_find_f(3,:,k) = [ sketch_1{k}(end,1) + sketch_center{k}(end-1,1) , sketch_1{k}(end,2) + sketch_center{k}(end-1,2)]/2 ;
    sketch_find_f(4,:,k) = [ sketch_2{k}(end,1) + sketch_center{k}(end-1,1) , sketch_2{k}(end,2) + sketch_center{k}(end-1,2)]/2 ;

    %  [7] : ODB提取路径点信息       
    sketch_path{k} = [ sketch_center{k}, ones(tar_n,1)*P1_AZ];
    
    %Bonding 数量处理：
    if bond(k) == -1
        sketch_ring{k}(2) = [];
        sketch_div{k}(3:4,:) = [];
        sketch_U{k,1}(2,:) = [];
        sketch_U{k,2}(2,:) = [];
    elseif bond(k) == 1
        sketch_ring{k}(1) = [];
        sketch_div{k}(1:2,:) = [];
        sketch_U{k,1}(1,:) = [];
        sketch_U{k,2}(1,:) = [];
    elseif bond(k) == 0
        sketch_ring{k} = [];
        sketch_div{k} = [];
        sketch_U{k,1} = [];
        sketch_U{k,2} = [];
    end    
    
end
sketch.w = w_ave;
sketch.tarline = sketch_tarline;
sketch.center = sketch_center;
sketch.line1 = sketch_1;
sketch.line2 = sketch_2;
sketch.line3 = sketch_3;
sketch.line4 = sketch_4;
sketch.path  = sketch_path;
sketch.div   = sketch_div ;
sketch.U     = sketch_U;
sketch.find_free = sketch_find_f;
sketch.find_bond = sketch_find_b;
sketch.ring = sketch_ring;
sketch.num_line = num_line;

sketch.Assemble_P1Z = P1_AZ;
sketch.Assemble_P2Z = P2_AZ ;
      

%%
%局部函数  -- 函数内
function Fun_TYPE_2
%处理带有截面扭转角度时的草图  
    for kkline = 1:num_line
        kk = line_list( kkline, 2) ;
        mas = line_list( kkline, 1) ;
        if JUDGE(2,kk) == 1  %如果启用角度参数
            K = sin(t_ave(:,kk) ) .*K1(:,kk);
            K1s   = spline( tar_s(:,kk), K);
            theta = - cal_InteSpline( tar_n, K1s);%Theta_2D=K1(S)对S求积分
            %根据转角方程Theta(S)求平面曲线方程r(S)
            theta_cos = cos(theta);
            theta_sin = sin(theta);
            theta_cos_s = spline( tar_s(:,kk), theta_cos);
            theta_sin_s = spline( tar_s(:,kk), theta_sin);
            %对θ-S正余弦插值
            %由转角-S求2D坐标
            X_2D = cal_InteSpline( tar_n, theta_cos_s);
            Y_2D = cal_InteSpline( tar_n, theta_sin_s);
            sketch_center{kk} = [ X_2D, Y_2D];
            %令2D构型的两个Bonding点均位于X轴上
            theta_bond1  = -atan2( Y_2D(end) , X_2D(end));
                matrix_rotation = [ cos(theta_bond1), sin(theta_bond1); -sin(theta_bond1), cos(theta_bond1)];
            sketch_center{kk} = sketch_center{kk} *matrix_rotation;
            Theta_2D_r = theta + theta_bond1;   
            ew = [ -sin(Theta_2D_r), cos(Theta_2D_r)];    %诱导公式，由原角度 +90°转化得到
            sketch_1{kk} = sketch_center{kk} + (w_ave(:,kk)).*ew; %上
            sketch_2{kk} = sketch_center{kk} - (w_ave(:,kk)).*ew; %下 
        else
            %[ sketch_center{kk}, sketch_1{kk},  sketch_2{kk}, theta_right(kk), theta_left(kk)] = Fun_TYPE_1( w_ave(:,kk), tar_s(:,kk), tar_n, theta_rotation(kk));
            X_2D = linspace( 0, tar_s(end,kk), tar_n)';
            sketch_1{kk} = [ X_2D , w_ave(:,kk)];
            sketch_2{kk} = [ X_2D , -w_ave(:,kk)]; 
            sketch_center{kk} = ( sketch_1{kk} + sketch_2{kk} )/2;
        end
            
            %定位 - Location
            %由于存在曲率，需要首先计算角度计算点的现有角度
            line_tar = sketch_center{kk};
            line_mas = sketch_center{mas};
            locmas1 = theta_loc(kk,2);
            locmas2 = theta_loc(kk,1);
            loctar1 = theta_loc(kk,4);
            loctar2 = theta_loc(kk,3);
            theta_now_mas = atan2( line_mas(locmas1,2) - line_mas(locmas2,2), line_mas(locmas1,1) - line_mas(locmas2,1) );
            theta_now_tar = atan2( line_tar(loctar1,2) - line_tar(loctar2,2), line_tar(loctar1,1) - line_tar(loctar2,1) );
            
            if kk ~= line_list(1,1) %若不为Location对应的线段，则需要将目标旋转角度(夹角)转化为方位角
                theta_rotation(kk) = theta_now_mas + theta_rotation(kk);
            end
            thetaR = theta_rotation(kk) - theta_now_tar;

                matrix_rotation = [ cos(thetaR), sin(thetaR); -sin(thetaR), cos(thetaR)];
            sketch_center{kk} = sketch_center{kk} * matrix_rotation;
            sketch_1{kk} = sketch_1{kk} * matrix_rotation;
            sketch_2{kk} = sketch_2{kk} * matrix_rotation;
            theta_left(kk) = atan2( sketch_center{kk}( 2, 2) - sketch_center{kk}( 1, 2) , sketch_center{kk}( 2, 1) - sketch_center{kk}( 1, 1) ) ;
            theta_right(kk) = atan2( sketch_center{kk}( end, 2) - sketch_center{kk}( end-1, 2) , sketch_center{kk}( end, 1) - sketch_center{kk}( end-1, 1) );
            

        
    end

end

%%


end%主函数结束

%%
%局部函数 - 函数外
function [ sketch_center, sketch_1,  sketch_2, theta_right, theta_left] = Fun_TYPE_1( w_ave, tar_s, tar_n, theta_rotation)
%数据提取
        X_2D = linspace( 0, tar_s(end), tar_n)';
        XY_1 = [ X_2D , w_ave];
        XY_2 = [ X_2D , -w_ave];
        sketch_center = [ X_2D*cos(theta_rotation ) , X_2D*sin(theta_rotation) ];
            matrix_rotation = [ cos(theta_rotation), sin(theta_rotation); -sin(theta_rotation), cos(theta_rotation)];
        sketch_1 = XY_1 * matrix_rotation;
        sketch_2 = XY_2 * matrix_rotation;
        theta_right = theta_rotation;
        theta_left  = theta_rotation;
end 


    function Inter=cal_InteSpline(N,pp)
    %Theta_2D=K1(S)对S求积分
        Inter=zeros(N,1);
        di=zeros(N-1,1);
        I_K=zeros(N-1,1);
        a=pp.coefs(:,1);
        b=pp.coefs(:,2);
        c=pp.coefs(:,3);
        d=pp.coefs(:,4);
        for i=1:N-1
            di(i)=pp.breaks(i+1)-pp.breaks(i);
            I_K(i)=0.25*a(i)*di(i)^4+(1/3)*b(i)*di(i)^3+0.5*c(i)*di(i)^2+d(i)*di(i);
            Inter(i+1)=Inter(i)+I_K(i);%假定Inter(1)=0,求i=2~N_inter-1的spline函数的积分
        end
    end
    
    function theta_loc = Fun_Point_judge( cross, tar_n ) 
        %判断交点位置
        % [ 靠后的点， 靠前的点]
if cross == 1
    tt1 = 1;
    tt2 = 2;
elseif cross == tar_n
    tt1 = tar_n - 1 ;
    tt2 = tar_n ;
else
    tt1 = cross - 1;
    tt2 = cross + 1;
end
theta_loc = [ tt1, tt2];
        
    end


