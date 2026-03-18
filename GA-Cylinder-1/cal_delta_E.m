function del = cal_delta_E(et,er,N)
%计算误差函数δC
up   = 0;   %分子
down = 0;   %分母
for i = 1 : N
    for j=1:3
            up   = up + ( et(i,j)- er(i,j))^2;
            down = down + et(i,j)^2;
    end
end
del = up^0.5 / down^0.5 ;
end