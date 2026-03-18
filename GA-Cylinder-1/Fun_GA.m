function NEWP = Fun_GA( DEL, PP)
%(NP,G_XY_2,DEL_value,G_Pc,G_Pm,G_L,R,Judge)
%执行遗传操作：复制、交叉、变异；返回十进制子代
%由于各个优化参数的二进制长度不一致，在复制操作仍然以十进制进行
%之后再通过循环依次转为二进制，进行交叉和编译
    
%数据提取
global tar GA
line= tar.line;
N = GA.N;
D = GA.D;
Pc = GA.Pc;
Pm = GA.Pm;
type = GA.TYPE;
BIT = GA.BIT;
%DEL : [ N, 1 ]

NEWP = PP;
% 构成：  （个体数，变量类型数） = cell（ N, type）
% cell构成： （条代数，对应优化变量数） = （ line， D(line) ）
%--------------------------------------------------------------------------基于轮盘赌的复制操作
%提取极值
    DEL_MIN = min( DEL );  %返回靠前的最小值
    DEL_MAX = max( DEL );           %返回最大值
%归一化适应度值
    if DEL_MIN == DEL_MAX
        DEL_A = ( DEL*0 + 1) /(N-1);   %预防所有DEL相等的情况导致分母为0
    else
        DEL_A   = ( DEL_MAX - DEL ) / ( DEL_MAX - DEL_MIN );
    end
%累计概率
    DEL_value1 = DEL_A ./ (sum(DEL_A));
    DEL_value  = cumsum(DEL_value1);     %求累加和
%将新一代的轮盘赌从小到大排列(复制N个个体)
    DEL_ms = sort( rand( N,1) );
    oldp = 1;
    newp = 1;
    while newp <= N
        if ( DEL_ms(newp) ) < DEL_value(oldp)
            NEWP( newp, :) = PP( oldp, :);
            newp = newp+1;
        else
            oldp = oldp+1;
        end
    end
    %%  %%  %%
    for t = 1: type
        unit1 = zeros( N, D(t) ); %预设存储每条带数据的数组
        unit2 = zeros( N, D(t) ); %预设存储更新后的数据
        for kl = 1:line
    %十进制转为二进制，首先进行数据整数化处理：
            %数据依次提取
            for kn = 1:N
                unit1( kn, :) = NEWP{ kn, t}( kl, :);  %第一次数据提取：每个个体的某一条带的数据，规整到一个数组中
            end
            
            for kd = 1:D(t)
                data = round( unit1( :, kd) * BIT);  %第二次数据提取：每个个体某一维度的数据，规整到一个数组里并取四位小数
                %十进制转为二进制
                BIN  = dec2bin( data );
                L    = size(BIN,2);
                %%
        %基于概率的 交叉 操作  (单点交叉)
                Parents = randperm(N) ;            %亲代组合序列
                for j = 1:2:N
                    if rand <= Pc
                        local = find(randperm( L )==1);%确定交叉位置,对此修改可改为双点交叉
                        temp  = BIN( Parents(j+1), local:end); 
                        BIN( Parents(j+1) , local:L ) = BIN( Parents(j), local:L);
                        BIN( Parents(j), local:L) = temp;
                    end
                end
                %%
        %基于概率的 变异 操作
                for j = 1 : N
                    if rand <= Pm   %由概率决定是否发生变异
                        local = find( randperm(L)==1 );%确定变异位置
                        if BIN( j, local) == '1'
                            BIN( j, local) = '0';
                        else
                            BIN( j, local) ='1';
                        end%单点变异
                    end
                end
        %将二进制转回十进制
                data = bin2dec(BIN) / BIT;
                %边界条件判定：[0,1]
                for j1 = 1:N
                    if data( j1) <0 || data( j1) >1
                        data( j1) = rand;
                    end
                end
                unit2(:,kd) = data ;  %对应第二次数据提取
            end
            
            %将新数据返回NEWP
            for kn = 1:N
                NEWP{ kn, t}( kl,:) = unit2( kn, :) ; %对应第一次数据提取
            end
            
        end
        
    end %结束一个维度的更新
    
end

