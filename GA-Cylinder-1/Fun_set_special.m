function PN = Fun_set_special( PP )
%对于特定情况才启用，否则为空白
global tar
PN = PP;
    %symmetry[:,3] == 11  ：单条带关于中点对称,第二个参数不使用
    %symmetry[:,3] == 12  ：单条带关于中点反对称
    %symmetry[:,3] == 21  ：左右对称，宽度复制，扭转角(若有)取反 【 若要上下反对称，在参数方程中令Y取反值 】
    %symmetry[:,3] == 22  ：左右对称，宽度复制，扭转角(若有)取反 【 若要上下反对称，在参数方程中令Y取反值 】
%例 ： 对称条带
if isfield(tar,'symmetry') == 1 %对称判断
    [ n, type] = size(PP);
    num_sym = size(tar.symmetry,1);
    line = tar.symmetry ;
    for i = 1:type %变量类型
        % 每种变量类型的控制点数量可能不同（D(1) 不一定等于 D(2)）。
        % 旧代码固定读取 PP{1,1} 的列数，会导致 D 不一致时对称映射错误。
        % 这里改为对每个 i 动态读取列数，保证对称约束与真实维度一致。
        [ ~, d ] = size(PP{1,i});
        for j = 1:n %个体总数
            for k = 1:num_sym %对称处理数
                
                if line(k,3) == 11 %单条带关于中点对称
                    for dd = 1:d/2
                        if i == 1
                            PN{ j, i}( line(k,1), d - dd + 1) = PN{ j, i}( line(k,1), dd);
                        elseif i == 2
                            PN{ j, i}( line(k,1), d - dd + 1) = 1 - PN{ j, i}( line(k,1), dd); %截面扭转角取反
                        end
                    end
                    
                elseif line(k,3) == 12   %单条带关于中点反对称
                    for dd = 1:d/2
                        if i == 1
                            PN{ j, i}( line(k,1), d - dd + 1) = PN{ j, i}( line(k,1), dd);
                        elseif i == 2
                            PN{ j, i}( line(k,1), d - dd + 1) = PN{ j, i}( line(k,1), dd);
                        end
                        if dd == d/2 %这一段其实没太懂 为什么要偶数项的时候除2
                            PN{ j, i}( line(k,1), dd) = PN{ j, i}( line(k,1), dd-1)/2;
                        end
                    end   
                
                elseif line(k,3) == 21 % 两条带左右对称
                    PN{ j, i}( line(k,2), :) = PN{ j, i}( line(k,1), :);
                    if i ==2
                        PN{ j, i}( line(k,2), :) = 1 - PN{ j, i}( line(k,2), :); %第二条带截面扭转角取反
                    end
                    
                elseif line(k,3) == 22 % 左右反对称
                    PN{ j, i}( line(k,2), :) = PN{ j, i}( line(k,1), :);
                    if i ==2
                        PN{ j, i}( line(k,2), :) = PN{ j, i}( line(k,2), :);
                    end
                
                 
                    
                end
                
            end %END : for k = 1:s %对称处理数
        end %END : for j = 1:n %个体总数
    end %END : for i = 1:type %变量类型
    
end %END : if isfield(tar,'symmetry') == 1 %对称判断

%============================ (新增) 第4类变量：lambda 作为圆柱全局属性 ============================
% 说明：GA内部仍然按“每条带一个值”存放PP{4}(line,1)，这里强制把它变成全局同一个值
if size(PN,2) >= 4
    [n, ~] = size(PN);                 % n = 个体数
    for j = 1:n
        lam = PN{j,4}(1,1);            % 取第1条带的lambda作为全局lambda
        PN{j,4}(:,:) = lam;            % 复制给所有条带
    end
%============================ (新增) 第5类变量：预扭转 theta_pre 作为全局属性 ============================
% 说明：同 lambda，一般作为圆柱/夹具的全局加载量；这里强制个体内所有条带相同
if size(PN,2) >= 5
    [n, ~] = size(PN);
    for j = 1:n
        thpre = PN{j,5}(1,1);
        PN{j,5}(:,:) = thpre;
    end
end

end

end

