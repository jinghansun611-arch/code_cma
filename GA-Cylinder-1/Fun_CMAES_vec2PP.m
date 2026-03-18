function PP = Fun_CMAES_vec2PP(x, map)
%FUN_CMAES_VEC2PP  将 CMA-ES 的向量变量 x 还原为原工程使用的 PP(cell)。
%
% 输入：
%   x   : 1×n 或 n×1，取值建议在 [0,1]
%   map : Fun_CMAES_build_map 输出
% 输出：
%   PP  : 1×GA.TYPE 的 cell，每个 cell 为 [tar.line, GA.D(type)]
%
% 说明：
%   - 非自由变量默认填 0.5（“中性”值），随后会被 Fun_set_special 修正对称/全局约束。
%   - type=3/4/5 若启用，则视作标量，默认会写入 PP{t}(1,1)，其余行列保持0.5。

global GA tar

x = x(:); % column

PP = map.template;

for k = 1:map.nvars
    t = map.type(k);
    l = map.line(k);
    c = map.col(k);
    PP{t}(l,c) = x(k);
end

% 数值精度对齐（保持你原来的 GA.BIT 逻辑：避免写入过长小数影响后续 dec2bin 等）
if isfield(GA,'BIT') && ~isempty(GA.BIT)
    for t = 1:GA.TYPE
        PP{t} = round(PP{t} * GA.BIT) / GA.BIT;
    end
end

% 保险：确保都在 [0,1]
for t = 1:GA.TYPE
    PP{t}(PP{t} < 0) = 0;
    PP{t}(PP{t} > 1) = 1;
end

% 第4类变量(lambda)是全局属性：让所有条带保持一致（也可依赖 Fun_set_special 再做一次）
if GA.TYPE >= 4 && numel(PP) >= 4
    lam = PP{4}(1,1);
    PP{4}(:,:) = lam;
end

% 第5类变量（预扭转 theta_pre）同理：统一为标量
if GA.TYPE >= 5 && numel(PP) >= 5
    thpre = PP{5}(1,1);
    PP{5}(:,:) = thpre;
end

% 第3类变量（单条带随机角）同理：统一为标量
if GA.TYPE >= 3 && numel(PP) >= 3
    th = PP{3}(1,1);
    PP{3}(:,:) = th;
end

end
