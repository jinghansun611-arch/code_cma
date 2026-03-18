function map = Fun_CMAES_build_map()
%FUN_CMAES_BUILD_MAP  构建“向量变量 x ↔ 细胞数组 PP”的映射关系。
%
% 你的工程里，优化变量以 cell 形式存放：
%   PP{type} 的尺寸通常为 [tar.line, GA.D(type)]，元素取值范围 [0,1]
%   Fun_cal_para_2D 会再用 tar.range 把它们映射到实际物理量范围。
%
% CMA-ES 的核心变量是一个连续向量 x ∈ R^n。
% 为了在不大改原有 Abaqus/Matlab 管线的前提下替换 GA，这里需要：
%   1) 定义哪些 PP 元素是“自由变量”（真正影响几何/边界条件）
%   2) 把这些自由变量按固定顺序拉直为向量 x
%   3) 由向量 x 反向重建 PP
%
% 关键点（考虑你当前代码里已经存在的规则）：
%   - tar.judge_type_line(type,line)=0 的参数不启用（对应条带该类变量固定）
%   - tar.symmetry 可能会把某些条带/控制点强制对称（在 Fun_set_special 里执行）
%       * 两条带左右对称(21/22)：第二条带是第一条带的复制/取反 → 第二条带的变量不应再作为自由变量
%       * 单条带关于中点对称/反对称(11/12)：后半控制点由前半生成 → 后半控制点不应再作为自由变量
%   - 第4/5类变量(lambda/pretwist)在你的工程里是“圆柱全局属性”：Fun_set_special 会强制所有条带相同
%       → 把 type=4 当成 1 个标量自由变量（而不是 tar.line 个）
%
% 输出 map 结构体：
%   map.nvars    : 自由变量总维度 n
%   map.type     : 每个自由变量对应的 type
%   map.line     : 每个自由变量对应的 line
%   map.col      : 每个自由变量对应的 col
%   map.template : 1×GA.TYPE 的 cell，每个元素是 [tar.line, GA.D(type)] 的初值(默认0.5)
%
% 说明：该函数依赖全局 tar / GA（与你现有工程风格保持一致）。

global tar GA

TYPE = GA.TYPE;
D    = GA.D;
L    = tar.line;

% -------- 1) 取 judge（允许 tar.judge_type_line 行数 < TYPE 的兼容） --------
if isfield(tar,'judge_type_line') && ~isempty(tar.judge_type_line)
    J = tar.judge_type_line;
else
    J = ones(TYPE, L);
end
if size(J,1) < TYPE
    J = [J; ones(TYPE-size(J,1), L)];
end

% -------- 2) 解析 symmetry：找出“依赖条带”和“依赖控制点” --------
dep_line = false(L,1);      % 依赖条带（由其它条带生成），对 type=1/2 生效
half_cols = false(TYPE,1);  % 是否存在单条带对称（11/12），需要只保留前半控制点
half_rule = zeros(L,1);     % 记录每条带的 11/12 规则（0=无，11/12=对应规则）

if isfield(tar,'symmetry') && ~isempty(tar.symmetry)
    sym = tar.symmetry;
    for k = 1:size(sym,1)
        base = sym(k,1);
        dep  = sym(k,2);
        code = sym(k,3);
        if code == 21 || code == 22
            if dep >= 1 && dep <= L
                dep_line(dep) = true;
            end
        elseif code == 11 || code == 12
            if base >= 1 && base <= L
                half_rule(base) = code;
            end
        end
    end
end

% -------- 3) 生成 template（默认0.5更“中性”） --------
template = cell(1,TYPE);
for t = 1:TYPE
    template{t} = 0.5 * ones(L, D(t));
end

% -------- 4) 构建自由变量索引表 --------
type_list = [];
line_list = [];
col_list  = [];

for t = 1:TYPE
    % type=3/4 通常是标量（D=1），且在你的工程里只取 PP{t}(1,1) 使用
    if (t == 3 || t == 4 || t == 5) && D(t) == 1
        if any(J(t,:) == 1)
            type_list(end+1,1) = t; %#ok<AGROW>
            line_list(end+1,1) = 1;
            col_list (end+1,1) = 1;
        end
        continue;
    end

    for l = 1:L
        if J(t,l) ~= 1
            continue; % 该条带该类变量不启用
        end

        % 两条带左右对称(21/22)：第二条带由第一条带生成 → 对 type=1/2 排除依赖条带
        if (t == 1 || t == 2) && dep_line(l)
            continue;
        end

        % 单条带中点对称/反对称(11/12)：只保留前半控制点作为自由变量（后半由 Fun_set_special 生成）
        if (t == 1 || t == 2) && (half_rule(l) == 11 || half_rule(l) == 12)
            keep_cols = 1:ceil(D(t)/2);
        else
            keep_cols = 1:D(t);
        end

        for c = keep_cols
            type_list(end+1,1) = t; %#ok<AGROW>
            line_list(end+1,1) = l;
            col_list (end+1,1) = c;
        end
    end
end

map.nvars    = numel(type_list);
map.type     = type_list;
map.line     = line_list;
map.col      = col_list;
map.template = template;

% 便于调试：打印一次维度信息
% fprintf('[CMA] decision dimension = %d\n', map.nvars);

end
