function x = Fun_CMAES_PP2vec(PP, map)
%FUN_CMAES_PP2VEC  将原工程的 PP(cell) 按 map 规则抽取为向量 x。
%
% 用途：
%   - 在 Fun_set_special 强制对称/全局约束之后，得到“可行域内”的 x
%     用于 CMA-ES 的 tell/update（让分布学习在可行域里发生）。
%
% 输入：
%   PP  : 1×GA.TYPE 的 cell（例如 P(i,:,T)）
%   map : Fun_CMAES_build_map 输出
% 输出：
%   x   : n×1

x = zeros(map.nvars, 1);
for k = 1:map.nvars
    t = map.type(k);
    l = map.line(k);
    c = map.col(k);
    x(k) = PP{t}(l,c);
end

end
