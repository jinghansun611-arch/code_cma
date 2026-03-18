function [X, cma] = Fun_CMAES_ask(cma)
%FUN_CMAES_ASK  根据当前分布采样 λ 个候选解。
%
% 输出：
%   X   : λ×nvars，且每个元素都在 [0,1]
%   cma : 更新后的 cma（主要是 counteval 等）

n = cma.nvars;
lambda = cma.lambda;

% 标准正态采样
arz = randn(n, lambda);
ary = cma.B * (cma.D .* arz);                    % ~ N(0, C)
arx = repmat(cma.m, 1, lambda) + cma.sigma * ary;

% 边界处理：[0,1]
arx = Fun_bound_reflect01(arx);

X = arx';

% 精英注入：把历史最优解放进当前种群（可减小噪声/偶发失败导致的退化）
if isfield(cma,'elite_inject') && cma.elite_inject && isfinite(cma.fbest)
    k = min(cma.elite_count, lambda);
    for i = 1:k
        X(i,:) = cma.xbest(:)';
    end
end

% 记录本代采样（调试/复现可用）
cma.arz = arz;
cma.ary = ary;
cma.arx = arx;

end
