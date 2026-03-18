function cma = Fun_CMAES_tell(cma, X, f)
%FUN_CMAES_TELL  用 (X,f) 更新 CMA-ES 分布（协方差/步长/均值）。
%
% 输入：
%   X : λ×nvars，候选解（建议是“可行域内”的值；你的工程里可在 Fun_set_special 之后再 PP2vec）
%   f : λ×1 或 1×λ，目标函数（越小越好）
%
% 输出：
%   cma : 更新后的 CMA-ES 状态

X = double(X);
f = double(f(:));

lambda = cma.lambda;
n = cma.nvars;

if size(X,1) ~= lambda
    error('Fun_CMAES_tell: X must have %d rows (lambda).', lambda);
end
if size(X,2) ~= n
    error('Fun_CMAES_tell: X must have %d columns (nvars).', n);
end
if numel(f) ~= lambda
    error('Fun_CMAES_tell: f must have length %d.', lambda);
end

% 计数/代数
cma.generation = cma.generation + 1;
cma.counteval  = cma.counteval + lambda;

% 记录历史最优
[fmin, imin] = min(f);
if fmin < cma.fbest
    cma.fbest = fmin;
    cma.xbest = X(imin,:).';
end

% 排序（最小化问题：f 小更好）
[~, idx] = sort(f, 'ascend');
X = X(idx, :);

% 重组：取前 mu 个
mu = cma.mu;
weights = cma.weights;

xold = cma.m;
xsel = X(1:mu, :).';            % n×mu

% 新均值
cma.m = xsel * weights;

% y = (x - xold) / sigma
y = (xsel - xold) / cma.sigma;  % n×mu

% 进化路径 ps
cma.ps = (1 - cma.cs) * cma.ps + ...
    sqrt(cma.cs * (2 - cma.cs) * cma.mueff) * (cma.invsqrtC * (cma.m - xold) / cma.sigma);

% hsig
ps_norm = norm(cma.ps);
hsig = ps_norm / sqrt(1 - (1 - cma.cs)^(2*cma.generation)) / cma.chiN < (1.4 + 2/(n+1));

% 进化路径 pc
cma.pc = (1 - cma.cc) * cma.pc + hsig * sqrt(cma.cc * (2 - cma.cc) * cma.mueff) * ((cma.m - xold) / cma.sigma);

% 协方差更新
% rank-one
C = (1 - cma.c1 - cma.cmu) * cma.C + cma.c1 * (cma.pc * cma.pc');
% 当 hsig 为 0 时的修正项
if ~hsig
    C = C + cma.c1 * cma.cc * (2 - cma.cc) * cma.C;
end

% rank-mu
% sum_i w_i * y_i * y_i'
for i = 1:mu
    C = C + cma.cmu * weights(i) * (y(:,i) * y(:,i)');
end

% 数值对称化
C = triu(C) + triu(C,1)';
cma.C = C;

% 步长更新
cma.sigma = cma.sigma * exp((cma.cs / cma.damps) * (ps_norm / cma.chiN - 1));

% 特征分解更新（不是每代都做，节省开销）
if cma.counteval - cma.eigeneval > (lambda / (cma.c1 + cma.cmu) / n / 10)
    cma.eigeneval = cma.counteval;
    % 防止数值不正定：加一个很小的 jitter
    C = (C + C')/2;
    jitter = 1e-12;
    [B, D] = eig(C + jitter*eye(n));
    D = real(diag(D));
    % 如果出现非正特征值，进一步抬升
    if any(D <= 0)
        D(D <= 0) = jitter;
    end
    cma.B = real(B);
    cma.D = sqrt(D);
    cma.invsqrtC = cma.B * diag(1./cma.D) * cma.B';
end

end
