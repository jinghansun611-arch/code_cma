function cma = Fun_CMAES_init(nvars, popsize, x0, sigma0)
%FUN_CMAES_INIT  初始化 CMA-ES 状态结构体。
%
% 该实现遵循经典 CMA-ES（Hansen 等）的标准参数设置，
% 但为了适配你工程的“变量范围 [0,1] + Abaqus 黑盒 + 高代价评价”，做了两点务实改动：
%   1) 采样后用 Fun_bound_reflect01 把候选解反射到 [0,1]
%   2) 保留 elitism 注入接口（可选），减少噪声/偶发失败导致的退化
%
% 输入：
%   nvars   : 自由变量维度 n
%   popsize : 种群规模 λ（建议直接用 GA.N）
%   x0      : 初始均值（n×1），可为空 -> 默认 0.5
%   sigma0  : 初始步长（标量），可为空 -> 默认 0.25
%
% 输出：
%   cma : CMA-ES 状态结构体（ask/tell 循环使用）

if nargin < 1 || isempty(nvars)
    error('Fun_CMAES_init: nvars is required');
end
if nargin < 2 || isempty(popsize)
    popsize = 4 + floor(3*log(nvars));
end
if nargin < 3 || isempty(x0)
    x0 = 0.5 * ones(nvars,1);
else
    x0 = x0(:);
end
if nargin < 4 || isempty(sigma0)
    sigma0 = 0.25;
end

lambda = popsize;
mu = floor(lambda/2);

% 对数权重（更稳）
weights = log(mu + 0.5) - log(1:mu)';
weights = weights / sum(weights);
mueff = 1 / sum(weights.^2);

% --- 策略参数（标准推荐值） ---
cc   = (4 + mueff/nvars) / (nvars + 4 + 2*mueff/nvars);
cs   = (mueff + 2) / (nvars + mueff + 5);
c1   = 2 / ((nvars + 1.3)^2 + mueff);
cmu  = min(1 - c1, 2 * (mueff - 2 + 1/mueff) / ((nvars + 2)^2 + mueff));
damps = 1 + 2*max(0, sqrt((mueff - 1)/(nvars + 1)) - 1) + cs;

% 期望 ||N(0,I)||
chiN = sqrt(nvars) * (1 - 1/(4*nvars) + 1/(21*nvars^2));

% 初始分布
cma.nvars  = nvars;
cma.lambda = lambda;
cma.mu     = mu;
cma.weights = weights;
cma.mueff   = mueff;

cma.m      = x0;                 % 均值
cma.sigma  = sigma0;             % 步长
cma.sigma0 = sigma0;

cma.C  = eye(nvars);             % 协方差
cma.B  = eye(nvars);
cma.D  = ones(nvars,1);
cma.invsqrtC = eye(nvars);

cma.pc = zeros(nvars,1);
cma.ps = zeros(nvars,1);

cma.cc = cc;
cma.cs = cs;
cma.c1 = c1;
cma.cmu = cmu;
cma.damps = damps;
cma.chiN = chiN;

cma.eigeneval = 0;
cma.counteval = 0;
cma.generation = 0;

% --- 可选：精英注入（默认开启 1 个） ---
cma.elite_inject = true;
cma.elite_count  = 1;
cma.xbest = x0;
cma.fbest = inf;

end
