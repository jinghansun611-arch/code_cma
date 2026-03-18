function seed = Fun_CMAES_multires_seed(bestfile_stage1, D_new, n_new, varargin)
%FUN_CMAES_MULTIRES_SEED  多分辨率精修：从粗分辨率最优解生成“高控制点 + 高插值点”的 CMA-ES 初值。
%
% 你现在的痛点：
%   - 用较少控制点(GA.D) + 较少插值点(tar.n) 粗搜能找到大体构型，但局部贴合不够；
%   - 直接把 GA.D、tar.n 拉高再从头跑，维度暴涨容易卡住、预算也不够。
%
% 本函数做的事：
%   1) 读取 stage1 保存的 best_stage1.mat（其中包含 BEST_P_G / GA / tar 等）；
%   2) 在不改变目标构型方程/点集的前提下，把 tar.n 提高到 n_new，并重算 tar.x/y/z/s/K1/location/cross；
%   3) 把粗搜得到的 BEST_P_G 在每个 type/line 上“重采样”为新的控制点数量 D_new；
%   4) 重新 build map，并把重采样后的 PP 转成 x0（用于 CMA-ES 的均值初始化）；
%   5) 给出一个 sigma0 的建议值（可选用）。
%
% 用法（推荐放到你的精修主脚本最前面，在 build_map 之前）：
%   A_set_tar_equation;
%   seed = Fun_CMAES_multires_seed(bestfile_stage1, [21,21,1,1,1], 101);
%   GA.D = seed.D_new;
%   tar.n = seed.n_new;
%   map  = seed.map;
%   GA.CMA.x0     = seed.x0;
%   GA.CMA.sigma0 = seed.sigma0_suggest;
%
% 输入：
%   bestfile_stage1 : stage1 的 best_stage1.mat 路径
%   D_new           : 1×GA.TYPE 的新控制点数，例如 [21,21,1,1,1]
%   n_new           : 新的插值点数 tar.n，例如 101
%
% 可选参数（name-value）：
%   'InterpMethod'  : 控制点重采样方法（默认 'pchip'，也可 'spline'/'linear'）
%   'Clamp01'       : 是否把 PP 限制到 [0,1]（默认 true）
%   'UseStage1Tar'  : 是否直接使用 stage1 文件里保存的 tar（默认 false；推荐 false，以保证与当前 A_set_tar_equation 一致）
%
% 输出 seed 结构体字段：
%   seed.PP0              : 1×TYPE 的 cell（新分辨率下的初始个体 PP）
%   seed.x0               : nvars×1 的向量（新 map 下的初始均值）
%   seed.map              : 新的 map（对应新 GA.D / tar.symmetry / judge 等）
%   seed.D_old / D_new    : 旧/新控制点数
%   seed.n_old / n_new    : 旧/新插值点数
%   seed.sigma0_suggest   : 建议的 sigma0（经验值，可自行调整）
%   seed.best_del_stage1  : stage1 最优误差
%
% 注意：本函数依赖你工程里的全局变量风格。

global GA tar

%% ----------------- 解析可选参数 -----------------
ip = inputParser;
ip.addParameter('InterpMethod', 'pchip', @(s)ischar(s) || isstring(s));
ip.addParameter('Clamp01', true, @(x)islogical(x) || isnumeric(x));
ip.addParameter('UseStage1Tar', false, @(x)islogical(x) || isnumeric(x));
ip.parse(varargin{:});
opt = ip.Results;
opt.InterpMethod = char(opt.InterpMethod);

%% ----------------- 读 stage1 最优 -----------------
if ~exist(bestfile_stage1,'file')
    error('找不到 bestfile_stage1：%s', bestfile_stage1);
end

S = load(bestfile_stage1);

if ~isfield(S,'BEST_P_G') || isempty(S.BEST_P_G)
    error('stage1 文件中缺少 BEST_P_G，无法进行多分辨率精修初始化。');
end
PP_old = S.BEST_P_G;   % 1×TYPE cell

% 旧分辨率
if isfield(S,'GA') && isfield(S.GA,'D')
    D_old = S.GA.D;
else
    % 兜底：从 PP_old 推断
    D_old = zeros(1, numel(PP_old));
    for t = 1:numel(PP_old)
        D_old(t) = size(PP_old{t},2);
    end
end
if isfield(S,'tar') && isfield(S.tar,'n')
    n_old = S.tar.n;
else
    n_old = NaN;
end

% stage1 最优误差
best_del_stage1 = NaN;
if isfield(S,'BEST_D_G')
    best_del_stage1 = S.BEST_D_G;
end

%% ----------------- 设置新分辨率（GA.D + tar.n） -----------------
TYPE = GA.TYPE;
if numel(D_new) ~= TYPE
    error('D_new 的长度必须等于 GA.TYPE=%d。当前 D_new 长度=%d。', TYPE, numel(D_new));
end

% 一些安全提示：对称12最好用偶数控制点（你现在默认是11，不受影响）
if isfield(tar,'symmetry') && ~isempty(tar.symmetry)
    sym = tar.symmetry;
    if any(sym(:,3)==12)
        if mod(D_new(1),2)~=0 || (TYPE>=2 && mod(D_new(2),2)~=0)
            warning('检测到 symmetry=12（反对称）。建议 type1/type2 的控制点数用偶数，否则 Fun_set_special 中 dd==d/2 逻辑可能不触发。');
        end
    end
end

GA.D = D_new;
tar.n = n_new;

%% ----------------- 目标曲线/点集重新离散（更新 tar.x/y/z/s/location/cross/center/K1） -----------------
% 你工程里两种目标定义方式：
%   - 方程（Fun_cal_para_tar_eqation）
%   - 点集（Fun_cal_para_tar_point_set）
% 这里按字段自动选择。

if opt.UseStage1Tar && isfield(S,'tar') && ~isempty(S.tar)
    % 用 stage1 的 tar 作为“方程/点集定义”，只把 n_new 覆盖后重算离散
    tar = S.tar;
    tar.n = n_new;
end

if isfield(tar,'linelist') && ~isempty(tar.linelist)
    Fun_cal_para_tar_point_set;
else
    Fun_cal_para_tar_eqation;
end

%% ----------------- 控制点重采样：PP_old -> PP0_new -----------------
PP0_new = cell(1, TYPE);
L = tar.line;

for t = 1:TYPE
    if t > numel(PP_old) || isempty(PP_old{t})
        % 没有该 type 的旧数据：用 0.5 初始化
        PP0_new{t} = 0.5 * ones(L, D_new(t));
        continue;
    end

    V_old = PP_old{t};
    % 兼容：有些 type 可能是 1×D 或 L×D
    if size(V_old,1) == 1 && L > 1
        V_old = repmat(V_old, L, 1);
    end

    d0 = size(V_old,2);
    d1 = D_new(t);

    if d0 == d1
        V_new = V_old;
    elseif d0 == 1
        % 旧的是标量，新的变成多控制点：直接复制（对 type=3/4/5 最常见）
        V_new = repmat(V_old(:,1), 1, d1);
    else
        u0 = linspace(0,1,d0);
        u1 = linspace(0,1,d1);
        V_new = zeros(size(V_old,1), d1);
        for l = 1:size(V_old,1)
            v = V_old(l,:);
            % 采用 pchip 更稳，不容易 overshoot；也允许用户切换 spline
            vv = interp1(u0, v, u1, opt.InterpMethod);
            % interp1 在极端情况下可能出 NaN（例如输入含 NaN）
            if any(isnan(vv))
                vv = interp1(u0, v, u1, 'linear', 'extrap');
            end
            V_new(l,:) = vv;
        end
    end

    if opt.Clamp01
        V_new = min(1, max(0, V_new));
    end
    PP0_new{t} = V_new;
end

% 重新施加对称/全局变量（lambda/theta_pre）等规则，保证是“可行初值”
PP0_new = Fun_set_special(PP0_new);

%% ----------------- 新 map + 新 x0 -----------------
map_new = Fun_CMAES_build_map();
x0_new  = Fun_CMAES_PP2vec(PP0_new, map_new);

%% ----------------- sigma0 建议值（经验：升维/加密后步长要更小） -----------------
sigma0_suggest = 0.10;
try
    if isfield(S,'GA') && isfield(S.GA,'CMA') && isfield(S.GA.CMA,'sigma0')
        sig1 = S.GA.CMA.sigma0;
        % 一个稳健的经验缩放：
        %   - 在粗搜基础上做局部精修：sigma 先减半
        %   - 升维后再按 sqrt(n_old/n_new) 缩放，避免“抖动过大”
        nvars_old = NaN;
        if isfield(S,'map') && isfield(S.map,'nvars')
            nvars_old = S.map.nvars;
        elseif isfield(S,'BEST_X_G') && ~isempty(S.BEST_X_G)
            nvars_old = numel(S.BEST_X_G);
        end
        if ~isnan(nvars_old)
            sigma0_suggest = 0.5 * sig1 * sqrt(nvars_old / map_new.nvars);
        else
            sigma0_suggest = 0.5 * sig1;
        end
        sigma0_suggest = max(0.02, min(0.20, sigma0_suggest));
    end
catch
    % 保持默认
end

%% ----------------- 输出 -----------------
seed = struct();
seed.bestfile = bestfile_stage1;
seed.best_del_stage1 = best_del_stage1;
seed.D_old = D_old;
seed.D_new = D_new;
seed.n_old = n_old;
seed.n_new = n_new;
seed.PP0   = PP0_new;
seed.x0    = x0_new(:);
seed.map   = map_new;
seed.sigma0_suggest = sigma0_suggest;

fprintf('[MultiRes] Stage1 DEL=%g,  D: %s -> %s,  n: %g -> %g,  new nvars=%d,  sigma0~%.3g\n', ...
    seed.best_del_stage1, mat2str(seed.D_old), mat2str(seed.D_new), seed.n_old, seed.n_new, seed.map.nvars, seed.sigma0_suggest);

end
