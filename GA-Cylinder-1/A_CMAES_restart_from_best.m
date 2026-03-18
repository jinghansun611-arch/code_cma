function A_CMAES_restart_from_best(bestfile, new_TMAX, new_sigma0)
% 从 best_stage1.mat 里的当前最优解出发，重新开一轮 CMA-ES 精修
%
% 用法示例：
%   A_CMAES_restart_from_best('E:\xxx\best_stage1.mat', 20, 0.12)
%
% 含义：
%   - 从 best_stage1.mat 读出 BEST_X_G
%   - 把 BEST_X_G 作为新的 CMA-ES 初始中心 x0
%   - 用更小的 sigma0 再跑 new_TMAX 代
%
% 注意：
%   这不是精确续跑第31代，而是“从第30代最优解继续精修”

if nargin < 1 || isempty(bestfile)
    error('请提供 best_stage1.mat 路径。');
end
if nargin < 2 || isempty(new_TMAX)
    new_TMAX = 20;
end
if nargin < 3 || isempty(new_sigma0)
    new_sigma0 = 0.12;
end

if exist(bestfile, 'file') ~= 2
    error('找不到文件：%s', bestfile);
end

S = load(bestfile);

need_vars = {'BEST_X_G','BEST_D_G','BEST_NUM','GA','tar'};
for i = 1:numel(need_vars)
    if ~isfield(S, need_vars{i})
        error('best_stage1.mat 缺少变量：%s', need_vars{i});
    end
end

fprintf('\n==================== 从 best 重新启动 CMA-ES ====================\n');
fprintf('best 文件        : %s\n', bestfile);
fprintf('原最优 DEL       : %.8g\n', S.BEST_D_G);
fprintf('原最优位置       : T=%d, Phe=%d\n', S.BEST_NUM(1), S.BEST_NUM(2));
fprintf('新的 T_MAX       : %d\n', new_TMAX);
fprintf('新的 sigma0      : %.6g\n', new_sigma0);
fprintf('===============================================================\n');

% 把 best 作为新一轮优化的初始中心
GA = S.GA;
GA.T_MAX = new_TMAX;
GA.CMA.x0 = S.BEST_X_G(:)';
GA.CMA.sigma0 = new_sigma0;

% 可以根据需要再加一点设置
% 例如精英注入始终开启
GA.CMA.elite_inject = true;
GA.CMA.elite_warmup = 0;

% 把修改后的 GA 传回 base workspace
assignin('base', 'GA_restart_seed', GA);
assignin('base', 'BEST_restart_seed', S.BEST_X_G(:)');

fprintf('\n已把新的 GA 参数放到 base workspace：\n');
fprintf('  GA_restart_seed\n');
fprintf('  BEST_restart_seed\n');

fprintf('\n下一步操作：\n');
fprintf('1) 打开你的主函数\n');
fprintf('2) 在参数设定区后面插入下面两句：\n\n');
fprintf('   if evalin(''base'',''exist(''''GA_restart_seed'''',''''var'''')'')\n');
fprintf('       GA = evalin(''base'',''GA_restart_seed'');\n');
fprintf('   end\n\n');
fprintf('3) 然后正常运行主函数即可\n');
fprintf('===============================================================\n');
end