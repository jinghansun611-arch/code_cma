function export_cae_one(T_target, Phe_target, workdir)
%EXPORT_CAE_ONE  导出指定代/个体的 CAE（只建模+保存，不求解）
%
% 放在 GA-Cylinder-1 下也没问题，因为函数内部会 cd 到 workdir。
%
% 用法：
%   export_cae_one(99, 9)
%   export_cae_one(99, 9,'D:\CMAES_inverse_design_pack\code_cma\A-CMAES-3')

    if nargin < 2
        error('用法: export_cae_one(T_target, Phe_target, [workdir])');
    end
    if nargin < 3 || isempty(workdir)
        workdir = 'D:\CMAES_inverse_design_pack\code_cma_withtwist\A-CMAES-8';
    end

    oldpwd = pwd;
    c = onCleanup(@() cd(oldpwd));
    cd(workdir);

    % --- 检查 ---
    if exist('mat-cae.mat', 'file') ~= 2
        error('找不到 mat-cae.mat：%s', fullfile(workdir,'mat-cae.mat'));
    end
    target_mat = sprintf('mat-%d-%d.mat', T_target, Phe_target);
    if exist(target_mat, 'file') ~= 2
        error('找不到目标 mat 文件：%s', fullfile(workdir, target_mat));
    end
    if exist('python_export_cae.py', 'file') ~= 2
        error('找不到 python_export_cae.py：%s', fullfile(workdir,'python_export_cae.py'));
    end

    % --- 备份 mat-cae.mat ---
    backup = 'mat-cae_backup_export.mat';
    copyfile('mat-cae.mat', backup, 'f');

    % --- 写入导出指令到 FEA（只写数值字段，scipy 读取最稳） ---
    S = load('mat-cae.mat');
    if ~isfield(S,'FEA')
        error('mat-cae.mat 中找不到 FEA 变量');
    end
    FEA = S.FEA;

    FEA.EXPORT_MODE = 1;
    FEA.EXPORT_T    = T_target;
    FEA.EXPORT_PHE  = Phe_target;

    save('mat-cae.mat', 'FEA');

    % --- 调用 Abaqus 执行导出脚本（noGUI） ---
    cmd = 'abaqus cae script=python_export_cae.py';
    status = system(cmd);

    % --- 恢复 mat-cae.mat ---
    copyfile(backup, 'mat-cae.mat', 'f');

    if status ~= 0
        error('导出失败（Abaqus返回码=%d）。请检查 workdir 下的日志输出。', status);
    end

    fprintf('✅ 导出完成：CAE-T%d-Phe%d.cae（目录：%s）\n', T_target, Phe_target, workdir);
end
