%% preview_target_centerlines.m
% 功能：批量/选择性预览“目标构型函数”生成的目标质心线三维图
% 用法：在 GA-Cylinder-1 根目录运行本脚本即可

clear; clc; close all;

addpath(genpath(pwd));   % 确保能找到工程内所有函数

% ====== 你只需要改这一个路径（如果你的文件夹不在根目录下）======
targetFolder = fullfile(pwd, '目标构型函数');   % 你截图里的那个文件夹
% ======================================================================

% 可选：是否自动保存图片
SAVE_PNG = true;
outFolder = fullfile(pwd, 'TargetPreviewPNGs');
if SAVE_PNG && ~exist(outFolder,'dir'), mkdir(outFolder); end

% 可选：是否每看完一个暂停（按任意键看下一个）
PAUSE_EACH = true;

% 扫描目标构型函数
files = dir(fullfile(targetFolder, '*.m'));
names = {files.name};

% 过滤掉不想预览的文件（按需增删）
skipList = {'A_set_tar_equation.m', 'Result_plot.m'};
names = names(~ismember(names, skipList));

if isempty(names)
    error('在 %s 下没找到 .m 文件，请检查 targetFolder 路径。', targetFolder);
end

% 让你选择要预览哪些目标函数
[idx, tf] = listdlg( ...
    'ListString', names, ...
    'SelectionMode', 'multiple', ...
    'PromptString', '选择要预览的目标构型函数（可多选）', ...
    'ListSize', [420 320] ...
);
if ~tf
    disp('已取消选择。');
    return;
end

% 逐个预览
for ii = 1:numel(idx)
    fname = names{idx(ii)};
    targetName = erase(fname, '.m');

    fprintf('\n=== 预览: %s ===\n', targetName);

    % 核心：执行目标构型函数 -> 生成 tar -> 画图
    preview_one_target(targetFolder, targetName, SAVE_PNG, outFolder);

    if PAUSE_EACH
        disp('按任意键继续预览下一个...');
        pause;
        close(gcf);
    end
end

disp('预览完成。');


%% ====================== 本文件的局部函数（不用单独存） ======================
function preview_one_target(targetFolder, targetName, SAVE_PNG, outFolder)
    % 说明：你工程里 tar 是 global，Fun_cal_para_tar_eqation 会用到它
    global tar GA FEA %#ok<NUSED>

    % 尽量“清空上一次目标”的残留，避免串台
    tar = struct();

    % ---- 1) 执行目标构型函数（兼容：函数/脚本两种写法）----
    % 首选 feval（如果它是 function），失败再用 run（如果它是 script 或函数名不匹配）
    try
        feval(targetName);
    catch
        run(fullfile(targetFolder, [targetName '.m']));
    end

    % ---- 2) 生成离散目标质心线 tar.x / tar.y / tar.z ----
    % 这一步用你工程原本的函数
    Fun_cal_para_tar_eqation;

    % ---- 3) 绘制 ----
    nLines = size(tar.x, 2);
    figure('Color','w'); hold on; grid on; axis equal;

    for k = 1:nLines
        x = tar.x(:,k); y = tar.y(:,k); z = tar.z(:,k);
        x = x(:); y = y(:); z = z(:);

        h = plot3(x, y, z, 'LineWidth', 2);

        % 起点/终点标记（可删）
        plot3(x(1),   y(1),   z(1),   'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', h.Color, 'MarkerEdgeColor', 'none');
        plot3(x(end), y(end), z(end), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', h.Color, 'MarkerEdgeColor', 'none');
    end

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Target centerline preview: %s', targetName), 'Interpreter','none');
    view(45, 25);
    rotate3d on;

    % 打印一下信息，方便你确认是不是“三条线/二维”
    zRange = max(tar.z(:)) - min(tar.z(:));
    fprintf('  -> 检测到线条数 nLines = %d\n', nLines);
    fprintf('  -> tar.z range = %.6g\n', zRange);

    % ---- 4) 可选：保存 PNG ----
    if SAVE_PNG
        outPng = fullfile(outFolder, [targetName '.png']);
        exportgraphics(gcf, outPng, 'Resolution', 300);
        fprintf('  -> 已保存: %s\n', outPng);
    end
end
