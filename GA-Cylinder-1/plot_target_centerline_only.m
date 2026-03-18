%% plot_target_centerline_only.m
% 只绘制目标构型质心线（三维）

clear; clc; close all;

addpath(genpath(pwd));   % 确保能找到工程函数

global tar

% 1) 生成目标构型（得到 tar.x / tar.y / tar.z）
A_set_tar_equation;
Fun_cal_para_tar_eqation;

% 2) 基本检查
if ~isfield(tar,'x') || isempty(tar.x)
    error('tar.x 为空：请检查 A_set_tar_equation / Fun_cal_para_tar_eqation 是否正确生成目标曲线。');
end

nLines = size(tar.x, 2);   % 目标线条数量（更可靠）
fprintf('检测到目标质心线条数 nLines = %d\n', nLines);

zmin = min(tar.z(:)); zmax = max(tar.z(:));
fprintf('tar.z 范围: [%.6g, %.6g], range = %.6g\n', zmin, zmax, zmax - zmin);
if abs(zmax - zmin) < 1e-6
    warning('tar.z 几乎不变：图看起来会像二维。这通常说明目标曲线本身在一个平面内或 z 定义没生效。');
end

% 3) 绘图
figure('Color','w'); hold on; grid on; axis equal;
for k = 1:nLines
    x = tar.x(:,k); y = tar.y(:,k); z = tar.z(:,k);
    x = x(:); y = y(:); z = z(:);   % 保证列向量

    h = plot3(x, y, z, 'LineWidth', 2);

    % 起点/终点标记（可删）
    plot3(x(1),   y(1),   z(1),   'o', 'MarkerSize', 6, ...
          'MarkerFaceColor', h.Color, 'MarkerEdgeColor', 'none');
    plot3(x(end), y(end), z(end), 's', 'MarkerSize', 6, ...
          'MarkerFaceColor', h.Color, 'MarkerEdgeColor', 'none');
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Target centerline(s)');
view(45, 25);     % 强制一个三维视角
rotate3d on;      % 鼠标可拖拽旋转
