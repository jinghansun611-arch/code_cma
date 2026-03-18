clear; clc; close all;
addpath(genpath(pwd));
global tar

A_set_tar_equation;
Fun_cal_para_tar_eqation;   % 生成 tar.x/y/z
Fun_cal_para_2D;            % 生成 sketch_center / sketch_U 等（你改URX的地方就在这里）

disp('--- check radius fields ---');
if isfield(tar,'r_simu'), disp(['tar.r_simu = ' num2str(tar.r_simu)]); end
if isfield(tar,'r_bond'), disp(['tar.r_bond = ' num2str(tar.r_bond)]); end
disp(['tar.r (legacy) = ' num2str(tar.r)]);

% 下面这两个变量通常在 Fun_cal_para_2D 里生成（若名称不同就只看你自己的输出）
whos sketch_center sketch_U
