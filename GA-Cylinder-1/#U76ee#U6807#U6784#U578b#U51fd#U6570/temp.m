clear; clc; close all;

%% ===== 参数 =====
r_simu = 50;
r = r_simu;

x  = 52;
rr = 8.5;

n = 401;
t = linspace(0,1,n)';

%% ===== 方程 =====
funr = r ...
     + rr * (1 - cos(2*pi*t)) ...
     .* (1 + 0.10*cos(4*pi*t));

funa = 0.18*pi*(t - 0.5) ...
     + 0.08*pi*sin(2*pi*t) ...
     + 0.04*pi*sin(4*pi*t);

funx = x * t;
funy = funr .* sin(funa);
funz = funr .* cos(funa) - r;

%% ===== 画图 =====
figure('Color','w','Position',[100 60 1350 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% --- 图1：三维构型 ---
nexttile(1);
hold on; grid on; axis equal;

xx = linspace(min(funx), max(funx), 80);
th = linspace(0,2*pi,80);
[TT,XX] = meshgrid(th, xx);
YY = r * sin(TT);
ZZ = r * cos(TT) - r;

surf(XX,YY,ZZ, ...
    'FaceAlpha',0.08, ...
    'EdgeAlpha',0.08, ...
    'FaceColor',[0.4 0.6 0.9], ...
    'EdgeColor',[0.4 0.6 0.9]);

plot3(funx, funy, funz, 'LineWidth', 2.5);
plot3(funx(1), funy(1), funz(1), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
plot3(funx(end), funy(end), funz(end), 's', 'MarkerSize', 8, 'LineWidth', 1.5);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D 目标曲线');
view(35,22);

% --- 图2：X-Y 投影 ---
nexttile(2);
plot(funx, funy, 'LineWidth', 2.0);
grid on; axis equal;
xlabel('X'); ylabel('Y');
title('X-Y 投影');

% --- 图3：X-Z 投影 ---
nexttile(3);
plot(funx, funz, 'LineWidth', 2.0);
grid on; axis equal;
xlabel('X'); ylabel('Z');
title('X-Z 投影');

% --- 图4：funr 和 funa ---
nexttile(4);
yyaxis left
plot(t, funr, 'LineWidth', 2.0);
ylabel('funr');

yyaxis right
plot(t, rad2deg(funa), 'LineWidth', 2.0);
ylabel('funa (deg)');

grid on;
xlabel('t');
title('径向抬升 funr 与周向角 funa');
legend('funr','funa','Location','best');

sgtitle('N3R / -chiralDoublePeakSoft 目标构型预览');