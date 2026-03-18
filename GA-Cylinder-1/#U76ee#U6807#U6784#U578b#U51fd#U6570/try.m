%% ===== 目标3（优化版）：低周向跨度的手性双峰条带 =====
GA.ORDER = 'N3R';
GA.Structure = '-chiralDoublePeakSoft';

n = 101;                      % 插值点数
r_simu = 50;
r = r_simu;
TYPE = 5;

% ---------------- FEA 参数（按你现有习惯可再微调） ----------------
FEA.load_disturb   = 5e-7;
FEA.load_pressall  = 5e-8;
FEA.load_pressbond = 5e-8;
FEA.mesh = [0.5; 0.5];
FEA.cpus = 2;
FEA.P1_AZ = r + 0.1;
FEA.P2_AZ = 0.2;

% ---------------- 几何主参数 ----------------
x  = 52;          % 轴向长度
rr = 8.5;         % 径向抬升幅值（比之前更温和）

t_lim = [0, 1];

% funr：双峰但更平滑，避免过尖、过激进
funr = r ...
     + rr * (1 - cos(2*pi*t)) ...
     .* (1 + 0.10*cos(4*pi*t));

% funa：显著缩小周向跨度
% 组成：
% 1) 一个较弱的整体手性偏置项
% 2) 一个温和的一次波动项
% 3) 一个更小的二次起伏项
%
% 整体目标：让条带有明显手性，但周向展开不要过大
funa = 0.18*pi*(t - 0.5) ...
     + 0.08*pi*sin(2*pi*t) ...
     + 0.04*pi*sin(4*pi*t);

funx = x * t;
funy = funr * sin(funa);
funz = funr * cos(funa) - r;

LINE = size(funx,1);

% ---------------- 定位/交点/键合 ----------------
location = [1, 0, 0, 1, 0, 0];
cross = [];
bond = [2];

% ---------------- 开关 ----------------
JUDGE = zeros(TYPE, LINE);
JUDGE(1,:) = 1;   % TYPE1 宽度分布
JUDGE(2,:) = 1;   % TYPE2 截面扭转分布
JUDGE(4,1) = 1;   % TYPE4 lambda_x
JUDGE(5,1) = 1;   % TYPE5 主动预扭转

% ---------------- 变量范围（你可按经验再调） ----------------
range(1,:) = [0.5, 1.5];        % 宽度 0.5~2.0
range(2,:) = [-pi, 2*pi];       % 截面扭转角
range(3,:) = [-pi/36, pi/18];   % TYPE3 当前不用，但保留
range(4,:) = [1.00, 0.5];       % lambda_x: 1.00~1.50
range(5,:) = [-pi/6, pi/3];     % TYPE5 主动扭转，建议先稍收窄一些

% ---------------- 写入 tar ----------------
tar.r = r_simu;
tar.r_simu = r_simu;

tar.n = n;
tar.t_lim = t_lim;
tar.location = location;
tar.cross = cross;
tar.bond = bond;
tar.funr = funr;
tar.funx = funx;
tar.funy = funy;
tar.funz = funz;
tar.line = LINE;
tar.range = range;
tar.judge_type_line = JUDGE;
GA.TYPE = TYPE;