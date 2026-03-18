function del = cal_delta_C_m_new(result)
% 改进版目标曲线误差函数（稳定修正版）
% 特点：
% 1) 按弧长重参数化，避免简单按节点序号比较
% 2) 加入切向误差，约束曲线走向
% 3) 加入主峰位置/高度误差，增强几何特征一致性
% 4) 对重复点/零弧长段进行容错，避免 interp1 报错
%
% author: ChatGPT

global tar

center = tar.center;
num_line = tar.line;

cx = center(1);
cy = center(2);
cz = center(3);

res_data = result.resdata;
tar_data = result.tardata;

% ---------------- 权重 ----------------
w_pos  = 0.65;   % 位置误差权重
w_tan  = 0.20;   % 切向误差权重
w_peak = 0.15;   % 主峰误差权重

% ---------------- 初始化 ----------------
E_pos_num = 0;
E_pos_den = 0;

E_tan_sum = 0;
N_tan_sum = 0;

E_peak_sum = 0;
N_peak_sum = 0;

for k = 1:num_line
    % ===== 读取单条曲线 =====
    xr = res_data{k,1}(:);
    yr = res_data{k,2}(:);
    zr = res_data{k,3}(:);

    xt = tar_data{k,1}(:);
    yt = tar_data{k,2}(:);
    zt = tar_data{k,3}(:);

    % 点数太少直接跳过
    if length(xr) < 2 || length(xt) < 2
        continue;
    end

    % ===== 计算弧长并去除重复采样点 =====
    [sr, xr, yr, zr] = local_arc_length_unique(xr, yr, zr);
    [st, xt, yt, zt] = local_arc_length_unique(xt, yt, zt);

    % 去重后点数仍然太少，跳过
    if length(sr) < 2 || length(st) < 2
        continue;
    end

    % ===== 弧长重参数化 =====
    n_sample = max([length(xr), length(xt), 200]);
    sq = linspace(0,1,n_sample)';

    xrq = interp1(sr, xr, sq, 'pchip');
    yrq = interp1(sr, yr, sq, 'pchip');
    zrq = interp1(sr, zr, sq, 'pchip');

    xtq = interp1(st, xt, sq, 'pchip');
    ytq = interp1(st, yt, sq, 'pchip');
    ztq = interp1(st, zt, sq, 'pchip');

    % ===== 1. 位置误差 =====
    E_pos_num = E_pos_num + sum((xrq-xtq).^2 + (yrq-ytq).^2 + (zrq-ztq).^2);
    E_pos_den = E_pos_den + sum((xtq-cx).^2 + (ytq-cy).^2 + (ztq-cz).^2);

    % ===== 2. 切向误差 =====
    trx = gradient(xrq);
    try_ = gradient(yrq);
    trz = gradient(zrq);

    ttx = gradient(xtq);
    tty = gradient(ytq);
    ttz = gradient(ztq);

    nr = sqrt(trx.^2 + try_.^2 + trz.^2) + 1e-12;
    nt = sqrt(ttx.^2 + tty.^2 + ttz.^2) + 1e-12;

    trx  = trx  ./ nr;
    try_ = try_ ./ nr;
    trz  = trz  ./ nr;

    ttx = ttx ./ nt;
    tty = tty ./ nt;
    ttz = ttz ./ nt;

    dotval = trx.*ttx + try_.*tty + trz.*ttz;
    dotval = max(min(dotval,1),-1);

    E_tan_sum = E_tan_sum + sum(1 - dotval);
    N_tan_sum = N_tan_sum + length(dotval);

    % ===== 3. 主峰误差 =====
    [~, idr] = max(zrq);
    [~, idt] = max(ztq);

    x_span = max(xtq) - min(xtq);
    z_span = max(ztq) - min(ztq);

    if x_span < 1e-12
        x_span = 1;
    end
    if z_span < 1e-12
        z_span = 1;
    end

    E_peak_k = abs(xrq(idr) - xtq(idt)) / x_span ...
             + abs(zrq(idr) - ztq(idt)) / z_span;

    E_peak_sum = E_peak_sum + E_peak_k;
    N_peak_sum = N_peak_sum + 1;
end

% ---------------- 汇总 ----------------
if E_pos_den < 1e-12
    E_pos = 1e6;
else
    E_pos = sqrt(E_pos_num) / sqrt(E_pos_den);
end

if N_tan_sum < 1
    E_tan = 1e6;
else
    E_tan = E_tan_sum / N_tan_sum;
end

if N_peak_sum < 1
    E_peak = 1e6;
else
    E_peak = E_peak_sum / N_peak_sum;
end

% ---------------- 总误差 ----------------
del = w_pos * E_pos + w_tan * E_tan + w_peak * E_peak;

end


function [s_unique, x_unique, y_unique, z_unique] = local_arc_length_unique(x, y, z)
% 计算归一化累计弧长，并去除重复弧长点
% 这样可避免 interp1 报“采样点必须唯一”

x = x(:);
y = y(:);
z = z(:);

% 删除 NaN / Inf
valid = isfinite(x) & isfinite(y) & isfinite(z);
x = x(valid);
y = y(valid);
z = z(valid);

n = length(x);

if n == 0
    s_unique = [];
    x_unique = [];
    y_unique = [];
    z_unique = [];
    return;
elseif n == 1
    s_unique = [0;1];
    x_unique = [x;x];
    y_unique = [y;y];
    z_unique = [z;z];
    return;
end

% 原始累计弧长
ds = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
s = [0; cumsum(ds)];

% 若整条线几乎完全重合
if s(end) < 1e-14
    s_unique = [0;1];
    x_unique = [x(1); x(end)];
    y_unique = [y(1); y(end)];
    z_unique = [z(1); z(end)];
    return;
end

% 归一化
s = s / s(end);

% 去除重复采样点：保留每个唯一弧长第一次出现的位置
tol = 1e-12;
keep = [true; diff(s) > tol];

s_unique = s(keep);
x_unique = x(keep);
y_unique = y(keep);
z_unique = z(keep);

% 如果去重后只剩一个点，再强制补成两个点
if length(s_unique) == 1
    s_unique = [0;1];
    x_unique = [x_unique; x_unique];
    y_unique = [y_unique; y_unique];
    z_unique = [z_unique; z_unique];
end

end