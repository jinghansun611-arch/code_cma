function Fun_showlogfile(path)
% 更清晰地显示 pylog.txt 摘要：
% 1) 如果有错误关键词，只打印关键异常行
% 2) 如果没有明显错误，只打印最后几行
% 3) 不再整份刷屏

matlabpath = pwd();
cd(path);

logfile = 'pylog.txt';
if exist(logfile, 'file') ~= 2
    cd(matlabpath);
    return;
end

txt = '';
try
    txt = fileread(logfile);
catch
    try
        txt = strjoin(cellstr(readlines(logfile)), newline);
    catch
        txt = '';
    end
end

if isempty(txt)
    cd(matlabpath);
    return;
end

lines = regexp(txt, '\r\n|\n|\r', 'split');
lines = string(lines);
lines = strtrim(lines);
lines = lines(strlength(lines) > 0);

if isempty(lines)
    cd(matlabpath);
    return;
end

% 过滤掉纯 license 噪声（不完全过滤，只去掉最常见的）
mask_noise = contains(upper(lines), 'ABAQUS LICENSE MANAGER CHECKED OUT') | ...
             contains(upper(lines), 'LICENSES REMAIN AVAILABLE');
lines_show = lines(~mask_noise);

% 错误关键词
mask_err = contains(upper(lines_show), 'ERROR') | ...
           contains(upper(lines_show), 'FAILED') | ...
           contains(upper(lines_show), 'EXCEPTION') | ...
           contains(upper(lines_show), 'TRACEBACK') | ...
           contains(upper(lines_show), 'ODBERROR');

if any(mask_err)
    fprintf('---- pylog 摘要（检测到异常关键词） ----\n');
    err_lines = lines_show(mask_err);
    for i = 1:numel(err_lines)
        fprintf('%s\n', err_lines(i));
    end
    fprintf('--------------------------------------\n');
else
    % 没有明显异常，只打印最后几行摘要
    tail_n = min(8, numel(lines_show));
    fprintf('---- pylog 摘要（最后 %d 行） ----\n', tail_n);
    for i = numel(lines_show)-tail_n+1 : numel(lines_show)
        fprintf('%s\n', lines_show(i));
    end
    fprintf('---------------------------------\n');
end

cd(matlabpath);
end