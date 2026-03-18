function [result, ok, msg] = abaqus_odb_m(path)
global tar T Phe
num_line = tar.line;
matlabpath = pwd();
cd(path);

ok = false;
msg = '';
result_x = cell(num_line,1);
result_y = cell(num_line,1);
result_z = cell(num_line,1);
result.rx = result_x;
result.ry = result_y;
result.rz = result_z;

runtime_dir = fullfile(path, '_abq_runtime');
tmp_dir     = fullfile(runtime_dir, 'tmp');
home_dir    = fullfile(runtime_dir, 'home');
if ~exist(runtime_dir, 'dir'); mkdir(runtime_dir); end
if ~exist(tmp_dir, 'dir'); mkdir(tmp_dir); end
if ~exist(home_dir, 'dir'); mkdir(home_dir); end
setenv('TMP',  tmp_dir);
setenv('TEMP', tmp_dir);
setenv('HOME', home_dir);
setenv('PYTHONIOENCODING', 'utf-8');

% 清理旧 TXT，避免误读上一轮
for k = 1:num_line
    fx = fullfile(path, ['TXT-RESULT-X-', num2str(k), '.txt']);
    fy = fullfile(path, ['TXT-RESULT-Y-', num2str(k), '.txt']);
    fz = fullfile(path, ['TXT-RESULT-Z-', num2str(k), '.txt']);
    if exist(fx,'file') == 2, delete(fx); end
    if exist(fy,'file') == 2, delete(fy); end
    if exist(fz,'file') == 2, delete(fz); end
end

save('mat-odb.mat','Phe','T');
ABAQUS_CMD = '"C:\SIMULIA\Commands\abq2024.bat"';
ret = system([ABAQUS_CMD ' cae noGUI=python_odb_M.py']);
if ret ~= 0
    msg = sprintf('python_odb_M.py 返回码非0：%d', ret);
    cd(matlabpath);
    return;
end

data = cell(3,1);
for k = 1:num_line
    fx = fullfile(path, ['TXT-RESULT-X-', num2str(k), '.txt']);
    fy = fullfile(path, ['TXT-RESULT-Y-', num2str(k), '.txt']);
    fz = fullfile(path, ['TXT-RESULT-Z-', num2str(k), '.txt']);
    if ~isfile(fx) || ~isfile(fy) || ~isfile(fz)
        msg = sprintf('缺少 TXT-RESULT 文件（line=%d）', k);
        cd(matlabpath);
        return;
    end
    try
        data{1} = readmatrix(fx);
        data{2} = readmatrix(fy);
        data{3} = readmatrix(fz);
    catch ME
        msg = sprintf('读取 TXT-RESULT 失败（line=%d）：%s', k, ME.message);
        cd(matlabpath);
        return;
    end

    if isempty(data{1})
        msg = sprintf('TXT-RESULT 为空（line=%d）', k);
        cd(matlabpath);
        return;
    end

    for j = 1:3
        range = max(data{j}(:,2)) - min(data{j}(:,2));
        num_point = size(data{j},1);
        if num_point <= 1
            continue;
        end
        cor = range / (num_point-1) * 5;
        for i = 2:num_point-1
            del = abs(data{j}(i,2) - data{j}(i-1,2));
            if del > cor
                data{j}(i,2) = (data{j}(i-1,2) + data{j}(i+1,2))/2;
            end
        end
    end

    result_x{k} = data{1};
    result_y{k} = data{2};
    result_z{k} = data{3};
end

cd(matlabpath);
Fun_showlogfile(path);
result.rx = result_x;
result.ry = result_y;
result.rz = result_z;
ok = true;
msg = 'OK';
end