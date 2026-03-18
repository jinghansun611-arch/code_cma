function result = abaqus_odb_m( path )
global tar T Phe
num_line = tar.line;
%运行python脚本，读取odb数据库的历史输出，并将参数写入Abaqus工作目录下的txt中
matlabpath = pwd() ;     %记录当前matlab目录，pwd()确认当前文件夹，返回其路径
cd(path);                %进入abaqus目录，即inp文件所在目录

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

save('mat-odb.mat','Phe','T');
%ABAQUS_CMD = '"C:\SIMULIA\Commands\abq2024.bat"'; %这个远程控制别的电脑新安装在C盘里的时候默认这个路径
ABAQUS_CMD = '"D:\Abaqus2024\SIMULIA\Commands\abq2024.bat"';
system([ABAQUS_CMD ' cae noGUI=python_odb_M.py']);
%system('abaqus cae noGUI=python_ODB_M.py');%调用python脚本
%system('abaqus cae script=Odb_path.py');
%读取结果
result_x = cell(num_line,1);
result_z = cell(num_line,1);
result_y = cell(num_line,1);
data = cell(3,1);
for k = 1:num_line
    data{1} = readmatrix( ['TXT-RESULT-X-',num2str(k),'.txt' ] );
    data{2} = readmatrix( ['TXT-RESULT-Y-',num2str(k),'.txt' ] );
    data{3} = readmatrix( ['TXT-RESULT-Z-',num2str(k),'.txt' ] );
    %增设：异常插值点修正
    %若某个点与前后两点的距离差距过大，将其与附近点进行线性插值
    if isempty(data{1}) == 1
        continue
    end
    for j = 1:3
        range = max( data{j}(:,2) ) - min( data{j}(:,2) );
        num_point = size(data{j},1);
        cor = range / (num_point-1) * 5 ;
        for i = 2:num_point-1
            del = abs( data{j}(i,2) - data{j}(i-1,2) );
            if del > cor
                data{j}(i,2) = ( data{j}(i-1,2) + data{j}(i+1,2) )/2;
            end
        end
    end
    result_x{k} = data{1};
    result_y{k} = data{2};
    result_z{k} = data{3};
        
        
end

cd(matlabpath)%返回matlab工作目录
Fun_showlogfile(path);%显示python运行的输出信息

result.rx = result_x;
result.ry = result_y;
result.rz = result_z;


end




