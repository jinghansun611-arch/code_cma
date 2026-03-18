function abaqus_odb_open( path, T, Phe )
%运行python脚本，读取odb数据库的历史输出，并将参数写入Abaqus工作目录下的txt中

matlabpath = pwd() ;     %记录当前matlab目录，pwd()确认当前文件夹，返回其路径
cd(path);                %进入abaqus目录，即inp文件所在目录

save('mat-odb.mat','Phe','T');
%system('abaqus cae script=python_ODB_open.py');%调用python脚本
%ABAQUS_CMD = '"C:\SIMULIA\Commands\abq2024.bat"'; %这个远程控制别的电脑新安装在C盘里的时候默认这个路径
ABAQUS_CMD = '"D:\Abaqus2024\SIMULIA\Commands\abq2024.bat"';
system([ABAQUS_CMD ' cae script=python_ODB_open.py']);

cd(matlabpath)%返回matlab工作目录
end

