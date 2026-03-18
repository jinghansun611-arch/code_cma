function [status, reason] = check_abaqus_job_status(path, jobname)
% 严格判断 Abaqus/Standard job 是否真正成功完成
% DONE    : .sta 明确写出分析已成功完成，且 .dat 无致命错误
% FAILED  : .sta/.dat/launcher 明确提示失败，或命令已结束但未出现成功标志
% RUNNING : .lck 存在，或 launcher 已启动但仍未结束
% PENDING : 尚未看到明确文件

status = "PENDING";
reason = "";

logf = fullfile(path, [jobname, '.log']);
staf = fullfile(path, [jobname, '.sta']);
datf = fullfile(path, [jobname, '.dat']);
lckf = fullfile(path, [jobname, '.lck']);
launcherf = fullfile(path, ['launcher_', jobname, '.out']);

if isfile(lckf)
    status = "RUNNING";
    reason = '.lck still exists';
    return;
end

sta_txt = safe_read_text(staf);
dat_txt = safe_read_text(datf);
launcher_txt = safe_read_text(launcherf);
log_txt = safe_read_text(logf);

usta = upper(sta_txt);
udat = upper(dat_txt);
ulaunch = upper(launcher_txt);
ulog = upper(log_txt);

% 1) 最严格的成功：sta 明确完成，且 dat 无 fatal error
if contains(usta, 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY') && ...
   ~contains(udat, '***ERROR') && ...
   ~contains(udat, 'A FATAL ERROR HAS OCCURRED') && ...
   ~contains(udat, 'THE ANALYSIS HAS NOT BEEN COMPLETED')
    status = "DONE";
    reason = 'sta says completed successfully';
    return;
end

% 2) 明确失败：sta / dat / launcher 有显式失败
if contains(usta, 'THE ANALYSIS HAS NOT BEEN COMPLETED') || ...
   contains(usta, 'A FATAL ERROR HAS OCCURRED') || ...
   contains(udat, '***ERROR') || ...
   contains(udat, 'A FATAL ERROR HAS OCCURRED') || ...
   contains(udat, 'THE ANALYSIS HAS NOT BEEN COMPLETED') || ...
   contains(udat, 'TOO MANY ATTEMPTS MADE FOR THIS INCREMENT') || ...
   contains(udat, 'THE TIME INCREMENT REQUIRED IS LESS THAN THE MINIMUM')
    status = "FAILED";
    reason = 'sta/dat reports failure';
    return;
end

% 3) launcher 已结束，但没有出现成功标志
if contains(ulaunch, 'EXITCODE=')
    if ~contains(ulaunch, 'END ') && ~contains(ulaunch, ' END ')
        status = "RUNNING";
        reason = 'launcher started but not ended';
        return;
    end

    if contains(ulaunch, 'EXITCODE=0')
        % 仅说明命令完成，不说明分析真正成功
        if ~isempty(sta_txt) || ~isempty(dat_txt)
            status = "FAILED";
            reason = 'launcher ended but no success flag in sta';
        else
            status = "FAILED";
            reason = 'launcher ended early and no sta/dat evidence';
        end
        return;
    else
        status = "FAILED";
        reason = 'launcher exitcode nonzero';
        return;
    end
end

% 4) 兼容某些环境存在 log 且写出 COMPLETED
if contains(ulog, 'COMPLETED')
    status = "DONE";
    reason = 'log says completed';
    return;
end

% 5) 有 sta/dat 但无明确成功失败，暂视为 RUNNING
if ~isempty(sta_txt) || ~isempty(dat_txt)
    status = "RUNNING";
    reason = 'sta/dat exists but no final marker yet';
    return;
end
end

function txt = safe_read_text(fname)
if ~isfile(fname)
    txt = '';
    return;
end
try
    txt = fileread(fname);
catch
    try
        txt = strjoin(cellstr(readlines(fname)), newline);
    catch
        txt = '';
    end
end
end