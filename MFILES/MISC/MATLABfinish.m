%MATLABfinish.m
%
% Tanya Maurer
% MBARI
% Dec 7, 2020
% Script to safely exit MATLAB when launched from batch file such that
% stale MATLAB sessions are not left hanging in the task manager!
%--------------------------------------------------------------------------

id  = feature('getpid');
if ispc
  cmd = sprintf('Taskkill /PID %d /F',id);
elseif (ismac || isunix)
  cmd = sprintf('kill -9 %d',id);
else
  disp('unknown operating system');
end
system(cmd);

% end