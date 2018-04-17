% Script to step through all floats in list for ARGO

% ************************************************************************
% FILTER EXPRESSIONS
% ************************************************************************
% %refine_expr  = ''; % REFINE FLOAT LIST FOR REGION
% refine_expr  = '^9\d+SoOcn'; % SOUTHERN OCEAN 9000 series
% %refine_expr  = 'SoOcn|Drake|Ross|SoAtl'; % REFINE FOR SOUTHERN OCEAN 
% exclude_expr = '(MTY)|(cor)|(Surface0)|(\d+\.txt)';

% ************************************************************************
% dirs  = A structure with directory strings where files are
%         located.
%
%   dirs.msg       = path to float message file directories
%   dirs.mat       = path to *.mat profile files for ARGO NETCDF
%   dirs.cal       = path to *.cal (nitrate) and cal####.mat (working float cal) files
%   dirs.config    = path to config.txt files
%   dirs.NO3config = path to ####.cal nitrate calibration files
%   dirs.temp      = path to temporary working dir
%   dirs.FV        = path to FloatViz files served to web
%   dirs.QCadj      = path to QC adjustment list for all floats
%   dirs.FVlocal   = path to Floatviz file made by matlab go here

% PATHS & NAMES - JOSH DEFAULTS
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    dirs.mat       = 'C:\Users\jplant\Documents\MATLAB\ARGO\DURA\';
    dirs.temp      = 'C:\temp\';
% ************************************************************************

Float_type_list2

%FLOAT_LIST = SOCCOM_list;
% FLOAT_LIST = NON_SOCCOM_list;
%FLOAT_LIST = [NON_SOCCOM_list;SOCCOM_list];

% jp = regexpi(NON_SOCCOM_list(:,1),'RosSea|SoAtl|SoPac|SoOcn|drake');
% t1 = cellfun(@isempty,jp);
% FLOAT_LIST = [SOCCOM_list; NON_SOCCOM_list(~t1,:)];
FLOAT_LIST = [SOCCOM_list; NON_SOCCOM_list];

[~, ind] = sort(FLOAT_LIST(:,2)); % sort by UW ID
FLOAT_LIST = FLOAT_LIST(ind,:);

% % ***************** SOCCOM ONLY *****************

clearvars -except SOCCOM_list NON_SOCCOM_list FLOAT_LIST update_str dirs


NO_GO_LIST = cell(1000,5);
NO_GO_ct   = 0;
for loop_ctr = 1: size(FLOAT_LIST,1)
    if FLOAT_LIST{loop_ctr,4} < 100 % All
        flt_str = FLOAT_LIST{loop_ctr,2};
        if strcmp('APEX',FLOAT_LIST{loop_ctr,5});
            Merge_dura_msgs(flt_str, dirs)
        end
    end
end










