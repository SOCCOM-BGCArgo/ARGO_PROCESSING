function Loop_ARGO_dura(float_num,float_type)

% Script to step through all floats in list for ARGO and create ODV-type
% txt file of dura diagnostics for use in E-Viz. http://www3.mbari.org/chemsensor/eviz.htm
%
% Created by Josh Plant, MBARI, May2020
%
% Updates: TM, June 2, 2020, turned into a function for calling within
% Loop_ARGO_float after routine float processing.  Only called for floats
% with new incoming msg files.
%
% INPUTS: float_num = MBARI float ID (ie '9095SOOCN')
%         float_type = ie 'APEX' or 'NAVIS'
%
%
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
dirs.temp  = 'C:\temp\';
dirs.msg   = '\\atlas\ChemWebData\floats\';
% dirs.cal   = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\';
% dirs.FV    = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
%     'DATA\FLOATVIZ\'];
dirs.cal   = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\';
dirs.FV    = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\';
dirs.save = ['C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\', ...
    'DATA\PH_DIAGNOSTICS\'];

% GET PH SENSOR SN & VERSION LIST, WILL ADD INFO TO ODV TEXT FILES
% Pass to merge_dura fuction via dirs structure
% dirs.ph_ver = get_ph_version; % THIS WAS CHANGED INTO A MAT FILE; WINDOWS
% TASK SCHEDULER HAVING TROUBLE WITH XLSREAD!!!
D = load([dirs.save,'MBARI_pHsensor_versions.mat']);
dirs.ph_ver = D.d;


% % LOAD MASTER FLOAT LIST & DO SOME FILTERING
% load([dirs.cal,'MBARI_float_list.mat'])
% tNAVIS = strcmp(list(:,4),'NAVIS');
% list(tNAVIS,:) =[];


% for ct = 1:rlist
%     MBARI_ID = list{ct,1};
if strcmp(float_type,'NAVIS') == 1 %dura diagnostic files do not exist for NAVIS floats
    return
end
MBARI_ID = float_num;
cal_fn = ['cal',MBARI_ID,'.mat'];
load([dirs.cal,cal_fn]);

if cal.info.pH_flag == 0
    disp(['No pH cal data found for ',float_type,' float ',MBARI_ID,'; moving to next float.']);
    return
end

disp(['Processing pH diagnostic data for ',float_type,' float ', MBARI_ID,'.'])
d = Merge_dura_msgs(MBARI_ID, dirs);
close all
% end

%% Move this to full copy batchfile, keep all copies in same place for now
% % COPY FILES TO NETWORK WITH BATCH FILE
% system(['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\batchfiles\', ...
%     'Copy_Eviz_to_network.bat'])