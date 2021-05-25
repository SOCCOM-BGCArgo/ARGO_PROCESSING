function Loop_ARGO_isus(float_num)

% Script to step through all floats in list for ARGO and create ODV-type
% txt file of isus diagnostics for use in N-Viz. http://www3.mbari.org/chemsensor/eviz.htm
%
% Created by Josh Plant, MBARI, May2020
%
% Updates:
%   TM, June 2, 2020, turned Loop_ARGO_dura into a function for calling within
%       Loop_ARGO_float after routine float processing.  Only called for floats
%       with new incoming msg files.
%   JP Aug 17,2020 JP converted Loop_ARGO_isus to a function using Tanya's pH adapatation
%      as a model
%  
%   TM, 8/17/20 finalized conversion to function for call within Loop_Argo_float.m
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

% PATHS & NAMES - 

user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir = [user_dir, '\Documents\MATLAB\'];

dirs.temp  = 'C:\temp\';
dirs.msg   = '\\atlas\ChemWebData\floats\';
% dirs.cal   = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\';
% dirs.FV    = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
%     'DATA\FLOATVIZ\'];
dirs.cal   = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\';
dirs.FV    = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\';
dirs.save = [user_dir,'\ARGO_PROCESSING\', ...
    'DATA\NO3_DIAGNOSTICS\'];

% % exclude_floats = '6966HAWAII|7622HAWAII';

% % % LOAD MASTER FLOAT LIST & DO SOME FILTERING
% % load([dirs.cal,'MBARI_float_list.mat'])
% %
% % % NOW STEP THROUGH FLOATS & PROCESS NO3 DIAGNOSTIC (*isus) FILES
% % rlist = size(list,1);
% %
% % for ct = 1:rlist
% %     MBARI_ID = list{ct,1};
% %
% %     if regexpi(MBARI_ID, exclude_floats,'once')
% %         disp([MBARI_ID, ' is on the exclusion list and will not be processed']);
% %         continue
% %     end
MBARI_ID = float_num;

cal_fn = ['cal',MBARI_ID,'.mat'];
load([dirs.cal,cal_fn]);

if cal.info.isus_flag == 0
    disp(['No NO3 cal data found for ',MBARI_ID,' moving to next float']);
    return
end

disp(['Processing isus  message files for ', MBARI_ID])
d = Merge_isus_msgs(MBARI_ID, dirs);
end


