function tf_float =Process_GUI_float_GLT(handles,dirs)
% Script to process one float at a time
% update_str = 'all';
% flt_str    = '9095';

tic;
% GET SOME INFO
tf_float     = 0;
float_IDs    = handles.float_IDs; %MBARI name, UW_ID, WMO#, type
MBARI_ID_str = handles.info.float_name;
MFLOATlist = load([dirs.cal,'MBARI_float_list.mat']);
t1           = strcmpi(MBARI_ID_str,MFLOATlist.list(:,1));
float_type   = MFLOATlist.list{t1,4};


% ************************************************************************
% dirs  = A structure with directory strings where files are
%         located.
%
%   dirs.mat       = path to *.mat profile files for ARGO NETCDF
%   dirs.cal       = path to *.cal (nitrate) and cal####.mat (working float cal) files
%   dirs.NO3config = path to ####.cal nitrate calibration files
%   dirs.FVlocal   = path to Floatviz file made by matlab go here
%   dirs.temp      = path to temporary working dir
%
%   dirs.msg       = path to float message file directories
%   dirs.config    = path to config.txt files
%   dirs.FV        = path to FloatViz files served to web
%   dirs.QCadj     = path to QC adjustment list for all floats


update_str = 'all';

% ************************************************************************

%Make FloatViz subdirectories if they don't exist
if ~exist([dirs.FVlocal,'QC',filesep], 'dir') 
    mkdir([dirs.FVlocal,'QC',filesep])
end
if ~exist([dirs.FVlocal,'HR',filesep], 'dir') 
    mkdir([dirs.FVlocal,'HR',filesep])
end
if ~exist([dirs.FVlocal,'HRQC',filesep], 'dir') 
    mkdir([dirs.FVlocal,'HRQC',filesep])
end

% set(handles.recumpute_text,'String', ...
%     {'Updating *.mat profile files for ', ...
%     [MBARI_ID_str, '(',float_type,') ...']})
% 
% drawnow

% ************************************************************************
% GET BAD SENSOR LIST
bad_sensor_list = parse_bad_sensor_list([dirs.cal,'bad_sensor_list.txt']);
iM   = find(strcmp('MBARI ID STR',bad_sensor_list.hdr) == 1);

% CHECK IF SPECIFC FLOAT HAS BAD SENSOR ISSUES
if ~isempty(bad_sensor_list.list)
    tSENSOR = strcmpi(MBARI_ID_str,bad_sensor_list.list(:,iM));
    if sum(tSENSOR) > 0
        disp([MBARI_ID_str,' found on the bad sensor list!'])
        BSL = bad_sensor_list;
        BSL.list = BSL.list(tSENSOR,:);
        dirs.BSL = BSL;
        clear BSL
    else
        dirs.BSL.hdr  = [];
        dirs.BSL.list = [];
    end
    clear tSENSOR
end

% ************************************************************************
% PROCESS MESSAGE FILES
if strcmp(float_type,'APEX')
    tf_float = Process_APEX_float(MBARI_ID_str, dirs,update_str);
    
elseif strcmp(float_type,'NAVIS')
    tf_float = Process_NAVIS_float(MBARI_ID_str, dirs,update_str);
else
    disp(['Unknown float type for ',MBARI_ID_str, ...
        '! processing next float'])
    return
end

% ***********************************************************************
% CREATE ODV COMPATIBLE TEXT FILES FOR FLOATVIZ AND SOCCOMVIZ
% set(handles.recumpute_text,'String', ...
%     {'Updating ODV compatible .txt files for ', ...
%     [MBARI_ID_str, '(',float_type,') ...']})
% drawnow

ODV_tf = argo2odv_LIAR(MBARI_ID_str, dirs, update_str);
if ODV_tf == 0
    disp(['A LIAR ODV FILE WAS NOT CREATED OR UPDATED FOR ',MBARI_ID_str])
else
    disp(['A NEW LIAR ODV FILE WAS CREATED FOR ',MBARI_ID_str])
end


elapsedTime = toc;
disp(['Processes finished in ',num2str((elapsedTime)/60, ...
    '%0.2f'), ' minutes'])










