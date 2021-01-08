function tf_float = Process_APEX_float(MBARI_ID_str, dirs, update_str)
% ************************************************************************
% PURPOSE:
%    This function processes raw message files for a given APEX float
%    (.msg, .isus, .dura), calculates useful values from the raw signals
%    (P, T, S, O2, NO3, pH, CHl, Bbp, CDOM) and then merges these results
%    for a given profile.
%
%    Each profile is saved as a *.mat file with the following format:
%    WMO_ID.PROFILE_#.mat (Ex 5904683.006.mat) in a directory using the
%    WMO_ID for its name. Each file contains 3 structures with variables
%    inside named using ARGO definitions:
%       LR   - Low resolution data
%       HR   - High resolution data
%       INFO - Some information about the profile and float
%
% USAGE:
%	tf_float = Process_APEX_float(MBARI_ID_str, dirs, update_str)
%
% INPUTS:
%   MBARI_ID_str  = MBARI Float ID as a string
%
%   update_str = "all" or "update"
%                 all    - to process all available msg files
%                 update - only process new msg files
%
%   dirs       = Either an empty variable or a structure with directory
%                strings where files are located. It must contain these
%                fields at a minimum or be empty. If dirs is empty default
%                paths will be used.
%
%   dirs.msg       = path to float message file directories
%   dirs.mat       = path to *.mat profile files for ARGO NETCDF
%   dirs.cal       = path to *.cal (nitrate) and cal####.mat (working float cal) files
%   dirs.config    = path to config.txt files
%   dirs.NO3config = path to ####.cal nitrate calibration files
%   dirs.temp      = path to temporary working dir
%   dirs.FV        = path to FloatViz files served to web
%   dirs.QCadj     = path to QC adjustment list for all floats
%   dirs.FVlocal   = path to Floatviz file made by matlab go here
%
% OUTPUTS:
%   tf_float =  1 if processing was a success, 0 otherwise
%
%   Files are saved to WMO dir. This output is intended for input to ARGO
%   NetCDF files via Annie Wong at Univ. of washington and for MBARI use
%   to build ODV compatible *.txt files for FloatViz and SOCCOMViz
%
% HELPER FUNCTIONS CALLED:
%   get_float_cals            get_FloatViz_data       calc_FLOAT_NO3
%   get_MBARI_WMO_list        calc_O2_4ARGO           theta
%   parseNO3cal               betasw_ZHH2009          density
%   get_QC_adjustments        apply_QC_corr
%   get_QCstep_dates          phcalc
%   get_last_cast             parse_pHmsg
%   get_msg_list              parse_NO3msg
%
% CHANGE LOG
% 02/23/2017 - If update string  = all, clear all *.mat files in dir 1st.
%       This fixed a potential issue with floats with duplicate UW_ID #'s
% 03/08/2017 - Added time delay to file processing. If msg file < 4 hrs old
%       don't process. This prevents a bunch of partial files from being
%       processed (~Line 180)
% 03/28/2017 - added range checks for IB_PH ansIK to flag bad pH from bad
%       IB. Added range checks for optode phase and optode T for bad O2
% 04/18/2017 - added code to create PSAL & TEMP QC and ADJUSTED variables.
%       This will be used to get S & T QC flags for the ODV files.
% 06/26/2017 - added code to each sensor processing section to querry bad
%       sensor list using the function, isbadsensor.m and set quality flags
%       for all profiles on the list to bad
% 06/27/2017 - added code to set unrealistic LR DOXY to "crazy_val"
% 08/01/2017 - fixed code to set QC for O2, NO3 & pH to bad if bad S or T
% 08/11/2017 - add code to catch deployed float with no msg files yet in
%              chemwebdata ~ line 211
% 08/28/17 added code to re-assign dirs.msg directory for special case
%   floats which include duplicate UW ID floats and floats with NO WMO
% 09/29/2017 -  Increased O2 max range from 450 to 550 to cover high O2
%               supersatuaration in artic ice edge blooms
% 10/02/2017 - tweaked code to wait for .isus and .dura files before
% processing.  Added this a couple of weeks ago, but made a minor tweak to
% keep these 'limbo' files included in the email list.
% 12/18/2017 - Added code to include ARGO QC flag == 5 for NPQ corrected data
%   & to process chl for in situ DC (gets added to cal file), NPQ corr.- j
% 02/5/2018 - Added code for including "<param>_DATA_MODE" variables in Argo mat files
%	This is used in identifying whether a cycle is real-time or delayed mode, for BRtransfer purposes.
% 04/02/2018 Added fix 9660 pH code. 9660 pH sensor exposed to abmbient
%            light - not in flow stream. Use IBb to adjust surface part
%            of daylight profile (50 m deep)
% 05/25/18 Add code so that CHL_ADJUSTED is always calculated if the raw
%           data exists. If an in situ DC can't be determined the factory
%           DC is used. - jp
% 05/31/2018 Added spiketest for DOXY, NITRATE, PH_IN_SITU_TOTAL (raw and adjusted)
% 08/09/2018, TM  Added bug fix for pH processing related to erroneous
%       flagging of good data as bad if pH data was processed from msg file
%       while dura was present but incomplete.
% 09/11/2018, TM changed calls to isbadsensor.m in support of adding QF='questionable' to bad sensor list capabilities.
% 01/16/19, TM, replaced phcalc_jp.m with phcalc.m
% 10/04/19, JP Added code to recheck _ADJ fill value flags for NO3 & pH
%              gps rollover bug is causing bad end date & no adj data values
%              oven though QC file exists. reset qc from 1 back to 99
% 10/08/19, JP Changed logical test for flbbmode. Used to be 1 or 0. For float
%               17534 flbbmode = 255
% 10/31/19, JP added code to deal with gps week number rollover bug (float
%       9634 for now) which adjusts the date if is determined to be bad
% 02/11/20, JP added PRES QC, PRES_ADJUSTED & PRES_ADJUSTED_QC as variables
%       to cacth bad pressure values as seen in 7593Hawaii cycle 32
% 04/07/20, JP Updated NO3 & pH QC flag assignment. Using a flag matix now.
%     find the max flag in each row for non fill value samples. This will not
%     work for CHL wherea 5 may be entered in the NPQ zone
%04/10/20, JP improvedQC lagging for pH incase dura diagnostics not there
%    or diagnostic data fails the file size test( i.e. 12884)
%6/11/20 TM - added the ".*1e9" to this line.  It was lacking in previous update from 4/7/20 
%    and causing erroneous flagging for certain cases, ie float 9634 surface samples of cycles 67 and 97.
%    However, this has pointed to the need to re-evaluate this code, and make it cleaner and more modular.  So, additional 
%    changes to how the flagging is handled may be forthcoming!
%9/30/20 TM - Added exclusion block to pH IK/IB diagnostic flagging for two recent EqPac floats (these floats are presenting railed diagnostics when pH sensor is actually working)
%10/27/20 TM, minor mods to accomadate change to parser (inclusion of sdn in gps vector)
% ************************************************************************

% *** FOR TESTING ***
% MBARI_ID_str = '9095SOOCN';
% update_str = 'all';
% dirs =[];


% ************************************************************************
% SET FORMATS, DEFAULT DIRS, PREDIMENSION STRUCTURE, BOOKKEEPING
% ************************************************************************
fclose all; % CLOSE ANY OPEN FILES
tf_float.status = 0; % Default flag for float processing success 0 = no good
tf_float.new_messages = {};
tf_float.bad_messages = {};

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    %dirs.NO3config = [user_dir,'CAL\'];
    dirs.FVlocal   = [user_dir,'FLOATVIZ\'];
    dirs.FV        = [user_dir,'FLOATVIZ\'];
    %dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    
    dirs.QCadj     = [user_dir,'CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    %dirs.msg       = '\\atlas\chemwebdata\floats\duplicate\';
    %dirs.msg       = 'C:\temp\';
    dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    
    dirs.log = [user_dir,'Processing_logs\'];
    
    bad_sensor_list = parse_bad_sensor_list([dirs.cal,'bad_sensor_list.txt']);
    iM   = find(strcmp('MBARI ID STR',bad_sensor_list.hdr) == 1);
    % CHECK IF SPECIFC FLOAT HAS BAD SENSOR ISSUES
    if ~isempty(bad_sensor_list.list)
        tSENSOR = strcmp(MBARI_ID_str,bad_sensor_list.list(:,iM));
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
    
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

% ************************************************************************
% SET MSG DIRECTORIES FOR SPECIAL CASE FLOATS. SO FAR THIS INCLUDES FLOATS
% WITH NO WMO #'s OR FLOATS WITH DUPLICATE UW ID #'s. THE MBARI ID IS THE
% ONLY ID THAT UNIQUELY IDENTIFIES ALL FLOATS
if strcmpi(MBARI_ID_str,'6966HAWAII') 
    dirs.msg       = '\\atlas\chemwebdata\floats\duplicate\';
end
if strcmpi(MBARI_ID_str,'8501CalCurrent')
    dirs.msg       = '\\atlas\chemwebdata\floats\duplicate\';
end
if strcmpi(MBARI_ID_str,'8514HAWAII')
    dirs.msg       = '\\atlas\chemwebdata\floats\duplicate\f8514_2\';
end
if strcmpi(MBARI_ID_str,'9254SOOCN')
    dirs.msg       = '\\atlas\chemwebdata\floats\duplicate\';
end
% ************************************************************************

% SET DATA FILL VALUES
fv.bio = 99999;
fv.QC  = 99;

% VALID BIO ARGO RANGE CHECKS - RAW DATA [MIN MAX]
RCR.P     = [0 12000];
RCR.S     = [26 38]; % from argo parameter list
RCR.T     = [-2.5 40]; % from argo parameter list
RCR.O     = [-5 550]; % from argo parameter list
RCR.OP    = [10 70]; % optode phase, from argo parameter list
RCR.OT    = [-2.5 40]; % optode temperature, from argo parameter list
RCR.CHL   = [-0.1 150]; % argoBGC QC manual 09July2016
RCR.BB700 = [-0.000025 0.1]; % argoBGC QC manual 09July2016
RCR.BB532 = [-0.000005 0.1]; % argoBGC QC manual 09July2016
RCR.CDOM = [-1000 1000]; % Place Holder
RCR.NO3  = [-15 65];
RCR.PH     = [7.0 8.8];
RCR.PHV  = [-1.2 -0.7]; % Range check on pH volts
RCR.IB   = [-100 100]; % Range check on pH Ib, nano amps
RCR.IK   = [-100 100]; % Range check on pH Ik, nano amps

% VALID BIO ARGO RANGE CHECKS - QC DATA [MIN MAX]
RC.P     = [0 12000];
RC.S     = [26 38]; % from argo parameter list
RC.T     = [-2.5 40]; % from argo parameter list
RC.O     = [-5 550]; % from argo parameter list
RC.OP    = [10 70]; % optode phase, from argo parameter list
RC.OT    = [-2.5 40]; % optode temperature, from argo parameter list
RC.CHL   = [-0.1 50]; % argoBGC QC manual 09July2016
RC.BB700 = [-0.000025 0.1]; % argoBGC QC manual 09July2016
RC.BB532 = [-0.000005 0.1]; % argoBGC QC manual 09July2016
RC.CDOM  = [-1000 1000]; % Place Holder
RC.NO3   = [-10 55];
RC.PH    = [7.3 8.5];
RC.PHV   = [-1.2 -0.7]; % Range check on pH volts
RC.IB    = [-100 100]; % Range check on pH Ib, nano amps
RC.IK    = [-100 100]; % Range check on pH Ik, nano amps

UW_ID_str = regexp(MBARI_ID_str,'^\d+', 'once','match');
crazy_val = 99990; %some ridiculous number that is less than the GDAC fill value of 99999

% ************************************************************************
% LOAD OR BUILD CALIBRATION DATA
% ************************************************************************
float_cal_path = [dirs.cal,'cal',MBARI_ID_str,'.mat'];
if exist(float_cal_path,'file')
    load(float_cal_path);
    disp(' ')
    disp(' ')
    disp(['FLOAT ',MBARI_ID_str, '(',cal.info.WMO_ID,')']);
    disp(['Existing calibration file loaded: ',float_cal_path])
%	if ~isempty(strfind(cal.info.WMO_ID,'NO_WMO'))
	if ~isempty(strfind(cal.info.WMO_ID,'NO_WMO')) || ~isempty(strfind(cal.info.WMO_ID,'NaN')) %TM 8/16/20
		old_WMO = cal.info.WMO_ID; %get old WMO
		delete(float_cal_path)
		disp(['Attempting to retrieve updated WMOnum by deleting and rebuilding cal file ',float_cal_path,'.'])
		cal = get_float_cals(MBARI_ID_str, dirs);
		if ~strcmp(cal.info.WMO_ID,old_WMO) % if different WMO found
			update_str = 'all';
		end
	end
else
    disp(' ')
    disp(' ')
    disp(['FLOAT ',MBARI_ID_str]);
    cal = get_float_cals(MBARI_ID_str, dirs);
end

if isempty(cal)
    disp(['NO CALIBRATION DATA FOR ',MBARI_ID_str])
    return
end

% CHECK FOR WL_OFFSET IN ISUS NITRATE CAL
% IF IT EXISTS AND IS NOT THE DEFAULT (210), COMMENT
if isfield(cal, 'N')
    if isfield(cal.N, 'WL_offset') && cal.N.WL_offset ~= 210
        disp(['NON STANDARD WAVE LENGTH OFFSET IN NITRATE CAL = ',...
            num2str(cal.N.WL_offset)]);
    end
end
clear dirs.config dirs.cal s1

% ************************************************************************
% GET QC ADJUSTMENT DATA
% ************************************************************************
QC = get_QC_adjustments(cal.info.name, dirs);

% ************************************************************************
% LOOK FOR MOST RECENT PROFILE FILE (*.mat) AND ANY MISSING CASTS
% ONLY WANT NEW OR MISSING CASTS IN THE COPY LIST
% YOU MUST DELETE A MAT FILE TO REDO IF IT IS PARTIAL FOR A GIVEN CAST
% ************************************************************************
[last_cast, missing_casts, first_sdn] = get_last_cast(dirs.mat,cal.info.WMO_ID);
%[last_cast, missing_casts] = get_last_cast(dirs.mat,cal.info.WMO_ID);
if isempty(last_cast)
    last_cast = 0; % set low for logic tests later
end

% GET FILE LISTS
mlist = get_msg_list(MBARI_ID_str,dirs,'msg');  % message file list
ilist = get_msg_list(MBARI_ID_str,dirs,'isus'); % isus file list
dlist = get_msg_list(MBARI_ID_str,dirs,'dura'); % isus file list

% CHECK FOR EMTPY MSG DIR - chemwebdata dir is empty - FLOAT HAS NOT SENT
% ANY MSG FILES YET
if isempty(mlist)
    disp([MBARI_ID_str, ' may be deployed but it has not sent any *.msg',...
        ' files yet!']);
    return
end

% GIVE FILE PROCESSING A TIME LAG - ONLY PROCESS FILES > 4 HRS
% SEEMS LIKE SEVERAL INTERATIONS OF INCOMPLETE FILES CAN BE TRANSMITTED
% DURING SURFACING

% CREATE ANONYMOUS FONCTION TO TEST FOR FILE AGE LESS THEN 4 HOURS OLD
% REMOVE THESE FILES FROM LIST
age_test = @(x) x > (now-4/24);
% age_test = @(x) x > (now);

% if strcmp(MBARI_ID_str,'12701SOOCN')==1
%     if ~isempty(mlist) && ~isempty(mlist.reg_list)
%         t1 = cellfun(age_test2, mlist.reg_sdn, 'UniformOutput', 1);
%         mlist.reg_list(t1) = [];
%         mlist.reg_sdn(t1)  = [];
%     end
% else
%     if ~isempty(mlist) && ~isempty(mlist.reg_list)
%         t1 = cellfun(age_test, mlist.reg_sdn, 'UniformOutput', 1);
%         mlist.reg_list(t1) = [];
%         mlist.reg_sdn(t1)  = [];
%     end
% end


if ~isempty(mlist) && ~isempty(mlist.reg_list)
    t1 = cellfun(age_test, mlist.reg_sdn, 'UniformOutput', 1);
    mlist.reg_list(t1) = [];
    mlist.reg_sdn(t1)  = [];
end
    
if ~isempty(mlist) && ~isempty(mlist.alt_list)
    t1 = cellfun(age_test, mlist.alt_sdn, 'UniformOutput', 1);
    mlist.alt_list(t1) = [];
    mlist.alt_sdn(t1)  = [];
end
if ~isempty(ilist) && ~isempty(ilist.reg_list)
    t1 = cellfun(age_test, ilist.reg_sdn, 'UniformOutput', 1);
    ilist.reg_list(t1) = [];
    ilist.reg_sdn(t1)  = [];
end
if ~isempty(ilist) && ~isempty(ilist.alt_list)
    t1 = cellfun(age_test, ilist.alt_sdn, 'UniformOutput', 1);
    ilist.alt_list(t1) = [];
    ilist.alt_sdn(t1)  = [];
end 
if ~isempty(dlist) && ~isempty(dlist.reg_list)
    t1 = cellfun(age_test, dlist.reg_sdn, 'UniformOutput', 1);
    dlist.reg_list(t1) = [];
    dlist.reg_sdn(t1)  = [];
end
if ~isempty(dlist) && ~isempty(dlist.alt_list)
    t1 = cellfun(age_test, dlist.alt_sdn, 'UniformOutput', 1);
    dlist.alt_list(t1) = [];
    dlist.alt_sdn(t1)  = [];
end  


% UPDATE MODE -- ONLY LOOK FOR NEW FILES
if strcmp(update_str, 'update')
    % Build anonymous function to get cast number from msg file name
    % reg exp pulls str between two periods & sscanf converts to #
    get_cast_num = @(str) sscanf(regexp(str,'(?<=\.)\d+(?=\.)','once', ...
        'match'),'%f');
 
    if ~isempty(last_cast) % ONLY KEEP NEW FILES IN LIST
        if ~isempty(mlist) && ~isempty(mlist.reg_list)
            casts = cellfun(get_cast_num, mlist.reg_list);
            mlist.reg_list(casts <= last_cast &  ...
                ~ismember(casts,missing_casts)) =[];
            mlist.reg_sdn(casts <= last_cast &  ...
                ~ismember(casts,missing_casts)) =[];            
        end
        if ~isempty(mlist) && ~isempty(mlist.alt_list)
            casts = cellfun(get_cast_num, mlist.alt_list);
            mlist.alt_list(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];
            mlist.alt_sdn(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];            
        end
       
        if ~isempty(ilist) && ~isempty(ilist.reg_list)
            casts = cellfun(get_cast_num, ilist.reg_list);
            ilist.reg_list(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];
            ilist.reg_sdn(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];            
        end
        if ~isempty(ilist) && ~isempty(ilist.alt_list)
            casts = cellfun(get_cast_num, ilist.alt_list);
            ilist.alt_list(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];
            ilist.alt_sdn(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];            
        end
        if ~isempty(dlist) && ~isempty(dlist.reg_list)
            casts = cellfun(get_cast_num, dlist.reg_list);
            dlist.reg_list(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];
            dlist.reg_sdn(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];            
        end
        if ~isempty(dlist) && ~isempty(dlist.alt_list)
            casts = cellfun(get_cast_num, dlist.alt_list);
            dlist.alt_list(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];
            dlist.alt_sdn(casts <= last_cast & ...
                ~ismember(casts,missing_casts)) =[];            
        end
    end
end

% ************************************************************************
% COPY FLOAT MESSAGE FILES TO LOCAL TEMPORARY DIRECTORY
% Three file groups *.msg, *.isus, *.dura
% ************************************************************************
if isempty(mlist.reg_list) && isempty(mlist.alt_list)
    disp(['No float message files found to process for float ',MBARI_ID_str]);
    disp(['Last processed cast: ',sprintf('%03.0f',last_cast)])
    return
else
    % CHECK FLOAT TYPE
    if regexp(mlist.reg_dir,'\\f\d{3}\d+','once') % APEX
        float_type = 'APEX';
        disp([cal.info.name,' is an APEX float'])
    elseif regexp(mlist.reg_dir,'\\n\d{3}\d+','once') % NAVIS
        float_type = 'NAVIS';
        disp(['Float ',cal.info.name, ' appears to be a NAVIS float.'])
        disp(['Processing should be done with Process_NAVIS_float.m',...
            ' instead'])
        return
    else
        float_type = 'unknown';
        disp('Unknown float type')
        return
    end
    
    disp(['Copying message files to ', dirs.temp, '  .......'])
    
    % COPY *.MSG files
    if ~isempty(ls([dirs.temp,'*.msg']))
        delete([dirs.temp,'*.msg']) % Clear any message files in temp dir
    end
    if ~isempty(mlist.reg_list)
        for i = 1:length(mlist.reg_list)
            copyfile([mlist.reg_dir, mlist.reg_list{1,i}],dirs.temp);
        end
    end
    if ~isempty(mlist.alt_list)
        for i = 1:length(mlist.alt_list)
            copyfile([mlist.alt_dir, mlist.alt_list{1,i}],dirs.temp);
        end
    end
    clear mlist i
    
    % COPY *.ISUS files
    if ~isempty(ls([dirs.temp,'*.isus']))
        delete([dirs.temp,'*.isus']) % Clear any message files in temp dir
    end
    if ~isempty(ilist)
        if ~isempty(ilist.reg_list)
            for i = 1:length(ilist.reg_list)
                copyfile([ilist.reg_dir, ilist.reg_list{1,i}],dirs.temp);
            end
        end
        if ~isempty(ilist.alt_list)
            for i = 1:length(ilist.alt_list)
                copyfile([ilist.alt_dir, ilist.alt_list{1,i}],dirs.temp);
            end
        end
        clear ilist i
    end
end

% COPY *.dura files
if ~isempty(ls([dirs.temp,'*.dura']))
    delete([dirs.temp,'*.dura']) % Clear any message files in temp dir
end
if ~isempty(dlist)
    if ~isempty(dlist.reg_list)
        for i = 1:length(dlist.reg_list)
            copyfile([dlist.reg_dir, dlist.reg_list{1,i}],dirs.temp);
        end
    end
    if ~isempty(dlist.alt_list)
        for i = 1:length(dlist.alt_list)
            copyfile([dlist.alt_dir, dlist.alt_list{1,i}],dirs.temp);
        end
    end
    clear dlist i
end

% ************************************************************************
% CHECK FOR WMO ID & DIR, IF NO WMO ID CREATE TEMPORARY
% IF DIR NOT THERE MAKE IT,
% IF UPDATE = ALL, CLEAR FILES
% ************************************************************************
% CHECK FOR WMO ID
WMO  = cal.info.WMO_ID;
if isempty(WMO) | isnan(WMO) %TM 8/16/20
    disp(['NO WMO# FOUND FOR FLOAT! CREATING TEMPORARY ', ...
        'DATA DIR FOR ', MBARI_ID_str])
    WMO = ['NO_WMO_',cal.info.UW_ID]; % CREATE TEMPORARY WMO NAME
end

% CHECK FOR EXISTING WMO DIR, CREATE IF NOT THERE
if ~exist([dirs.mat,WMO,'\'],'dir')
    status = mkdir([dirs.mat,WMO,'\']);
    if status == 0
        disp(['Directory could not be created at: ', ...
            [dirs.mat,WMO,'\']]);
        tf_float.status = 0;
        return
    end
end

% IF UPDATE STR = ALL, CLEAR MAT FILES IN WMO DIR
if strcmp(update_str,'all')
    file_chk = ls([dirs.mat,WMO,'\*.mat']);
    if ~isempty(file_chk)
        disp(['Updating all files, clearing existing *.mat files from ',...
            dirs.mat,WMO,'\ first!'])
        delete([dirs.mat,WMO,'\*.mat']);
    end
end

% ************************************************************************
% PROCESS APEX MESSAGE FILES FOR GIVEN FLOAT
% ************************************************************************
if strcmp(float_type,'NAVIS')
    disp(['Float ',cal.info.name, 'appears to be a NAVIS float.'])
    disp('moving to next float')
    return
end

% ***********************************************************************
% GET DIR LIST AS STRUCTURE - THIS WILL BE USED TO SEND AN EMAIL OF
% ARRIVING NEW MSG FILES
mdir = dir([dirs.temp,'*.msg']);
[rr,~] = size(mdir);
if rr > 0
    new_msgs = cell(rr,1);
    ct = 0;
    for m_ct = 1:rr
        tmp_cycle = str2double(regexp(mdir(m_ct).name,'\d+(?=\.msg)', ...
                    'once','match')); % profile # from name
        t1 = missing_casts == tmp_cycle;% missing or new? only add new
        if sum(t1) == 0
            ct = ct+1;
            str = sprintf('%s\t%s\t%0.1f', mdir(m_ct).name, ...
                datestr(mdir(m_ct).datenum,'mm/dd/yy HH:MM'), ...
                mdir(m_ct).bytes/1000);
            new_msgs{ct} = str;
            tf_float.status = 1;
        end
    end
    tf_float.new_messages = new_msgs(1:ct); % list of new msgs from float
    %clear mdir rr new_msgs str
end

% ***********************************************************************


msg_list   = ls([dirs.temp,'*.msg']); % get list of file names to process

lr_ind_chk = 0; % indice toggle - find indices once per float
hr_ind_chk = 0; % toggle

master_FLBB = 0; % some floats (i.e.7558) change mode, start 1 , always 1

% GET BAD SENSOR LIST FOR FLOAT IF ANY
BSL = dirs.BSL;

disp(['Processing ARGO float ' cal.info.name, '.........'])

for msg_ct = 1:size(msg_list,1)
    clear LR HR INFO
    msg_file = strtrim(msg_list(msg_ct,:));
    disp(['PROCESSING MSG FILE ', msg_file])
    NO3_file = regexprep(msg_file,'msg','isus');
    pH_file  = regexprep(msg_file,'msg','dura');
    % find block of numbers then look ahead to see if '.msg' follows
    cast_num = regexp(msg_file,'\d+(?=\.msg)','once','match');
    %         if rem(msg_ct,30) ~= 0
    %             fprintf('%s ',cast_num)
    %         else
    %             fprintf('\r\n')
    %         end
    
    d = parse_APEXmsg4ARGO([dirs.temp,msg_file]);
    
    % IF MSG FILE EXIST BUT NO LR DATA REMOVE FROM NEW LIST AND
    % PUT IN BAD LIST
    if isempty(d.lr_d)
        msg_size = size(msg_file,2);
        t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
        if sum(t1>0)
            tf_float.bad_messages = [tf_float.bad_messages; ...
                tf_float.new_messages(t1)];
            tf_float.new_messages(t1) =[];
            tf_float.status = 0;
        end
    end
    
    fileinfos = dir([dirs.temp,msg_file]);
    timestamps = fileinfos.date;
	timediff = now-datenum(timestamps); % was using 'd.sdn', but for floats just coming up from under ice that doesn't make sense.
    %make sure all 3 files are present if float includes all 3 sensors,
    %otherwise, end processing for this cycle (but if msg file is > 20 days old, and still no isus or dura, then process as usual)
%     timediff = 50; % for manual processing override.
    %disp(['MBARI_ID_str ===>  ',MBARI_ID_str])
    if (strcmp(MBARI_ID_str,'6091SOOCN')==1) || ...
            (strcmp(MBARI_ID_str,'0569SOOCN')==1) || ...
            (strcmp(MBARI_ID_str,'7622HAWAII')==1) || ...
            (strcmp(MBARI_ID_str,'18340CalCurrent')==1)%these 4 will never have .isus files
    else
        if (exist([dirs.temp, NO3_file],'file')==0 && isfield(cal,'N') && ...
                timediff<=20) || (exist([dirs.temp, pH_file],'file')==0 && ...
                isfield(cal,'pH') && timediff<=20) 
            disp('*******************************************************')
            disp(['WARNING: .isus OR .dura FILE IS MISSING FOR: ',msg_file])
            disp(['PROCESSING HALTED FOR ',msg_file,' UNTIL ALL FILES ARE PRESENT.'])
            disp('*******************************************************')
            msg_size = size(msg_file,2);
            t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
            if sum(t1>0)
                tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' waiting for .isus/.dura']}; %mark on list as in limbo
            end       
            continue
        end
    end
	
	if exist([dirs.temp,pH_file],'file')==1
	    dura = parse_pHmsg([dirs.temp,pH_file]);
		if dura.EOT~=1 && timediff<=20
			disp('*******************************************************')
            disp('WARNING: .dura FILE IS MISSING <EOT>')
            disp(['PROCESSING HALTED FOR ',msg_file,' UNTIL ALL FILES ARE PRESENT IN FULL.'])
            disp('*******************************************************')
            msg_size = size(msg_file,2);
            t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
            if sum(t1>0)
                tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' waiting for complete .dura file']}; %mark on list as in limbo
            end       
            continue
		end
	end

    % GATHER INFO FOR BUILDING ODV FILE LATER AND JUST TO MAKE LIFE
    % EASIER
    INFO.CpActivationP = d.CpActivationP;
    INFO.FwRev    = d.FwRev;
    INFO.CTDtype  = d.CTDtype;
    INFO.CTDsn    = d.CTDsn;
    INFO.FlbbMode = d.FlbbMode;
    if isfinite(INFO.FlbbMode)  && INFO.FlbbMode >0 
        master_FLBB = 1;
    end
    
    INFO.sdn      = d.sdn;
    
    
    % CHECK FOR BAD GPS TIME 10/30/19 -jp
    % WE ARE GETTING TIME FROM PROFILE TERMINATION TIME NOT GPS FIX TIME!!!
    if ~isempty(first_sdn) && ~isempty(INFO.sdn)  
        if INFO.sdn > first_sdn + 365*20 && dvec(1) == 2099 % 20 yrs from start?
            disp(['GPS time for this profile is > 20 years past start ', ...
                '- gps week day number bug?!!'])
            dvec = datevec(INFO.sdn);
            dvec(1)  = 1999; % per aoml suggestion don't quite understand jump to 2099
            INFO.sdn = datenum(dvec) + 1024*7;
        elseif INFO.sdn < first_sdn % bad gps time fix 10/30/19
            disp('GPS time for this profile is unreasonable - gps week day number bug!!')
            disp(['days since first profile = ',num2str((INFO.sdn - ...
                first_sdn),'%0.0f')]);
            INFO.sdn = INFO.sdn + 1024*7;
        end
    end

    INFO.cast     = d.cast;
    INFO.gps      = d.gps;
    
    INFO.UW_ID  = cal.info.UW_ID;
    INFO.name   = cal.info.name;
    INFO.WMO_ID = cal.info.WMO_ID;
    INFO.float_type = float_type;
    
    INFO.EOT    = d.EOT;
    
    % ****************************************************************
    % DEAL WITH LOW RESOLUTION DATA FIRST
    % ****************************************************************
    if isempty(d.lr_d) || size(d.lr_d,1)<=3 % CHECK FOR DATA; sometimes incomplete msg files will come through with minimal LR data (ie 12633.127.msg). 
        disp(['No low res data in message file for ', ...
            strtrim(msg_list(msg_ct,:))])
        continue
    else % if data also header variable
        % Low resolution data
        lr_rows = size(d.lr_d,1);
        IX      = (lr_rows: -1 : 1)'; % build reverse indices - ARGO WANTS SHALLOW TO DEEP
        lr_d    = d.lr_d(IX,:); % ARGO WANTS SHALLOW TO DEEP
        clear B IX lr_rows
    end
    
    % DO A BUNCH OF ONE TIME TASKS
    % GET VARIABLE INDICES - ONLY NEED TO DO THIS ONCE PER FLOAT
    if lr_ind_chk == 0
        lr_ind_chk = 1;
        
        iP   = find(strcmp('p',      d.lr_hdr) == 1); % CTD P
        iT   = find(strcmp('t',      d.lr_hdr) == 1); % CTD T
        iS   = find(strcmp('s',      d.lr_hdr) == 1); % CTD S
        iTo  = find(strcmp('Topt',   d.lr_hdr) == 1); % optode temp
        iTPh = find(strcmp('TPhase', d.lr_hdr) == 1); % T Phase 4330
        iRPh = find(strcmp('RPhase', d.lr_hdr) == 1); % R Phase 4330
        iBPh = find(strcmp('BPhase', d.lr_hdr) == 1); % B Phase 3830
        ipH  = find(strcmp('pH(V)',  d.lr_hdr) == 1); % pH
        iChl = find(strcmp('FSig',   d.lr_hdr) == 1); % CHL fluor
        iBb  = find(strcmp('BbSig',  d.lr_hdr) == 1); % Backscatter
        iNO3 = find(strcmp('no3',    d.lr_hdr) == 1); % Nitrate
        
        % CDOM PLACEHOLDER FOR NOW....
        iCdm = find(strcmp('Cdm',   d.lr_hdr) == 1); % CDOM, NAVIS
        
        % CHECK AANDERAA OPTODE INDICE TYPE & SET TO COMMON NAME
        if strcmp(cal.O.type,'3830')
            iPhase = iBPh;
        elseif strcmp(cal.O.type,'4330')
            iPhase = iTPh;
        else
            iPhase = [];
        end
        
        % GET FLOATVIZ DATA - REG AND QC & GET HR DATA TOO
        % WILL BE USED TO EXTRACT QF DATA FLAGS LATER
        % DATA IS BEING PULLED FROM SIROCCO
        FVQC_flag = 1;
        FV_data   = get_FloatViz_data(INFO.name);
        FV_QCdata = get_FloatViz_data([INFO.name,'QC']);
        FV_HRdata   = get_FloatViz_data([INFO.name,'_HR']);
        FV_HRQCdata = get_FloatViz_data([INFO.name,'_HRQC']);
        if isempty(FV_QCdata)
            FV_QCdata = FV_data;
            FV_HRQCdata = FV_HRdata;
            FVQC_flag = 0;
        end
    end
    
    %GET SOME ARGO CORE PARAMTERS
    LR.PRES      = lr_d(:,iP);
    
    % MAKE AN ARRAY OF ZEROS FOR FILLING ARRAYS LATER
    fill0  = ones(size(LR.PRES))* 0;
    
    LR.PRES_QC          = fill0 + fv.QC;
    LR.PRES_ADJUSTED    = LR.PRES;
    LR.PRES_ADJUSTED_QC = fill0 + fv.QC;

    LR.PSAL             = lr_d(:,iS);
    LR.PSAL_QC          = fill0 + fv.QC;
    LR.PSAL_ADJUSTED    = LR.PSAL;
    LR.PSAL_ADJUSTED_QC = fill0 + fv.QC;
    
    LR.TEMP             = lr_d(:,iT);
    LR.TEMP_QC          = fill0 + fv.QC;
    LR.TEMP_ADJUSTED    = LR.TEMP;
    LR.TEMP_ADJUSTED_QC = fill0 + fv.QC;
    
    % CHECK FOR BAD PRESS VALUES
    LRQF_P = LR.PRES < RCR.P(1) | LR.PRES > RCR.P(2);
    LR.PRES_QC(LRQF_P)  = 4;  % BAD
    LR.PRES_QC(~LRQF_P & LR.PRES ~= fv.bio) = 1; % GOOD    
    
    LR.PRES_ADJUSTED_QC(LRQF_P)  = 4;  % BAD
    LR.PRES_ADJUSTED_QC(~LRQF_P & LR.PRES_ADJUSTED ~= fv.bio) = 1; % GOOD     
    
    % FIND BAD SALINITY & TEMP QF BECAUSE BAD S PERCOLATES TO
    % O, N and pH, density
    [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'S');
    LRQF_S = LR.PSAL < RCR.S(1) | LR.PSAL > RCR.S(2);
    t_bio  = LR.PSAL ~= fv.bio;
    
    LR.PSAL_QC(LRQF_S)  = 4;  % BAD
    LR.PSAL_QC(~LRQF_S & LR.PSAL ~= fv.bio) = 1; % GOOD
    LR.PSAL_QC(t_bio) = LR.PSAL_QC(t_bio) * ~BSLflag + BSLflag*theflag;
    LR.PSAL_ADJUSTED_QC(LRQF_S)  = 4;  % BAD
    LR.PSAL_ADJUSTED_QC(~LRQF_S & LR.PSAL_ADJUSTED ~= fv.bio) = 1; % GOOD
    LR.PSAL_ADJUSTED_QC(t_bio) = LR.PSAL_ADJUSTED_QC(t_bio) * ...
        ~BSLflag + BSLflag*theflag;

    [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'T');
    LRQF_T  = LR.TEMP < RCR.T(1) | LR.TEMP > RCR.T(2);
    t_bio   = LR.TEMP ~= fv.bio;
    
    LR.TEMP_QC(LRQF_T)  = 4;  % BAD
    LR.TEMP_QC(~LRQF_T & LR.TEMP ~= fv.bio) = 1; % GOOD
    LR.TEMP_QC(t_bio) = LR.TEMP_QC(t_bio) * ~BSLflag + BSLflag*theflag;
    LR.TEMP_ADJUSTED_QC(LRQF_T)  = 4;  % BAD
    LR.TEMP_ADJUSTED_QC(~LRQF_T & LR.TEMP_ADJUSTED ~= fv.bio) = 1; % GOOD
    LR.TEMP_ADJUSTED_QC(t_bio) = LR.TEMP_ADJUSTED_QC(t_bio) * ...
        ~BSLflag + BSLflag*theflag;
    
    % CALCULATE POTENTIAL DENSITY (STILL NEED TO UPDATE TO TEOS 10 !)
    potT   = theta(LR.PRES, LR.TEMP, LR.PSAL,0); % potential temp, p=0
    % USING POTENTIAL DENSITY ( per: Processing Argo OXYGEN data at the
    % DAC level Version 2.0 October 1st 2015, Chap 5)
    lr_den = density(LR.PSAL, potT); % potential density, kg /m^3
    
    % ****************************************************************
    % CALCULATE OXYGEN CONCENTRATION
    % ****************************************************************
    if ~isempty(iPhase) % Bphase or Tphase exists so O2 data should exist
        t_nan = isnan(lr_d(:,iPhase));
        O2_chk = lr_d(:,iPhase) < RCR.OP(1) | ...
            lr_d(:,iPhase) > RCR.OP(2) | ...
            lr_d(:,iTo) < RCR.OT(1) | ...
            lr_d(:,iTo) > RCR.OT(2); % Flag phase, Temp out of range
        if any(O2_chk)
            disp(['Out of range phase or optode T detected for ', ...
                msg_file])
        end
        
        if strcmp(cal.O.type,'3830')
            LR.BPHASE_DOXY            = fill0 + fv.bio; % predim w fill val
            LR.BPHASE_DOXY_QC         = fill0 + fv.QC;
            LR.BPHASE_DOXY(~t_nan)    = lr_d(~t_nan,iBPh);
            LR.BPHASE_DOXY_QC(~t_nan) = fv.QC;
        elseif strcmp(cal.O.type,'4330')
            LR.TPHASE_DOXY            = fill0 + fv.bio;
            LR.TPHASE_DOXY_QC         = fill0 + fv.QC;
            LR.TPHASE_DOXY(~t_nan)    = lr_d(~t_nan,iTPh);
            LR.TPHASE_DOXY_QC(~t_nan) =  fv.QC;
        end
        LR.DOXY                = fill0 + fv.bio;
        LR.DOXY_QC             = fill0 + fv.QC;
        LR.TEMP_DOXY           = fill0 + fv.bio;
        LR.TEMP_DOXY_QC        = fill0 + fv.QC;
        LR.DOXY_ADJUSTED       = fill0 + fv.bio;
        LR.DOXY_ADJUSTED_QC    = fill0 + fv.QC;
        LR.DOXY_ADJUSTED_ERROR = fill0 + fv.bio;
        INFO.DOXY_SCI_CAL_EQU  = 'not applicable';
        INFO.DOXY_SCI_CAL_COEF = 'not applicable';
        INFO.DOXY_SCI_CAL_COM  = 'not applicable';
        INFO.DOXY_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
        
        %[S, T, P, Phase, CalPhase, [O2], O2sol, pO2]
        O2 = calc_O2_4ARGO(LR.PSAL(~t_nan), LR.TEMP(~t_nan), ...
            LR.PRES(~t_nan),lr_d((~t_nan),iPhase), cal.O); % O2 in µmol/L + more
        
        LR.DOXY(~t_nan) = O2(:,6) ./ lr_den(~t_nan) .* 1000; % µmol/kg
        
        tDOXY = abs(LR.DOXY) > crazy_val & ~t_nan; % Unrealistic bad value
        LR.DOXY(tDOXY) = crazy_val; % SET TO crazy bad value
        LR.DOXY_QC(~t_nan)      = 3;
        % JP: QC assignment below not needed - crazy value greater than range limits
		%LR.DOXY_QC(tDOXY) = 4; % set crazy bad value QF to 4 
        LR.TEMP_DOXY(~t_nan)    = lr_d(~t_nan,iTo);
        LR.TEMP_DOXY_QC(~t_nan) = fv.QC;
        % Save O2sol, pO2?
        clear O2
        
        if isfield(QC,'O')
            % !! ONE TIME GAIN CORRECTION ONLY !!
%             LR.DOXY_ADJUSTED(~t_nan)  = LR.DOXY(~t_nan) .* QC.O.steps(1,3);
            QCD = [LR.PRES(~t_nan), LR.TEMP(~t_nan), LR.PSAL(~t_nan), LR.DOXY(~t_nan)];
            LR.DOXY_ADJUSTED(~t_nan) = ...
                apply_QC_corr(QCD, INFO.sdn, QC.O);
            tDOXY_ADJ = abs(LR.DOXY_ADJUSTED) > crazy_val & ~t_nan; % Unrealistic bad value
            LR.DOXY_ADJUSTED(tDOXY_ADJ) = crazy_val; % SET TO crazy bad value
            LR.DOXY_ADJUSTED_QC(~t_nan) = 1; % set=1 9/27/16 vs 2 = probably good
            % JP: QC assignment below not needed - crazy value greater than range limits
			%LR.DOXY_ADJUSTED_QC(tDOXY_ADJ) = 4; %set crazy val QF to bad
            LR.DOXY_ADJUSTED_ERROR(~t_nan) = LR.DOXY_ADJUSTED(~t_nan) * 0.01;
            INFO.DOXY_SCI_CAL_EQU  = 'DOXY_ADJUSTED=DOXY*G';
            INFO.DOXY_SCI_CAL_COEF = ['G=', ...
                num2str(QC.O.steps(3),'%0.4f')];
            
            if isfield(cal.O,'SVUFoilCoef')
                O2_cal_str = 'SVU Foil calibration coeficients were used. ';
            else
                O2_cal_str = 'Polynomial calibration coeficients were used. ';
            end
            if ~isempty(d.air)
                INFO.DOXY_SCI_CAL_COM  = [O2_cal_str,'G determined from ',...
                    'float  measurements in air. See Johnson et al.,2015,', ...
                    'doi:10.1175/JTECH-D-15-0101.1'];
            else
                INFO.DOXY_SCI_CAL_COM  = [O2_cal_str,'G determined by surface' ...
                    ' measurement comparison to World Ocean Atlas 2009.', ...
                    'See Takeshita et al.2013,doi:10.1002/jgrc.20399'];
            end
        end
        
        if ~isempty(d.air) & sum(d.air(:)) ~= 0 % 0 means ice detection on
            zero_fill = d.air(:,1) * 0; % make array of zeros
            O2 = calc_O2_4ARGO(zero_fill, d.air(:,1), zero_fill, ...
                d.air(:,2), cal.O);
            AIR_O2 = O2; % µmol/L still NOT USED YET
            clear O2
        end
        
        % NEW AIR CAL MEASUREMEMTS 4 just below surf and 8 above
        % 0 for phase means ice detection on
        % aircal = [sdn, bladder(?), pres, temperature, Tphase, Rphase]
        if ~isempty(d.aircal) & sum(d.aircal(:,5)) ~= 0
            zero_fill = d.aircal(:,1) * 0; % make array of zeros
            tpress    = d.aircal(:,3) > 0; % find poss press
            tzmin     = lr_d(:,iP) == min(lr_d(:,iP)); % shallowest value
            S0 = zero_fill + nanmean(lr_d(tzmin,iS)); % surface salinity
            O2 = calc_O2_4ARGO(S0.*tpress, d.aircal(:,4), ...
                d.aircal(:,3).*tpress, d.aircal(:,5), cal.O);
            AIRCAL_O2 = [d.aircal(:,1:3),O2]; % µmol/L still NOT USED YET
            clear O2
        end
        
        
        % DO SOME FINAL CHECK ON VALUES, IF BAD SET QF = 4
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'O');
        t_bio = LR.DOXY ~= fv.bio; % Non fill value samples
        %tST   = LR.PSAL_QC == 4 | LR.TEMP_QC == 4; % Bad S or T will affect O2
        tST   = LR.PSAL_QC == 4 | LR.TEMP_QC == 4 | LR.PRES_QC == 4; % Bad S or T will affect O2
        t_chk = LR.DOXY < RCR.O(1)| LR.DOXY > RCR.O(2); % range check
        t_chk = t_chk | O2_chk | tST; % BAD O2 T or phase
        
        t_chk = t_chk & t_bio; % & not MVI either
        
%         % TESTING
%          if INFO.cast == 24
%              BSLflag
%              LR.PSAL_QC
%              tST
%              t_chk
%              pause
%          end
        
        if isfield(LR,'BPHASE_DOXY')
            LR.BPHASE_DOXY_QC(t_bio) = LR.BPHASE_DOXY_QC(t_bio) * ...
                ~BSLflag + BSLflag*theflag;
            LR.BPHASE_DOXY_QC(t_chk) = 4;
        elseif isfield(LR,'TPHASE_DOXY')
            LR.TPHASE_DOXY_QC(t_bio) = LR.TPHASE_DOXY_QC(t_bio) *  ...
                ~BSLflag + BSLflag*theflag;
            LR.TPHASE_DOXY_QC(t_chk) = 4;
        end
        
        % VERY BAD PSAL & TEMP = BAD O2 TOO
        LR.DOXY_QC(t_bio) = LR.DOXY_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        LR.DOXY_QC(t_chk) = 4;
        
        if isfield(QC,'O')
            t_bio = LR.DOXY_ADJUSTED ~= fv.bio;
            %tST   = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4; % Bad S or T will affect O2
            tST   = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4 ...
                | LR.PRES_ADJUSTED_QC == 4;
            t_chk = LR.DOXY_ADJUSTED < RC.O(1)| ...
                LR.DOXY_ADJUSTED > RC.O(2) | O2_chk;
            t_chk = (t_chk | tST) & t_bio;

            LR.DOXY_ADJUSTED_QC(t_bio) = LR.DOXY_ADJUSTED_QC(t_bio) * ...
                ~BSLflag + BSLflag*theflag;
            LR.DOXY_ADJUSTED_QC(t_chk) = 4;
        end
        
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
        % DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
        % (BAD).  BGC_spiketest screens for nans and fill values and pre-identified bad values.
        %
        % RUN TEST ON RAW DOXY
		% however....if Arctic float, do not perform O2 spiketest (very large gradient near surface...does not work well in this region!)
		if strcmp(MBARI_ID_str,'7596ARCTIC') ~=1 && strcmp(MBARI_ID_str,'7564ARCTIC') ~=1
			QCscreen_O = LR.DOXY_QC == 4; % screen for BAD data already assigned.
			[spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[LR.PRES LR.DOXY],'O2',dirs.cal,fv.bio,QCscreen_O);
			if ~isempty(spike_inds)
				LR.DOXY_QC(spike_inds) = quality_flags;
				disp(['LR.DOXY QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
	%         else
	%             disp(['NO LR.DOXY SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
			end
			%
			% RUN TEST ON QC DOXY_ADJUSTED
			QCscreen_Oadj = LR.DOXY_ADJUSTED_QC == 4; % screen for BAD data already assigned.
			[spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[LR.PRES LR.DOXY_ADJUSTED],'O2',dirs.cal,fv.bio,QCscreen_Oadj);
			if ~isempty(spike_inds)
				LR.DOXY_ADJUSTED_QC(spike_inds) = quality_flags;
				disp(['LR.DOXY_ADJUSTED QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
	%         else
	%             disp(['NO LR.DOXY_ADJUSTED SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
			end
			% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
			clear QCscreen_O QCscreenOadj spike_inds quality_flags
		end
    end
    
        % ****************************************************************
    % CALCULATE CHLOROPHYLL CONCENTRATION (µg/L or mg/m^3)
    % ****************************************************************
    if (~isempty(iChl) && master_FLBB ~= 0) || (~isempty(iChl) && (strcmp(MBARI_ID_str,'12542SOOCN')==1)) || (~isempty(iChl) && (strcmp(MBARI_ID_str,'18169SOOCN')==1))
        % 5/11/20, add exception for 12542, flbb dying and sampling turned
        % off by Dana Swift on cycle 117.  But, still need to create BR file
        % fields as this float has an flbb sensor so these variables should
        % be present.  Not sure why master_FLBB variable part of the
        % statement, seems having just ~istempty(iChl) would suffice?  Add
        % special case for this float for now.
        %disp(INFO.FlbbMode)
        % Predim variables with fill values then adjust as needed
        LR.FLUORESCENCE_CHLA    = fill0 + fv.bio;
        LR.FLUORESCENCE_CHLA_QC = fill0 + fv.QC;
        LR.CHLA                 = fill0 + fv.bio;
        LR.CHLA_QC              = fill0 + fv.QC;
        LR.CHLA_ADJUSTED        = fill0 + fv.bio;
        LR.CHLA_ADJUSTED_QC     = fill0 + fv.QC;
        LR.CHLA_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CHLA_SCI_CAL_EQU   = 'not applicable';
        INFO.CHLA_SCI_CAL_COEF  = 'not applicable';
        INFO.CHLA_SCI_CAL_COM   = 'not applicable';
        INFO.CHLA_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

        t_nan = isnan(lr_d(:,iChl)); % NaN's in data if any
        
        LR.FLUORESCENCE_CHLA(~t_nan)    = lr_d(~t_nan,iChl);
        LR.FLUORESCENCE_CHLA_QC(~t_nan) = fv.QC;
        
        if isfield(cal,'CHL') % Sensor could be bad so maybe no cal info
            LR.CHLA(~t_nan) = (lr_d(~t_nan,iChl) - cal.CHL.ChlDC) .* ...
                cal.CHL.ChlScale;
            LR.CHLA_QC(~t_nan) =  3; % 3 do not use w/o adjusting
            
            % ADJUSTED DATA BASED ON ADMT18 CONCENSUS -WILL BE UPDATED
            % WITHIN THE YEAR - jp 12/13/2017
            
            % 1st CHECK FOR IN SITU DARKS & TRY AND GET THEM IF NOT THERE
            if ~isfield(cal.CHL, 'SWDC')
                SWDC = get_CHLdark(MBARI_ID_str, dirs, 5); % 1st 5 good profiles
                if ~isempty(SWDC)
                    cal.CHL.SWDC = SWDC; % add structure to cal file
                    save(float_cal_path,'cal') % update cal file
                    disp(['Cal file updated: ',float_cal_path]);
                end
            end
            
            % FIGURE OUT WHICH DC TO USE BUT ALWAYS CREATE ADJUSTED CHL
            % DATA IF RAW EXISTS
            if isfield(cal.CHL, 'SWDC')
                CHL_DC = cal.CHL.SWDC.DC;
            else
                CHL_DC = cal.CHL.ChlDC; % NO IN SITU DARKS YET OR EVER
            end
            
            
            LR.CHLA_ADJUSTED(~t_nan) = (lr_d(~t_nan,iChl) - ...
                CHL_DC) .* cal.CHL.ChlScale ./ 2;
            LR.CHLA_ADJUSTED_QC(~t_nan) =  2;
            % NPQ NEXT
            NPQ_CHL = LR.CHLA_ADJUSTED;
            NPQ_CHL(t_nan) = NaN; % fill back to NaN
%                 NPQ = get_NPQcorr([INFO.sdn, nanmean(INFO.gps,1)], ...
%                 %INFO.gps now includes sdn, TM 10/27/20
                NPQ = get_NPQcorr(nanmean(INFO.gps,1), ...
                    [lr_d(:,[iP,iT,iS]),NPQ_CHL], dirs);
            NPQ.data(t_nan,2:end) = fv.bio; % nan back to fill
            
            if ~isempty(NPQ.data)
                iXing   = find(strcmp('Xing_MLD',NPQ.hdr) == 1);
                iSPIKE  = find(strcmp('CHLspike',NPQ.hdr) == 1);
                tNPQ = lr_d(:,iP) <= NPQ.XMLDZ;
                LR.CHLA_ADJUSTED(tNPQ & ~t_nan) = ...
                    NPQ.data(tNPQ & ~t_nan,iXing) + ...
                    NPQ.data(tNPQ & ~t_nan,iSPIKE);
                LR.CHLA_ADJUSTED_QC(tNPQ & ~t_nan) =  5;
                
%                 figure(10) % TESTING
%                 plot(NPQ.data(:,2),NPQ.data(:,1),'go-','MarkerFaceColor','g')
%                 hold on
%                 plot(LR.CHLA_ADJUSTED,NPQ.data(:,1),'b*-', ...
%                     NPQ.data(tNPQ,iXing)+NPQ.data(tNPQ,iSPIKE),NPQ.data(tNPQ,1),'ko-', ...
%                     NPQ.data(:,3),NPQ.data(:,1),'r-')
%                 plot(xlim,xlim*0+NPQ.XMLDZ,'y','LineWidth',2)
%                 set(gca,'Ydir','reverse','Ylim',[0 100])
%                 legend('SWDC & *0.5', 'Corrected', 'Xing + spike', ...
%                     '3pt median','NPQ start')
%                 hold off
%                 pause
            end
            LR.CHLA_ADJUSTED_ERROR(~t_nan) = ...
                abs(LR.CHLA_ADJUSTED(~t_nan) * 2);
            
            INFO.CHLA_SCI_CAL_EQU  = ['CHLA_ADJUSTED=CHLA/A, '...
                'NPQ corrected (Xing et al., 2012), spike profile ', ...
                'added back in'];
            INFO.CHLA_SCI_CAL_COEF = 'A=2';
            INFO.CHLA_SCI_CAL_COM  =['A is best estimate ', ...
                'from Roesler et al., 2017, doi: 10.1002/lom3.10185'];
            
            
            
            % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CHL');
            t_bio   = LR.CHLA ~= fv.bio;
            LR.FLUORESCENCE_CHLA_QC(t_bio) = LR.FLUORESCENCE_CHLA_QC(t_bio) ...
                * ~BSLflag + BSLflag*theflag;
            LR.CHLA_QC(t_bio) = LR.CHLA_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            
            t_chk = t_bio & (LR.CHLA < RCR.CHL(1)|LR.CHLA > RCR.CHL(2));
            LR.CHLA_QC(t_chk) = 4;
            LR.FLUORESCENCE_CHLA_QC(t_chk) = 4;
            
            
            if isfield(cal.CHL, 'SWDC')
                t_bio   = LR.CHLA_ADJUSTED ~= fv.bio;
                LR.CHLA_ADJUSTED_QC(t_bio) = LR.CHLA_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                t_chk = t_bio & ...
                    (LR.CHLA_ADJUSTED < RC.CHL(1)|LR.CHLA_ADJUSTED > RC.CHL(2));
                LR.CHLA_ADJUSTED_QC(t_chk) = 4;
            end
        end
    end

    % ****************************************************************
    % CALCULATE PARTICLE BACKSCATTER COEFFICIENT FROM VOLUME
    % SCATTERING FUNCTION (VSF) (m^-1)
    % APEX FLBB
    % ****************************************************************
    if (~isempty(iBb) && master_FLBB ~= 0) || (~isempty(iBb) && (strcmp(MBARI_ID_str,'12542SOOCN')==1)) || (~isempty(iChl) && (strcmp(MBARI_ID_str,'18169SOOCN')==1))
        %     if ~isempty(iBb) && master_FLBB ~= 0
        VSF                          = fill0 + fv.bio;
        BETA_SW                      = fill0 + fv.bio;
        LR.BETA_BACKSCATTERING700    = fill0 + fv.bio;
        LR.BETA_BACKSCATTERING700_QC = fill0 + fv.QC;
        LR.BBP700                    = fill0 + fv.bio;
        LR.BBP700_QC                 = fill0 + fv.QC;
        LR.BBP700_ADJUSTED           = fill0 + fv.bio;
        LR.BBP700_ADJUSTED_QC        = fill0 + fv.QC;
        LR.BBP700_ADJUSTED_ERROR     = fill0 + fv.bio;
        INFO.BBP700_SCI_CAL_EQU   = 'not applicable';
        INFO.BBP700_SCI_CAL_COEF  = 'not applicable';
        INFO.BBP700_SCI_CAL_COM   = 'not applicable';
        INFO.BBP700_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
        
        t_nan = isnan(lr_d(:,iBb)); % NaN's in data if any
        
        if isfield(cal,'BB') % Sensor could be bad so maybe no cal info
            % VOLUME SCATTERING FUNCTION (VSF) (m^-1 sr^-1)
            VSF(~t_nan) = (lr_d(~t_nan, iBb) - cal.BB.BetabDC) ...
                .* cal.BB.BetabScale;
            % NOW CALCULATE PARTICLE BACKSCATTER COEFFICIENT
            %X       = 1.13*2*pi; % (Barnes and Antoine, 2014)
            %X       = 1.17*2*pi; % (email from E. Boss 24 May 2016)
            
            % (Proc. Bio-Argo particle backscattering at the DAC level
            % Version 1.2, July 21th 2016
            X      = 1.097*2*pi; % FLBB,BBP processing doc July 2016
            LAMBDA = 700;
            DELTA  = 0.039;      % depolarization ratio
            THETA  = 142;        % FLBBAP2 , DAC manual says 120
            BETA_SW_ind = find(t_nan == 0);
            if ~isempty(BETA_SW_ind)
                for ct = 1:size(BETA_SW_ind,1)
                    [BETA_SW(BETA_SW_ind(ct)), b90sw, bsw] = ...
                        betasw_ZHH2009(LAMBDA,lr_d(BETA_SW_ind(ct),iT), ...
                        THETA, lr_d(BETA_SW_ind(ct),iS), DELTA);
                end
            end
            %LR.BETA_BACKSCATTERING700 = VSF; % (DOUBLE CHECK SPEC say counts???)
            LR.BETA_BACKSCATTERING700(~t_nan) = lr_d(~t_nan,iBb); % counts
            LR.BETA_BACKSCATTERING700_QC(~t_nan) = fv.QC; % counts
            LR.BBP700(~t_nan)    = (VSF(~t_nan) - BETA_SW(~t_nan)) * X; %b_bp m^-1
            LR.BBP700_QC(~t_nan) = 3; % 3 do not use w/o adjusting
            
            % CALCULATE ADJUSTED DATA _ NO ADJUSTMENTS AT THIS TIME
            if isfield(QC,'BB')
                QCD = [LR.PRES, LR.TEMP, LR.PSAL, LR.BBP700];
                LR.BBP700_ADJUSTED(~t_nan) = ...
                    apply_QC_corr(QCD(~t_nan,:), INFO.sdn, QC.BB);
                LR.BBP700_ADJUSTED_QC(~t_nan) =  2;
                LR.BBP700_ADJUSTED_ERROR(~t_nan) = fv.bio; % PLACE HOLDER FOR NOW
                INFO.BBP700_SCI_CAL_EQU  = 'BBP700_ADJUSTED=BBP700*A-B';
                INFO.BBP700_SCI_CAL_COEF = ['A=', ...
                    num2str(QC.BB.steps(3),'%0.4f'),',B=',...
                    num2str(QC.BB.steps(4),'%0.4f')];
                INFO.BBP700_SCI_CAL_COM  =['A and B determined by comparison', ...
                    ' to discrete samples from post deployment calibration',...
                    ' rosette cast'];
            end
        end
        
        % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'BBP');
        t_bio   = LR.BBP700 ~= fv.bio;
        LR.BETA_BACKSCATTERING700_QC(t_bio) = ...
            LR.BETA_BACKSCATTERING700_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        LR.BBP700_QC(t_bio) = LR.BBP700_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        
        
        t_chk = t_bio & ...
            (LR.BBP700 < RCR.BB700(1)| LR.BBP700 > RCR.BB700(2));
        LR.BBP700_QC(t_chk) = 4;
        LR.BETA_BACKSCATTERING700_QC(t_chk) = 4;
        
        t_bio = LR.BBP700_ADJUSTED ~= fv.bio;
        LR.BBP700_ADJUSTED_QC(t_bio) = LR.BBP700_ADJUSTED_QC(t_bio) ...
            * ~BSLflag + BSLflag*theflag;
        t_chk = t_bio & (LR.BBP700_ADJUSTED < RC.BB700(1)| ...
            LR.BBP700_ADJUSTED > RC.BB700(2));
        LR.BBP700_ADJUSTED_QC(t_chk) = 4;
        
        clear BETA_SW X VSF ct b90sw bsw
    end
    
    % ****************************************************************
    % CALCULATE CDOM (ppb)
    % ****************************************************************
    if ~isempty(iCdm)
        LR.FLUORESCENCE_CDOM    = fill0 + fv.bio;
        LR.FLUORESCENCE_CDOM_QC = fill0 + fv.QC;
        LR.CDOM                 = fill0 + fv.bio;
        LR.CDOM_QC              = fill0 + fv.QC;
        LR.CDOM_ADJUSTED        = fill0 + fv.bio;
        LR.CDOM_ADJUSTED_QC     = fill0 + fv.QC;
        LR.CDOM_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CDOM_SCI_CAL_EQU   = 'not applicable';
        INFO.CDOM_SCI_CAL_COEF  = 'not applicable';
        INFO.CDOM_SCI_CAL_COM   = 'not applicable';
        INFO.CDOM_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
        
        t_nan = isnan(lr_d(:,iCdm)); % NaN's in data if any
        if isfield(cal,'CDOM') % Sensor could be bad so maybe no cal info
            
            LR.FLUORESCENCE_CDOM(~t_nan)    = lr_d(~t_nan,iCdm);
            LR.FLUORESCENCE_CDOM_QC(~t_nan) = fv.QC;
            LR.CDOM(~t_nan)  = (lr_d(~t_nan,iCdm) - cal.CDOM.CDOMDC) ...
                .* cal.CDOM.CDOMScale;
            LR.CDOM_QC(~t_nan) =  3; % 3 do not use w/o adjusting
            
            % NO ADJUSTED DATA YET !!!! THESE ARE PLACE HOLDERS FOR NOW !!!!
            if isfield(QC,'CDOM')
                LR.CDOM_ADJUSTED(~t_nan)        = fv.bio;
                LR.CDOM_ADJUSTED_QC(~t_nan)     = fv.QC;
                LR.CDOM_ADJUSTED_ERROR(~t_nan) =  fv.bio; % PLACE HOLDER FOR NOW
            end
        end
        
        % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CDOM');
        t_bio = LR.CDOM ~= fv.bio;
        LR.FLUORESCENCE_CDOM_QC(t_bio) = FLUORESCENCE_CDOM_QC(t_bio) ...
            * ~BSLflag + BSLflag*theflag;        
        LR.CDOM_QC(t_bio) = LR.CDOM_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        
        t_chk = t_bio & (LR.CDOM < RCR.CDOM(1)|LR.CDOM > RCR.CDOM(2));
        LR.CDOM_QC(t_chk) = 4;
        LR.FLUORESCENCE_CDOM_QC(t_chk) = 4;
        
        t_bio = LR.CDOM_ADJUSTED ~= fv.bio;
        LR.CDOM_ADJUSTED_QC(t_bio) = LR.CDOM_ADJUSTED_QC(t_bio)  ...
            * ~BSLflag + BSLflag*theflag;
        t_chk = t_bio & ...
            (LR.CDOM_ADJUSTED < RC.CDOM(1)|LR.CDOM_ADJUSTED > RC.CDOM(2));
        LR.CDOM_ADJUSTED_QC(t_chk) = 4;
        
    end
    
    % ****************************************************************
    % CALCULATE pH (µmol / kg scale)
    % ****************************************************************
    if ~isempty(ipH)
        LR.VRS_PH                          = fill0 + fv.bio;
        LR.VRS_PH_QC                       = fill0 + fv.QC;
        LR.PH_IN_SITU_FREE                 = fill0 + fv.bio;
        LR.PH_IN_SITU_FREE_QC              = fill0 + fv.QC;
        LR.PH_IN_SITU_TOTAL                = fill0 + fv.bio;
        LR.PH_IN_SITU_TOTAL_QC             = fill0 + fv.QC;
        LR.IB_PH                           = fill0 + fv.bio;
        LR.IB_PH_QC                        = fill0 + fv.QC;
        IK                                 = fill0 + fv.bio;
        LR.PH_IN_SITU_TOTAL_ADJUSTED       = fill0 + fv.bio;
        LR.PH_IN_SITU_TOTAL_ADJUSTED_QC    = fill0 + fv.QC;
        LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR = fill0 + fv.bio;
        INFO.PH_SCI_CAL_EQU  = 'not applicable';
        INFO.PH_SCI_CAL_COEF = 'not applicable';
        INFO.PH_SCI_CAL_COM  = 'not applicable';
        INFO.PH_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
        
        t_nan = isnan(lr_d(:,ipH)); %NaN's?
        
        LR.VRS_PH(~t_nan)      = lr_d(~t_nan,ipH); % I param
        LR.VRS_PH_QC(~t_nan)   = fv.QC;
        if isfield(cal,'pH')
            [phfree,phtot] = phcalc(LR.VRS_PH(~t_nan), ...
                LR.PRES(~t_nan), LR.TEMP(~t_nan), LR.PSAL(~t_nan), ...
                cal.pH.k0, cal.pH.k2, cal.pH.pcoefs);
            LR.PH_IN_SITU_FREE(~t_nan)    = phfree; % I param
            LR.PH_IN_SITU_FREE_QC(~t_nan) = fv.QC;
            
            LR.PH_IN_SITU_TOTAL(~t_nan)    = phtot;
            LR.PH_IN_SITU_TOTAL_QC(~t_nan) = 3;
            
            LR_inf = isinf(LR.PH_IN_SITU_FREE); % happens if S = 0
            LR.PH_IN_SITU_FREE(LR_inf)    = 20.1; %UNREAL #
            LR.PH_IN_SITU_FREE_QC(LR_inf) = 4;
            LR.PH_IN_SITU_TOTAL(LR_inf)    = 20.1; %UNREAL #
            LR.PH_IN_SITU_TOTAL_QC(LR_inf) = 4;
        else % data returned  but no cal
            disp(['pH data exists but no cal file for ',UW_ID_str])
        end
        
        % GET  DATA FROM *.dura SO I CAN GRAB Ib & Ik
        dura = parse_pHmsg([dirs.temp,pH_file]);
        myINFO = dir([dirs.temp,pH_file]);
        if ~isempty(myINFO)
            ph_filesize = myINFO.bytes;
        else 
            ph_filesize = 0;
        end
        if (isempty(dura.data) || ph_filesize< 10000) % usually > 10 Kb for full profile
            disp(['Not enough pH data returned for ',pH_file])
        else
            pH_data = dura.data;
            pH_hdr  = dura.hdr;
            % GET some pH column indices
            ipH_p  = find(strcmp('CTD Pres', pH_hdr) == 1);
            ipH_t  = find(strcmp('CTD Temp', pH_hdr) == 1);
            ipH_s  = find(strcmp('CTD Sal', pH_hdr) == 1);
            ipH_Ib = find(strcmp('Ib', pH_hdr) == 1); % Base current
            ipH_Ik = find(strcmp('Ik', pH_hdr) == 1); % Base current
            
            ph_rows = size(pH_data,1);
            IX      = (ph_rows:-1:1)';
            %[B,IX]  = sort(pH_data(:,ipH_p));
            pH_data = pH_data(IX,:);
            clear B IX ph_rows
            
            % Mismatch in sample lengths
            if size(LR.PRES,1) ~= size(pH_data(:,1),1)
                rP       = size(LR.PRES,1); % # of CTD samples (*.msg)
                [rN, cN] = size(pH_data);       % size pH data
                disp(['*.dura(',num2str(rN),') & *.msg(',num2str(rP), ...
                    ') sample counts not equal for ',msg_file])
                
                pHtmp    = ones(rP,cN)* fv.bio; % start with fill values
                % fill msg file data then replace with pH data
                pHtmp(:,[ipH_p,ipH_t,ipH_s]) = [LR.PRES, LR.TEMP, LR.PSAL];
                
                for i = 1: rP
                    ind  = find(pH_data(:,ipH_p) == LR.PRES(i),1,'first');
                    if ~isempty(ind)
                        pHtmp(i,:)  = pH_data(ind,:);
                    end
                end
                pH_data = pHtmp;
                clear rP rN cN NO3tmp i pHtmp
            end
            LR.IB_PH(~t_nan)    = pH_data(~t_nan,ipH_Ib) * 1e9; % nano amps
            LR.IB_PH_QC(~t_nan) = fv.QC;
            IK(~t_nan) = pH_data(~t_nan, ipH_Ik) * 1e9; % nano amps, local not for ARGO
            clear ipH_p ipH_t ipH_t ipH_Ib pH_data pH_hdr
        end
        
        if isfield(cal,'pH') && isfield(QC,'pH')
            QCD = [LR.PRES(~t_nan), LR.TEMP(~t_nan), LR.PSAL(~t_nan), phtot];
            LR.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan) = ...
                apply_QC_corr(QCD, INFO.sdn, QC.pH);
            LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_nan) = 1; % set=1 9/27/16 vs
            LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR(~t_nan) = 0.01;
            
            LR.PH_IN_SITU_TOTAL_ADJUSTED(LR_inf)    = 20.1; %UNREAL #
            LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR_inf) = 4;
            
            INFO.PH_SCI_CAL_EQU  = ['PH_ADJUSTED=PH+[PUMP_OFFSET', ...
                '-SUM(OFFSET(S)+DRIFT(S))]* TCOR'];
            INFO.PH_SCI_CAL_COEF = ['PUMP_OFFSET=Pump may add ', ...
                'interfernece in CP mode,OFFSET(S) and DRIFT(S) from ',...
                'climatology comparisons at 1000m or 1500m,','TCOR=', ...
                '(TREF+273.15)./(T+273.15).  TREF = TEMP at 1500m.'];
            INFO.PH_SCI_CAL_COM  =['Contact Tanya Maurer ',...
                '(tmaurer@mbari.org) or Josh Plant (jplant@mbari.org) ',...
                'for more information'];
            
            % TEMPORARY ADJUSTED pH FIX 08/02/2016
            % FLOATVIZ pH CALCULATED WITH OLDER FUNCTION. QC STEPS
            % DETERMINED WITH OLD pH VALUES, BUT A CONSTANT OFFSET
            % JP pH - FV pH = 0.0167)
            % Commented out 9/28/16 doing Qc on JP files now
            %                 disp(['!!!! APPLYING TEMPORARY pH CORRECTION TO ADJUSTED',...
            %                     ' VALUES (adj_pH = adj_pH - 0.0167)']);
            %                 LR.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan) = ...
            %                     LR.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan) - 0.0167;
            
            % FIX 9660 04/02/18 -jp
            if strcmpi(MBARI_ID_str, '9660SOOCN')
                % pH sensor is not in pumped flow stream. Sensor is exposed to
                % ambient light. Use Ib to correct for light sensitivity
                %d(pH)/d(iB), daylight, <50m, Model I R*R > 0.8, = 1.488e6 
                gps9660 = mean(INFO.gps,1);
                [~,El] = SolarAzElq(gps9660(1), gps9660(3),gps9660(2), 0);% sdn lat lon al
                if El > 0
                    dz1 = LR.PRES <= 50;
                    dz2 = LR.PRES <= 120 & LR.PRES >= 70;
                    IBo = median(LR.IB_PH(dz2),'omitnan');
                    LR.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan & dz1) = ...
                        LR.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan & dz1) + ...
                        (LR.IB_PH(~t_nan & dz1) - IBo)./1e9 * 1.488e6;
                    LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_nan & dz1) = 2;
                    clear dz1 dz2 IBo gps9660
                end
            end    
        end
        
%         % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
%         % RAW DATA
%         %tST     = LR.PSAL_QC == 4 | LR.TEMP_QC == 4; % Bad S or T will affect pH
%         tST     = LR.PSAL_QC == 4 | LR.TEMP_QC == 4 | LR.PRES_QC == 4; 
%         t_bio   = LR.PH_IN_SITU_TOTAL ~= fv.bio;        
%         
%         t_chk1 = t_bio & (LR.PH_IN_SITU_TOTAL < RCR.PH(1)| ...
%                  LR.PH_IN_SITU_TOTAL > RCR.PH(2) | tST); % RANGE
%         t_chk2 = t_bio & LR.IB_PH ~= fv.bio.*1e9 & LR.IB_PH ~= fv.bio ...
%             & (LR.IB_PH < RCR.IB(1) | LR.IB_PH > RCR.IB(2) | ...
%             IK < RCR.IK(1) | IK > RCR.IK(2)); % DIAGNOSTIC
%         %t_chk3 = t_bio & (LRQF_S | LRQF_T);
%         t_chk4 = LR.VRS_PH ~= fv.bio &  (LR.VRS_PH < RCR.PHV(1) | ...
%             LR.VRS_PH > RCR.PHV(2));
%         t_chk5 = t_bio & LR.IB_PH ~= fv.bio.*1e9 & LR.IB_PH ~= fv.bio  ...
%             & (LR.IB_PH < RCR.IB(1).*25 | LR.IB_PH > RCR.IB(2).*25 | ...
%             IK < RCR.IK(1).*25 | IK > RCR.IK(2).*25); % DIAGNOSTIC
% 
%         %LR.PH_IN_SITU_TOTAL_QC(t_chk1 | t_chk2 | t_chk3) = 4;
%         LR.PH_IN_SITU_TOTAL_QC(t_chk2) = 3; % set pH with out of range Ik/Ib to questionable, unless pH is also out of range (then set to bad)
%         LR.PH_IN_SITU_TOTAL_QC(t_chk5) = 4;
% 		LR.PH_IN_SITU_TOTAL_QC(t_chk1) = 4;
%         LR.VRS_PH_QC(t_chk4) = 4;
%         LR.IB_PH_QC(t_chk2)  = 3;
%         LR.IB_PH_QC(t_chk5)  = 4;        
%         
% 		[BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'PH');
%         LR.VRS_PH_QC(t_bio) = LR.VRS_PH_QC(t_bio) * ~BSLflag + BSLflag*theflag;
%         LR.IB_PH_QC(t_bio)  = LR.IB_PH_QC(t_bio) * ~BSLflag + BSLflag*theflag;
%         LR.PH_IN_SITU_FREE_QC(t_bio) = LR.PH_IN_SITU_FREE_QC(t_bio) ...
%             * ~BSLflag + BSLflag*theflag;
%         LR.PH_IN_SITU_TOTAL_QC(t_bio) = LR.PH_IN_SITU_TOTAL_QC(t_bio) ...
%             * ~BSLflag + BSLflag*theflag;
        
       
        % DO A FINAL QC CHECK
        % 04/07/2020 JP - Now that the BSL flag can be 3 or 4 it has the
        % potential to overwrite a 4 with a 3 so don't let this happen
        t_bio = LR.PH_IN_SITU_TOTAL ~= fv.bio;
        
        %6/11/20 TM need this incase no diagnostic data, 0 if no data.  
        % added the ".*1e9" to this line.  It was lacking in Josh's update 
        % from 4/7/20 and causing erroneous flagging for certain cases, ie 
        % float 9634 surface samples of cycles 67 and 97.
        %t_diag = LR.IB_PH ~= fv.bio.*1e9; 
        
        % JP 12/17/20 slight modification to TM's fix. If the dura file is
        % too small no Ib & Ik tests happen & LR.IB_PH, LR.IK_PH is never
        % reassigned so you need to check for fv.bio*1e9 & fv.bio as fill
        % values
        t_diag = ~(LR.IB_PH == fv.bio.*1e9 | LR.IB_PH == fv.bio); % JP 12/17/20
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'PH');
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (LR.PSAL_QC == 4 | LR.TEMP_QC == 4 | LR.PRES_QC == 4)*4; % Bad STP will affect nitrate
        tRC   = (LR.PH_IN_SITU_TOTAL < RCR.PH(1)|LR.PH_IN_SITU_TOTAL > RCR.PH(2))*4; % Range check
        tVRS  = (LR.VRS_PH < RCR.PHV(1) | LR.VRS_PH > RCR.PHV(2))*4;
        tIB3  = (LR.IB_PH < RCR.IB(1) | LR.IB_PH > RCR.IB(2)).*t_diag*3; % DIAGNOSTIC
        tIK3  = (IK < RCR.IK(1) | IK > RCR.IK(2)).*t_diag*3; % DIAGNOSTIC
        tIB4  = (LR.IB_PH < RCR.IB(1)*25 | LR.IB_PH > RCR.IB(2)*25).*t_diag*4; % DIAGNOSTIC
        tIK4  = (IK < RCR.IK(1)*25 | IK > RCR.IK(2)*25).*t_diag*4; % DIAGNOSTIC
		
		% 09/30/20 TM - EXCLUSION BLOCK FOR NEWER EQPAC FLOATS EXHIBITING THE ERRONEOUS RAILED DIAG VALUES
		if strcmp(MBARI_ID_str,'17534EqPacE') == 1 || strcmp(MBARI_ID_str,'18601EqPacE')==1 || strcmp(MBARI_ID_str,'18114EqPacE')==1 || strcmp(MBARI_ID_str,'17350EqPacE')==1
		%keep syntax the same to ensure proper function!  Set flags in this block to 0 (always good)
			disp('Excluding float from pH sensor diagnostic checks (erroneous railed values!)!')
		    tIB3  = (LR.IB_PH < RCR.IB(1) | LR.IB_PH > RCR.IB(2)).*t_diag*0; % DIAGNOSTIC
			tIK3  = (IK < RCR.IK(1) | IK > RCR.IK(2)).*t_diag*0; % DIAGNOSTIC
			tIB4  = (LR.IB_PH < RCR.IB(1)*25 | LR.IB_PH > RCR.IB(2)*25).*t_diag*0; % DIAGNOSTIC
			tIK4  = (IK < RCR.IK(1)*25 | IK > RCR.IK(2)*25).*t_diag*0; % DIAGNOSTIC
		end
		   
        tALL  = max([LR.PH_IN_SITU_TOTAL_QC, tBSL, tSTP, tRC, ...
            tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        LR.PH_IN_SITU_TOTAL_QC(t_bio) = tALL(t_bio);
        LR.PH_IN_SITU_FREE_QC(t_bio)  = tALL(t_bio);
        
        tALL  = max([LR.VRS_PH_QC, tBSL, tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag           
        LR.VRS_PH_QC(t_bio) = tALL(t_bio);

        tALL  = max([LR.IB_PH_QC, tBSL, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag       
        LR.IB_PH_QC(t_bio)  = tALL(t_bio);

       
		
%         % *****   ADJUSTED DATA   *****
%         t_bio = LR.PH_IN_SITU_TOTAL_ADJUSTED ~= fv.bio;
%         %tST     = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4; % Bad S or T will affect pH
%         tST     = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4 |...
%                   LR.PRES_ADJUSTED_QC == 4; 
%         
%         t_chk1 = t_bio & (LR.PH_IN_SITU_TOTAL_ADJUSTED < RC.PH(1)| ...
%             LR.PH_IN_SITU_TOTAL_ADJUSTED > RC.PH(2) | tST); 
%         t_chk2 = t_bio & LR.IB_PH ~= fv.bio.*1e9 & LR.IB_PH ~= fv.bio & ...
%             (LR.IB_PH < RC.IB(1) | LR.IB_PH > RC.IB(2) | ...
%             IK < RC.IK(1) | IK > RC.IK(2));
%         %t_chk3 = t_bio & (LRQF_S | LRQF_T);
%         t_chk5 = t_bio & LR.IB_PH ~= fv.bio.*1e9 & LR.IB_PH ~= fv.bio ...
%             & (LR.IB_PH < RCR.IB(1).*25 | LR.IB_PH > RCR.IB(2).*25 | ...
%             IK < RCR.IK(1).*25 | IK > RCR.IK(2).*25); % DIAGNOSTIC
%         
%         %LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_chk1 | t_chk2 |t_chk3) = 4;
%         LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_chk2) = 3;  % set pH with out of range Ik/Ib to questionable, unless pH is also out of range (then set to bad)
% 		LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_chk5) = 4;
%         LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_chk1) = 4;
% 		
% 		% double check fill value QC flags - 9634 returns fv's from
% 		qc_adjustmet   
% 		% due to bad gps date bug - jp 10/04/19
% 		LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_bio)  = fv.QC;
%  
% 		LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_bio) =  ...
%             LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        
        
        % 04/07/2020 JP - Now that the BSL flag can be 3 or 4 it has the
        % potential to overwrite a 4 with a 3 so don't let this happen
        t_bio = LR.PH_IN_SITU_TOTAL_ADJUSTED ~= fv.bio;  
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4 | ....
                LR.PRES_ADJUSTED_QC == 4)*4; % Bad STP will affect nitrate
        tRC   = (LR.PH_IN_SITU_TOTAL_ADJUSTED < RC.PH(1) | ...
                LR.PH_IN_SITU_TOTAL_ADJUSTED > RC.PH(2))*4; % Range check
        tVRS  = (LR.VRS_PH < RC.PHV(1) | LR.VRS_PH > RC.PHV(2))*4;

        % double check fill value QC flags - 9634 returns fv's from qc_adjustmet
		% due to bad gps date bug - jp 10/04/19
        LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_bio)  = fv.QC;
        
         tALL  = max([LR.PH_IN_SITU_TOTAL_ADJUSTED_QC, tBSL, tSTP, tRC, ...
                 tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_bio) = tALL(t_bio);               
			
					% 09/30/20 TM - EXCLUSION BLOCK FOR NEWER EQPAC FLOATS EXHIBITING THE ERRONEOUS RAILED DIAG VALUES
		if strcmp(MBARI_ID_str,'17534EqPacE') ~= 1 && strcmp(MBARI_ID_str,'18601EqPacE') ~=1 && strcmp(MBARI_ID_str,'18114EqPacE')~=1
        if ~isempty(dura.data) && ph_filesize> 10000 && ... % usually > 10000 bytes for full profile...
                any(LR.IB_PH < RC.IB(1) | LR.IB_PH > RC.IB(2))
            disp([pH_file,' has out of range pH Ib values. pH ', ...
                'QF''s for these samples will be set to bad or questionable'])
        end
        if ~isempty(dura.data) && ph_filesize> 10000 && ... % usually > 10000 bytes for full profile...
                any(IK < RC.IK(1) | IK > RC.IK(2))
            disp([pH_file,' has out of range pH Ik values. pH ', ...
                'QF''s for these samples will be set to bad or questionable'])
        end
		end
        
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
        % FINALLY DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
        % (BAD).  BGC_spiketest screens for nans and fill values.
        %
        % RUN TEST ON RAW PH_IN_SITU_TOTAL
        QCscreen_phT = LR.PH_IN_SITU_TOTAL_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[LR.PRES LR.PH_IN_SITU_TOTAL],'PH',dirs.cal,fv.bio,QCscreen_phT);
        if ~isempty(spike_inds)
            % Could add optional code to compare qf already assigned (ie in range
            % checks above), and qf assigned during spiketest.  Keep
            % whichever is greater.  (At this point, I'm pre-screening the vector using
            % current quality flags.  Screening occurs within Function.  Currently assigning 4 (bad) to spikes, although may
            % be modified in the future).  5/22/18. TM.
            %spikeQF = LR.PH_INSITU_FREE_QC(spike_inds) < quality_flags;
            %LR.PH_INSITU_FREE_QC(spike_inds(spikeQF)) = quality_flags(spikeQF);
            LR.PH_IN_SITU_TOTAL_QC(spike_inds) = quality_flags;
            disp(['LR.PH_IN_SITU_TOTAL QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
%         else
%             disp(['NO LR.PH_IN_SITU_TOTAL SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        %
        % RUN TEST ON QC PH_IN_SITU_TOTAL_ADJUSTED
        QCscreen_phTadj = LR.PH_IN_SITU_TOTAL_ADJUSTED_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[LR.PRES LR.PH_IN_SITU_TOTAL_ADJUSTED],'PH',dirs.cal,fv.bio,QCscreen_phTadj);
        if ~isempty(spike_inds)
            LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(spike_inds) = quality_flags;
            disp(['LR.PH_IN_SITU_TOTAL_ADJUSTED QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
%         else
%             disp(['NO LR.PH_IN_SITU_TOTAL_ADJUSTED SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  

        clear phfree phtot QCD dura t_chk1 t_chk2 t_chk3 QCscreen_phF QCscreen_phT QCscreen_phTadj QCscreen_O QCscreenOadj spike_inds quality_flags

    end
    
    % ****************************************************************
    % CALCULATE NITRATE (µmol / kg scale)
    % DO DEPTH CORRECTION 1st
    % CONVERT TO µmol/kg
    % ****************************************************************
    if ~isempty(iNO3)
        LR.NITRATE                      = fill0 + fv.bio;
        LR.NITRATE_QC                   = fill0 + fv.QC;
        LR.UV_INTENSITY_DARK_NITRATE    = fill0 + fv.bio;
        LR.UV_INTENSITY_DARK_NITRATE_QC = fill0 + fv.QC;
        LR.UV_INTENSITY_NITRATE         = fill0 + fv.bio;
        LR.UV_INTENSITY_NITRATE_QC      = fill0 + fv.QC;
        
        LR.NITRATE_ADJUSTED             = fill0 + fv.bio;
        LR.NITRATE_ADJUSTED_QC          = fill0 + fv.QC;
        LR.NITRATE_ADJUSTED_ERROR       = fill0 + fv.bio;
        INFO.NITRATE_SCI_CAL_EQU  = 'not applicable';
        INFO.NITRATE_SCI_CAL_COEF = 'not applicable';
        INFO.NITRATE_SCI_CAL_COM  = 'not applicable';
        INFO.NITRATE_DATA_MODE = 'not applicable';
        
        % IF ISUS FILE & CAL INFO EXIST, TRY AND PROCESS IT
        if exist([dirs.temp, NO3_file],'file') && isfield(cal,'N')
            spec = parse_NO3msg([dirs.temp,NO3_file]); % return struct
            UV_INTEN = spec.UV_INTEN;
            if ~isempty(UV_INTEN)
                % [SDN, DarkCur, Pres, Temp, Sal, NO3, BL_int,BL_slope,
                %  RMS_ER, Wl~240, ABS~240] !!! NITRATE STILL µmol/L !!!
                
                NO3  = calc_FLOAT_NO3(spec, cal.N, 1); % ESW P corr
                %NO3  = calc_APEX_NO3_JP(spec, cal.N, 0); % NO P corr
                
                IX = (size(NO3,1):-1:1)'; % FLIP SHALLOW TO DEEP FOR ARGO
                %[B,IX]   = sort(NO3(:,3)); % SORT SHALLOW TO DEEP FOR ARGO
                NO3      = NO3(IX,:); % µmol/L
                UV_INTEN = UV_INTEN(IX,:);
                clear B IX
            end
            
        elseif ~isfield(cal,'N')
            disp('No calibration info to process nitrate with')
        else
            disp([dirs.temp, NO3_file,' not found!!'])
        end
        
        if exist('NO3', 'var') == 1
            % MATCH NITRATE DEPTHS TO CTD DEPTHS IF LENGTHS DIFFERENT
            % MISSING SPECTRA  AND NO3 CONC. SET TO NaN
            % could also  replace missing NO3 from *.isus rwith *.msg
            
            if size(LR.PRES,1) ~= size(NO3,1) % Mismatch in sample lengths
                cS       = size(spec.UV_INTEN,2); % cols of WL's
                rP       = size(LR.PRES,1); % # of CTD samples (*.msg)
                [rN, cN] = size(NO3); % RN = # of isus samples (*.isus)
                disp(['*.isus(',num2str(rN),') & *.msg(',num2str(rP), ...
                    ') sample counts not equal for ',msg_file])
                NO3tmp   = ones(rP,cN)* NaN;
                spectmp  = ones(rP,cS)* fv.bio;
                % fill msg file data then replace with NO3 data
                %NO3tmp(:,3:6) = [LR.PRES, LR.TEMP, LR.PSAL lr_d(:,iNO3)];
                NO3tmp(:,3:5) = [LR.PRES, LR.TEMP, LR.PSAL];
                
                for i = 1: rP
                    ind  = find(NO3(:,3) == LR.PRES(i),1,'first');
                    if ~isempty(ind)
                        NO3tmp(i,:)  = NO3(ind,:);
                        spectmp(i,:) = UV_INTEN(ind,:);
                    end
                end
                NO3      = NO3tmp;
                
                UV_INTEN = spectmp; % n x m matrix
                clear rP rN cN NO3tmp i
            end
            
            t_nan  = isnan(NO3(:,6));
            % CHECK FIT ERROR AND BASELINE ABSORBANCE
            tABS08 = NO3(:,11) > 0.8; %ABS240 > 0.8 (QF=3)
            tABS11 = NO3(:,9) > 0.003 | NO3(:,11) > 1.1; % RMS & ABS240 (QF=4)
            
            LR.UV_INTENSITY_DARK_NITRATE    = NO3(:,2);
            LR.UV_INTENSITY_DARK_NITRATE(t_nan) = fv.bio;
            LR.UV_INTENSITY_DARK_NITRATE_QC = fill0 + fv.QC;
            
            LR.UV_INTENSITY_NITRATE    =  UV_INTEN;
            LR.UV_INTENSITY_NITRATE_QC =  UV_INTEN * 0 + fv.QC;
            
            NpotT  = theta(NO3(:,3), NO3(:,4), NO3(:,5),0);
            N_den  = density(NO3(:,5), NpotT);
            %NO3_kg = NO3(:,6) ./ N_den * 1000;
            
            
            % ********************************************************
            % CALCULATE NITRATE IN umol/L FOR PRESSURE = 0
            % (vs in situ p). ISUS calibrated at p = 0 but
            % ISUS MEASURES IN SITU CONCENTRATION SO ISUS NO3 SLIGHTLY
            % GREATER THAN BOTTLE or SURFACE NO3 (p=0) AT DEPTH
            % (more moles NO3 per liter at depth)
            % NO3(p=0) = N03(p=z) *[density(S,T,p=0) /density(S,T,p@z)]
            
            NO3_p0 = NO3(:,6) .* sw_dens(NO3(:,5),NpotT,0) ./ ...
                sw_dens(NO3(:,5),NO3(:,4),NO3(:,3));
            
            NO3_p0_kg         = NO3_p0 ./ N_den * 1000; % µmol/kg
            LR.NITRATE        = NO3_p0_kg;
            LR.NITRATE(t_nan) = fv.bio;
            
            % QC RAW
            LR.NITRATE_QC = fill0 + 3; % questionable to start
            LR.NITRATE_QC(t_nan) = fv.QC;
            LR.NITRATE_QC(~t_nan & tABS11) = 4;
            
            % ********************************************************
            % APPLY QC CORRECTIONS
            % CORRECTIONS DETERMINED ON µmol/L scale so adjust on that
            % scale and then convert
            if isfield(cal,'N') && isfield(QC,'N')
                QCD = [LR.PRES, LR.TEMP, LR.PSAL,NO3_p0_kg];
                LR.NITRATE_ADJUSTED = fill0 + fv.bio;
                LR.NITRATE_ADJUSTED(~t_nan) = ...
                    apply_QC_corr(QCD(~t_nan,:), INFO.sdn, QC.N);
                
                LR.NITRATE_ADJUSTED_QC = fill0 + 1; % set=1 9/27/16
                LR.NITRATE_ADJUSTED_QC(t_nan) = fv.QC; % nan's to fill value
                
                LR.NITRATE_ADJUSTED_QC(~t_nan & tABS08) = 3; %Bad ABS240
                if sum(tABS11) > 0
                    disp(['Nitrate ABS@240 > 1.1 or RMS > 0.003 ',...
                        'detected. Setting QF = 4 for flagged values'])
                    LR.NITRATE_ADJUSTED_QC(~t_nan & tABS11) = 4; % Really bad abs or RMS
                end
                
                LR.NITRATE_ADJUSTED_ERROR = (abs(LR.NITRATE - ...
                    LR.NITRATE_ADJUSTED)) * 0.1 + 0.5;
                LR.NITRATE_ADJUSTED_ERROR(t_nan) = fv.bio;
               
                
                INFO.NITRATE_SCI_CAL_EQU  = ['NITRATE_ADJUSTED=', ...
                    '[NITRATE-SUM(OFFSET(S)+DRIFT(S))]/GAIN'];
                INFO.NITRATE_SCI_CAL_COEF = ['OFFSET(S) and DRIFT(S) ', ...
                    'from climatology comparisons at 1000m or 1500m. GAIN ',...
                    'from surface/deep comparison where surface values ',...
                    'are known'];
                INFO.NITRATE_SCI_CAL_COM  =['Contact Tanya Maurer ',...
                    '(tmaurer@mbari.org) or Josh Plant (jplant@mbari.org) ',...
                    'for more information'];
            end
            clear QCD UV_INTEN
        end
        
%         % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
%         % RAW
%         [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'N');
%         t_bio   = LR.NITRATE ~= fv.bio;
%         %tST     = LR.PSAL_QC == 4 | LR.TEMP_QC == 4; % Bad S or T will affect nitrate
%         tST     = LR.PSAL_QC == 4 | LR.TEMP_QC == 4 | LR.PRES_QC == 4; % Bad STP will affect nitrate
%         LR.NITRATE_QC(t_bio) = LR.NITRATE_QC(t_bio) * ~BSLflag + BSLflag*theflag;
%         LR.UV_INTENSITY_DARK_NITRATE_QC(t_bio) = ...
%             LR.UV_INTENSITY_DARK_NITRATE_QC(t_bio) * ~BSLflag + BSLflag*theflag;
%         
%         t_chk1 = t_bio &(LR.NITRATE < RCR.NO3(1)| LR.NITRATE > RCR.NO3(2)); % RANGE
%         t_chk2 = t_bio & tST; % SALT
%         LR.NITRATE_QC(t_chk1 | t_chk2) = 4;
%         LR.UV_INTENSITY_DARK_NITRATE_QC(t_chk1)  = 4;
        
        % DO A FINAL RANGE CHECK BASED ON BSL & OTHER DEPENDANT PARAMETERS
        % MODIFIED 04/07/2020 BY JP BECAUSE BSL FLAG OF 3 COULD OVERWRITE A 4
        t_bio = LR.NITRATE ~= fv.bio; % non fill values
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'N');
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (LR.PSAL_QC == 4 | LR.TEMP_QC == 4 | LR.PRES_QC == 4)*4; % Bad STP will affect nitrate
        tRC   = (LR.NITRATE < RCR.NO3(1)| LR.NITRATE > RCR.NO3(2))*4; % RAnge check
        
        tALL  = max([LR.NITRATE_QC, tBSL, tSTP, tRC],[],2); % get highest flag
        LR.NITRATE_QC(t_bio) = tALL(t_bio);
        
        tALL  = max([LR.UV_INTENSITY_DARK_NITRATE_QC, tBSL, tRC],[],2); % get highest flag
        LR.UV_INTENSITY_DARK_NITRATE_QC(t_bio) = tALL(t_bio);
        
        
        
        % *****  ADJUSTED  *****
%         t_bio   = LR.NITRATE_ADJUSTED ~= fv.bio;
%         %tST     = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4; % Bad S or T will affect nitrate
%         tST     = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4 | ...
%                   LR.PRES_ADJUSTED_QC == 4;         
%         LR.NITRATE_ADJUSTED_QC(t_bio) = LR.NITRATE_ADJUSTED_QC(t_bio) ...
%             * ~BSLflag + BSLflag*theflag;
%         
%         t_chk1 = t_bio & (LR.NITRATE_ADJUSTED < RC.NO3(1)| ...
%                  LR.NITRATE_ADJUSTED > RC.NO3(2));
%         t_chk2 = t_bio & tST;
%         LR.NITRATE_ADJUSTED_QC(t_chk1 | t_chk2) = 4;
        
        t_bio = LR.NITRATE_ADJUSTED ~= fv.bio;
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4 | ...
                  LR.PRES_ADJUSTED_QC == 4)*4; 
        tRC   = (LR.NITRATE_ADJUSTED < RC.NO3(1)| LR.NITRATE_ADJUSTED > RC.NO3(2))*4;
        
        tALL  = max([LR.NITRATE_ADJUSTED_QC, tBSL, tSTP, tRC],[],2); % get highest flag
        LR.NITRATE_ADJUSTED_QC(t_bio) = tALL(t_bio);
              
        % double check fill value QC flags - 9634 returns fv's from qc_adjustmet
        % due to bad gps date bug - jp 10/04/19
        LR.NITRATE_ADJUSTED_QC(~t_bio)  = fv.QC;
        
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
        % DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
        % (BAD).  BGC_spiketest screens for nans and fill values.
        %
        % RUN TEST ON RAW NITRATE
        QCscreen_N = LR.NITRATE_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[LR.PRES LR.NITRATE],'NO3',dirs.cal,fv.bio,QCscreen_N);
        if ~isempty(spike_inds)
            LR.NITRATE_QC(spike_inds) = quality_flags;
            disp(['LR.NITRATE QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
%         else
%             disp(['NO LR.NITRATE SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        %
        % RUN TEST ON QC NITRATE_ADJUSTED
        QCscreen_Nadj = LR.NITRATE_ADJUSTED_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[LR.PRES LR.NITRATE_ADJUSTED],'NO3',dirs.cal,fv.bio,QCscreen_Nadj);
        if ~isempty(spike_inds)
            LR.NITRATE_ADJUSTED_QC(spike_inds) = quality_flags;
            disp(['LR.NITRATE_ADJUSTED QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
%         else
%             disp(['NO LR.NITRATE_ADJUSTED SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
% -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  
        clear tchk tABS08 tABS11 NO3 QCscreen_N QCscreen_Nadj spike_inds quality_flags
    end
    
    % ********************************************************************
    % INCORPORATE QUALITY FLAGS FROM EXISTING ODV FILES IF THEY ARE GREATER
    % THAN ZERO. BRUTE FORCE LOOK UP. USE PRESURE AND TEMPERATURE TO MAKE
    % MATCH AND KEEP TOLERENCES TIGHT
    %              !!! EVENTUALLY MAKE THIS A FUNCTION !!!
    % ********************************************************************
    % BUILD LIST OF ODV VARIABLES MATCHING ARGO VARIABLES
    
    %     QCvars(1,:) = {'Temperature[°C]'    'TEMP'};
    %     QCvars(2,:) = {'Salinity'           'PSAL'};
    %         QCvars(1,:) = {'Oxygen[µM]'         'DOXY'};
    %         QCvars(2,:) = {'Nitrate[µM]'        'NITRATE'};
    %         QCvars(3,:) = {'Chlorophyll[µg/l]'  'CHLA'};
    %         QCvars(4,:) = {'BackScatter[/m/sr]' 'BBP700'};
    %         QCvars(5,:) = {'CDOM[ppb]'          'CDOM'};
    %         QCvars(6,:) = {'pHinsitu[Total]'    'PH_IN_SITU_TOTAL'};
    
    QCvars(1,:) = {'Temperature[°C]'    'TEMP'}; % NEW FLOATVIZ
    QCvars(2,:) = {'Salinity[pss]'      'PSAL'};
    QCvars(3,:) = {'Oxygen[µmol/kg]'    'DOXY'};
    QCvars(4,:) = {'Nitrate[µmol/kg]'   'NITRATE'};
    QCvars(5,:) = {'Chl_a[mg/m^3]'      'CHLA'};
    QCvars(6,:) = {'b_bp700[1/m]'       'BBP700'};
    QCvars(7,:) = {'pHinsitu[Total]'    'PH_IN_SITU_TOTAL'};
    QCvars(8,:) = {'b_bp532[1/m]'       'BBP532'};
    QCvars(9,:) = {'CDOM[ppb]'          'CDOM'};
    
    % DO UNADJUSTED QF's FIRST
    if ~isempty(FV_data)
        iQF = find(strcmp(FV_data.hdr,'QF') == 1); % QF column indices
        tFV         = FV_data.data(:,2) == INFO.cast;
        FV_cast     = FV_data.data(tFV,:);   % get individual FloatViz cast
        FV_QF_sum   = sum(FV_cast(:,iQF),1); % sum of QF columns
        if sum(FV_QF_sum) > 0             % ANY ODV QF FLAGS GREATER THAN ZERO?
            indQF   = iQF(FV_QF_sum > 0); % ODV QF colums w/ non zero flags
            for QF_ct = 1 : size(indQF,2)
                ind = strcmp(FV_data.hdr{indQF(QF_ct)-1},QCvars(:,1));
                if ~isfield(LR, [QCvars{ind,2},'_QC']) % N col in floatviz but really no data
                    continue
                end
                % FLOAVIZ VAR MATCHES LIST, GET MATCHING QF's
                if sum(ind) > 0 && isfield(LR,QCvars{ind,2})
                    ODV_QF  = [FV_cast(:,6),FV_cast(:,8), ...
                        FV_cast(:,indQF(QF_ct))]; % P, T & QC
                    ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                    ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                    
                    
                    ARGO_QF = [LR.PRES, LR.TEMP, LR.([QCvars{ind,2},'_QC'])];
                    ct = 0;
                    for i = 1:size(ARGO_QF,1) % line by line just in case
                        dP = abs(ODV_QF(:,1) - ARGO_QF(i,1)); % press
                        dT = abs(ODV_QF(:,2) - ARGO_QF(i,2)); % temp
                        min_dP = min(dP); % this should be 0 but def < 1m
                        min_dT = min(dT); % this should be 0 but def < 1m
                        if min_dP > 1 % poor match
                            disp(['Poor QF match for ARGO PRES = ', ...
                                num2str(ARGO_QF(i,1))])
                            continue
                        end
                        t1 = dP == min_dP & dT == min_dT; % there should only be one
                        if ODV_QF(t1,3) > ARGO_QF(i,3) % ODV QF worse than ARGO
                            ARGO_QF(i,3) = ODV_QF(t1,3); % replace argo w/ ODV
                            ct =ct+1;
                        end
                    end
                    LR.([QCvars{ind,2},'_QC']) = ARGO_QF(:,3); %update structure
                    if ct > 0
                        disp([num2str(ct), ' QC flags added from ODV file for ', ...
                            QCvars{ind,2}, '_QC'])
                    end
                end
            end
            clear tFV FV_cast FV_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i dP t1
            clear min_dP min_dT dT
        end
    end
    
    % NOW DO ADJUSTED DATA QF's
    if ~isempty(FV_QCdata) && FVQC_flag == 1
        iQF = find(strcmp(FV_QCdata.hdr,'QF') == 1); % QF column indices
        tFVQC       = FV_QCdata.data(:,2) == INFO.cast;
        FVQC_cast   = FV_QCdata.data(tFVQC,:); % get individual QC FloatViz cast
        FVQC_QF_sum = sum(FVQC_cast(:,iQF),1); % sum of columns
        
        if sum(FVQC_QF_sum) > 0             % ANY ODV QF FLAGS GREATER THAN ZERO?
            indQF   = iQF(FVQC_QF_sum > 0); % ODV QF colums w/ non zero flags
            
            for QF_ct = 1 : size(indQF,2)
                ind = strcmp(FV_QCdata.hdr{indQF(QF_ct)-1},QCvars(:,1));
                ind_ct = find(ind == 1,1);
                
                % FLOAVIZ VAR MATCHES LIST, GET MATCHING QF's
                if sum(ind) > 0 && (isfield(LR,[QCvars{ind,2},'_ADJUSTED']))
                    
                    ODV_QF  = [FVQC_cast(:,6), FVQC_cast(:,8), ...
                        FVQC_cast(:,indQF(QF_ct))];% P&T&QC
                    ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                    ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                    
                    ARGO_QF = [LR.PRES, LR.TEMP, ...
                        LR.([QCvars{ind,2},'_ADJUSTED_QC'])];
                    
                    ct = 0;
                    for i = 1:size(ARGO_QF,1) % line by line just in case
                        dP = abs(ODV_QF(:,1) - ARGO_QF(i,1));
                        dT = abs(ODV_QF(:,2) - ARGO_QF(i,2));
                        min_dP = min(dP); % this should be 0 but def < 1m
                        min_dT = min(dT); % this should be 0 but def < 1m
                        if min_dP > 1 % poor match
                            disp(['Poor QCQF match for ARGO PRES = ', ...
                                num2str(ARGO_QF(i,1))])
                            continue
                        end
                        t1 = dP == min_dP & dT == min_dT; % there should only be one
                        if ODV_QF(t1,3) > ARGO_QF(i,3) % ODV QF worse than ARGO
                            ARGO_QF(i,3) = ODV_QF(t1,3); % replace argo w/ ODV
                            ct =ct+1;
                        end
                    end
                    
                    LR.([QCvars{ind,2},'_ADJUSTED_QC']) = ARGO_QF(:,3);
                    if ct > 0
                        disp([num2str(ct), ' QC flags added from ODV file for ', ...
                            QCvars{ind,2}, '_ADJUSTED_QC'])
                    end
                end
            end
            clear tFV FVQC_cast FVQC_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i
            clear dP t1 min_dP
        end
    end
    
    
    % ****************************************************************
    % NOW DEAL WITH HIGH RESOLUTION (CP) DATA
    % ****************************************************************
    if isempty(d.hr_d) % CHECK FOR DATA
        disp(['No high res data in message file for ', ...
            strtrim(msg_list(msg_ct,:))])
        HR.PRES = [];
		HR.PRES_ADJUSTED = [];
        HR.PSAL = [];
        HR.PSAL_ADJUSTED = [];
        HR.TEMP = [];
        HR.TEMP_ADJUSTED = [];
        HR.NBIN_CTD = [];
    else % if data also header variable
        hr_d = d.hr_d; % Low resolution data
        % GET VARIABLE INDICES - ONLY NEED TO DO THIS ONCE
        if hr_ind_chk == 0
            hr_ind_chk = 1;
            
            iPP   = find(strcmp('p', d.hr_hdr) == 1); % CTD P
            iTT   = find(strcmp('t', d.hr_hdr) == 1); % CTD T
            iSS   = find(strcmp('s', d.hr_hdr) == 1); % CTD S
            iBIN   = find(strcmp('nbin ctd', d.hr_hdr) == 1); % CTD S
        end
        HR.PRES = hr_d(:,iPP);
        HR.PRES_ADJUSTED = HR.PRES;
        HR.PSAL = hr_d(:,iSS);
        HR.PSAL_ADJUSTED = HR.PSAL; % NO CORRECTIONS DONE FOR S OR T
        HR.TEMP = hr_d(:,iTT);
        HR.TEMP_ADJUSTED = HR.TEMP;
        HR.NBIN_CTD = hr_d(:,iBIN);
        
        HR.PSAL_QC = HR.PRES * 0 + fv.QC; % Predimmension QF's
        HR.PSAL_ADJUSTED_QC = HR.PSAL_QC;
        HR.TEMP_QC = HR.PSAL_QC;
        HR.TEMP_ADJUSTED_QC = HR.PSAL_QC;
        HR.PRES_QC = HR.PSAL_QC;
        HR.PRES_ADJUSTED_QC = HR.PSAL_QC;
        
        HRQF_P  = HR.PRES < RC.P(1) | HR.PRES > RC.P(2); % BAD HR S
        t_bio  = HR.PRES ~= fv.bio;
        HR.PRES_QC(HRQF_P)  = 4;  % BAD
        HR.PRES_QC(~HRQF_P & HR.PRES ~= fv.bio) = 1; % GOOD
        HR.PRES_ADJUSTED_QC(HRQF_P)  = 4;  % BAD
        HR.PRES_ADJUSTED_QC(~HRQF_P & HR.PRES_ADJUSTED ~= fv.bio) = 1; % GOOD
        
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'S');
        HRQF_S  = HR.PSAL < RC.S(1) | HR.PSAL > RC.S(2); % BAD HR S
        t_bio  = HR.PSAL ~= fv.bio;
        HR.PSAL_QC(HRQF_S)  = 4;  % BAD
        HR.PSAL_QC(~HRQF_S & HR.PSAL ~= fv.bio) = 1; % GOOD
        HR.PSAL_QC(t_bio) = HR.PSAL_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        HR.PSAL_ADJUSTED_QC(HRQF_S)  = 4;  % BAD
        HR.PSAL_ADJUSTED_QC(~HRQF_S & HR.PSAL_ADJUSTED ~= fv.bio) = 1; % GOOD
        HR.PSAL_ADJUSTED_QC(t_bio) = HR.PSAL_ADJUSTED_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        
        
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'T');
        HRQF_T  = HR.TEMP < RC.T(1) | HR.TEMP > RC.T(2); % BAD HR T
        t_bio   = HR.TEMP ~= fv.bio;
        HR.TEMP_QC(HRQF_T)  = 4;  % BAD
        HR.TEMP_QC(~HRQF_T & HR.TEMP ~= fv.bio) = 1; % GOOD
        HR.TEMP_QC(t_bio) = HR.TEMP_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        HR.TEMP_ADJUSTED_QC(HRQF_T)  = 4;  % BAD
        HR.TEMP_ADJUSTED_QC(~HRQF_T & HR.TEMP_ADJUSTED ~= fv.bio) = 1; % GOOD
        HR.TEMP_ADJUSTED_QC(t_bio) = HR.TEMP_ADJUSTED_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        
        % CHECK FOR ANY HIGH RESOLUTION RAW QUALITY FLAGS TO ADD
        if ~isempty(FV_HRdata)
            
            iQF = find(strcmp(FV_HRdata.hdr,'QF') == 1); % QF indices
            tFV         = FV_HRdata.data(:,2) == INFO.cast;
            FV_cast     = FV_HRdata.data(tFV,:);   % get FV cast
            FV_QF_sum   = sum(FV_cast(:,iQF),1); % sum of QF columns
            if sum(FV_QF_sum) > 0 % ANY ODV QF FLAGS GREATER THAN ZERO?
                indQF   = iQF(FV_QF_sum > 0); % ODV QF colums w/ non zero flags
                for QF_ct = 1 : size(indQF,2)
                    ind = strcmp(FV_HRdata.hdr{indQF(QF_ct)-1}, ...
                        QCvars(:,1));
                    % N col in floatviz but really no data
                    if ~isfield(HR, [QCvars{ind,2},'_QC'])
                        continue
                    end
                    
                    % FLOAVIZ VAR MATCHES LIST, GET MATCHING QF's
                    if sum(ind) > 0 && isfield(HR,QCvars{ind,2})
                        ODV_QF  = [FV_cast(:,6),FV_cast(:,8), ...
                            FV_cast(:,indQF(QF_ct))]; % P, T & QC
                        ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                        ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                        
                        ARGO_QF = [HR.PRES, HR.TEMP, ...
                            HR.([QCvars{ind,2},'_QC'])];
                        ct = 0;
                        for i = 1:size(ARGO_QF,1) % line by line just in case
                            dP = abs(ODV_QF(:,1) - ARGO_QF(i,1)); % press
                            dT = abs(ODV_QF(:,2) - ARGO_QF(i,2)); % temp
                            min_dP = min(dP); % this should be 0 but def < 1m
                            min_dT = min(dT); % this should be 0 but def < 1m
                            if min_dP > 1 % poor match
                                disp(['Poor HR QF match for ARGO PRES = ', ...
                                    num2str(ARGO_QF(i,1))])
                                continue
                            end
                            t1 = dP == min_dP & dT == min_dT; % there should only be one
                            if ODV_QF(t1,3) > ARGO_QF(i,3) % ODV QF worse than ARGO
                                ARGO_QF(i,3) = ODV_QF(t1,3); % replace argo w/ ODV
                                ct =ct+1;
                            end
                        end
                        HR.([QCvars{ind,2},'_QC']) = ARGO_QF(:,3); %update structure
                        if ct > 0
                            disp([num2str(ct), ' QC flags added from HR ODV file for ', ...
                                QCvars{ind,2}, '_QC'])
                        end
                    end
                end
                clear tFV FV_cast FV_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i dP t1
                clear min_dP min_dT dT
            end
        end
        
        % CHECK FOR ANY HIGH RESOLUTION ADJUSTED QUALITY FLAGS TO ADD
        if ~isempty(FV_HRQCdata)
            
            iQF = find(strcmp(FV_HRQCdata.hdr,'QF') == 1); % QF indices
            tFV         = FV_HRQCdata.data(:,2) == INFO.cast;
            FV_cast     = FV_HRQCdata.data(tFV,:);   % get FV cast
            FV_QF_sum   = sum(FV_cast(:,iQF),1); % sum of QF columns
            if sum(FV_QF_sum) > 0 % ANY ODV QF FLAGS GREATER THAN ZERO?
                indQF   = iQF(FV_QF_sum > 0); % ODV QF colums w/ non zero flags
                for QF_ct = 1 : size(indQF,2)
                    ind = strcmp(FV_HRQCdata.hdr{indQF(QF_ct)-1}, ...
                        QCvars(:,1));
                    % N col in floatviz but really no data
                    if ~isfield(HR, [QCvars{ind,2},'_ADJUSTED_QC'])
                        continue
                    end
                    
                    % FLOAVIZ VAR MATCHES LIST, GET MATCHING QF's
                    if sum(ind) > 0 && isfield(HR,[QCvars{ind,2},'_ADJUSTED'])
                        ODV_QF  = [FV_cast(:,6),FV_cast(:,8), ...
                            FV_cast(:,indQF(QF_ct))]; % P, T & QC
                        ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                        ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                        
                        ARGO_QF = [HR.PRES, HR.TEMP, ...
                            HR.([QCvars{ind,2},'_ADJUSTED_QC'])];
                        ct = 0;
                        for i = 1:size(ARGO_QF,1) % line by line just in case
                            dP = abs(ODV_QF(:,1) - ARGO_QF(i,1)); % press
                            dT = abs(ODV_QF(:,2) - ARGO_QF(i,2)); % temp
                            min_dP = min(dP); % this should be 0 but def < 1m
                            min_dT = min(dT); % this should be 0 but def < 1m
                            if min_dP > 1 % poor match
                                disp(['Poor HR QF match for ARGO PRES = ', ...
                                    num2str(ARGO_QF(i,1))])
                                continue
                            end
                            t1 = dP == min_dP & dT == min_dT; % there should only be one
                            if ODV_QF(t1,3) > ARGO_QF(i,3) % ODV QF worse than ARGO
                                ARGO_QF(i,3) = ODV_QF(t1,3); % replace argo w/ ODV
                                ct =ct+1;
                            end
                        end
                        HR.([QCvars{ind,2},'_ADJUSTED_QC']) = ARGO_QF(:,3); %update structure
                        if ct > 0
                            disp([num2str(ct), ' QC flags added from HR ODV file for ', ...
                                QCvars{ind,2}, '_ADJUSTED_QC'])
                        end
                    end
                end
                clear tFV FV_cast FV_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i dP t1
                clear min_dP min_dT dT
            end
        end
    end
    
    %-----------------------------------------------------------
    % Add LR.parameter_DATA_MODE for each parameter.  12/21/17

%     cycdate = INFO.sdn;
    cycdate = datenum(timestamps); %use date when file came in, not internal cycle date

    % OXYGEN
    if isfield(cal,'O') 
    XEMPT_O = find(LR.DOXY ~= 99999,1); % if empty, then cycle has no data
        if isempty(QC) || isempty(XEMPT_O) % no adjustments have been made yet, or all data is 99999
            INFO.DOXY_DATA_MODE = 'R';
        elseif ~isempty(QC) && isfield(QC,'O')
            if cycdate > QC.date
                INFO.DOXY_DATA_MODE = 'A';
            else
                INFO.DOXY_DATA_MODE = 'D';
            end
        end
    end
    % NITRATE
    if isfield(cal,'N') 
        XEMPT_N = find(LR.NITRATE ~= 99999,1); % if empty, then cycle has no data
        if isempty(QC) || isempty(XEMPT_N) % no adjustments have been made yet, or all data is 99999
            INFO.NITRATE_DATA_MODE = 'R';
        elseif ~isempty(QC) && isfield(QC,'N')
            if cycdate > QC.date
                INFO.NITRATE_DATA_MODE = 'A';
            else
                INFO.NITRATE_DATA_MODE = 'D';
            end
        end
    end
    % PH
    if isfield(cal,'pH') 
         XEMPT_PH = find(LR.PH_IN_SITU_TOTAL ~= 99999,1); % if empty, then cycle has no data       
        if isempty(QC) || isempty(XEMPT_PH) % no adjustments have been made yet, or all data is 99999
            INFO.PH_DATA_MODE = 'R';
        elseif ~isempty(QC) && isfield(QC,'pH')
            if cycdate > QC.date
                INFO.PH_DATA_MODE = 'A';
            else
                INFO.PH_DATA_MODE = 'D';
            end
        end
    end
    % BBP700
    if isfield(cal,'BB') 
        INFO.BBP700_DATA_MODE = 'R';
    end
    % CDOM
    if isfield(cal,'CDOM') 
        INFO.CDOM_DATA_MODE = 'R';
    end
    % CHL
    if isfield(cal,'CHL')
%         if isfield(cal.CHL,'SWDC') && isfield(cal.CHL.SWDC,'DC') && sum(LR.CHLA_ADJUSTED<99999)>0  % median dark count has been quantified (and there is data for that profile) --> adjustment has been made
        if sum(LR.CHLA_ADJUSTED<99999)>0  % median dark count has been quantified (and there is data for that profile) --> adjustment has been made
            INFO.CHLA_DATA_MODE = 'A';
        else
            INFO.CHLA_DATA_MODE = 'R';
        end
    end
    %-----------------------------------------------------------

    
    
    % *********************************************************************
    % SAVE THE PROFILE AS WMO_ID#.PROFILE.mat
    % THE *.mat file will contain 3 structures:
    %   LR      for low resolution data
    %   HR      for high resolution data (constant profiling)
    %   info	float info that may be of use for ARGO data stream, also
    %           used to build ODV compatible text file from 8.mat files
    % *********************************************************************
    if exist('AIR_O2', 'var')
        INFO.AIR_O2   = AIR_O2; % [p t s phase cor_phase uM uMsat pO2]
    end
    if exist('AIRCAL_O2', 'var')
        % [sdn bladderPress P p t s phase cor_phase uM uMsat pO2]
        INFO.AIRCAL_O2   = AIRCAL_O2;
    end
    
    %         WMO  = INFO.WMO_ID; % string
    %         if isempty(WMO)
    %             disp(['NO WMO# FOUND FOR FLOAT! CREATING TEMPORARY ', ...
    %                   'DATA DIR FOR ', INFO.UW_ID])
    %             WMO = ['NO_WMO_',INFO.UW_ID]; % CREATE TEMPORARY WMO NAME
    %         end
    
    save_str = [dirs.mat, WMO,'\', WMO,'.', cast_num,'.mat'];
    save(save_str,'LR','HR','INFO');
    if msg_ct == 1
        copyfile(float_cal_path, [dirs.mat, WMO,'\']); % copy over cal file
    end
    
    %         % CHECK FOR EXISTING WMO DIR, CREATE IF NOT THERE
    %         if exist([dirs.mat,WMO,'\'],'dir')
    %             save(save_str,'LR','HR','INFO');
    %             if msg_ct == 1
    %                 copyfile(float_cal_path, [dirs.mat, WMO,'\'])
    %             end
    %             tf_float.status = 1;
    %         else
    %             status = mkdir([dirs.mat,WMO,'\']);
    %             if status
    %                 save(save_str,'LR','HR','INFO');
    %                 if msg_ct == 1
    %                     copyfile(float_cal_path, [dirs.mat, WMO,'\'])
    %                 end
    %                 tf_float.status = 1;
    %             else
    %                 disp(['Directory could not be created at: ', ...
    %                     [dirs.mat,WMO,'\']]);
    %                 tf_float.status = 0;
    %             end
    %         end
end
    %fprintf('\r\n')
    
    % *********************************************************************
    % CLEAN UP
    if ~isempty(ls([dirs.temp,'*.msg']))
        delete([dirs.temp,'*.msg']);
    end
    if ~isempty(ls([dirs.temp,'*.isus']))
        delete([dirs.temp,'*.isus']);
    end
    if ~isempty(ls([dirs.temp,'*.dura']))
        delete([dirs.temp,'*.dura']);
    end
%     tf_float.status = 1;
%end












