function [tf_float,FULLreprocess] = Process_APEX_float(MBARI_ID_str, dirs, update_str, ISDEAD)
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
%   dirs.temp      = path to temporary working dir
%   dirs.FV        = path to FloatViz files served to web
%   dirs.QCadj     = path to QC adjustment list for all floats
%   ISDEAD         = binary, 1 = inactive (0 = active).  carrying along
%                       because we don't need to trigger a reprocess on old floats that are no
%                       longer reporting (but have partials coming in that are divisible by 5,
%                       but no data to process)
%
% OUTPUTS:
%   tf_float =  1 if processing was a success, 0 otherwise
%   FULLreprocess = 1 if every 5th cycle (trigger a full reprocess in  Loop_Argo_float)
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
% 03/28/2017 - added range checks for IB_PH and IK to flag bad pH from bad
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
% 09/11/2018, TM changed calls to isbadsensor.m in support of adding
%             QF='questionable' to bad sensor list capabilities.
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
%6/11/20 TM - added the ".*1e9" to this line.  It was lacking in previous
%             update from 4/7/20 and causing erroneous flagging for certain
%             cases, ie float 9634 surface samples of cycles 67 and 97.
%             However, this has pointed to the need to re-evaluate this
%             code, and make it cleaner and more modular.  So, additional
%             changes to how the flagging is handled may be forthcoming!
%9/30/20 TM - Added exclusion block to pH IK/IB diagnostic flagging for two
%             recent EqPac floats (these floats are presenting railed
%             diagnostics when pH sensor is actually working)
%10/27/20 TM, minor mods to accomadate change to parser (inclusion of sdn in gps vector)
% 12/17/20 JP, minor adjustment to TM 6/11/20 fix. Fill value can be
%              99999*1e9 or 99999(if no dura file)
% 03/11/2021 JP, modified for GOBGC file name change and updates to some
%             other functions. Still could use a thourough clean up when time permits
% 6/3/21 TM, Modifications in DOXY section in support of new SBE83 sensor
% 6/21/21 JP, Modifications for weekday rollover bug to fix values for
%             INFO.sdn & INFO.gps(:,1). Fixed small bug I found for dealing
%             with year = 2099. ua8482 & ua9634 are the floats affected
% 6/22/21 TM, Added code to transfer BBP700 to BBP700_ADJUSTED in real time
%       (same for BBP532).  Also included are slight modifications to the
%       specification of DOXY SCIENTIFIC_CALIB fields.
% 07/13/21 JP, Addded ua18431 to the "bad_chl_filter" string
% 8/25/21 TM, Added code to write pH and NO3 adjustment coeffs to INFO
%           structure (for inclusion in Argo Bfiles)  Also added code to
%           deal with QC for floats with failed optodes. This requires
%           updates to <param>_ADJUSTED_ERROR (to account for any increases
%           in uncertainty from the LIReqn8 without O2)
% 9/28/21 TM, Added code to write VK_PH and IK_PH to the matfiles (for
%           population of Bfiles).  Up until now, only IB_PH was getting
%           propagated to the Bfiles(?).  Navis pH diagnostics are still
%           outstanding.
% 11/15/21 TM, Added ua12540 to the "bad_chl_filter" string
% 12/15/21 TM, Added ua17328 to the "skip_ph_diag" exclusion block.  Data
%           looks ok so far, but only 34 cycles in.  Keep an eye on this one.
% 12/16/21 JP, Added ua17350 to the "skip_ph_diag" exclusion block& 1- QC=3
%           in BSL for non bad profiles
% 1/26/22 TM, Added code to implement OCR RTQC (range checks only from
%           https://doi.org/10.13155/62466)
% 02/24/2022 JP - Major modification to O2 processing chain to include
%            2X02 APEX flots (i.e. ua19298 & ua19843). pull in air O2 from
%            "*.srf" files & loop through O2 sensors using dynamic field
%            names. Changes to OCR processing to accept s001 style OCR config
% 05/17/2022 TM - Added capability to utilize the "bad_sample_list" to
%               replace manual flagging.  The check on sirocco flags is still
%               in placefor now, but can be deprecated in the near future.
% 7/25/2022 EC - Added capabilities for processing the park depth BGC data.
% 07/26/2022 TM, Added assinment of PPOX_DOXY to the internal matlab
%                structures (for use internally at MBARI)
% 12/02/2022 TM, added check on OCR data structure (if empty, continue)
%                triggered by error in 19314 cycle 127
% 03/8/23 TM, Fixed small data-mode assignment discrepancy (affects failed
%           sensors returning all nans)
% 04/19/23 TM, Added code enhancements for inflating parameter-adjusted
%              error, and sci-cal-comment for special case floats. These
%              exceptions are outlined in the file "Define_ArgoSpecs_SPECIALCASES"
% 05/10/23 TM, Amended code for bad-sample-list (was in error for
%              multi-line entries of same float, same cycle, different pres levels)
% 05/21/23 TM, Amended code to process APEX 6 sensor OCR with OCR data in
%              the msg file LR data block
% 05/31/23 JP, Amended code to process CHL 435 for profile & park data
% 06/14/23 TM, Small bug fix to indexing on the telemetry cycle air
%               measurement calculation (not the air-series.  That is ok).
% 08/20/23 EC, Park QC flagging, pulling in BSL to park depth section
%
% 11/16/23 TM, Removed processing lag for cycles when isus or dura file is
%                missing

% 01/04/24, TM Added ua21291 to bad_chl_filter; flbb dropouts and Dana
%           switched FLBBmode to 0 starting cycle 60 
% 01/18/24, TM modified calls to Calc_SBE63_O2, which now requires float-type as an input.
% 02/12/24, TM, new variable added for helping with triggering full
%           reprocess on every fifth cycle.
% 4/2/24,  TM, added iridium position fix to INFO structure storage
% 4/30/24, TM, added INFO.ice_flag (T/F) for ice detection logical
% 5/14/24, TM, added CHLA_FLUORESCENCE variable (per Argo documentation) to LR output
% 6/10/24 TM, fixed some old hard-wired indices that were mucking up the
%           sirocco flag grab (post addition of the ice column)!
%
% 9/3/2024 LG Fixed issue with OCR test floats not having data written into
% ODV files. Edits are made on lines 3081 and 3101. OCR data structure no
% longer contained headers that were used to retreive data from .msg files,
% now it gets header indices from correct location and data can be found.
% ************************************************************************

% ************************************************************************
% ************************************************************************
% *** FOR TESTING ***

% MBARI_ID_str = 'ua18739';
% MBARI_ID_str = 'ua19719';
% MBARI_ID_str = 'ua17350';
%MBARI_ID_str = 'ua19314'; % OCR 5906446	APEX
%MBARI_ID_str = 'ua19191'; % OCR '5906320'; % ua19191
%MBARI_ID_str = 'ua19298'; % 2XO2
%MBARI_ID_str = 'ua19843'; % 2XO2
% MBARI_ID_str = 'ua12733'; % 	12733	5905131	APEX
%
% % % % MBARI_ID_str = 'ua20532'; % EMILY TEST FLBB
% % % MBARI_ID_str = 'ua12363'; % EMILY TEST FLBB
% % % % % %
% % % update_str = 'all';
% % % % update_str = 'update';
% % % dirs =[];

% MBARI_ID_str = 'ua22218';
% update_str   = 'all';
% dirs         = [];

% ************************************************************************
% ************************************************************************
get_existing_QC = 1;
FULLreprocess = 0;  %This gets changed to '1' for every 5th cycle coming in in 'update' mode.  TM 2/12/24
springchicken = 30;

if get_existing_QC == 0
    disp('get_existing_QC = 0!! NO QC flags will be pulled from the existing file')
    str = input('Would you like to continue processing [Y/N]','s');
    if isempty(regexpi(str,'^Y','once'))
        return
    end
end

% ************************************************************************
% SET FORMATS, DEFAULT DIRS, PREDIMENSION STRUCTURE, FILTERS, BOOKKEEPING
% ************************************************************************
fclose all; % CLOSE ANY OPEN FILES
tf_float.status = 0; % Default flag for float processing success 0 = no good
tf_float.new_messages = {};
tf_float.bad_messages = {};


% THIS IS AN EXCEPTION STRING FOR GDF TEST FLOAT ua19719 (5906441). The
% optode TPhase appears mostly OK while optode T & RPhase do not so use CTD
% T instead of optode T & skip Optode T range check
gdf_wierd_O2 = 'ua19719';

% THIS IS AN EXCEPTION FOR FLOAT 19314 and 19191, TEST APEX WITH OCR - OCR DATA IS
% RETURNED ON SEPARATE PRES AXIS (ONLY THIS FLOAT!)
ocr_testflt = 'ua19314|ua19191';

% THIS IS A LIST OF FLOATS WITH BAD NO3 FROM THE START OR NO3 LISTED IN THE
% MSG HEADER BUT NO SENSOR ON BOARD
bad_no3_filter = 'un0691|un0569|ua7622|ua18340';

% THIS IS A LIST OF FLOATS WITH WITH STRONG O2 GRADIENTS AT THE SURFACE
% NO SPIKE TESTS PERFORMED ON THESE FLOATS (CURRENTLY JUST ARCTIC FLOATS)
no_o2_spike = 'un0691|ua7564|';

% THIS IS FOR CHL PROCESSING EXCEPTIONS THAT ARE NOT CAUGHT BY OTHER MEANS
% 5/11/20, add exception for 12542, flbb dying and sampling turned
% off by Dana Swift on cycle 117.  But, still need to create BR file
% fields as this float has an flbb sensor so these variables should
% be present.  Not sure why master_FLBB variable part of the
% statement, seems having just ~istempty(iChl) would suffice?  Add
% special case for this float for now.  Same for 19875?
bad_chl_filter = 'ua12542|ua18169|ua19875|ua18431|ua12540|ua21291';

% THIS IS FOR pH PROCESSING EXCEPTIONS THAT ARE NOT CAUGHT BY OTHER MEANS
% 09/30/20 TM - EXCLUSION BLOCK FOR NEWER EQPAC FLOATS EXHIBITING THE
% ERRONEOUS RAILED DIAG VALUES
skip_ph_diag = 'ua17534|ua18601|ua18114|ua17328|ua17350|ua19018|ua18081|ua18829';

% THIS IS FOR FLOATS WITH FAILED OPTODES THAT ARE BEING QC'D USING
% LI(PH/N)R EQN 8 (T/S/LOC INPUTS ONLY) - TM 8/25/21
% NOTE THAT AT THIS TIME (8/25/21) THE FOLLOWING APEX FLOAT ALSO HAS A
% FAILED OPTODE, BUT PH AND NO3 QC IS NOT YET AFFECTED: 19142 (MAY BE ADDED
% IN FUTURE AT NEXT DMQC ASSESSMENT)
bad_O2_filter = 'ua9274|ua9642|ua9659|ua12551|ua12573|ua18082|ua19327|ua19412|ua19142';

% THIS IS A LIST OF FLOATS WITH WITH RAPIDLY DRIFTING PSAL OR FAILED PSAL
% FOR WHICH THE PSAL FROM THE ARMOR3D PRODUCT WILL BE USED AS PROXY
% format for multiple float entris will be as such:
%psal_proxy_flts = {'un0949',99;'un0999',100}
% psal_proxy_flts = {'ua9018',121;'ua9630',24;'ua9631',79;'ua9652',88;'ua12369',101;'ua12379',24;'ua12380',47;'ua12542',122;'ua12702',1};
psal_proxy_flts = {'ua9630',24; 'ua9018', 121; 'ua9652', 88; 'ua12369', 101; 'ua12379', 49; 'ua12380', 46; 'ua12702',0};

% THIS IS A LIST OF FLOATS WITH BAD PH PUMP-OFFSET FOR WHICH WE'VE QC'D TO
% A DEPTH OF ABOVE 985m (PUMP-ON LEVEL).  THUS DATA DEEPER THAN 980m MUST
% BE MARKED QUESTIONABLE.  ADDITIONALLY ADJUSTED_ERROR ON THESE PROFILES
% WILL BE INFLATED.  NOTE THAT THIS PROTOCOL IS TEMPORARY UNTIL A MORE
% ROBUST APPROACH FOR PUMP-OFFSET CORRECTION IS IMPLEMENTED.
pH_pumpoffset_980


% Load 'special-case' Argo exceptions for error inflation & special comment
Define_ArgoSpecs_SPECIALCASES

yesBSAML = 0;

% ************************************************************************
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    % user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    dirs.FV        = [user_dir,'FLOATVIZ\'];

    dirs.QCadj     = [user_dir,'CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\seaecho.shore.mbari.org\floats\';

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

% PARSE BAD SAMPLE LIST-----------------------------------------------
bad_sample_list = parse_bad_sample_list([dirs.cal,'bad_sample_list.txt']);
iM   = find(strcmp('MBARI ID STR',bad_sample_list.hdr) == 1);
ibsSENS   = find(strcmp('SENSOR',bad_sample_list.hdr) == 1);
ibsCYC   = find(strcmp('CYCLE',bad_sample_list.hdr) == 1);
ibsD   = find(strcmp('DEPTH',bad_sample_list.hdr) == 1);
ibsDB   = find(strcmp('DEPTH BLOCKS',bad_sample_list.hdr) == 1);
ibsFL   = find(strcmp('FLAG',bad_sample_list.hdr) == 1);

% CHECK IF SPECIFC FLOAT HAS BAD SAMPLES NOT CAUGHT BY RT TESTS
if ~isempty(bad_sample_list.list)
    tSENSOR = strcmp(MBARI_ID_str,bad_sample_list.list(:,iM));
    if sum(tSENSOR) > 0
        disp([MBARI_ID_str,' found on the bad samples list!'])
        BSAML = bad_sample_list;
        dirs.BSAML = BSAML;
        dirs.BSAML.list = BSAML.list(tSENSOR,:);
        yesBSAML = 1;
        %             clear BSL SbsIND
    else
        dirs.BSAML.hdr  = [];
        dirs.BSAML.list = [];
        yesBSAML = 0;
    end
    clear tSENSOR
end
% ------------------------------------------------------------------------



% ************************************************************************
% SET DATA FILL VALUES
fv.bio = 99999;
fv.QC  = 99;

% VALID BIO ARGO RANGE CHECKS - RAW DATA [MIN MAX]
RCR.P     = [0 12000];
RCR.S     = [26 38]; % from argo parameter list
RCR.T     = [-2.5 40]; % from argo parameter list
RCR.O     = [-5 550]; % from argo parameter list
RCR.PO2   = [-5 5000]; % from argo parameter list; added to this code 5/25/21, TM: is this range right?  Should be max 500??
RCR.OP    = [10 70]; % optode phase, from argo parameter list (the range is the same for RPhase and TPhase)
RCR.RP    = [0 15]; % RPhase limits should not be the same as TPhase limits - this is just a guess!!
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
RCR.OCR380 = [-1 1.7]; %from  https://doi.org/10.13155/62466 ;
RCR.OCR412 = [-1 2.9];
RCR.OCR443 = [-1 3.2]; 
RCR.OCR490 = [-1 3.4];
RCR.OCRPAR = [-1 4672];

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

INST_ID_str = regexp(MBARI_ID_str,'^\d+', 'once','match');
crazy_val = 99990; %some ridiculous number that is less than the GDAC fill value of 99999

% ************************************************************************
% LOAD OR BUILD CALIBRATION DATA
% ************************************************************************
fp_cal = [dirs.cal,'cal',MBARI_ID_str,'.mat'];
if exist(fp_cal,'file') == 2
    load(fp_cal);
    disp(' ')
    disp(['FLOAT ',MBARI_ID_str, '(',cal.info.WMO_ID,')']);
    disp(['Existing calibration file loaded: ',fp_cal])

    % IF TEMPORARY WMO STILL - TRY & UPDATE WITH ACTUAL WMO

    if strncmp(cal.info.WMO_ID,'NO', 2) && cal.info.tf_bfile ~= 0
        disp(['Attempting to retrieve updated WMO num by deleting ', ...
            'and rebuilding cal file ',fp_cal])
        tmp_WMO = cal.info.WMO_ID; % need this to wipe old files if new WMO found
        delete(fp_cal)
        cal = get_float_cals(MBARI_ID_str, dirs); %
        if regexp(cal.info.WMO_ID, '^\d{7}', 'once') % WMO found - update successful!
            update_str = 'all';

            %             % WIPE OLD FILES (sirocco, chem, ftp & local
            %             disp(['Temporary WMO text & mat files (',tmp_WMO,') for ', ...
            %                 cal.info.WMO_ID,' have been removed']);
            %             % *.MAT FILES FIRST
            %             mat_dirs = {dirs.mat,'\\atlas\chem\ARGO_PROCESSING\DATA\FLOATS\', ...
            %                 '\\atlas\ftp\pub\ARGO_DATA\'};
            %             for dct = 1:size(mat_dirs,2)% step through *.mat profile dirs
            %                 fp_mat = [mat_dirs{dct}, tmp_WMO,'\'];
            %                 if ~isempty(ls([fp_mat,'*.mat']))
            %                     delete([fp_mat,'*.mat']);
            %                     rmdir(fp_mat)
            %                 end
            %             end
            %             % NOW WIPE TEMPORARY TXT FILES
            %             TXT_dirs = {dirs.FV, ...
            %                 '\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATVIZ\', ...
            %                 '\\sirocco\wwwroot\lobo\Data\FloatVizData\'};
            %             TXT_subdirs ={'' 'QC\' 'HR\' 'HRQC\' 'CANYON_QC' 'MLR_QC'};
            %             for dct = 1:size(TXT_dirs,2)% step through *.TXT file dirs
            %                 fp_txt = TXT_dirs{dct};
            %                 for fct = 1:size(TXT_subdirs,2)
            %                     fp_txt2 = [fp_txt, TXT_subdirs{fct}, tmp_WMO,'.TXT'];
            %                     if isfile(fp_txt2)
            %                         delete(fp_txt2);
            %                     end
            %                 end
            %             end
            %             clear tmp_WMO mat_dirs fp_mat dct
            %             clear TXT_dirs TXT_subdirs fct fp_txt fp_txt2
        end
    end
else % No cal file built so build it
    fprintf('\n\n');
    disp(['FLOAT ',MBARI_ID_str]);
    cal = get_float_cals(MBARI_ID_str, dirs);
end

if isempty(cal)
    disp(['NO CALIBRATION DATA FOR ',MBARI_ID_str])
    return
end

% CHECK FOR FLOAT TYPE, IF NOT APEX: EXIT
float_type = cal.info.float_type;
if strcmp(float_type,'APEX') % APEX
    disp([cal.info.name,' is an APEX float'])
elseif strcmp(float_type,'NAVIS') % NAVIS
    disp(['Float ',cal.info.name, ' appears to be a NAVIS float.'])
    disp(['Processing should be done with Process_NAVIS_float.m',...
        ' instead'])
    return
elseif strcmp(float_type,'SOLO') % SOLO
    disp(['Float ',cal.info.name, ' appears to be a SOLO float.'])
    disp(['Processing should be done with Process_SOLO_float.m',...
        ' instead'])
    return
else
    disp('Unknown float type')
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

% ************************************************************************
% GET QC ADJUSTMENT DATA
% ************************************************************************
QC = get_QC_adjustments(cal.info.WMO_ID, dirs);

% ************************************************************************
% GET MSG, ISUS, AND DURA FILE LISTS
% call returns a structure of hdr & list{file name, file dir, & file date}
mlist = get_msg_list(cal.info, 'msg');
ilist = get_msg_list(cal.info, 'isus');
dlist = get_msg_list(cal.info, 'dura');
slist = get_msg_list(cal.info, 'srf');

if cal.info.O2_flag > 1
    fprintf('%0.0f oxygen senors detected for %s - looking for *.srf files', ...
        cal.info.O2_flag, MBARI_ID_str);
    %    slist = get_msg_list(cal.info, 'srf');
end

% CHECK FOR EMTPY MSG DIR - if msgdir is empty - FLOAT HAS NOT SENT
% ANY MSG FILES YET
if isempty(mlist.list)
    disp([MBARI_ID_str, ' may be deployed but it has not sent any *.msg',...
        ' files yet!']);
    return
end

% GIVE FILE PROCESSING A TIME LAG - ONLY PROCESS FILES > 4 HRS
% SEEMS LIKE SEVERAL INTERATIONS OF INCOMPLETE FILES CAN BE TRANSMITTED
% DURING SURFACING

% TEST FOR FILE AGE LESS THEN 4 HOURS OLD
% REMOVE THESE FILES FROM LIST
%age_limit = now - 4/24;
age_limit = now - 2/24;
%age_limit = now; % no lag if you want an imediate process
if ~isempty(mlist.list)
    t1 = cell2mat(mlist.list(:,3)) < age_limit;
    mlist.list = mlist.list(t1,:);
end

if ~isempty(ilist.list)
    t1 = cell2mat(ilist.list(:,3)) < age_limit;
    ilist.list = ilist.list(t1,:);
end

if ~isempty(dlist.list)
    t1 = cell2mat(dlist.list(:,3)) < age_limit;
    dlist.list = dlist.list(t1,:);
end

if cal.info.O2_flag > 1 && ~isempty(slist.list)
    t1 = cell2mat(slist.list(:,3)) < age_limit;
    slist.list = slist.list(t1,:);
end
clear t1
% ************************************************************************
% LOOK FOR MOST RECENT PROFILE FILE (*.mat) AND ANY MISSING CASTS
% ONLY WANT NEW OR MISSING CASTS IN THE COPY LIST
% YOU MUST DELETE A MAT FILE TO REDO IF IT IS PARTIAL FOR A GIVEN CAST
% ************************************************************************
% NOTE:: TM,JP 9/27/22, this call to get_last_cast is somewhat circular!!  If a bug creeps in to the matfiles,
%        and you want to reprocess, this call references the bad matfile sdn prior to deleting the matfiles for reprocess....
% NOTE:: TM, JP 12/3/23, this issue cropped up again during cross-talk of
%          processing runs!  A bad date crept into the matfiles.  Note for future!!:if a bad date is persisting, clear all .mat files manually prior to reprocessing!
[last_cast, missing_casts, first_sdn] = get_last_cast(dirs.mat, ...
    cal.info.WMO_ID);
if isempty(last_cast)
    last_cast = 0; % set low for logic tests later
end
% UPDATE MODE -- ONLY LOOK FOR NEW FILES
if strcmp(update_str, 'update') && last_cast > 0
    if ~isempty(mlist.list)
        tmp   = (regexp(mlist.list(:,1),'(?<=\d+\.)\d{3}(?=\.msg)','match','once'));
        casts = str2double(tmp);
        t1    = casts > last_cast; % new cycles in file list
        t2    = ismember(casts, missing_casts); %loc of missing in file list
        mlist.list = mlist.list(t1|t2,:);
    end
    
    timechk = cell2mat(mlist.list(:,end)); %timestamp of receipt for 'new' files identified
%     springchicken = 30;
    xxtimechk = find(timechk>datenum(now-springchicken)); %find timestamps for 'new' files that are less than 'springchicken' days old
    if ISDEAD && isempty(xxtimechk) % Remember, we are still in UPDATE MODE ONLY!  None of this is relevant for an 'all' full reprocess of a dead float.
        FULLreprocess = 0;
        disp(['Stale files found for dead float ',MBARI_ID_str,'! Continuing on to next float...'])
        return
    end

    if ~isempty(ilist.list)
        tmp   = (regexp(ilist.list(:,1),'(?<=\d+\.)\d{3}(?=\.isus)','match','once'));
        casts = str2double(tmp);
        t1    = casts > last_cast; % new cycles in file list
        t2    = ismember(casts, missing_casts); %loc of missing in file list
        ilist.list = ilist.list(t1|t2,:);
    end

    if ~isempty(dlist.list)
        tmp   = (regexp(dlist.list(:,1),'(?<=\d+\.)\d{3}(?=\.dura)','match','once'));
        casts = str2double(tmp);
        t1    = casts > last_cast; % new cycles in file list
        t2    = ismember(casts, missing_casts); %loc of missing in file list
        dlist.list = dlist.list(t1|t2,:);
    end

    if ~isempty(slist.list) % *.srf files, multile 4330 O2 floats only
        tmp   = (regexp(slist.list(:,1),'(?<=\d+\.)\d{3}(?=\.srf)','match','once'));
        casts = str2double(tmp);
        t1    = casts > last_cast; % new cycles in file list
        t2    = ismember(casts, missing_casts); %loc of missing in file list
        slist.list = slist.list(t1|t2,:);
    end
end
clear tmp t1 t2 casts

% ************************************************************************
% COPY FLOAT MESSAGE FILES TO LOCAL TEMPORARY DIRECTORY
% Three file groups *.msg, *.isus, *.dura
% ************************************************************************
if isempty(mlist.list)
    disp(['No float message files found to process for float ',MBARI_ID_str]);
    disp(['Last processed cast: ',sprintf('%03.0f',last_cast)])
    return
end

disp(['Copying message files to ', dirs.temp, '  .......'])

% COPY *.MSG files
if ~isempty(ls([dirs.temp,'*.msg']))
    delete([dirs.temp,'*.msg']) % Clear any message files in temp dir
end
if ~isempty(mlist.list)
    for i = 1:size(mlist.list,1)
        fp = fullfile(mlist.list{i,2}, mlist.list{i,1});
        copyfile(fp, dirs.temp);
    end
end

% COPY *.ISUS files
if ~isempty(ls([dirs.temp,'*.isus']))
    delete([dirs.temp,'*.isus']) % Clear any message files in temp dir
end
if ~isempty(ilist.list)
    for i = 1:size(ilist.list,1)
        fp = fullfile(ilist.list{i,2}, ilist.list{i,1});
        copyfile(fp, dirs.temp);
    end
end

% COPY *.DURA files
if ~isempty(ls([dirs.temp,'*.dura']))
    delete([dirs.temp,'*.dura']) % Clear any message files in temp dir
end
if ~isempty(dlist.list)
    for i = 1:size(dlist.list,1)
        fp = fullfile(dlist.list{i,2}, dlist.list{i,1});
        copyfile(fp, dirs.temp);
    end
end

% COPY *.SRF files
if cal.info.O2_flag > 1 &&  ~isempty(ls([dirs.temp,'*.srf']))
    delete([dirs.temp,'*.serf']) % Clear any message files in temp dir
end
if cal.info.O2_flag > 1 && ~isempty(slist.list)
    for i = 1:size(slist.list,1)
        fp = fullfile(slist.list{i,2}, slist.list{i,1});
        copyfile(fp, dirs.temp);
    end
end

% ************************************************************************
% CHECK FOR WMO ID & DIR, IF NO WMO ASSIGNED YET NOTIFY
% IF WMO DIR NOT THERE MAKE IT,
% IF UPDATE = ALL, CLEAR FILES
% ************************************************************************
% CHECK FOR WMO ID
WMO  = cal.info.WMO_ID;

if strncmp(WMO,'^NO_WMO',6) && cal.info.tf_bfile == 0
    fprintf('WMO will never be assigned to %s\n', cal.info.name);
elseif strncmp(WMO,'^NO_WMO',6)
    fprintf('No WMO number assigned to %s yet\n', cal.info.name);
end

% if strncmp(WMO,'^NO_WMO',6) && cal.info.tf_bfile == 1
%     fprintf('No WMO number assigned to %s yet\n', cal.info.name);
% elseif strncmp(WMO,'^NO_WMO',6) && cal.info.tf_bfile == 0
%     fprintf('WMO will never be assigned to %s\n', cal.info.name);
% end

% CHECK FOR EXISTING WMO DIR, CREATE IF NOT THERE
if exist([dirs.mat,WMO,'\'],'dir') ~= 7
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

% ***********************************************************************
% GET DIR LIST AS STRUCTURE - THIS WILL BE USED TO SEND AN EMAIL OF
% ARRIVING NEW MSG FILES
mdir = dir([dirs.temp,'*.msg']);
% keyboard
[rr,~] = size(mdir);
if rr > 0
    new_msgs = cell(rr,1);
    ct = 0;
    for m_ct = 1:rr
        tmp_cycle = str2double(regexp(mdir(m_ct).name,'\d+(?=\.msg)', ...
            'once','match')); % profile # from name
        %TM 3/2/21.  Why did we exclude the missing casts from the
        %final reporting?  This is useful info for floats in
        %catchup mode.  Not sure what the original reasoning was?
        %Perhaps I'm missing something.  I think we want to modify
        %this moving forward, but maybe still distinguish between
        %the two (["new"] VS ["previously missing" or "catch-up"]).
        t1 = missing_casts == tmp_cycle;% missing or new? only add new
        if sum(t1) == 0
            ct = ct+1;
            str = sprintf('%s\t%s\t%s\t%s\t%0.1f','new     :', WMO, mdir(m_ct).name, ...
                datestr(mdir(m_ct).datenum,'mm/dd/yy HH:MM'), ...
                mdir(m_ct).bytes/1000);
            new_msgs{ct} = str;
            tf_float.status = 1;
        else
            ct = ct+1;
            str = sprintf('%s\t%s\t%s\t%s\t%0.1f','catch-up:', WMO, mdir(m_ct).name, ...
                datestr(mdir(m_ct).datenum,'mm/dd/yy HH:MM'), ...
                mdir(m_ct).bytes/1000);
            new_msgs{ct} = str;
            tf_float.status = 1;
        end
        REMAINDER = rem(tmp_cycle,5); % If the remainder is 0 then store a flag to trigger a full reprocess in loop-argo-float (but only if in update mode!!).
        if REMAINDER == 0 && strcmp(update_str, 'update')
            if mdir(m_ct).datenum > now-springchicken %new files
                FULLreprocess = 1;
                %         return %if we are doing a full reprocess...then don't need to progress further as the full float will get reprocessed in the next line of Loop_Argo_float (so break out here).
            end
        end

    end
    tf_float.new_messages = new_msgs(1:ct); % list of new msgs from float
    %clear mdir rr new_msgs str
% % % %     REMAINDER = rem(tmp_cycle,5); % If the remainder is 0 then store a flag to trigger a full reprocess in loop-argo-float (but only if in update mode!!).
% % % %     if REMAINDER == 0 && strcmp(update_str, 'update')
% % % %         if mdir(ct).datenum < springchicken
% % % %             FULLreprocess = 1;
% % % %             %         return %if we are doing a full reprocess...then don't need to progress further as the full float will get reprocessed in the next line of Loop_Argo_float (so break out here).
% % % %         end
% % % %     end
end

% ***********************************************************************

msg_list   = ls([dirs.temp,'*.msg']); % get list of file names to process

lr_ind_chk = 0; % indice toggle - find indices once per float
hr_ind_chk = 0; % toggle
pk_ind_chk = 0; % park toggle

master_FLBB = 0; % some floats (i.e.7558) change mode, start 1 , always 1

% GET BAD SENSOR LIST FOR FLOAT IF ANY
BSL = dirs.BSL;

disp(['Processing ARGO float ' cal.info.name, '.........'])
for msg_ct = 1:size(msg_list,1)

    clear LR HR INFO

    if exist('TRAJ','var')
        clear TRAJ
    end
    if exist('OCR','var')
        clear OCR
    end

    msg_file = strtrim(msg_list(msg_ct,:));
    disp(['PROCESSING MSG FILE ', msg_file])
    NO3_file = regexprep(msg_file,'msg','isus');
    pH_file  = regexprep(msg_file,'msg','dura');
    srf_file  = regexprep(msg_file,'msg','srf'); % for 2XO2 floats only
    % find block of numbers then look ahead to see if '.msg' follows
    cast_num = regexp(msg_file,'\d+(?=\.msg)','once','match');
    if regexp(MBARI_ID_str, ocr_testflt, 'once')
        %     if strcmp(MBARI_ID_str,ocr_testflt)==1
        d = parse_APEXmsg4ARGO_OCR([dirs.temp,msg_file]);
        if isempty(d.lr_d) && isempty(d.ocr_d)
            continue
        end
    else

        d = parse_APEXmsg4ARGO([dirs.temp,msg_file]); % PARSER
        if cal.info.O2_flag > 1 && isfile([dirs.temp,srf_file])
            d.srf_air = parse_APEX_srf([dirs.temp,srf_file]); % add surface air
            delete([dirs.temp,srf_file]);
        else
            d.srf_air = [];
        end
    end
    
    % Add ice detection true/false
    INFO.ice_flag = 0;
    if isfield(d,'ice_flag')
        if d.ice_flag
            INFO.ice_flag = 1;
        end
    end


    %     % IF MSG FILE EXIST BUT NO LR DATA REMOVE FROM NEW LIST AND
    %     % PUT IN BAD LIST
    %     if isempty(d.lr_d) || size(d.lr_d,1)<=3
    %         %         msg_size = size(msg_file,2);
    %         %         t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
    %         Indext1 = strfind(tf_float.new_messages,msg_file);
    %         t1 = find(not(cellfun('isempty',Indext1)));
    %         if ~isempty(t1)
    %             tf_float.bad_messages = [tf_float.bad_messages; ...
    %                 tf_float.new_messages(t1)];
    %             tf_float.new_messages(t1) =[];
    %             tf_float.status = 0;
    %         end
    %     end

    % IF MSG FILE EXIST BUT NO LR DATA REMOVE FROM NEW LIST AND
    % PUT IN BAD LIST. IF NOTHING GOOD LEFT, EXIT
    %     if ~strcmp(MBARI_ID_str,ocr_testflt)==1
    if isempty(regexp(MBARI_ID_str, ocr_testflt, 'once'))
        if isempty(d.lr_d) || size(d.lr_d,1)<=3
            t1 = ~cellfun(@isempty,regexp(tf_float.new_messages, msg_file,'once'));
            tf_float.bad_messages = [tf_float.bad_messages; ...
                tf_float.new_messages(t1)];
            tf_float.new_messages(t1) =[];
            if isempty(tf_float.new_messages)
                tf_float.status = 0; % NO VALID MESSAGES LEFT TO PROCESS
                disp(['No complete float message files found to process', ...
                    'for float ',MBARI_ID_str]);
                disp(['Last processed cast: ',sprintf('%03.0f',last_cast)])
                return
            end
        end
    end

    fileinfos = dir([dirs.temp,msg_file]);
    timestamps = fileinfos.date;
    timediff = now-datenum(timestamps); % was using 'd.sdn', but for floats just coming up from under ice that doesn't make sense.
    %make sure all 3 files are present if float includes all 3 sensors,
    %otherwise, end processing for this cycle (but if msg file is > 20 days old, and still no isus or dura, then process as usual)
    %     timediff = 50; % for manual processing override.
	TIMELAG = 0;
    if regexp(MBARI_ID_str, bad_no3_filter, 'once') % EXCEPTIONS
        disp([MBARI_ID_str,' Nitrate sensor died from the start & ', ...
            'never any isus files to  process'])
    else
        %         exist([dirs.temp, NO3_file],'file')
        %         isfield(cal,'N')
        %         timediff
        %         exist([dirs.temp, pH_file],'file')
        %         isfield(cal,'pH')
        %         pause
%         if (exist([dirs.temp, NO3_file],'file')==0 && isfield(cal,'N') && ...
%                 timediff<=TIMELAG)
        if (exist([dirs.temp, NO3_file],'file')==0 && isfield(cal,'N')) %TM, remove timelag, always process, but add message note if a file is missing
            disp('*******************************************************')
            disp(['WARNING: .isus FILE IS MISSING FOR: ',msg_file])
%             disp(['PROCESSING HALTED FOR ',msg_file,' UNTIL ALL FILES ARE PRESENT.'])
            disp(['PROCESSING ALL AVAILABLE DATA FOR ',msg_file,'.'])
            disp('*******************************************************')
            %             msg_size = size(msg_file,2);
            %             t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
            %             if sum(t1>0)
            Indext1 = strfind(tf_float.new_messages,msg_file);
            t1 = find(not(cellfun('isempty',Indext1)));
            if ~isempty(t1)
                tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' .isus file missing!']}; %mark on list as in limbo
            end
%             continue %TM Nov2023 don't skip processing!  but keep the
%             message
        end

%         if (exist([dirs.temp, pH_file],'file')==0 && ...
%                 isfield(cal,'pH') && timediff<=TIMELAG)
        if (exist([dirs.temp, pH_file],'file')==0 && ...
                isfield(cal,'pH'))
            disp('*******************************************************')
            disp(['WARNING: .dura FILE IS MISSING FOR: ',msg_file])
%             disp(['PROCESSING HALTED FOR ',msg_file,' UNTIL ALL FILES ARE PRESENT.'])
            disp(['PROCESSING AVAILABLE DATA FOR ',msg_file,'.'])
            disp('*******************************************************')
            %             msg_size = size(msg_file,2);
            %             t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
            %             if sum(t1>0)
            Indext1 = strfind(tf_float.new_messages,msg_file);
            t1 = find(not(cellfun('isempty',Indext1)));
            if ~isempty(t1)
                tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' .dura file missing!']}; %mark on list as in limbo
            end
%             continue
        end
    end

    if exist([dirs.temp,pH_file],'file')==1
        dura = parse_pHmsg([dirs.temp,pH_file]);
        if dura.EOT~=1 && timediff<=20
            disp('*******************************************************')
            disp('WARNING: .dura FILE IS MISSING <EOT>')
            disp(['PROCESSING HALTED FOR ',msg_file,' UNTIL ALL FILES ARE PRESENT IN FULL.'])
            disp('*******************************************************')
            %             msg_size = size(msg_file,2);
            %             t1 = strncmp(tf_float.new_messages, msg_file, msg_size);
            %             if sum(t1>0)
            Indext1 = strfind(tf_float.new_messages,msg_file);
            t1 = find(not(cellfun('isempty',Indext1)));
            if ~isempty(t1)
                tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' .dura size']}; %mark on list as in limbo
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
%     if d.cast == 1
%         keyboard
%     end
    if ~isempty(first_sdn) && ~isempty(INFO.sdn)
        dvec = datevec(INFO.sdn);
        if INFO.sdn > first_sdn + 365*20 && dvec(1) == 2099 % 20 yrs from start?
            disp(['GPS time for this profile is > 20 years past start ', ...
                '- gps week day number bug?!!'])
            %dvec = datevec(INFO.sdn); % I think this should be moved above 2nd IF? jp 06/20/21
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
    if isfield(d,'irid')
        INFO.iridium  = d.irid;
    else
        INFO.iridium = [];
    end

    yesBSAMLcyc=0;
    if yesBSAML==1
        SbsIND = find(strcmp(num2str(INFO.cast),dirs.BSAML.list(:,ibsCYC)));
        if ~isempty(SbsIND)
            tmpBSAML = dirs.BSAML;
            tmpBSAML.list = tmpBSAML.list(SbsIND,:);
            yesBSAMLcyc = 1;
        else
            yesBSAMLcyc = 0;
        end
    end

    % NEED TO CORRECT GPS SDN FOR WEEK DAY ROLLOVER BUG TOO (ua9634 & ua8482)
    if ~isempty(first_sdn) && ~isempty(INFO.gps)
        dvec = datevec(INFO.gps(:,1));
        if INFO.gps(1,1) > first_sdn + 365*20 && dvec(1,1) == 2099 % 20 yrs from start?
            disp(['GPS time for this profile is > 20 years past start ', ...
                '- gps week day number bug?!!'])
            dvec(:,1)  = 1999; % per aoml suggestion don't quite understand jump to 2099
            INFO.gps(:,1)= datenum(dvec) + 1024*7;
        elseif INFO.gps(1,1) < first_sdn % bad gps time fix 10/30/19
            disp('GPS time for this profile is unreasonable - gps week day number bug!!')
            disp(['days since first profile = ',num2str((INFO.gps(1,1) - ...
                first_sdn),'%0.0f')]);
            INFO.gps(:,1) = INFO.gps(:,1) + 1024*7;
        end

    end

    INFO.INST_ID    = cal.info.INST_ID;
    INFO.name       = cal.info.name;
    INFO.WMO_ID     = cal.info.WMO_ID;
    INFO.float_type = float_type;
    INFO.EOT        = d.EOT;

    % ****************************************************************
    % DEAL WITH LOW RESOLUTION DATA FIRST
    % ****************************************************************

    if isempty(regexp(MBARI_ID_str, ocr_testflt, 'once')) && (isempty(d.lr_d)) %|| size(d.lr_d,1)<=3) % CHECK FOR DATA; sometimes incomplete msg files will come through with minimal LR data (ie 12363.127.msg).
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

    if lr_ind_chk <2
        lr_ind_chk = lr_ind_chk+1;  %TM try check ind twice, in case faulty msg ie 000 skips index chk

        iP   = find(strcmp('p',      d.lr_hdr) == 1); % CTD P
        iT   = find(strcmp('t',      d.lr_hdr) == 1); % CTD T
        iS   = find(strcmp('s',      d.lr_hdr) == 1); % CTD S
        iPh83  = find(strcmp('Phase',   d.lr_hdr) == 1); % SBE83 O2 phase delay
        iT83 = find(strcmp('T83', d.lr_hdr) == 1); % SBE83 Temp (degC)
        iTo  = find(strcmp('Topt',   d.lr_hdr) == 1); % optode temp
        iTPh = find(strcmp('TPhase', d.lr_hdr) == 1); % T Phase 4330

        iTo1  = find(strcmp('Topt(1)',   d.lr_hdr) == 1); % optode temp
        iTPh1 = find(strcmp('TPhase(1)', d.lr_hdr) == 1); % T Phase 4330
        iTo2  = find(strcmp('Topt(2)',   d.lr_hdr) == 1); % optode temp
        iTPh2 = find(strcmp('TPhase(2)', d.lr_hdr) == 1); % T Phase 4330

        iRPh = find(strcmp('RPhase', d.lr_hdr) == 1); % R Phase 4330 % not used
        iBPh = find(strcmp('BPhase', d.lr_hdr) == 1); % B Phase 3830
        ipH  = find(strcmp('pH(V)',  d.lr_hdr) == 1); % pH
        iChl = find(strcmp('FSig',   d.lr_hdr) == 1); % CHL fluor
        if isempty(iChl)
            iChl = find(strcmp('FSig[0]',   d.lr_hdr) == 1); % CHL fluor
        end
        iChl435 = find(strcmp('FSig[1]',   d.lr_hdr) == 1); % CHL fluor

        iBb  = find(strcmp('BbSig',  d.lr_hdr) == 1); % Backscatter
        iNO3 = find(strcmp('no3',    d.lr_hdr) == 1); % Nitrate

        % CDOM PLACEHOLDER FOR NOW....
        iCdm = find(strcmp('Cdm',   d.lr_hdr) == 1); % CDOM, NAVIS

        % GET FLOATVIZ DATA - REG AND QC & GET HR DATA TOO
        % WILL BE USED TO EXTRACT QF DATA FLAGS LATER
        % DATA IS BEING PULLED FROM SIROCCO
        FVQC_flag   = 1;
        %         mymockSirocco = 'C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING_siroccoMock\DATA\FLOATVIZ\';
        %         FV_data     = get_FloatViz_data([mymockSirocco,INFO.WMO_ID,'.TXT']);
        %         FV_QCdata   = get_FloatViz_data([mymockSirocco,'QC\',INFO.WMO_ID,'QC.TXT']);
        %         FV_HRdata   = get_FloatViz_data([mymockSirocco,'HR\',INFO.WMO_ID,'_HR.TXT']);
        %         FV_HRQCdata = get_FloatViz_data([mymockSirocco,'HRQC\',INFO.WMO_ID,'_HRQC.TXT']);
        FV_data     = get_FloatViz_data(INFO.WMO_ID);
        FV_QCdata   = get_FloatViz_data([INFO.WMO_ID,'QC']);
        FV_HRdata   = get_FloatViz_data([INFO.WMO_ID,'_HR']);
        FV_HRQCdata = get_FloatViz_data([INFO.WMO_ID,'_HRQC']);

        if isempty(FV_QCdata)
            FV_QCdata = FV_data;
            FV_HRQCdata = FV_HRdata;
            FVQC_flag = 0;
        end
    end

    %GET SOME ARGO CORE PARAMETERS
    LR.PRES      = lr_d(:,iP);

    % MAKE AN ARRAY OF ZEROS FOR FILLING ARRAYS LATER
    fill0  = ones(size(LR.PRES))* 0;

    LR.PRES_QC          = fill0 + fv.QC;
    LR.PRES_ADJUSTED    = LR.PRES;
    LR.PRES_ADJUSTED_QC = fill0 + fv.QC;

    % PSAL proxy??
    %     if regexp(MBARI_ID_str, psal_proxy_flts, 'once')
    Indexp = strfind(psal_proxy_flts(:,1),MBARI_ID_str);
    Indexpsp = find(not(cellfun('isempty',Indexp)));
    if ~isempty(Indexpsp) && str2num(cast_num) > psal_proxy_flts{Indexpsp,2}
        INFO.PSAL_PROXY_USED = 1;
        proxydirname = ['\\atlas\Chem\ARGO_PROCESSING\DATA\ARMOR3D_PSAL_PROXY\',WMO,'\', WMO,'.', cast_num,'.mat'];
        load(proxydirname)
        if sum(LRproxy.PRES - LRproxy.PRES)>0
            disp('WARNING, PSAL proxy PRES MISMATCH!!')
            %             pause
        end
        LR.PSAL = LRproxy.PSAL;
        LR.PSAL_PROXY = LRproxy.PSAL; %include within the matfile storage as well.
        lr_potT   = theta(lr_d(:,iP), lr_d(:,iT), LR.PSAL,0);
        lr_den    = density(LR.PSAL, lr_potT); % kg/ m^3, pot den
    else
        LR.PSAL             = lr_d(:,iS); % else use data direct from msg file parser
        INFO.PSAL_PROXY_USED = 0;
    end

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
    % I don't think any change required to isbadsensor.m -JP 03/10/2021
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

    if yesBSAML == 1 && yesBSAMLcyc==1
        SbsIND = find(strcmp('S',tmpBSAML.list(:,ibsSENS)));
        if ~isempty(SbsIND)
            TMPsbs = tmpBSAML.list(SbsIND,:);
            singleBADs = TMPsbs(:,ibsD);
            singleBADSflags = TMPsbs(:,ibsFL);
            rangeBADs = TMPsbs(:,ibsDB);
            rangeBADsflags = TMPsbs(:,ibsFL);

            if ~isempty(singleBADs{1})
                for i2 = 1:length(singleBADs)
                    xxtmp = find(LR.PRES == singleBADs{i2});
                    if isempty(xxtmp)
                        disp('WARNING!! PSAL PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                    else
                        if LR.PSAL~=fv.bio
                            LR.PSAL_QC(xxtmp) = str2double(singleBADSflags{i2});
                        end
                        if LR.PSAL_ADJUSTED~=fv.bio
                            LR.PSAL_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                        end
                    end
                end
                clear i2
            end
            for i3 = 1:length(rangeBADs)
                LR.PSAL_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.PSAL~=fv.bio) = str2double(rangeBADsflags{i3});
                LR.PSAL_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.PSAL_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
            end
            clear i3
        end
    end


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
    O2_flag = cal.info.O2_flag;
    if O2_flag > 0 && isfield(cal,'O') % At least 1 oxygen calibration
        doxy_fields = {'DOXY' 'O' 'RAW'; ... % Data field, cal field names
            'DOXY2' 'O2' 'RAW2'};

        for O2ct = 1:O2_flag % Loop through sensors, usually just 1
            % SET UP STRINGS TO BUILD FIELD NAMES
            O2fn      = doxy_fields{O2ct,1};
            O2fnQC    = [O2fn,'_QC'];
            O2fnadj   = [O2fn,'_ADJUSTED'];
            O2fnadjQC = [O2fnadj,'_QC'];
            O2cal     = doxy_fields{O2ct,2};
            air_raw   = doxy_fields{O2ct,3};

            % CHECK OPTODE TYPE & SET INDICE TYPE TO COMMON INDICE
            if O2_flag == 1 && strcmp(cal.(O2cal).type,'4330')
                iPhase = iTPh;
            elseif O2_flag == 1 && strcmp(cal.(O2cal).type,'3830')
                iPhase = iBPh;
            elseif O2_flag == 1 && strcmp(cal.(O2cal).type,'SBE83')
                iPhase = iPh83;
                iTo    = iT83;
            elseif O2_flag > 1 && strcmp(cal.(O2cal).type,'4330') ...
                    && O2ct == 1
                iPhase = iTPh1;
                iTo    = iTo1;
            elseif O2_flag > 1 && strcmp(cal.(O2cal).type,'4330') ...
                    && O2ct == 2
                iPhase = iTPh2;
                iTo    = iTo2;
            end

            % CHECK FOR GDF FLOAT & set optode T to CTD T, JP 11/08/2021
            if regexp(MBARI_ID_str,  gdf_wierd_O2, 'once')
                fprintf(['%s has bad T optode but good TPhase. Using CTD T', ...
                    'instead\n'], MBARI_ID_str)
                iTo = iT;
            end

            t_nan = isnan(lr_d(:,iPhase));
            O2_chk = lr_d(:,iPhase) < RCR.OP(1) | ...
                lr_d(:,iPhase) > RCR.OP(2) | ...
                lr_d(:,iTo) < RCR.OT(1) | ...
                lr_d(:,iTo) > RCR.OT(2); % Flag phase, Temp out of range
            if any(O2_chk)
                disp(['Out of range phase or optode T detected for ', ...
                    msg_file])
            end

            if strcmp(cal.(O2cal).type,'3830')
                LR.(['BPHASE_',O2fn])     = fill0 + fv.bio; % predim w fill val
                LR.(['BPHASE_',O2fnQC])         = fill0 + fv.QC;
                LR.(['BPHASE_',O2fn])(~t_nan)    = lr_d(~t_nan,iPhase);
                LR.(['BPHASE_',O2fnQC])(~t_nan) = fv.QC;
            elseif strcmp(cal.(O2cal).type,'4330')
                LR.(['TPHASE_',O2fn])            = fill0 + fv.bio;
                LR.(['TPHASE_',O2fnQC])         = fill0 + fv.QC;
                LR.(['TPHASE_',O2fn])(~t_nan)    = lr_d(~t_nan,iPhase);
                LR.(['TPHASE_',O2fnQC])(~t_nan) =  fv.QC;
            elseif strcmp(cal.(O2cal).type,'SBE83')
                LR.(['PHASE_DELAY_',O2fn])            = fill0 + fv.bio;
                LR.(['PHASE_DELAY_',O2fnQC])         = fill0 + fv.QC;
                LR.(['PHASE_DELAY_',O2fn])(~t_nan)    = lr_d(~t_nan,iPhase);
                LR.(['PHASE_DELAY_',O2fnQC])(~t_nan) =  fv.QC;
            end
            LR.(['PPOX_',O2fn])               = fill0 + fv.bio;
            LR.(O2fn)               = fill0 + fv.bio;
            LR.(O2fnQC)             = fill0 + fv.QC;
            LR.(['TEMP_',O2fn])            = fill0 + fv.bio;
            LR.(['TEMP_',O2fnQC])        = fill0 + fv.QC;
            LR.(O2fnadj)       = fill0 + fv.bio;
            LR.(O2fnadjQC)    = fill0 + fv.QC;
            LR.([O2fnadj,'_ERROR']) = fill0 + fv.bio;
            INFO.([O2fn,'_SCI_CAL_EQU'])  = 'not applicable';
            INFO.([O2fn,'_SCI_CAL_COEF']) = 'not applicable';
            INFO.([O2fn,'_SCI_CAL_COM'])  = 'not applicable';
            INFO.([O2fn,'_DATA_MODE'])  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

            if strcmp(cal.(O2cal).type,'SBE83')
                myOdata = [LR.PRES(~t_nan) LR.TEMP(~t_nan) LR.PSAL(~t_nan) lr_d((~t_nan),iPhase),lr_d((~t_nan),iTo)];
                [ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myOdata, cal.(O2cal),0);
                LR.(O2fn)(~t_nan) = O2_uM ./ lr_den(~t_nan) .* 1000; % µmol/kg
                LR.(['PPOX_',O2fn])(~t_nan) = ppoxdoxy; %TM PPOX_DOXY added for storage internally
            else
                %[S, T, P, Phase, CalPhase, [O2], O2sol, pO2]
                O2 = calc_O2_4ARGO(LR.PSAL(~t_nan), LR.TEMP(~t_nan), ...
                    LR.PRES(~t_nan),lr_d((~t_nan),iPhase),lr_d((~t_nan),iTo), cal.(O2cal)); % O2 in µmol/L + more
                LR.(O2fn)(~t_nan) = O2(:,6) ./ lr_den(~t_nan) .* 1000; % µmol/kg
                LR.(['PPOX_',O2fn])(~t_nan) = O2(:,end); %TM PPOX_DOXY added for storage internally
            end
            tDOXY = abs(LR.(O2fn)) > crazy_val & ~t_nan; % Unrealistic bad value
            LR.(O2fn)(tDOXY) = crazy_val; % SET TO crazy bad value
            LR.(O2fnQC)(~t_nan)      = 3;
            % JP: QC assignment below not needed - crazy value greater than range limits
            %LR.DOXY_QC(tDOXY) = 4; % set crazy bad value QF to 4
            LR.(['TEMP_',O2fn])(~t_nan)    = lr_d(~t_nan,iTo);
            LR.(['TEMP_',O2fnQC])(~t_nan) = fv.QC;
            % Save O2sol, pO2?
            clear O2

            if isfield(QC,O2cal)
                % !! ONE TIME GAIN CORRECTION ONLY !!
                %             LR.DOXY_ADJUSTED(~t_nan)  = LR.DOXY(~t_nan) .* QC.O.steps(1,3);
                QCD = [LR.PRES(~t_nan), LR.TEMP(~t_nan), LR.PSAL(~t_nan), LR.(O2fn)(~t_nan)];
                LR.(O2fnadj)(~t_nan) = apply_QC_corr(QCD, INFO.sdn, QC.(O2cal));
                tDOXY_ADJ = abs(LR.(O2fnadj)) > crazy_val & ~t_nan; % Unrealistic bad value
                LR.(O2fnadj)(tDOXY_ADJ) = crazy_val; % SET TO crazy bad value
                LR.(O2fnadjQC)(~t_nan) = 1; % set=1 9/27/16 vs 2 = probably good
                % JP: QC assignment below not needed - crazy value greater than range limits
                %LR.DOXY_ADJUSTED_QC(tDOXY_ADJ) = 4; %set crazy val QF to bad
                LR.([O2fnadj,'_ERROR'])(~t_nan) = LR.(O2fnadj)(~t_nan) * 0.01;
                INFO.([O2fn,'_SCI_CAL_EQU'])  = [O2fnadj,'=',O2fn, ...
                    '*G; G = G_INIT + G_DRIFT*(JULD_PROF - JULD_INIT)/365'];
                % QC matrix entry relevant to current cycle.
                steptmp = find(QC.(O2cal).steps(:,2)<=INFO.cast,1,'last');
                juld_prof = INFO.sdn-datenum(1950,01,01); %convert to JULD
                juld_init = QC.(O2cal).steps(steptmp,1)-datenum(1950,01,01); %convert to JULD
                juld_end = QC.date -datenum(1950,01,01); %date at last DMQC, converted to JULD

                % Always define LAST_DOXY_ERROR first using the static mbar
                % error spec, then if the profile time is later than date
                % of last dmqc, inflate error.
                if ~isempty(d.air)  % 2mbar error spec for air-calibrated
                    LR.([O2fnadj,'_ERROR'])(~t_nan) = ...
                        convert_O2mb_error_to_conc(LR.TEMP(~t_nan),LR.PSAL(~t_nan),2);
                    %LAST_DOXY_ERROR = LR.([O2fnadj,'_ERROR']);
                else % 4-6mbar error spec for NON air-calibrated
                    LR.([O2fnadj,'_ERROR'])(~t_nan) = ...
                        convert_O2mb_error_to_conc(LR.TEMP(~t_nan),LR.PSAL(~t_nan),5);
                    %LAST_DOXY_ERROR = LR.([O2fnadj,'_ERROR']);
                end

                if juld_prof>juld_end %1 mb per year error inflation per Argo rec
                    extra_error_ppox = 1.*(juld_prof-juld_end)./365;
                    O2error = convert_O2mb_error_to_conc(LR.TEMP(~t_nan),...
                        LR.PSAL(~t_nan),extra_error_ppox);
                    %LR.DOXY_ADJUSTED_ERROR(~t_nan) = LAST_DOXY_ERROR(~t_nan)+O2error;
                    LR.([O2fnadj,'_ERROR'])(~t_nan) = LR.([O2fnadj,'_ERROR'])(~t_nan)+O2error;
                end

                INFO.([O2fn,'_SCI_CAL_COEF']) = ['G_INIT = ', ...
                    num2str(QC.(O2cal).steps(steptmp,3),'%0.4f'),...
                    '; G_DRIFT = ',num2str(QC.(O2cal).steps(steptmp,5),'%0.4f'),...
                    '; JULD_PROF = ',num2str(juld_prof,'%9.4f'),...
                    '; JULD_INIT = ',num2str(juld_init,'%9.4f')];

                if isfield(cal.(O2cal),'SVUFoilCoef') && ~strcmp(cal.O.type,'SBE83')
                    O2_cal_str = 'SVU Foil calibration coeficients were used. ';
                else
                    O2_cal_str = 'Polynomial calibration coeficients were used. ';
                end
                if ~isempty(d.air) || ~isempty(d.srf_air)
                    INFO.([O2fn,'_SCI_CAL_COM'])  = [O2_cal_str,'G determined from ',...
                        'float  measurements in air. See Johnson et al.,2015,', ...
                        'doi:10.1175/JTECH-D-15-0101.1'];
                else
                    INFO.([O2fn,'_SCI_CAL_COM'])   = [O2_cal_str,'G determined by surface' ...
                        ' measurement comparison to World Ocean Atlas 2018.', ...
                        'See Takeshita et al.2013,doi:10.1002/jgrc.20399'];
                end
            end

            % IN AIR MEASUREMENTS ASSOCIATED WITH THE TELEMETRY CYCLE. THESE
            % ARE THE OLDER IN AIR MEASUREMENTS.
            if ~isempty(d.air) && sum(d.air(:)) ~= 0 % 0 means ice detection on

                % CALCULATE pO2
                zero_fill = d.air(:,1) * 0; % make array of zeros
                if strcmp(cal.(O2cal).type,'SBE83')
                    myadata = [zero_fill d.air(:,1) zero_fill d.air(:,2) d.air(:,1)]; %use optode T here for CTD T input?
                    [ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myadata, cal.(O2cal),0);
                else
                    O2 = calc_O2_4ARGO(zero_fill, d.air(:,1), zero_fill, ...
                        d.air(:,2), d.air(:,1), cal.(O2cal));
                    ppoxdoxy = O2(:,end); %last column is pO2
                end

                SurfaceObs.(air_raw)          = d.air;
                SurfaceObs.(['TEMP_',O2fn])   = d.air(:,1);
                SurfaceObs.(['TPHASE_',O2fn]) = d.air(:,2);
                %tmmp = d.air;
                if size(d.air,2)>2
                    SurfaceObs.(['RPHASE_',O2fn])   = d.air(:,3);
                    SurfaceObs.(['RPHASE_',O2fnQC]) = zero_fill + fv.bio;
                    rphaseBAD = SurfaceObs.(['RPHASE_',O2fn]) < RCR.RP(1) | ...
                        SurfaceObs.(['RPHASE_',O2fn]) > RCR.RP(2);
                    SurfaceObs.(['RPHASE_',O2fnQC])(rphaseBAD)  = 4;
                    SurfaceObs.(['RPHASE_',O2fnQC])(~rphaseBAD) = 1;
                    clear rphaseBAD
                end
                SurfaceObs.(['PPOX_',O2fn])   = ppoxdoxy;

                % initialize qc variables and perform qc:
                SurfaceObs.(['PPOX_',O2fnQC])   = zero_fill + fv.bio;
                SurfaceObs.(['TEMP_',O2fnQC])   = zero_fill + fv.bio;
                SurfaceObs.(['TPHASE_',O2fnQC]) = zero_fill + fv.bio;

                po2BAD    = SurfaceObs.(['PPOX_',O2fn]) < RCR.PO2(1) | ...
                    SurfaceObs.(['PPOX_',O2fn]) > RCR.PO2(2);
                tempBAD   = SurfaceObs.(['TEMP_',O2fn]) < RCR.OT(1) | ...
                    SurfaceObs.(['TEMP_',O2fn]) > RCR.OT(2);
                tphaseBAD = SurfaceObs.(['TPHASE_',O2fn]) < RCR.OP(1) | ...
                    SurfaceObs.(['TPHASE_',O2fn]) > RCR.OP(2);

                SurfaceObs.(['PPOX_',O2fnQC])(po2BAD)       = 4;
                SurfaceObs.(['PPOX_',O2fnQC])(~po2BAD)      = 1;
                SurfaceObs.(['TEMP_',O2fnQC])(tempBAD)      = 4;
                SurfaceObs.(['TEMP_',O2fnQC])(~tempBAD)     = 1;
                SurfaceObs.(['TPHASE_',O2fnQC])(tphaseBAD)  = 4;
                SurfaceObs.(['TPHASE_',O2fnQC])(~tphaseBAD) = 1;
                %            SurfaceObs.RPHASE_DOXY_QC(rphaseBAD)=4;
                %            SurfaceObs.RPHASE_DOXY_QC(~rphaseBAD)=1;
                clear O2 po2BAD tempBAD tphaseBAD
            end

            % NEW AIR CAL MEASUREMEMTS: 4 JUST BELOW THE SURFACE & 8 IN AIR.
            % THE IN AIR MEASUREMENTS CAN BE ADVERSELY AFFECTED IF THE FLOAT
            % EXPERIENCES AIR BLADDER INFLATION PROBLEMS!
            % aircal = [sdn, bladder(?), pres, temperature, Tphase, Rphase]
            if ~isempty(d.aircal) % Not sure why the second part of this was needed?  When would all phase meas be zero??  & sum(d.aircal(:,5)) ~= 0
                zero_fill = d.aircal(:,1) * 0; % make array of zeros
                tzmin     = lr_d(:,iP) == min(lr_d(:,iP)); % shallowest value
                S0 = zero_fill + mean(lr_d(tzmin,iS),'omitnan'); % surface salinity.  Needed in O2 calc for near sfc in-air sequence samples
                % TM 5/25/21; reorganize the in-air data structures for Annie
                % to access.  The in-air PPOX_DOXY will be used for populating
                % the Dtraj files.  The variables will require QC as well.
                %AIRCAL_O2 = [d.aircal(:,1:3),O2]; % µmol/L still NOT USED YET.

                %run the check on nearsfc & in-air based on pneumatic pressure
                diffPneuPres = abs(diff(d.aircal(:,2)));
                Xmax = find(diffPneuPres == max(diffPneuPres,[],'omitnan'));
                %MAXind = Xmax+1; %index of end of near sfc data  Should reflect a large shift in pneumatic pressure.
                %Bdef_s = Xmax; %indices of air-cal surface samples
                %Bdef_a = Xmax+1:size(d.aircal,1);
                %inair_bin = zeros(size(d.aircal,1),1); %inair_bin is the last column in the raw_inair_series and will serve as a logical indicator of true "in-air" part of the series (1=inair; 0=nearsurface)
                %inair_bin(Bdef_a) = 1;
                inair_bin             = zero_fill;
                inair_bin(Xmax+1:end) = 1;
                raw_inair_series      = [d.aircal inair_bin];
                if strcmp(cal.(O2cal).type,'SBE83')
                    myadata = [zero_fill d.aircal(:,4) S0.*~logical(inair_bin) d.aircal(:,5) d.aircal(:,4)]; %use optode T here for CTD T input?
                    [ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myadata, cal.(O2cal),0);
                else
                    O2 = calc_O2_4ARGO(zero_fill, d.aircal(:,4), zero_fill, ...
                        d.aircal(:,5), d.aircal(:,4), cal.(O2cal));
                    ppoxdoxy = O2(:,end);
                end
                OptodeAirCal.(air_raw) = raw_inair_series; %store the raw data;
                OptodeAirCal.JULD      = raw_inair_series(:,1)-datenum(1950,01,01); %convert to JULD
                OptodeAirCal.PRES      = raw_inair_series(:,3);
                OptodeAirCal.(['TEMP_',O2fn]) = raw_inair_series(:,4);
                if strcmp(cal.(O2cal).type,'SBE83')
                    OptodeAirCal.(['PHASE_DELAY_',O2fn]) = raw_inair_series(:,5);
                else
                    OptodeAirCal.(['TPHASE_',O2fn]) = raw_inair_series(:,5);
                    OptodeAirCal.(['RPHASE_',O2fn]) = raw_inair_series(:,6);
                end

                OptodeAirCal.IN_AIR_LOGICAL   = raw_inair_series(:,end);
                OptodeAirCal.(['PPOX_',O2fn]) = ppoxdoxy; %This is pO2 returned from calc_O2_4ARGO (or Calc_SBE63_O2 for the SBE83s)
                %initialize qc vars and check qc
                OptodeAirCal.(['PPOX_',O2fnQC]) = zero_fill + fv.bio;
                OptodeAirCal.(['TEMP_',O2fnQC]) = zero_fill + fv.bio;
                if strcmp(cal.(O2cal).type,'SBE83')
                    OptodeAirCal.(['PHASE_DELAY_',O2fnQC]) = zero_fill + fv.bio;
                else
                    OptodeAirCal.(['TPHASE_',O2fnQC]) = zero_fill + fv.bio;
                    OptodeAirCal.(['RPHASE_',O2fnQC]) = zero_fill + fv.bio;
                end
                po2BAD  = OptodeAirCal.(['PPOX_',O2fn]) < RCR.PO2(1) | ...
                    OptodeAirCal.(['PPOX_',O2fn]) > RCR.PO2(2);
                tempBAD = OptodeAirCal.(['TEMP_',O2fn]) < RCR.OT(1) | ...
                    OptodeAirCal.(['TEMP_',O2fn]) > RCR.OT(2);
                if strcmp(cal.(O2cal).type,'SBE83')
                    tphaseBAD = OptodeAirCal.(['PHASE_DELAY_',O2fn]) < RCR.OP(1) | ...
                        OptodeAirCal.(['PHASE_DELAY_',O2fn]) > RCR.OP(2);
                else
                    tphaseBAD = OptodeAirCal.(['TPHASE_',O2fn]) < RCR.OP(1) | ...
                        OptodeAirCal.(['TPHASE_',O2fn]) > RCR.OP(2);
                    rphaseBAD = OptodeAirCal.(['RPHASE_',O2fn]) < RCR.RP(1) | ...
                        OptodeAirCal.(['RPHASE_',O2fn]) > RCR.RP(2);
                end
                %             save('tanyatemp.mat','AIRCAL_O2')
                % cludgy for now until we get more guidance from Dana...Apex
                % starting at 18000 series (Apf11) have valid pneumatic pres
                % range for in-air data as
                if str2num(INFO.INST_ID)>18000
                    pneupresBAD = (OptodeAirCal.(air_raw)(:,end)==1 & OptodeAirCal.(air_raw)(:,2)>160) ...
                        | (OptodeAirCal.(air_raw)(:,end)==1 & OptodeAirCal.(air_raw)(:,2)<150); %column 2 is pneumpres.  Only due this check on the "inair" part of the series.
                else
                    %pneupresBAD = false(size(OptodeAirCal.PPOX_DOXY));
                    pneupresBAD = logical(zero_fill);
                end
                OptodeAirCal.(['PPOX_',O2fnQC])(po2BAD)=4;
                OptodeAirCal.(['PPOX_',O2fnQC])(~po2BAD)=1;
                OptodeAirCal.(['PPOX_',O2fnQC])(pneupresBAD)=4;
                OptodeAirCal.(['TEMP_',O2fnQC])(tempBAD) = 4;
                OptodeAirCal.(['TEMP_',O2fnQC])(~tempBAD) = 1;
                if strcmp(cal.O.type,'SBE83')
                    OptodeAirCal.(['PHASE_DELAY_',O2fnQC])(tphaseBAD) = 4;
                    OptodeAirCal.(['PHASE_DELAY_',O2fnQC])(~tphaseBAD) = 1;
                else
                    OptodeAirCal.(['TPHASE_',O2fnQC])(tphaseBAD) = 4;
                    OptodeAirCal.(['TPHASE_',O2fnQC])(~tphaseBAD) = 1;
                    OptodeAirCal.(['RPHASE_',O2fnQC])(rphaseBAD) = 4;
                    OptodeAirCal.(['RPHASE_',O2fnQC])(~rphaseBAD) = 1;
                end
                clear O2 diffPneuPres inair_bin raw_inair_series po2BAD
                clear pneupresBAD Xmax MAXind Bdef_s Bdef_a tempBAD tphaseBAD rphaseBAD
            end

            % DUAL O2 FLOATS IN AIR MEASUREMENTS IN SEPERATE *.srf FILES WHICH
            % CONTAINING IN AIR MEASUREMENTS FOR BOTH SENSORS (4330's only so far)
            if  ~isempty(d.srf_air)
                SurfaceObs.(air_raw)          = d.srf_air.data;
                SurfaceObs.([air_raw,'_hdr']) = d.srf_air.hdr;

                % GET INDICES FOR PROPER T, S & PHASE EXTRACTION
                srfT  = strcmp(d.srf_air.hdr, sprintf('T%0.0f(C)',   O2ct));
                srfTP = strcmp(d.srf_air.hdr, sprintf('TPhase%0.0f', O2ct));
                srfRP = strcmp(d.srf_air.hdr, sprintf('RPhase%0.0f', O2ct));

                % CALCULATE pO2 % DUAL 4330 O2 FLOATS
                zero_fill = d.srf_air.data(:,1) * 0; % make array of zeros

                O2 = calc_O2_4ARGO(zero_fill, d.srf_air.data(:,srfT), ...
                    zero_fill, d.srf_air.data(:,srfTP), d.srf_air.data(:,srfT), cal.(O2cal));
                ppoxdoxy = O2(:,end); %last column is pO2

                SurfaceObs.(['TEMP_',O2fn])   = d.srf_air.data(:,srfT);
                SurfaceObs.(['TPHASE_',O2fn]) = d.srf_air.data(:,srfTP);
                SurfaceObs.(['RPHASE_',O2fn]) = d.srf_air.data(:,srfRP);
                SurfaceObs.(['PPOX_',O2fn])   = ppoxdoxy;

                % initialize qc variables and perform qc:
                SurfaceObs.(['PPOX_',O2fnQC])   = zero_fill + fv.bio;
                SurfaceObs.(['TEMP_',O2fnQC])   = zero_fill + fv.bio;
                SurfaceObs.(['TPHASE_',O2fnQC]) = zero_fill + fv.bio;
                SurfaceObs.(['RPHASE_',O2fnQC]) = zero_fill + fv.bio;

                po2BAD    = SurfaceObs.(['PPOX_',O2fn]) < RCR.PO2(1) | ...
                    SurfaceObs.(['PPOX_',O2fn]) > RCR.PO2(2);
                tempBAD   = SurfaceObs.(['TEMP_',O2fn]) < RCR.OT(1) | ...
                    SurfaceObs.(['TEMP_',O2fn]) > RCR.OT(2);
                tphaseBAD = SurfaceObs.(['TPHASE_',O2fn]) < RCR.OP(1) | ...
                    SurfaceObs.(['TPHASE_',O2fn]) > RCR.OP(2);
                rphaseBAD = SurfaceObs.(['RPHASE_',O2fn]) < RCR.RP(1) | ...
                    SurfaceObs.(['RPHASE_',O2fn]) > RCR.RP(2);

                SurfaceObs.(['PPOX_',O2fnQC])(po2BAD)       = 4;
                SurfaceObs.(['PPOX_',O2fnQC])(~po2BAD)      = 1;
                SurfaceObs.(['TEMP_',O2fnQC])(tempBAD)      = 4;
                SurfaceObs.(['TEMP_',O2fnQC])(~tempBAD)     = 1;
                SurfaceObs.(['TPHASE_',O2fnQC])(tphaseBAD)  = 4;
                SurfaceObs.(['TPHASE_',O2fnQC])(~tphaseBAD) = 1;
                SurfaceObs.(['RPHASE_',O2fnQC])(rphaseBAD)  = 4;
                SurfaceObs.(['RPHASE_',O2fnQC])(~rphaseBAD) = 1;

                clear O2 po2BAD tempBAD tphaseBAD rphaseBAD
                clear hdr data sO2_ct zero_fill ppoxdoxy O2 sT po2BAD
            end


            % DO SOME FINAL CHECK ON VALUES, IF BAD SET QF = 4
            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, O2cal);
            t_bio = LR.(O2fn) ~= fv.bio; % Non fill value samples
            %tST   = LR.PSAL_QC == 4 | LR.TEMP_QC == 4; % Bad S or T will affect O2
            tST   = LR.PSAL_QC == 4 | LR.TEMP_QC == 4 | LR.PRES_QC == 4; % Bad S or T will affect O2
            t_chk = LR.(O2fn) < RCR.O(1)| LR.(O2fn) > RCR.O(2); % range check
            t_chk = t_chk | O2_chk | tST; % BAD O2 T or phase

            t_chk = t_chk & t_bio; % & not MVI either

            if isfield(LR,['BPHASE_',O2fn])
                LR.(['BPHASE_',O2fnQC])(t_bio) = LR.(['BPHASE_',O2fnQC])(t_bio) * ...
                    ~BSLflag + BSLflag*theflag;
                LR.(['BPHASE_',O2fnQC])(t_chk) = 4;
            elseif isfield(LR,['TPHASE_',O2fn])
                LR.(['TPHASE_',O2fnQC])(t_bio) = LR.(['TPHASE_',O2fnQC])(t_bio) *  ...
                    ~BSLflag + BSLflag*theflag;
                LR.(['TPHASE_',O2fnQC])(t_chk) = 4;
            end

            % VERY BAD PSAL & TEMP = BAD O2 TOO
            LR.(O2fnQC)(t_bio) = LR.(O2fnQC)(t_bio) * ~BSLflag + BSLflag*theflag;
            LR.(O2fnQC)(t_chk) = 4;

            if isfield(QC, O2cal)
                t_bio = LR.(O2fnadj) ~= fv.bio;
                %tST   = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4; % Bad S or T will affect O2
                tST   = LR.PSAL_ADJUSTED_QC == 4 | LR.TEMP_ADJUSTED_QC == 4 ...
                    | LR.PRES_ADJUSTED_QC == 4;
                t_chk = LR.(O2fnadj) < RC.O(1)| ...
                    LR.(O2fnadj) > RC.O(2) | O2_chk;
                t_chk = (t_chk | tST) & t_bio;

                LR.(O2fnadjQC)(t_bio) = LR.(O2fnadjQC)(t_bio) * ...
                    ~BSLflag + BSLflag*theflag;
                LR.(O2fnadjQC)(t_chk) = 4;
            end

            % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            % DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
            % (BAD).  BGC_spiketest screens for nans and fill values and pre-identified bad values.
            %
            % RUN TEST ON RAW DOXY
            % however....if Arctic float, do not perform O2 spiketest
            %(very large gradient near surface...does not work well in this region!)
            if regexp(MBARI_ID_str ,no_o2_spike,'once')
                disp([MBARI_ID_str,' is on the no DOXY spike test list'])
            else
                QCscreen_O = LR.(O2fnQC) == 4; % screen for BAD data already assigned.
                [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str, ...
                    str2double(cast_num),[LR.PRES LR.(O2fn)],'O2',dirs.cal,fv.bio,QCscreen_O);
                if ~isempty(spike_inds)
                    LR.(O2fnQC)(spike_inds) = quality_flags;
                    disp(['LR.',O2fn,' QUALITY FLAGS ADJUSTED FOR IDENTIFIED ', ...
                        'SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
                end
                %
                % RUN TEST ON QC DOXY_ADJUSTED
                QCscreen_Oadj = LR.(O2fnadjQC) == 4; % screen for BAD data already assigned.
                [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str, ...
                    str2double(cast_num),[LR.PRES LR.(O2fnadj)],'O2',...
                    dirs.cal,fv.bio,QCscreen_Oadj);
                if ~isempty(spike_inds)
                    LR.(O2fnadjQC)(spike_inds) = quality_flags;
                    disp(['LR.',O2fnadj,'DOXY_ADJUSTED QUALITY FLAGS ADJUSTED FOR ', ...
                        'IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
                end
                % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
                clear QCscreen_O QCscreenOadj spike_inds quality_flags
            end

            %------------------------------------------------------------------------------------
            %FLAG ANY BAD SAMPLES FROM 'BAD-SAMPLE LIST':
            if yesBSAML == 1 && yesBSAMLcyc==1
                SbsIND = find(strcmp('O',tmpBSAML.list(:,ibsSENS)));
                if ~isempty(SbsIND)
                    TMPsbs = tmpBSAML.list(SbsIND,:);
                    singleBADs = TMPsbs(:,ibsD);
                    singleBADSflags = TMPsbs(:,ibsFL);
                    rangeBADs = TMPsbs(:,ibsDB);
                    rangeBADsflags = TMPsbs(:,ibsFL);
                    if ~isempty(singleBADs{1})
                        for i2 = 1:length(singleBADs)
                            xxtmp = find(LR.PRES == singleBADs{i2});
                            if isempty(xxtmp)
                                disp('WARNING!! DOXY PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                            else
                                if LR.DOXY~=fv.bio
                                    LR.DOXY_QC(xxtmp) = str2double(singleBADSflags{i2});
                                end
                                if LR.DOXY_ADJUSTED~=fv.bio
                                    LR.DOXY_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                                end
                            end
                        end
                        clear i2
                    end
                    for i3 = 1:length(rangeBADs)
                        LR.DOXY_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.DOXY~=fv.bio) = str2double(rangeBADsflags{i3});
                        LR.DOXY_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.DOXY_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                    end
                    clear i3
                end
            end



        end
    end


    % ****************************************************************
    % CALCULATE CHLOROPHYLL CONCENTRATION (µg/L or mg/m^3)
    % ****************************************************************

    if (~isempty(iChl) && master_FLBB ~= 0) || ...
            (~isempty(iChl) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))

        % Not sure why master_FLBB variable part of the statement, seems having
        % just ~istempty(iChl) would suffice?
        % JP 03/10/20201 - There are a handful of earlier floats where the msg
        % file header indicates that an FLBB exists but there isn't one on the
        % float. For these floats the FLBB_mode = 0 from the start maybe there
        % is a better way...

        % Predim variables with fill values then adjust as needed
        LR.FLUORESCENCE_CHLA    = fill0 + fv.bio;
        LR.FLUORESCENCE_CHLA_QC = fill0 + fv.QC;
        LR.CHLA                 = fill0 + fv.bio;
        LR.CHLA_QC              = fill0 + fv.QC;
        LR.CHLA_ADJUSTED        = fill0 + fv.bio;
        LR.CHLA_ADJUSTED_QC     = fill0 + fv.QC;
        LR.CHLA_ADJUSTED_ERROR  = fill0 + fv.bio;
		LR.CHLA_FLUORESCENCE                 = fill0 + fv.bio;
        LR.CHLA_FLUORESCENCE_QC              = fill0 + fv.QC;
        LR.CHLA_FLUORESCENCE_ADJUSTED        = fill0 + fv.bio;
        LR.CHLA_FLUORESCENCE_ADJUSTED_QC     = fill0 + fv.QC;
        LR.CHLA_FLUORESCENCE_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CHLA_SCI_CAL_EQU   = 'not applicable';
        INFO.CHLA_SCI_CAL_COEF  = 'not applicable';
        INFO.CHLA_SCI_CAL_COM   = 'not applicable';
        INFO.CHLA_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
        INFO.CHLA_FLUORESCENCE_SCI_CAL_EQU   = 'not applicable';
        INFO.CHLA_FLUORESCENCE_SCI_CAL_COEF  = 'not applicable';
        INFO.CHLA_FLUORESCENCE_SCI_CAL_COM   = 'not applicable';
        INFO.CHLA_FLUORESCENCE_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
		
        t_nan = isnan(lr_d(:,iChl)); % NaN's in data if any

        LR.FLUORESCENCE_CHLA(~t_nan)    = lr_d(~t_nan,iChl);
        LR.FLUORESCENCE_CHLA_QC(~t_nan) = fv.QC;

        if isfield(cal,'CHL') % Sensor could be bad so maybe no cal info
            LR.CHLA(~t_nan) = (lr_d(~t_nan,iChl) - cal.CHL.ChlDC) .* ...
                cal.CHL.ChlScale;
            LR.CHLA_QC(~t_nan) =  3; % 3 do not use w/o adjusting
			% same for CHLA_FLUORESCENCE (just different units!) TM 5/14/24
            LR.CHLA_FLUORESCENCE(~t_nan) = (lr_d(~t_nan,iChl) - cal.CHL.ChlDC) .* ...
                cal.CHL.ChlScale;
            LR.CHLA_FLUORESCENCE_QC(~t_nan) =  3; % 3 do not use w/o adjusting

            % ADJUSTED DATA BASED ON ADMT18 CONCENSUS -WILL BE UPDATED
            % WITHIN THE YEAR - jp 12/13/2017

            % 1st CHECK FOR IN SITU DARKS & TRY AND GET THEM IF NOT THERE
            if ~isfield(cal.CHL, 'SWDC')
                SWDC = get_CHLdark(cal, dirs, 5); % 1st 5 good profiles JP 07/12/23
                if ~isempty(SWDC)
                    cal.CHL.SWDC = SWDC; % add structure to cal file
                    save(fp_cal,'cal') % update cal file
                    disp(['Cal file updated: ',fp_cal]);
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
            LR.CHLA_ADJUSTED_QC(~t_nan) =  1;
			% CHLA_FLUORESCENCE_ADJUSTED gets the in situ dark correction and factory scale, but no factor of '2' update!
            LR.CHLA_FLUORESCENCE_ADJUSTED(~t_nan) = (lr_d(~t_nan,iChl) - ...
                CHL_DC) .* cal.CHL.ChlScale;
            LR.CHLA_FLUORESCENCE_ADJUSTED_QC(~t_nan) =  1;
			
			% NPQ NEXT (but not for CHLA_FLUORESCENCE!!)
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
            end
            LR.CHLA_ADJUSTED_ERROR(~t_nan) = ...
                abs(LR.CHLA_ADJUSTED(~t_nan) * 2);

            INFO.CHLA_SCI_CAL_EQU  = ['CHLA_ADJUSTED=CHLA/A, '...
                'NPQ corrected (Xing et al., 2012), spike profile ', ...
                'added back in'];
            INFO.CHLA_SCI_CAL_COEF = 'A=2';
            INFO.CHLA_SCI_CAL_COM  =['A is best estimate ', ...
                'from Roesler et al., 2017, doi: 10.1002/lom3.10185'];
			% NOW ADD SCI-CAL-STUFF FOR CHLA_FLUORESCENCE_ADJUSTED
			INFO.CHLA_FLUORESCENCE_SCI_CAL_EQU  = ['CHLA_FLUORESCENCE_ADJUSTED = ((FLUORESCENCE_CHLA-DARK"_CHLA)*SCALE_CHLA)'];
            INFO.CHLA_FLUORESCENCE_SCI_CAL_COEF = ['DARK"_CHLA = ',num2str(CHL_DC),', SCALE_CHLA = ',num2str(cal.CHL.ChlScale)];
            INFO.CHLA_FLUORESCENCE_SCI_CAL_COM  =['CHLA_FLUORESCENCE RT adj specified in http://dx.doi.org/10.13155/35385'];


            % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CHL');
            t_bio   = LR.CHLA ~= fv.bio;
            LR.FLUORESCENCE_CHLA_QC(t_bio) = LR.FLUORESCENCE_CHLA_QC(t_bio) ...
                * ~BSLflag + BSLflag*theflag;
            LR.CHLA_QC(t_bio) = LR.CHLA_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            LR.CHLA_FLUORESCENCE_QC(t_bio) = LR.CHLA_FLUORESCENCE_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            t_chk = t_bio & (LR.CHLA < RCR.CHL(1)|LR.CHLA > RCR.CHL(2)); % same range chk for CHLA_FLUORESCENCE (just different units)
            LR.CHLA_QC(t_chk) = 4;
            LR.FLUORESCENCE_CHLA_QC(t_chk) = 4;


            if isfield(cal.CHL, 'SWDC')
                t_bio   = LR.CHLA_ADJUSTED ~= fv.bio;
                LR.CHLA_ADJUSTED_QC(t_bio) = LR.CHLA_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                LR.CHLA_FLUORESCENCE_ADJUSTED_QC(t_bio) = LR.CHLA_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;                
                t_chk = t_bio & ...
                    (LR.CHLA_ADJUSTED < RC.CHL(1)|LR.CHLA_ADJUSTED > RC.CHL(2));
                LR.CHLA_ADJUSTED_QC(t_chk) = 4;
                LR.CHLA_FLUORESCENCE_ADJUSTED_QC(t_chk) = 4;
            end
            %-------------------------------------------------------------------------
            % CHECK FOR BAD SAMPLES ON 'BAD SAMPLE LIST'

            if yesBSAML == 1 && yesBSAMLcyc==1
                SbsIND = find(strcmp('CHL',tmpBSAML.list(:,ibsSENS)));
                if ~isempty(SbsIND)
                    TMPsbs = tmpBSAML.list(SbsIND,:);
                    singleBADs = TMPsbs(:,ibsD);
                    singleBADSflags = TMPsbs(:,ibsFL);
                    rangeBADs = TMPsbs(:,ibsDB);
                    rangeBADsflags = TMPsbs(:,ibsFL);
                    if ~isempty(singleBADs{1})
                        for i2 = 1:length(singleBADs)
                            xxtmp = find(LR.PRES == singleBADs{i2});
                            if isempty(xxtmp)
                                disp('WARNING!! CHLA PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                            else
                                if LR.CHLA~=fv.bio
                                    LR.CHLA_QC(xxtmp) = str2double(singleBADSflags{i2});
                                    LR.CHLA_QC(xxtmp) = str2double(singleBADSflags{i2});                                    
                                end
                                if LR.CHLA_ADJUSTED~=fv.bio
                                    LR.CHLA_FLUORESCENCE_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                                    LR.CHLA_FLUORESCENCE_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});                                   
                                end
                            end
                        end
                        clear i2
                    end
                    for i3 = 1:length(rangeBADs)
                        LR.CHLA_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.CHLA~=fv.bio) = str2double(rangeBADsflags{i3});
                        LR.CHLA_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.CHLA_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                        LR.CHLA_FLUORESCENCE_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.CHLA_FLUORESCENCE~=fv.bio) = str2double(rangeBADsflags{i3});
                        LR.CHLA_FLUORESCENCE_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.CHLA_FLUORESCENCE_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                    end
                    clear i3
                end
            end
        end
    end

    % ****************************************************************
    % CALCULATE CHLOROPHYLL 435 CONCENTRATION (µg/L or mg/m^3)
    %  **** 435 ... 435 ... 435 ... 435 .. 435 ... 435 ... 435 *******
    % ****************************************************************

    if (~isempty(iChl435) && master_FLBB ~= 0) || ...
            (~isempty(iChl435) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))

        % Not sure why master_FLBB variable part of the statement, seems having
        % just ~istempty(iChl) would suffice?
        % JP 03/10/20201 - There are a handful of earlier floats where the msg
        % file header indicates that an FLBB exists but there isn't one on the
        % float. For these floats the FLBB_mode = 0 from the start maybe there
        % is a better way...

        % Predim variables with fill values then adjust as needed
        LR.FLUORESCENCE_CHLA435    = fill0 + fv.bio;
        LR.FLUORESCENCE_CHLA435_QC = fill0 + fv.QC;
        LR.CHLA435                 = fill0 + fv.bio;
        LR.CHLA435_QC              = fill0 + fv.QC;
        LR.CHLA435_ADJUSTED        = fill0 + fv.bio;
        LR.CHLA435_ADJUSTED_QC     = fill0 + fv.QC;
        LR.CHLA435_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CHLA435_SCI_CAL_EQU   = 'not applicable';
        INFO.CHLA435_SCI_CAL_COEF  = 'not applicable';
        INFO.CHLA435_SCI_CAL_COM   = 'not applicable';
        INFO.CHLA435_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

        t_nan = isnan(lr_d(:,iChl435)); % NaN's in data if any

        LR.FLUORESCENCE_CHLA435(~t_nan)    = lr_d(~t_nan,iChl435);
        LR.FLUORESCENCE_CHLA435_QC(~t_nan) = fv.QC;

        if isfield(cal,'CHL435') % Sensor could be bad so maybe no cal info
            LR.CHLA435(~t_nan) = (lr_d(~t_nan,iChl435) - cal.CHL435.ChlDC) .* ...
                cal.CHL435.ChlScale;
            LR.CHLA435_QC(~t_nan) =  3; % 3 do not use w/o adjusting

            % ADJUSTED DATA BASED ON ADMT18 CONCENSUS -WILL BE UPDATED
            % WITHIN THE YEAR - jp 12/13/2017

            % 1st CHECK FOR IN SITU DARKS & TRY AND GET THEM IF NOT THERE
            if ~isfield(cal.CHL435, 'SWDC')
                %SWDC = get_CHLdark(MBARI_ID_str, dirs, 5); % 1st 5 good profiles
                %SWDC = get_CHLdark(cal.info, dirs, 5); % 1st 5 good profiles
                SWDC = get_CHLdark(cal, dirs, 5); % already extract but do it again for now
                if ~isempty(SWDC)
                    cal.CHL435.SWDC = SWDC; % add structure to cal file
                    save(fp_cal,'cal') % update cal file
                    disp(['Cal file updated: ',fp_cal]);
                end
            end

            % FIGURE OUT WHICH DC TO USE BUT ALWAYS CREATE ADJUSTED CHL
            % DATA IF RAW EXISTS
            if isfield(cal.CHL435, 'SWDC')
                CHL435_DC = cal.CHL435.SWDC.DC435;
            else
                CHL435_DC = cal.CHL435.ChlDC; % NO IN SITU DARKS YET OR EVER
            end

            LR.CHLA435_ADJUSTED(~t_nan) = (lr_d(~t_nan,iChl435) - ...
                CHL435_DC) .* cal.CHL435.ChlScale ./ 2;
            LR.CHLA435_ADJUSTED_QC(~t_nan) =  1;
            % NPQ NEXT
            NPQ_CHL = LR.CHLA435_ADJUSTED;
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
                LR.CHLA435_ADJUSTED(tNPQ & ~t_nan) = ...
                    NPQ.data(tNPQ & ~t_nan,iXing) + ...
                    NPQ.data(tNPQ & ~t_nan,iSPIKE);
                LR.CHLA435_ADJUSTED_QC(tNPQ & ~t_nan) =  5;

            end
            LR.CHLA435_ADJUSTED_ERROR(~t_nan) = ...
                abs(LR.CHLA435_ADJUSTED(~t_nan) * 2);

            INFO.CHLA435_SCI_CAL_EQU  = ['CHLA435_ADJUSTED=CHLA435/A, '...
                'NPQ corrected (Xing et al., 2012), spike profile ', ...
                'added back in'];
            INFO.CHLA435_SCI_CAL_COEF = 'A=2';
            INFO.CHLA435_SCI_CAL_COM  =['A is best estimate ', ...
                'from Roesler et al., 2017, doi: 10.1002/lom3.10185'];



            % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CHL435');
            t_bio   = LR.CHLA ~= fv.bio;
            LR.FLUORESCENCE_CHLA435_QC(t_bio) = LR.FLUORESCENCE_CHLA435_QC(t_bio) ...
                * ~BSLflag + BSLflag*theflag;
            LR.CHLA435_QC(t_bio) = LR.CHLA435_QC(t_bio) * ~BSLflag + BSLflag*theflag;

            t_chk = t_bio & (LR.CHLA435 < RCR.CHL(1)|LR.CHLA435 > RCR.CHL(2));
            LR.CHLA435_QC(t_chk) = 4;
            LR.FLUORESCENCE_CHLA435_QC(t_chk) = 4;


            if isfield(cal.CHL435, 'SWDC')
                t_bio   = LR.CHLA435_ADJUSTED ~= fv.bio;
                LR.CHLA435_ADJUSTED_QC(t_bio) = LR.CHLA435_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                t_chk = t_bio & ...
                    (LR.CHLA435_ADJUSTED < RC.CHL(1)|LR.CHLA435_ADJUSTED > RC.CHL(2));
                LR.CHLA435_ADJUSTED_QC(t_chk) = 4;
            end
            %-------------------------------------------------------------------------
            % CHECK FOR BAD SAMPLES ON 'BAD SAMPLE LIST'

            if yesBSAML == 1 && yesBSAMLcyc==1
                SbsIND = find(strcmp('CHL435',tmpBSAML.list(:,ibsSENS)));
                if ~isempty(SbsIND)
                    TMPsbs = tmpBSAML.list(SbsIND,:);
                    singleBADs = TMPsbs(:,ibsD);
                    singleBADSflags = TMPsbs(:,ibsFL);
                    rangeBADs = TMPsbs(:,ibsDB);
                    rangeBADsflags = TMPsbs(:,ibsFL);
                    if ~isempty(singleBADs{1})
                        for i2 = 1:length(singleBADs)
                            xxtmp = find(LR.PRES == singleBADs{i2});
                            if isempty(xxtmp)
                                disp('WARNING!! CHLA435 PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                            else
                                if LR.CHLA435~=fv.bio
                                    LR.CHLA435_QC(xxtmp) = str2double(singleBADSflags{i2});
                                end
                                if LR.CHLA435_ADJUSTED~=fv.bio
                                    LR.CHLA435_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                                end
                            end
                        end
                        clear i2
                    end
                    for i3 = 1:length(rangeBADs)
                        LR.CHLA435_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.CHLA435~=fv.bio) = str2double(rangeBADsflags{i3});
                        LR.CHLA435_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.CHLA435_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                    end
                    clear i3
                end
            end
        end
    end

    % ****************************************************************
    % CALCULATE PARTICLE BACKSCATTER COEFFICIENT FROM VOLUME
    % SCATTERING FUNCTION (VSF) (m^-1)
    % APEX FLBB
    % ****************************************************************
    if (~isempty(iBb) && master_FLBB ~= 0) || ....
            (~isempty(iBb) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))
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
                        THETA, LR.PSAL(BETA_SW_ind(ct)), DELTA);
                end
            end
            %LR.BETA_BACKSCATTERING700 = VSF; % (DOUBLE CHECK SPEC say counts???)
            LR.BETA_BACKSCATTERING700(~t_nan) = lr_d(~t_nan,iBb); % counts
            LR.BETA_BACKSCATTERING700_QC(~t_nan) = fv.QC; % counts
            LR.BBP700(~t_nan)    = (VSF(~t_nan) - BETA_SW(~t_nan)) * X; %b_bp m^-1
            LR.BBP700_QC(~t_nan) = 2; % 3 do not use w/o adjusting ... 6/10/21 modify qcraw flag from 3 to 2.

            % CALCULATE ADJUSTED DATA _ NO ADJUSTMENTS AT THIS TIME
            % BUT...6/10/21 START POPULATING BBP700_ADJUSTED WITH BBP
            % DIRECTLY!!!  KEEP ORIGINAL BLOCK IN CASE AN ACTUAL ADJUSTMENT
            % GETS IMPLEMENTED/STORED IN THE QC MATRIX.
            LR.BBP700_ADJUSTED(~t_nan) = LR.BBP700(~t_nan);
            LR.BBP700_ADJUSTED_QC(~t_nan) = 1;
            LR.BBP700_ADJUSTED_ERROR(~t_nan) = fv.bio; % PLACE HOLDER FOR NOW
            INFO.BBP700_SCI_CAL_EQU  = 'BBP700_ADJUSTED=BBP700';
            INFO.BBP700_SCI_CAL_COEF = [''];
            INFO.BBP700_SCI_CAL_COM  =['BBP700_ADJUSTED is being filled with BBP700 directly in real time.  Adjustment method may be enhanced in the future.'];

            %             if isfield(QC,'BB')
            %                 QCD = [LR.PRES, LR.TEMP, LR.PSAL, LR.BBP700];
            %                 LR.BBP700_ADJUSTED(~t_nan) = ...
            %                     apply_QC_corr(QCD(~t_nan,:), INFO.sdn, QC.BB);
            %                 LR.BBP700_ADJUSTED_QC(~t_nan) =  2;
            %                 LR.BBP700_ADJUSTED_ERROR(~t_nan) = fv.bio; % PLACE HOLDER FOR NOW
            %                 INFO.BBP700_SCI_CAL_EQU  = 'BBP700_ADJUSTED=BBP700*A-B';
            %                 INFO.BBP700_SCI_CAL_COEF = ['A=', ...
            %                     num2str(QC.BB.steps(3),'%0.4f'),',B=',...
            %                     num2str(QC.BB.steps(4),'%0.4f')];
            %                 INFO.BBP700_SCI_CAL_COM  =['A and B determined by comparison', ...
            %                     ' to discrete samples from post deployment calibration',...
            %                     ' rosette cast'];
            %             end
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

        clear BETA_SW VSF ct b90sw bsw

        %----------------------------------------------------------------
        % CHECK FOR BAD SAMPLES ON THE 'BAD-SAMPLE LIST'
        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('BBP',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2 = 1:length(singleBADs)
                        xxtmp = find(LR.PRES == singleBADs{i2});
                        if isempty(xxtmp)
                            disp('WARNING!! BBP700 PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if LR.BBP700~=fv.bio
                                LR.BBP700_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if LR.BBP700_ADJUSTED~=fv.bio
                                LR.BBP700_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3 = 1:length(rangeBADs)
                    LR.BBP700_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.BBP700~=fv.bio) = str2double(rangeBADsflags{i3});
                    LR.BBP700_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.BBP700_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end
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

        % CHECK FOR BAD SAMPLES ON THE 'BAD-SAMPLE LIST'
        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('CDOM',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2 = 1:length(singleBADs)
                        xxtmp = find(LR.PRES == singleBADs{i2});
                        if isempty(xxtmp)
                            disp('WARNING!! CDOM PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if LR.CDOM~=fv.bio
                                LR.CDOM_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if LR.CDOM_ADJUSTED~=fv.bio
                                LR.CDOM_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3 = 1:length(rangeBADs)
                    LR.CDOM_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.CDOM~=fv.bio) = str2double(rangeBADsflags{i3});
                    LR.CDOM_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.CDOM_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end

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
        LR.IK_PH                           = fill0 + fv.bio;
        LR.IK_PH_QC                        = fill0 + fv.QC;
        LR.VK_PH                           = fill0 + fv.bio;
        LR.VK_PH_QC                        = fill0 + fv.QC;
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
            disp(['pH data exists but no cal file for ',INST_ID_str])
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
            ipH_Ik = find(strcmp('Ik', pH_hdr) == 1); % Counter electrode current
            ipH_Vk = find(strcmp('Vk', pH_hdr) == 1); % Counter electrode V

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
            LR.IK_PH(~t_nan) = pH_data(~t_nan, ipH_Ik) * 1e9; % nano amps
            LR.IK_PH_QC(~t_nan) = fv.QC;
            LR.VK_PH(~t_nan) = pH_data(~t_nan, ipH_Vk); % nano amps
            LR.VK_PH_QC(~t_nan) = fv.QC;

            clear ipH_p ipH_t ipH_t ipH_Ib pH_data pH_hdr
        end

        if isfield(cal,'pH') && isfield(QC,'pH') && ~isempty(phtot)
            QCD = [LR.PRES(~t_nan), LR.TEMP(~t_nan), LR.PSAL(~t_nan), phtot];
            LR.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan) = ...
                apply_QC_corr(QCD, INFO.sdn, QC.pH);
            LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_nan) = 1; % set=1 9/27/16 vs
            isPOF = find(str2double(WMO)==pH_pumpoffset_980_floats);
            
            %             if ~isempty(isPOF) % JP commented out 09/08/23
            %                 try
            %                 LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR.PRES>980 & ~t_nan) = 3;
            %                 catch
            %                     keyboard
            %                 end
            %             end
            if ~isempty(isPOF) % jp 09/08/23 removed try catch test block
                LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR.PRES>980 & ~t_nan) = 3;
            end

            %LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR(~t_nan) = 0.01;
            LR.PH_IN_SITU_TOTAL_ADJUSTED(LR_inf)    = 20.1; %UNREAL #
            LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR_inf) = 4;
            step_tmpPH = find(QC.pH.steps(:,2)<=INFO.cast,1,'last');
            juld_prof = INFO.sdn-datenum(1950,01,01); %convert to JULD
            juld_init = QC.pH.steps(step_tmpPH,1)-datenum(1950,01,01); %convert to JULD
            juld_end = QC.date - datenum(1950,01,01); %date at last DMQC, converted to JULD
            LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR = 0.01 + LR.DOXY_ADJUSTED_ERROR.*0.0016;
            if juld_prof>juld_end
                LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR = LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR + 0.03.*(juld_prof-juld_end)./365;
            end
            LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR(t_nan) = fv.bio;


            INFO.PH_SCI_CAL_EQU  = ['PH_IN_SITU_TOTAL_ADJUSTED=', ...
                '[PH_IN_SITU_TOTAL+[PUMP_OFFSET - [OFFSET + DRIFT(JULD-JULD_PIVOT)/365]*TCOR]]/GAIN;',...
                'TCOR=(TREF+273.15)./(TEMP+273.15);  TREF = TEMP at 1500m.'];
            INFO.PH_SCI_CAL_COEF = ['PUMP_OFFSET = ',num2str(QC.pH.pHpumpoffset),...
                '; OFFSET = ',num2str(QC.pH.steps(step_tmpPH,4),'%6.4f'),...
                '; DRIFT = ',num2str(QC.pH.steps(step_tmpPH,5),'%6.4f'),...
                '; GAIN = ',num2str(QC.pH.steps(step_tmpPH,3),'%6.4f'),...
                '; JULD = ',num2str(juld_prof,'%9.4f'),...
                '; JULD_PIVOT = ',num2str(juld_init,'%9.4f')];
            INFO.PH_SCI_CAL_COM  =['DMQC follows '...
                'Maurer et al., 2021 (https://doi.org/10.3389/fmars.2021.683207).'];
            if regexp(MBARI_ID_str, bad_O2_filter, 'once')
                [LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR, INFO.PH_SCI_CAL_COM,~,~] = Reassign_ArgoSpecs_LIReqn8(MBARI_ID_str,INFO.cast,LR.PH_IN_SITU_TOTAL_ADJUSTED_ERROR,INFO.PH_SCI_CAL_COM,0,0);
            end

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
            if strcmpi(MBARI_ID_str, 'ua9660')
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
        tIK3  = (LR.IK_PH < RCR.IK(1) | LR.IK_PH > RCR.IK(2)).*t_diag*3; % DIAGNOSTIC
        tIB4  = (LR.IB_PH < RCR.IB(1)*25 | LR.IB_PH > RCR.IB(2)*25).*t_diag*4; % DIAGNOSTIC
        tIK4  = (LR.IK_PH < RCR.IK(1)*25 | LR.IK_PH > RCR.IK(2)*25).*t_diag*4; % DIAGNOSTIC
        % 09/30/20 TM - EXCLUSION BLOCK FOR NEWER EQPAC FLOATS EXHIBITING
        % THE ERRONEOUS RAILED DIAG VALUES
        if regexp(MBARI_ID_str, skip_ph_diag, 'once')
            % keep syntax the same to ensure proper function!
            % Set flags in this block to 0 (always good)
            disp('Excluding float from pH sensor diagnostic checks (erroneous railed values!)!')
            tIB3  = (LR.IB_PH < RCR.IB(1) | LR.IB_PH > RCR.IB(2)).*t_diag*0; % DIAGNOSTIC
            tIK3  = (LR.IK_PH < RCR.IK(1) | LR.IK_PH > RCR.IK(2)).*t_diag*0; % DIAGNOSTIC
            tIB4  = (LR.IB_PH < RCR.IB(1)*25 | LR.IB_PH > RCR.IB(2)*25).*t_diag*0; % DIAGNOSTIC
            tIK4  = (LR.IK_PH < RCR.IK(1)*25 | LR.IK_PH > RCR.IK(2)*25).*t_diag*0; % DIAGNOSTIC
        end

        tALL  = max([LR.PH_IN_SITU_TOTAL_QC, tBSL, tSTP, tRC, ...
            tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        LR.PH_IN_SITU_TOTAL_QC(t_bio) = tALL(t_bio);
        LR.PH_IN_SITU_FREE_QC(t_bio)  = tALL(t_bio);

        tALL  = max([LR.VRS_PH_QC, tBSL, tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        LR.VRS_PH_QC(t_bio) = tALL(t_bio);

        tALL  = max([LR.IB_PH_QC, tBSL, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        LR.IB_PH_QC(t_bio)  = tALL(t_bio);
        LR.IK_PH_QC(t_bio)  = tALL(t_bio);
        LR.VK_PH_QC(t_bio)  = tALL(t_bio);


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

        % 09/30/20 TM - EXCLUSION BLOCK FOR NEWER EQPAC FLOATS EXHIBITING
        % THE ERRONEOUS RAILED DIAG VALUES
        if regexp(MBARI_ID_str, skip_ph_diag, 'once')
            if ~isempty(dura.data) && ph_filesize> 10000 && ... % usually > 10000 bytes for full profile...
                    any(LR.IB_PH < RC.IB(1) | LR.IB_PH > RC.IB(2))
                disp([pH_file,' has out of range pH Ib values. pH ', ...
                    'QF''s for these samples will be set to bad or questionable'])
            end
            if ~isempty(dura.data) && ph_filesize> 10000 && ... % usually > 10000 bytes for full profile...
                    any(LR.IK_PH < RC.IK(1) | LR.IK_PH > RC.IK(2))
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

        %if INFO.cast == 30,pause,end % TESTING

        %------------------------------------------------------------------
        % OK NOW MARK BAD SAMPLES AS LISTED ON THE BAD SAMPLE LIST!
        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('PH',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2 = 1:length(singleBADs)
                        xxtmp = find(LR.PRES == singleBADs{i2});
                        if isempty(xxtmp)
                            disp('WARNING!! PH PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if LR.PH_IN_SITU_TOTAL~=fv.bio
                                LR.PH_IN_SITU_TOTAL_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if LR.PH_IN_SITU_TOTAL_ADJUSTED~=fv.bio
                                LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    LR.PH_IN_SITU_TOTAL_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.PH_IN_SITU_TOTAL~=fv.bio) = str2double(rangeBADsflags{i3});
                    LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.PH_IN_SITU_TOTAL_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end
        %------------------------------------------------------------
        %OK, HERE IS WHERE WE WILL ADD THE OVERRIDE OF DEEP DATA QCFLAGS
        %FOR PH PUMP OFFSET FLOATS THAT HAVE BEEN QC'D TO 980.  ALL QC
        %BELOW 980M SHOULD BE SET TO 'QUESTIONABLE' FOR THESE EXCEPTION
        %FLOATS.
        % if regexp(MBARI_ID_str, ph_pumpoffset_flts, 'once')
        % LR.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR.PRES>=980) = 3;
        % end

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
        INFO.NITRATE_DATA_MODE = 'R';

        % IF ISUS FILE & CAL INFO EXIST, TRY AND PROCESS IT
        if exist([dirs.temp, NO3_file],'file') && isfield(cal,'N')
            spec = parse_NO3msg([dirs.temp,NO3_file]); % return struct
            if ~isempty(Indexpsp) && str2num(cast_num) > psal_proxy_flts{Indexpsp,2}
                [Cc,Aa,Bb] = intersect(LR.PRES,spec.P); %if a no3 file is truncated, it will have different pres (and psal) length.  ie 12702.005.isus...!
                spec.S = flipud(LR.PSAL(Aa));
            end
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
            %if sum(tABS11) >0 & updatemode. Tanya
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
            if isfield(cal,'N') && isfield(QC,'N') && ~isempty(NO3_p0_kg)
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

                %LR.NITRATE_ADJUSTED_ERROR = (abs(LR.NITRATE - ...
                %    LR.NITRATE_ADJUSTED)) * 0.1 + 0.5;
                %LR.NITRATE_ADJUSTED_ERROR(t_nan) = fv.bio;
                step_tmpN = find(QC.N.steps(:,2)<=INFO.cast,1,'last');
                juld_prof = INFO.sdn-datenum(1950,01,01); %convert to JULD
                juld_init = QC.N.steps(step_tmpN,1)-datenum(1950,01,01); %convert to JULD
                juld_end = QC.date - datenum(1950,01,01); %date at last DMQC
                LR.NITRATE_ADJUSTED_ERROR = 1 + LR.DOXY_ADJUSTED_ERROR./10;
                if juld_prof>juld_end
                    LR.NITRATE_ADJUSTED_ERROR = LR.NITRATE_ADJUSTED_ERROR + 1.*(juld_prof-juld_end)./365;
                end
                LR.NITRATE_ADJUSTED_ERROR(t_nan) = fv.bio;

                INFO.NITRATE_SCI_CAL_EQU  = ['NITRATE_ADJUSTED=', ...
                    '[NITRATE-[OFFSET + DRIFT(JULD-JULD_PIVOT)/365]]/GAIN'];
                INFO.NITRATE_SCI_CAL_COEF = ['OFFSET = ', ...
                    num2str(QC.N.steps(step_tmpN,4),'%6.4f'),...
                    '; DRIFT = ',num2str(QC.N.steps(step_tmpN,5),'%6.4f'),...
                    '; GAIN = ',num2str(QC.N.steps(step_tmpN,3),'%6.4f'),...
                    '; JULD = ',num2str(juld_prof,'%9.4f'),...
                    '; JULD_PIVOT = ',num2str(juld_init,'%9.4f')];
                INFO.NITRATE_SCI_CAL_COM  =['DMQC follows '...
                'Maurer et al., 2021 (https://doi.org/10.3389/fmars.2021.683207).'];
                if regexp(MBARI_ID_str, bad_O2_filter, 'once')
                    [~, ~, LR.NITRATE_ADJUSTED_ERROR, INFO.NITRATE_SCI_CAL_COM] = Reassign_ArgoSpecs_LIReqn8(MBARI_ID_str,INFO.cast,0,0,LR.NITRATE_ADJUSTED_ERROR,INFO.NITRATE_SCI_CAL_COM);
                end
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

        %------------------------------------------------------------------
        % NOW MARK BAD SAMPLES AS LISTED ON THE BAD-SAMPLE-LIST
        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('N',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2 = 1:length(singleBADs)
                        xxtmp = find(LR.PRES == singleBADs{i2});
                        if isempty(xxtmp)
                            disp('WARNING!! PH PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if LR.NITRATE~=fv.bio
                                LR.NITRATE_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if LR.NITRATE_ADJUSTED~=fv.bio
                                LR.NITRATE_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    LR.NITRATE_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2) & LR.NITRATE~=fv.bio) = str2double(rangeBADsflags{i3});
                    LR.NITRATE_ADJUSTED_QC(LR.PRES>=rangeBADs{i3}(1) & LR.PRES<=rangeBADs{i3}(2)& LR.NITRATE_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end
    end


    % ********************************************************************
    % INCORPORATE QUALITY FLAGS FROM EXISTING ODV FILES IF THEY ARE GREATER
    % THAN ZERO. BRUTE FORCE LOOK UP. USE PRESSURE AND TEMPERATURE TO MAKE
    % MATCH AND KEEP TOLERENCES TIGHT
    %              !!! EVENTUALLY MAKE THIS A FUNCTION !!!
    % ********************************************************************
    % BUILD LIST OF ODV VARIABLES MATCHING ARGO VARIABLES

    QCvars(1,:)  = {'Temperature[°C]'    'TEMP'}; % NEW FLOATVIZ
    QCvars(2,:)  = {'Salinity[pss]'      'PSAL'};
    QCvars(3,:)  = {'Oxygen[µmol/kg]'    'DOXY'};
    QCvars(4,:)  = {'Nitrate[µmol/kg]'   'NITRATE'};
    QCvars(5,:)  = {'Chl_a[mg/m^3]'      'CHLA'};
    QCvars(6,:)  = {'b_bp700[1/m]'       'BBP700'};
    QCvars(7,:)  = {'pHinsitu[Total]'    'PH_IN_SITU_TOTAL'};
    QCvars(8,:)  = {'b_bp532[1/m]'       'BBP532'};
    QCvars(9,:)  = {'CDOM[ppb]'          'CDOM'};
    QCvars(10,:) = {'DOWN_IRRAD380[W/m^2/nm]' 'DOWN_IRRADIANCE380'};
    QCvars(11,:) = {'DOWN_IRRAD412[W/m^2/nm]' 'DOWN_IRRADIANCE412'};
    QCvars(12,:) = {'DOWN_IRRAD490[W/m^2/nm]' 'DOWN_IRRADIANCE490'};
    QCvars(13,:) = {'DOWNWELL_PAR[µmol Quanta/m^2/sec]' 'DOWNWELLING_PAR'};
    QCvars(14,:) = {'Oxygen2[µmol/kg]'   'DOXY2'};
    QCvars(15,:) = {'Chl_a_435[mg/m^3]'   'CHLA435'};


    % DO UNADJUSTED QF's FIRST
    if ~isempty(FV_data) && get_existing_QC == 1
        ifvT    = find(strcmp('Temperature[°C]',FV_data.hdr)   == 1);
        ifvP    = find(strcmp('Pressure[dbar]',FV_data.hdr)   == 1);
        ifvSTN    = find(strcmp('Station',FV_data.hdr)   == 1);
        iQF = find(strcmp(FV_data.hdr,'QF') == 1); % QF column indices
        tFV         = FV_data.data(:,ifvSTN) == INFO.cast;
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
                    ODV_QF  = [FV_cast(:,ifvP),FV_cast(:,ifvT), ...
                        FV_cast(:,indQF(QF_ct))]; % P, T & QC
                    if strcmp(QCvars{ind,2},'BBP700') == 1
                        ODV_QF(ODV_QF(:,3) == 4,3) = 2;
                        ODV_QF(ODV_QF(:,3) == 8,3) = 4;
                    else
                        ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                        ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                    end


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
            clear ifvT ifvP ifvSTN tFV FV_cast FV_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i dP t1
            clear min_dP min_dT dT
        end
    end

    % NOW DO ADJUSTED DATA QF's
    if ~isempty(FV_QCdata) && FVQC_flag == 1
        ifvT    = find(strcmp('Temperature[°C]',FV_QCdata.hdr)   == 1);
        ifvP    = find(strcmp('Pressure[dbar]',FV_QCdata.hdr)   == 1);
        ifvSTN    = find(strcmp('Station',FV_QCdata.hdr)   == 1);

        iQF = find(strcmp(FV_QCdata.hdr,'QF') == 1); % QF column indices
        tFVQC       = FV_QCdata.data(:,ifvSTN) == INFO.cast;
        FVQC_cast   = FV_QCdata.data(tFVQC,:); % get individual QC FloatViz cast
        FVQC_QF_sum = sum(FVQC_cast(:,iQF),1); % sum of columns

        if sum(FVQC_QF_sum) > 0             % ANY ODV QF FLAGS GREATER THAN ZERO?
            indQF   = iQF(FVQC_QF_sum > 0); % ODV QF colums w/ non zero flags

            for QF_ct = 1 : size(indQF,2)
                ind = strcmp(FV_QCdata.hdr{indQF(QF_ct)-1},QCvars(:,1));
                ind_ct = find(ind == 1,1);

                % FLOAVIZ VAR MATCHES LIST, GET MATCHING QF's
                if sum(ind) > 0 && (isfield(LR,[QCvars{ind,2},'_ADJUSTED']))

                    ODV_QF  = [FVQC_cast(:,ifvP), FVQC_cast(:,ifvT), ...
                        FVQC_cast(:,indQF(QF_ct))];% P&T&QC
                    if strcmp(QCvars{ind,2},'BBP700') == 1
                        ODV_QF(ODV_QF(:,3) == 4,3) = 2;
                        ODV_QF(ODV_QF(:,3) == 8,3) = 4;
                    else
                        ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                        ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                    end

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
            clear ifvT ifvP ifvSTN tFV FVQC_cast FVQC_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i
            clear dP t1 min_dP
        end
    end

    % ****************************************************************
    % NOW DEAL WITH APEX OCR FLOAT DATA EXCEPTION - SEPARATE PRES AXIS
    % ** This should probably go above the ultimate flagging check block? **
    % ****************************************************************
    if d.OCRMode == 1
        if ~isempty(regexp(MBARI_ID_str, ocr_testflt, 'once'))

            % SET UP OCR PRES
            iPocr        = find(strcmp('p', d.hr_hdr) == 1); % P
            if isempty(d.ocr_d) % CHECK FOR DATA & DEAL WITH PRES FIRST; issue instance with test float....
                disp(['No ocr axis data in message file for ', strtrim(msg_list(msg_ct,:))])
                OCR.PRES             = [];
                OCR.PRES_QC          = [];
                OCR.PRES_ADJUSTED    = [];
                OCR.PRES_ADJUSTED_QC = [];
                continue %TM dec2,2022
            else
                OCR.PRES = d.ocr_d(:,iPocr);
                tbio     = OCR.PRES ~= fv.bio; % non fill value rows
                tchk     = OCR.PRES < RC.P(1) & OCR.PRES > RC.P(2); % out of range
                fill0    = ones(size(OCR.PRES))* 0;
                OCR.PRES_QC               = fill0 + fv.QC;
                OCR.PRES_QC(tchk & tbio)  = 4;
                OCR.PRES_QC(~tchk & tbio) = 1;
                OCR.PRES_ADJUSTED         = fill0 + fv.bio;
                OCR.PRES_ADJUSTED_QC      = fill0 + fv.QC;
            end
            OCRdata = d.ocr_d; %OK pres block is done that is unique to OCR-only APEX test floats.  Now move on to process the parameters

            %%%%%%%% LG 9/3/2024 added ocr hdr which appeared to missing for current data structure format
            OCRhdr = d.ocr_hdr;

        else
            OCRdata = d.lr_d;
            OCRhdr = d.lr_hdr;
        end

        % NOW LOOP THROUGH CHANNELS & BUILD DATA FIELDS
        ocr_cal_flds = fieldnames(cal.OCR);
        % Get the CHANNEL fields
        ocr_cal_flds = ocr_cal_flds(~cellfun(@isempty,regexp(ocr_cal_flds,'^CH','once')));
        for ch_ct = 1: size(ocr_cal_flds,1)
            wl_str = cal.OCR.(ocr_cal_flds{ch_ct}).WL;
            if ~isempty(regexp(MBARI_ID_str, ocr_testflt, 'once')) %need to bring this portion in line...initial OCR test float data parser has the wavelengths hard coded in the headers...
                % iCH    = find(strcmp(['ocr',wl_str], OCRdata) == 1);
                % LG 9/3/2024 ocr data matrix was being searched for string characters that didn't exist, so the index was empty, now it retrieves correct index
                % from header (SEE ABOVE HEADER COMMENT ON L3081)
                iCH = find(strcmp(['ocr',wl_str], OCRhdr) == 1);
            else
                iCH    = find(strcmp(['Ocr[',num2str(ch_ct-1),']'], OCRhdr) == 1);
            end
            t_nan  = isnan(OCRdata(:,iCH));
            a0     = cal.OCR.(ocr_cal_flds{ch_ct}).a0; % get cal ceofs
            a1     = cal.OCR.(ocr_cal_flds{ch_ct}).a1;
            im     = cal.OCR.(ocr_cal_flds{ch_ct}).im;
            unit_conv = 0.01;

            param_str = 'DOWN_IRRADIANCE';
            raw_param_str = 'RAW_DOWNWELLING_IRRADIANCE';
            if strcmp(wl_str, 'PAR')
                param_str = 'DOWNWELLING_';
                raw_param_str = 'RAW_DOWNWELLING_';
                unit_conv = 1;
            end

            % META DATA PARAMS
            INFO.([param_str, wl_str,'_SCI_CAL_EQU'])  = 'not applicable';
            INFO.([param_str, wl_str,'_SCI_CAL_COEF']) = 'not applicable';
            INFO.([param_str, wl_str,'_SCI_CAL_COM'])  = 'not applicable';
            INFO.([param_str, wl_str,'_DATA_MODE'])    = 'R';

            % DATA PARAMS
            if isempty(OCRdata)
                OCR.([param_str, wl_str])                   = [];
                OCR.([param_str, wl_str,'_ADJUSTED'])       = [];
                OCR.([param_str, wl_str,'_ADJUSTED_ERROR']) = [];
            else
                % PREDIM WITH FILL VALUES
                OCR.([raw_param_str, wl_str])               = fill0 + fv.bio;
                OCR.([raw_param_str, wl_str,'_QC'])         = fill0 + fv.QC;
                OCR.([param_str, wl_str])                   = fill0 + fv.bio;
                OCR.([param_str, wl_str,'_QC'])             = fill0 + fv.QC;
                OCR.([param_str, wl_str,'_ADJUSTED'])       = fill0 + fv.bio;
                OCR.([param_str, wl_str,'_ADJUSTED_QC'])    = fill0 + fv.QC;
                OCR.([param_str, wl_str,'_ADJUSTED_ERROR']) = fill0 + fv.bio;

                % NOW FILL WITH REAL DATA
                raw_irr = OCRdata(~t_nan,iCH); % raw counts for param
                OCR.([raw_param_str, wl_str])(~t_nan) = raw_irr;
                OCR.([param_str, wl_str])(~t_nan) =  ...
                    (raw_irr - a0) *a1 .* im * unit_conv; %units W/m2/nm or PAR units
                t_bio = OCR.([param_str, wl_str]) ~= fv.bio; % Non fill value samples
                t_chk = OCR.([param_str, wl_str]) < RCR.(['OCR',wl_str])(1)| ...
                    OCR.([param_str, wl_str]) > RCR.(['OCR',wl_str])(2); % range check
                % set non NaN, non bad & non fill QC = 2
                OCR.([param_str, wl_str,'_QC'])(~t_nan & ~t_chk & t_bio) = 2;
                % set non NaN, and BAD & non fill QC = 4
                OCR.([param_str, wl_str,'_QC'])(~t_nan & t_chk & t_bio) = 4;
            end
        end
        %         keyboard
        % OK, now for our 6-sensor APEX, we need to insert the OCR into the LR
        % structure.  A bit cludgy since we are working within the code block
        % originally written for the OCR-only APEX (which had an odd msg file
        % format)
        if isempty(regexp(MBARI_ID_str, ocr_testflt, 'once')) %NOT the OCR-only test float
            OCRfn = fieldnames(OCR);
            for ifn = 1:length(OCRfn)
                LR.(OCRfn{ifn}) = flipud(OCR.(OCRfn{ifn})); %shallow to deep
            end
        end
        %end
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

        %PSAL PROXY?
        %         if regexp(MBARI_ID_str, psal_proxy_flts, 'once')
        if ~isempty(Indexpsp) && str2num(cast_num) > psal_proxy_flts{Indexpsp,2}
            proxydirname = ['\\atlas\Chem\ARGO_PROCESSING\DATA\ARMOR3D_PSAL_PROXY\',WMO,'\', WMO,'.', cast_num,'.mat'];
            load(proxydirname)
            HR.PSAL = HRproxy.PSAL;
            HR.PSAL_PROXY = HRproxy.PSAL; %store within the internal matfiles as well.
            HR.PSAL_ADJUSTED = HR.PSAL;
        else
            HR.PSAL = hr_d(:,iSS);
            HR.PSAL_ADJUSTED = HR.PSAL;
        end

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

        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('S',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(HR.PRES == singleBADs{i2});
                        if isempty(xxtmp)
                            disp('WARNING!! PSAL PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if HR.PSAL ~= fv.bio
                                HR.PSAL_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if HR.PSAL_ADJUSTED ~=fv.bio
                                HR.PSAL_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    HR.PSAL_QC(HR.PRES>=rangeBADs{i3}(1) & HR.PRES<=rangeBADs{i3}(2) & HR.PSAL~=fv.bio) = str2double(rangeBADsflags{i3});
                    HR.PSAL_ADJUSTED_QC(HR.PRES>=rangeBADs{i3}(1) & HR.PRES<=rangeBADs{i3}(2) & HR.PSAL_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end


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
        if ~isempty(FV_HRdata) && get_existing_QC == 1
            ifvT    = find(strcmp('Temperature[°C]',FV_HRdata.hdr)   == 1);
            ifvP    = find(strcmp('Pressure[dbar]',FV_HRdata.hdr)   == 1);
            ifvSTN    = find(strcmp('Station',FV_HRdata.hdr)   == 1);
            iQF = find(strcmp(FV_HRdata.hdr,'QF') == 1); % QF indices
            tFV         = FV_HRdata.data(:,ifvSTN) == INFO.cast;
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
                        ODV_QF  = [FV_cast(:,ifvP),FV_cast(:,ifvT), ...
                            FV_cast(:,indQF(QF_ct))]; % P, T & QC
                        if strcmp(QCvars{ind,2},'BBP700') == 1
                            ODV_QF(ODV_QF(:,3) == 4,3) = 2;
                            ODV_QF(ODV_QF(:,3) == 8,3) = 4;
                        else
                            ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                            ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                        end

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
                clear ifvT ifvP ifvSTN tFV FV_cast FV_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i dP t1
                clear min_dP min_dT dT
            end
        end

        % CHECK FOR ANY HIGH RESOLUTION ADJUSTED QUALITY FLAGS TO ADD
        if ~isempty(FV_HRQCdata) && get_existing_QC == 1
            ifvT    = find(strcmp('Temperature[°C]',FV_HRQCdata.hdr)   == 1);
            ifvP    = find(strcmp('Pressure[dbar]',FV_HRQCdata.hdr)   == 1);
            ifvSTN    = find(strcmp('Station',FV_HRQCdata.hdr)   == 1);

            iQF = find(strcmp(FV_HRQCdata.hdr,'QF') == 1); % QF indices
            tFV         = FV_HRQCdata.data(:,ifvSTN) == INFO.cast;
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
                        ODV_QF  = [FV_cast(:,ifvP),FV_cast(:,ifvT), ...
                            FV_cast(:,indQF(QF_ct))]; % P, T & QC
                        if strcmp(QCvars{ind,2},'BBP700') == 1
                            ODV_QF(ODV_QF(:,3) == 4,3) = 2;
                            ODV_QF(ODV_QF(:,3) == 8,3) = 4;
                        else
                            ODV_QF(ODV_QF(:,3) == 4,3) = 3; % CONVERT TO ARGO VALUES
                            ODV_QF(ODV_QF(:,3) == 8,3) = 4; % CONVERT TO ARGO VALUES
                        end

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
                clear ifvP ifvT ifvSTN tFV FV_cast FV_QF_sum indQF QF_ct ODV_QF ARGO_QF ct i dP t1
                clear min_dP min_dT dT
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%% END OF HIGH RES %%%%%%%%%%%%%%%%%%%%%%%

    % ****************************************************************
    %    Emily's Park Depth Cave Baby :) <3
    % ****************************************************************

    % pk_ind_chk = 1 % Emily hasn't decided if there's a need to implement a
    % shortcut with this var yet

    if isempty(d.pk_d) % CHECK FOR DATA
        disp(['No park data in message file for ', ...
            strtrim(msg_list(msg_ct,:))])
    else % if data also header variable
        pk_d = d.pk_d;

        %     if pk_ind_chk <2
        %         pk_ind_chk = pk_ind_chk+1;  % toggle on/"'if it is msg #1 or #2, otherwise:"
        iPT   = find(strcmp('t',      d.pk_hdr) == 1); % pk T

        % JP 05/31/2023 should "Fsig" & Bbsig be "FSig" & "BbSig" spot
        % check says yes
        % if two flavors then maybe strcmpi vs strcmp - I did that
        iPChl = find(strcmpi('Fsig',   d.pk_hdr) == 1); % pk CHL fluor
        if isempty(iPChl)
            iPChl = find(strcmpi('FSig[0]',   d.pk_hdr) == 1); % pk CHL fluor
        end
        iPChl435 = find(strcmpi('FSig[1]',   d.pk_hdr) == 1); % pk CHL fluor
        iPkPres  = find(strcmpi('p',   d.pk_hdr) == 1); % pk pressure
        iPBb     = find(strcmpi('Bbsig',   d.pk_hdr) == 1); % pk Bb
        iPsdn    = find(strcmpi('Date',   d.pk_hdr) == 1); % pk sdn

        %%%%%%%%%%%% SET UP GENERAL VARS
        %GET SOME ARGO CORE PARAMTERS
        TRAJ.PARK.PRES      = pk_d(:,iPkPres);
        TRAJ.PARK.SDN       = pk_d(:,iPsdn); % TM: Matlab data format.  Annie may require JULD for the Argo Traj file but we will leave that conversion to her for now for ease-of-use at MBARI.
        % MAKE AN ARRAY OF ZEROS FOR FILLING ARRAYS LATER
        pkfill0 = zeros(size(TRAJ.PARK.PRES));

        TRAJ.PARK.PRES_QC          = pkfill0 + fv.QC;
        TRAJ.PARK.PRES_ADJUSTED    = TRAJ.PARK.PRES;
        TRAJ.PARK.PRES_ADJUSTED_QC = pkfill0 + fv.QC;

        TRAJ.PARK.TEMP             = d.pk_d(:,iPT);
        TRAJ.PARK.TEMP_QC          = pkfill0 + fv.QC;
        TRAJ.PARK.TEMP_ADJUSTED    = TRAJ.PARK.TEMP;
        TRAJ.PARK.TEMP_ADJUSTED_QC = pkfill0 + fv.QC;

        % CHECK FOR BAD PRESS VALUES
        PKQF_P = TRAJ.PARK.PRES < RCR.P(1) | TRAJ.PARK.PRES > RCR.P(2);
        TRAJ.PARK.PRES_QC(PKQF_P)  = 4;  % BAD
        TRAJ.PARK.PRES_QC(~PKQF_P & TRAJ.PARK.PRES ~= fv.bio) = 1; % GOOD
        TRAJ.PARK.PRES_ADJUSTED_QC(PKQF_P)  = 4;  % BAD
        TRAJ.PARK.PRES_ADJUSTED_QC(~PKQF_P & TRAJ.PARK.PRES_ADJUSTED ~= fv.bio) = 1; % GOOD

        % CHECK FOR BAD TEMP VALUES
        %         [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'T');
        PKQF_T  = TRAJ.PARK.TEMP < RCR.T(1) | TRAJ.PARK.TEMP > RCR.T(2);
        t_bio   = TRAJ.PARK.TEMP ~= fv.bio;

        TRAJ.PARK.TEMP_QC(PKQF_T)  = 4;  % BAD
        TRAJ.PARK.TEMP_QC(~PKQF_T & TRAJ.PARK.TEMP ~= fv.bio) = 1; % GOOD
        %         TRAJ.PARK.TEMP_QC(t_bio) = TRAJ.PARK.TEMP_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        TRAJ.PARK.TEMP_ADJUSTED_QC(PKQF_T)  = 4;  % BAD
        TRAJ.PARK.TEMP_ADJUSTED_QC(~PKQF_T & TRAJ.PARK.TEMP_ADJUSTED ~= fv.bio) = 1; % GOOD
        %         TRAJ.PARK.TEMP_ADJUSTED_QC(t_bio) = TRAJ.PARK.TEMP_ADJUSTED_QC(t_bio) * ...
        %             ~BSLflag + BSLflag*theflag;
        % ******************************************************************************
        % CALCULATE park <3 CHLOROPHYLL CONCENTRATION (µg/L or mg/m^3)
        % ******************************************************************************
        if (~isempty(iPChl) && master_FLBB ~= 0) || ...
                (~isempty(iPChl) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))

            % EMILY: master_FLBB is a flag that starts as 0 and goes to 1
            % Predim variables with fill values then adjust as needed
            TRAJ.PARK.FLUORESCENCE_CHLA    = pkfill0 + fv.bio;
            TRAJ.PARK.FLUORESCENCE_CHLA_QC = pkfill0 + fv.QC;
            TRAJ.PARK.CHLA                 = pkfill0 + fv.bio;
            TRAJ.PARK.CHLA_QC              = pkfill0 + fv.QC;
            TRAJ.PARK.CHLA_ADJUSTED        = pkfill0 + fv.bio;
            TRAJ.PARK.CHLA_ADJUSTED_QC     = pkfill0 + fv.QC;
            TRAJ.PARK.CHLA_ADJUSTED_ERROR  = pkfill0 + fv.bio;
            t_nan = isnan(pk_d(:,iPChl)); % NaN's in data if any

            TRAJ.PARK.FLUORESCENCE_CHLA(~t_nan)    = pk_d(~t_nan,iPChl);
            TRAJ.PARK.FLUORESCENCE_CHLA_QC(~t_nan) = fv.QC;

            if isfield(cal,'CHL') % Sensor could be bad so maybe no cal info
                TRAJ.PARK.CHLA(~t_nan) = (pk_d(~t_nan,iPChl) - cal.CHL.ChlDC) .* ...
                    cal.CHL.ChlScale;
                %                 TRAJ.PARK.CHLA_QC(~t_nan) =  3; % 3 do not use w/o adjusting

                TRAJ.PARK.CHLA_ADJUSTED(~t_nan) = (pk_d(~t_nan,iPChl) - ...
                    CHL_DC) .* cal.CHL.ChlScale ./ 2;
                TRAJ.PARK.CHLA_ADJUSTED_QC(~t_nan) =  1; % process chl data

                TRAJ.PARK.CHLA_ADJUSTED_ERROR(~t_nan) = ...
                    abs(TRAJ.PARK.CHLA_ADJUSTED(~t_nan) * 2);

                % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
                [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CHL'); % current cycle
                [BSLflag1, theflag1] = isbadsensor(BSL, MBARI_ID_str, INFO.cast-1, 'CHL'); % previous cycle

                flagsum = BSLflag + BSLflag1;

                t_bio   = TRAJ.PARK.CHLA ~= fv.bio; % get rid of empty rows

                if flagsum == 1 % if one of the bounding message files is not bad, mark 3!
                    theflag = 3;
                    BSLflag = 1;
                end

                TRAJ.PARK.FLUORESCENCE_CHLA_QC(t_bio) = TRAJ.PARK.FLUORESCENCE_CHLA_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                TRAJ.PARK.CHLA_QC(t_bio) = TRAJ.PARK.CHLA_QC(t_bio) * ~BSLflag + BSLflag*theflag;

                t_chk = t_bio & (TRAJ.PARK.CHLA < RCR.CHL(1)|TRAJ.PARK.CHLA > RCR.CHL(2));
                TRAJ.PARK.CHLA_QC(t_chk) = 4;
                TRAJ.PARK.FLUORESCENCE_CHLA_QC(t_chk) = 4;

                t_bio   = TRAJ.PARK.CHLA_ADJUSTED ~= fv.bio;

                TRAJ.PARK.CHLA_ADJUSTED_QC(t_bio) = TRAJ.PARK.CHLA_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;

                t_chk = t_bio & ...
                    (TRAJ.PARK.CHLA_ADJUSTED < RC.CHL(1)|TRAJ.PARK.CHLA_ADJUSTED > RC.CHL(2));
                TRAJ.PARK.CHLA_ADJUSTED_QC(t_chk) = 4;

            end % end of 'if there is cal info'
        end % end of 'if there is pk chl data present'


        % ******************************************************************************
        % CALCULATE park <3 CHLOROPHYLL 435 CONCENTRATION (µg/L or mg/m^3)
        % 435 ... 435 .... 435 ... 435 ... 435 .... 435 ... 435 ... 435 .... 435 ...
        % ******************************************************************************
        if (~isempty(iPChl435) && master_FLBB ~= 0) || ...
                (~isempty(iPChl435) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))

            % EMILY: master_FLBB is a flag that starts as 0 and goes to 1
            % Predim variables with fill values then adjust as needed
            TRAJ.PARK.FLUORESCENCE_CHLA435    = pkfill0 + fv.bio;
            TRAJ.PARK.FLUORESCENCE_CHLA435_QC = pkfill0 + fv.QC;
            TRAJ.PARK.CHLA435                 = pkfill0 + fv.bio;
            TRAJ.PARK.CHLA435_QC              = pkfill0 + fv.QC;
            TRAJ.PARK.CHLA435_ADJUSTED        = pkfill0 + fv.bio;
            TRAJ.PARK.CHLA435_ADJUSTED_QC     = pkfill0 + fv.QC;
            TRAJ.PARK.CHLA435_ADJUSTED_ERROR  = pkfill0 + fv.bio;
            t_nan = isnan(pk_d(:,iPChl435)); % NaN's in data if any

            TRAJ.PARK.FLUORESCENCE_CHLA435(~t_nan)    = pk_d(~t_nan,iPChl435);
            TRAJ.PARK.FLUORESCENCE_CHLA435_QC(~t_nan) = fv.QC;

            if isfield(cal,'CHL435') % Sensor could be bad so maybe no cal info
                TRAJ.PARK.CHLA435(~t_nan) = (pk_d(~t_nan,iPChl435) - cal.CHL435.ChlDC) .* ...
                    cal.CHL435.ChlScale;
                TRAJ.PARK.CHLA435_QC(~t_nan) =  3; % 3 do not use w/o adjusting

                TRAJ.PARK.CHLA435_ADJUSTED(~t_nan) = (pk_d(~t_nan,iPChl435) - ...
                    CHL435_DC) .* cal.CHL435.ChlScale ./ 2;
                TRAJ.PARK.CHLA435_ADJUSTED_QC(~t_nan) =  1; % process chl data

                TRAJ.PARK.CHLA435_ADJUSTED_ERROR(~t_nan) = ...
                    abs(TRAJ.PARK.CHLA435_ADJUSTED(~t_nan) * 2);

                % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
                [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CHL435');

                t_bio   = TRAJ.PARK.CHLA435 ~= fv.bio; % get rid of empty rows

                flagsum = BSLflag + BSLflag1;

                if flagsum == 1 % if one of the bounding message files is not bad, mark 3!
                    theflag = 3;
                    BSLflag = 1;
                end

                TRAJ.PARK.FLUORESCENCE_CHLA_QC(t_bio) = TRAJ.PARK.FLUORESCENCE_CHLA_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                TRAJ.PARK.CHLA_QC(t_bio) = TRAJ.PARK.CHLA_QC(t_bio) * ~BSLflag + BSLflag*theflag;

                t_chk = t_bio & (TRAJ.PARK.CHLA435 < RCR.CHL(1)|TRAJ.PARK.CHLA435 > RCR.CHL(2));
                TRAJ.PARK.CHLA435_QC(t_chk) = 4;
                TRAJ.PARK.FLUORESCENCE_CHLA435_QC(t_chk) = 4;

                t_bio   = TRAJ.PARK.CHLA435_ADJUSTED ~= fv.bio;
                TRAJ.PARK.CHLA435_ADJUSTED_QC(t_bio) = TRAJ.PARK.CHLA435_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                t_chk = t_bio & ...
                    (TRAJ.PARK.CHLA435_ADJUSTED < RC.CHL(1)|TRAJ.PARK.CHLA435_ADJUSTED > RC.CHL(2));
                TRAJ.PARK.CHLA435_ADJUSTED_QC(t_chk) = 4;

            end % end of 'if there is cal info'
        end % end of 'if there is pk chl 435 data present'


        % ****************************************************************
        % CALCULATE park <3 PARTICLE BACKSCATTER COEFFICIENT FROM VOLUME
        % SCATTERING FUNCTION (VSF) (m^-1)
        % APEX FLBB
        % ****************************************************************
        if (~isempty(iPBb) && master_FLBB ~= 0) || ....
                (~isempty(iPBb) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))
            %     if ~isempty(iBb) && master_FLBB ~= 0
            pVSF                          = pkfill0 + fv.bio;
            BETA_SW                      = pkfill0 + fv.bio;
            TRAJ.PARK.BETA_BACKSCATTERING700    = pkfill0 + fv.bio;
            TRAJ.PARK.BETA_BACKSCATTERING700_QC = pkfill0 + fv.QC;
            TRAJ.PARK.BBP700                    = pkfill0 + fv.bio;
            TRAJ.PARK.BBP700_QC                 = pkfill0 + fv.QC;
            TRAJ.PARK.BBP700_ADJUSTED           = pkfill0 + fv.bio;
            TRAJ.PARK.BBP700_ADJUSTED_QC        = pkfill0 + fv.QC;
            TRAJ.PARK.BBP700_ADJUSTED_ERROR     = pkfill0 + fv.bio;

            t_nan = isnan(pk_d(:,iPBb)); % NaN's in data if any

            if isfield(cal,'BB') % Sensor could be bad so maybe no cal info
                % VOLUME SCATTERING FUNCTION (VSF) (m^-1 sr^-1)
                pVSF(~t_nan) = (pk_d(~t_nan, iPBb) - cal.BB.BetabDC) ...
                    .* cal.BB.BetabScale;
                % NOW CALCULATE PARTICLE BACKSCATTER COEFFICIENT
                %X       = 1.13*2*pi; % (Barnes and Antoine, 2014)
                %X       = 1.17*2*pi; % (email from E. Boss 24 May 2016)

                % (Proc. Bio-Argo particle backscattering at the DAC level
                % Version 1.2, July 21th 2016
                % see LR for generation of constant vars

                SALest = 33.5; % ESTIMATED SALINITY, USED IN BSW CALCULATION ONLY
                % A 4.0 difference (35.5-31.5) in salinity yields a beta SW range of only 0.175x10^-5
                disp('CAUTION: Using estimated sal value of 33.5 in BSW calculation');
                BETA_SW_ind = find(t_nan == 0); % SEAWATER, dependent on salinity, only looking for particles

                if ~isempty(BETA_SW_ind)
                    for ct = 1:size(BETA_SW_ind,1)
                        [BETA_SW(BETA_SW_ind(ct)), b90sw, bsw] = ...
                            betasw_ZHH2009(LAMBDA,pk_d(BETA_SW_ind(ct),iPT), ...
                            THETA, SALest, DELTA);
                    end
                end

                TRAJ.PARK.BETA_BACKSCATTERING700(~t_nan) = pk_d(~t_nan,iPBb); % counts
                TRAJ.PARK.BETA_BACKSCATTERING700_QC(~t_nan) = fv.QC; % counts
                TRAJ.PARK.BBP700(~t_nan)    = (pVSF(~t_nan) - BETA_SW(~t_nan)) * X; %b_bp m^-1
                TRAJ.PARK.BBP700_QC(~t_nan) = 2; % 3 do not use w/o adjusting ... 6/10/21 modify qcraw flag from 3 to 2.

                TRAJ.PARK.BBP700_ADJUSTED(~t_nan) = TRAJ.PARK.BBP700(~t_nan);
                TRAJ.PARK.BBP700_ADJUSTED_QC(~t_nan) = 1;
                TRAJ.PARK.BBP700_ADJUSTED_ERROR(~t_nan) = fv.bio; % PLACEHOLDER FOR NOW

            end

            % BSL AND RANGE CHECKS  **new**

            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'BBP'); % does this cast have a flag
            [BSLflag1, theflag1] = isbadsensor(BSL, MBARI_ID_str, INFO.cast-1, 'BBP'); % previous cycle

            t_bio   = TRAJ.PARK.BBP700 ~= fv.bio; % logic: is there park BBP data present

            flagsum = BSLflag + BSLflag1;

            if flagsum == 1 % if one of the bounding message files is not bad, mark 3!
                BSLflag = 1;
                theflag = 3;
            end

            TRAJ.PARK.BETA_BACKSCATTERING700_QC(t_bio) = ...
                TRAJ.PARK.BETA_BACKSCATTERING700_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            TRAJ.PARK.BBP700_QC(t_bio) = TRAJ.PARK.BBP700_QC(t_bio) * ~BSLflag + BSLflag*theflag;

            % Range checks
            t_chk = t_bio & ...
                (TRAJ.PARK.BBP700 < RCR.BB700(1)| TRAJ.PARK.BBP700 > RCR.BB700(2));
            TRAJ.PARK.BBP700_QC(t_chk) = 4;
            TRAJ.PARK.BETA_BACKSCATTERING700_QC(t_chk) = 4;

            t_bio = TRAJ.PARK.BBP700_ADJUSTED ~= fv.bio; % logic: is there adjusted park BBP data present

            TRAJ.PARK.BBP700_ADJUSTED_QC(t_bio) = TRAJ.PARK.BBP700_ADJUSTED_QC(t_bio) ...
                * ~BSLflag + BSLflag*theflag;

            % Range checks
            t_chk = t_bio & (TRAJ.PARK.BBP700_ADJUSTED < RC.BB700(1)| ...
                TRAJ.PARK.BBP700_ADJUSTED > RC.BB700(2));
            TRAJ.PARK.BBP700_ADJUSTED_QC(t_chk) = 4;

            clear BETA_SW X pVSF ct b90sw bsw
        end
        %     end % end of msg loop for park
    end
    %%%%%%%%%%%%%%%%%%%%% END OF PARK DATA %%%%%%%%%%%%%%%%%%%%%%%
    %-----------------------------------------------------------
    % TM 4/10/23; for floats with bad PSAL, now that we've used the ARMOR3D
    % to reprocess the float...replace PSAL with original (we want to keep
    % the original in the files...and not propagate the proxy to the user
    % files). NOTE: COMMENT THESE LINES OUT WHEN AN INTERNAL FILE FILLED WITH PROXY DATA IS DESIRED!!!!!

        LR.PSAL =  lr_d(:,iS);
        LR.PSAL_ADJUSTED = LR.PSAL;
        %[r_hr, c_hr] = size(hr_d); % high res data dimensions
		if ~isempty(d.hr_d)
        %if r_hr > 0
            HR.PSAL = hr_d(:,iSS);
            HR.PSAL_ADJUSTED = HR.PSAL;
        end
        if ~isempty(Indexpsp) && str2num(cast_num) > psal_proxy_flts{Indexpsp,2}
            %if r_hr > 0
			if ~isempty(d.hr_d)
                HR.PSAL_QC = repmat(4,size(hr_d,1),1);
                HR.PSAL_ADJUSTED_QC = repmat(4,size(hr_d,1),1);
            end
            LR.PSAL_QC = repmat(4,size(lr_d,1),1);
            LR.PSAL_ADJUSTED_QC = repmat(4,size(lr_d,1),1);
        end
        %clear r_hr
    %----------------------------------------------------------------%
    %%%  NOW ADD FINAL MODIFICATION TO PARAM_ADJ_ERROR AND COMMENT %%%
    %%%  THIS IS FOR SPECIAL-CASE FLOATS, IE ROSS-SEA FLOATS WITH  %%%
    %%%  SHALLOW QC ASSESSMENT.  NOTE O2 FAILURE                   %%%
    %%%  FLOATS THAT REQUIRE ERROR INFLATION ARE DONE THROUGH A    %%%
    %%%  DIFFERENT ROUTINE.                                        %%%
    %----------------------------------------------------------------%
%     Define_ArgoSpecs_SPECIALCASES
    [LR, HR, INFO] = Reassign_ArgoSpecs_SPECIALCASES(LR,HR,INFO,FLOATS);
    clear indexpsp

    %-------------------------------------------------------------------------------------------------------------------------------------------
    % Add LR.parameter_DATA_MODE for each parameter.  12/21/17

    %     cycdate = INFO.sdn;
    cycdate = datenum(timestamps); % date when file came in, not internal cycle date

    % OXYGEN
    if isfield(cal,'O')
        %         XEMPT_O = find(LR.DOXY ~= 99999,1); % if empty, then cycle has no data
        %         XEMPT_O = find(LR.DOXY ~= 99999 & ~isnan(LR.DOXY)); % if empty, then cycle has no data
        %         if isempty(QC) || isempty(XEMPT_O) % no adjustments have been made yet, or all data is 99999
        %             INFO.DOXY_DATA_MODE = 'R';
        if ~isempty(QC) && isfield(QC,'O')
            if cycdate > QC.date
                INFO.DOXY_DATA_MODE = 'A';
            else
                INFO.DOXY_DATA_MODE = 'D';
            end
        end
    end
    % NITRATE
    if isfield(cal,'N')
        %         XEMPT_N = find(LR.NITRATE ~= 99999,1); % if empty, then cycle has no data
        %         if isempty(QC) || isempty(XEMPT_N) % no adjustments have been made yet, or all data is 99999
        %             INFO.NITRATE_DATA_MODE = 'R';
        if ~isempty(QC) && isfield(QC,'N')
            if cycdate > QC.date
                INFO.NITRATE_DATA_MODE = 'A';
            else
                INFO.NITRATE_DATA_MODE = 'D';
            end
        end
    end
    % PH
    if isfield(cal,'pH')
        %         XEMPT_PH = find(LR.PH_IN_SITU_TOTAL ~= 99999,1); % if empty, then cycle has no data
        %         if isempty(QC) || isempty(XEMPT_PH) % no adjustments have been made yet, or all data is 99999
        %             INFO.PH_DATA_MODE = 'R';
        if ~isempty(QC) && isfield(QC,'pH')
            if cycdate > QC.date
                INFO.PH_DATA_MODE = 'A';
            else
                INFO.PH_DATA_MODE = 'D';
            end
        end
    end
    % BBP700
    if isfield(cal,'BB')
        if sum(LR.BBP700_ADJUSTED<99999)>0  % there is data for that profile --> adjustment has been made
            INFO.BBP700_DATA_MODE = 'A';
        else
            INFO.BBP700_DATA_MODE = 'R';
        end
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
            INFO.CHLA_FLUORESCENCE_DATA_MODE = 'A';
        else
            INFO.CHLA_DATA_MODE = 'R';
            INFO.CHLA_FLUORESCENCE_DATA_MODE = 'R';
        end
    end

    % CHL435
    if isfield(cal,'CHL435')
        %         if isfield(cal.CHL,'SWDC') && isfield(cal.CHL.SWDC,'DC') && sum(LR.CHLA_ADJUSTED<99999)>0  % median dark count has been quantified (and there is data for that profile) --> adjustment has been made
        if sum(LR.CHLA435_ADJUSTED<99999)>0  % median dark count has been quantified (and there is data for that profile) --> adjustment has been made
            INFO.CHLA435_DATA_MODE = 'A';
        else
            INFO.CHLA435_DATA_MODE = 'R';
        end
    end

    % ^ for the creation of B-files, add later? (for park)
    % EMILY are we adding park to the above checks? Adjusted park data?
    % getting QC from profiles? separate?
    %------------------------------------------------------------------------------------------------------------------------------
    %-------------------------------------------------------------------------------------------------------------------------------------------

    % *********************************************************************
    % SAVE THE PROFILE AS WMO_ID#.PROFILE.mat
    % THE *.mat file will contain 3 structures:
    %   LR      for low resolution data
    %   HR      for high resolution data (constant profiling)
    %   info	float info that may be of use for ARGO data stream, also
    %           used to build ODV compatible text file from 8.mat files
    % *********************************************************************

    if exist('SurfaceObs', 'var')
        TRAJ.SurfaceObs   = SurfaceObs; % [p t s phase cor_phase uM uMsat pO2]
    end
    if exist('OptodeAirCal', 'var')
        % [sdn bladderPress P p t s phase cor_phase uM uMsat pO2]
        TRAJ.OptodeAirCal   = OptodeAirCal;
    end
    %         WMO  = INFO.WMO_ID; % string
    %         if isempty(WMO)
    %             disp(['NO WMO# FOUND FOR FLOAT! CREATING TEMPORARY ', ...
    %                   'DATA DIR FOR ', INFO.UW_ID])
    %             WMO = ['NO_WMO_',INFO.UW_ID]; % CREATE TEMPORARY WMO NAME
    %         end

    save_str = [dirs.mat, WMO,'\', WMO,'.', cast_num,'.mat'];
    if regexp(MBARI_ID_str, ocr_testflt, 'once')
        save(save_str,'OCR','LR','HR','INFO');
    else
        if exist('TRAJ','var')
            save(save_str,'LR','HR','INFO','TRAJ');
        else
            save(save_str,'LR','HR','INFO');
        end
    end
    %     if msg_ct == 1
    %         copyfile(fp_cal, [dirs.mat, WMO,'\']); % copy over cal file
    %     end

end
copyfile(fp_cal, [dirs.mat, WMO,'\']); % copy over cal file
% %fprintf('\r\n')
%
% % *********************************************************************
% % CLEAN UP
if ~isempty(ls([dirs.temp,'*.msg']))
    delete([dirs.temp,'*.msg']);
end
if ~isempty(ls([dirs.temp,'*.isus']))
    delete([dirs.temp,'*.isus']);
end
if ~isempty(ls([dirs.temp,'*.dura']))
    delete([dirs.temp,'*.dura']);
end
%    tf_float.status = 1;
%end

