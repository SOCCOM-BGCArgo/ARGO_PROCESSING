function tf_float = Process_SOLO_float(MBARI_ID_str, dirs, update_str)
% ************************************************************************
% PURPOSE:
%    This function processes raw bgc data for a given SOLO float.
%    Data is ingested via pre-processed sensor-specific txt files
%    delivered from SIO (jgilson@ucsd.edu) to MBARI.  Parsers were developed
%    in collaboration with John Gilson to reach joint-desired formats.
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
%	tf_float = Process_SOLO_float(MBARI_ID_str, dirs, update_str)
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
%   get_msg_list              parse_NO3sio
%
% CHANGE LOG
% 01/25/2022 TM, Code initialization, built off of Process_APEX_float template.
% 02/16/2022 TM, Some changes to the TRAJ structure format for holding
%                inAir data.
% 03/28/2022 TM, Added in-air logical per John's code modifications
%
% 04/20/2022 TM, Added final checks for missing values on computed b-params
%           (ones affected by missing salinity from potential faulty output from the
%           interp routine.  All of this was spurred by missing psal at 10.75 db
%           level on ss0003 cycle 1.)
%
% 06/15/2022 JP, Updated msg file deleteion in temp dir before bringing new
%                new files over locally (was deleting in file in temp dir
%                & trying to delete directories)
% 07/18/2022 TM, Updated code for processing files for which no ALK file
%                exists.
%
% 8/11/22     TM: Modification to keep missing samples (if shallower than max-pres)
%
% 10/20/2022 TM, Updated code to better handle NaN and associated flagging on decimated pH diagnostics (ie cycles 1-12 on ss0001)
%
% 11/21/2022 TM, Finished brute force code for OCR implementation of bad
%                  sample list.
%
% 3/6/2023 TM, Added error trapping for case where all pressure levels are
%               missing in phy file.... 4903026.181...)
% 05/11/23 TM, Amended code for bad-sample-list (was in error for multi-line entries of same float, same cycle, different pres levels)
%
% 11/06/23 TM  -- missing no3 spectra, so wait; should get returned with
%           the next cycle.  ie solo 002.103
% 02/01/2024 TM, modifications in support of ss4003.
% 05/16/2024 TM, integrated CHLA_FLUORESCENCE variable into processing (per Argo documentation).
% ************************************************************************

% ************************************************************************
% ************************************************************************
% *** FOR TESTING ***

% MBARI_ID_str = 'ss0001';
% update_str = 'all';
% dirs =[];
%
% MBARI_ID_str = 'ss0001';
% dirs = [];
% update_str = 'all';
% ************************************************************************
% ************************************************************************


% ************************************************************************
% SET FORMATS, DEFAULT DIRS, PREDIMENSION STRUCTURE, FILTERS, BOOKKEEPING
% ************************************************************************
fclose all; % CLOSE ANY OPEN FILES
tf_float.status = 0; % Default flag for float processing success 0 = no good
tf_float.new_messages = {};
tf_float.bad_messages = {};

% THIS IS AN EXCEPTION FOR FLOAT XXXX
% exception_float = '';

% pH SENSOR THAT HAS BEEN TURNED OFF (DUE TO IE FAILURE):
bad_ALK_filter = 'ss0001|ss0004';


% FAILED OPTODE, BUT PH AND NO3 QC IS NOT YET AFFECTED: 19142 (MAY BE ADDED
% IN FUTURE AT NEXT DMQC ASSESSMENT)
bad_O2_filter = 'xxx';

% ************************************************************************
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    % user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.bfile = [user_dir,'\Documents\MATLAB\ARGO\ARGO_MBARI2AOML\'];
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    dirs.FV        = [user_dir,'FLOATVIZ\'];

    dirs.QCadj     = [user_dir,'CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\seaecho.shore.mbari.org\floats\';

    dirs.log = [user_dir,'Processing_logs\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE BAD SENSOR LIST-----------------------------------------------
% SET DEFAULTS
dirs.BSL.hdr  = [];
dirs.BSL.list = [];
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
    end
    clear tSENSOR
end
% PARSE BAD SAMPLE LIST-----------------------------------------------
bad_sample_list = parse_bad_sample_list([dirs.cal,'bad_sample_list.txt']);
iM   = find(strcmp('MBARI ID STR',bad_sample_list.hdr) == 1);
ibsSENS   = find(strcmp('SENSOR',bad_sample_list.hdr) == 1);
ibsCYC   = find(strcmp('CYCLE',bad_sample_list.hdr) == 1);
ibsD   = find(strcmp('DEPTH',bad_sample_list.hdr) == 1);
ibsDB   = find(strcmp('DEPTH BLOCKS',bad_sample_list.hdr) == 1);
ibsFL   = find(strcmp('FLAG',bad_sample_list.hdr) == 1);

% SET DEFAULTS FOR BSAML
dirs.BSAML.hdr  = [];
dirs.BSAML.list = [];
yesBSAML = 0;
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
    end
    clear tSENSOR
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ************************************************************************
% SET DATA FILL VALUES
fv.bio = 99999;
fv.QC  = 99;

% CTD_Nax = 2; % Two CTD pressure axes.
% VALID BIO ARGO RANGE CHECKS - RAW DATA [MIN MAX]
RCR.P     = [0 12000];
RCR.S     = [26 38]; % from argo parameter list
RCR.T     = [-2.5 40]; % from argo parameter list
RCR.O     = [-5 550]; % from argo parameter list
RCR.PO2   = [-5 5000]; % from argo parameter list; added to this code 5/25/21, TM: is this range right?  Should be max 500??
RCR.OP    = [10 70]; % optode phase, from argo parameter list (the range is the same for RPhase and TPhase)
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
RCR.OCR490 = [-1 3.4];
RCR.OCR443 = [-1 3.2]; 
RCR.OCR555 = [-1 100]; % crude placeholder for now
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
        cal = get_float_cals(MBARI_ID_str, dirs);
        if regexp(cal.info.WMO_ID, '^\d{7}', 'once') % WMO found - update successful!
            update_str = 'all';

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

% CHECK FOR FLOAT TYPE IF NOT SOLO EXIT
float_type = cal.info.float_type;
if strcmp(float_type,'SOLO') % SOLO
    disp([cal.info.name,' is an SOLO float'])
elseif strcmp(float_type,'NAVIS') % NAVIS
    disp(['Float ',cal.info.name, ' appears to be a NAVIS float.'])
    disp(['Processing should be done with Process_NAVIS_float.m',...
        ' instead'])
    return
elseif strcmp(float_type,'APEX') % APEX
    disp(['Float ',cal.info.name, ' appears to be a APEX float.'])
    disp(['Processing should be done with Process_APEX_float.m',...
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
% GET FILE LISTS.
% For B-SOLOs, JG provides separate txt files of raw data for each sensor type.
% get_bsolo_file_list call returns a structure of hdr & list{file name, file dir, & file date} for specified file type
% % % cal.info.msg_dir = 'C:\Users\tmaurer\Documents\GO_BGC\SOLO\OCEAN_TEST_0001\ss0001\'; %TESTING.  temporary hack.
% % % cal.info.msg_dir = 'C:\Users\tmaurer\Documents\GO_BGC\SOLO_code_dev\OCEAN_TEST_0001\ss0001\';
FileList{1} = get_bsolo_file_list(cal.info, 'ALK'); %sio using 'alk' as pH identifier
FileList{2} = get_bsolo_file_list(cal.info, 'CTD');
FileList{3} = get_bsolo_file_list(cal.info, 'DOX');
FileList{4} = get_bsolo_file_list(cal.info, 'ECO');
FileList{5} = get_bsolo_file_list(cal.info, 'NO3');
FileList{6} = get_bsolo_file_list(cal.info, 'OCR');


% CHECK FOR EMTPY CTD dir - if CTD dir is empty - CANNOT CONTINUE WITH PROCESSING.  FLOAT HASN'T SENT ANYTHING YET??

if isempty(FileList{2}.list)
    disp([MBARI_ID_str, ' may be deployed but it has not sent any *.phy',...
        ' files yet!']);
    return
end

% GIVE FILE PROCESSING A TIME LAG - ONLY PROCESS FILES > 1 HR
% SEEMS LIKE SEVERAL INTERATIONS OF INCOMPLETE FILES CAN BE TRANSMITTED
% DURING SURFACING

% TEST FOR FILE AGE LESS THEN XX HOURS OLD
% REMOVE THESE FILES FROM LIST
XXhrs = 4; %start with 4hr for solo until we get a better handle on timing on JG's side.  For UW, we finally reduced this lag to one hr.
age_limit = now - XXhrs/24;
age_limit = now; % no lag if you want an imediate process


for ifl = 1:length(FileList)
    if ~isempty(FileList{ifl}.list)
        t1 = cell2mat(FileList{ifl}.list(:,3)) < age_limit;
        FileList{ifl}.list = FileList{ifl}.list(t1,:);
    end
end

clear t1
% ************************************************************************
% LOOK FOR MOST RECENT PROFILE FILE (*.mat) AND ANY MISSING CASTS
% ONLY WANT NEW OR MISSING CASTS IN THE COPY LIST
% YOU MUST DELETE A MAT FILE TO REDO IF IT IS PARTIAL FOR A GIVEN CAST
% ************************************************************************
[last_cast, missing_casts, first_sdn] = get_last_cast(dirs.mat, ...
    cal.info.WMO_ID);
if isempty(last_cast)
    last_cast = 0; % set low for logic tests later
end
% UPDATE MODE -- ONLY LOOK FOR NEW FILES
if strcmp(update_str, 'update') && last_cast > 0
    for ifl2 = 1:length(FileList)
        if ~isempty(FileList{ifl2}.list)
            tmp   = (regexp(FileList{ifl2}.list(:,1),'(?<=\d+\_)\d{4}(?=\.)','match','once'));
            casts = str2double(tmp);
            t1    = casts > last_cast; % new cycles in file list
            t2    = ismember(casts, missing_casts); %loc of missing in file list
            FileList{ifl2}.list = FileList{ifl2}.list(t1|t2,:);
        end
    end
end
clear tmp t1 t2 casts

% ************************************************************************
% COPY FLOAT MESSAGE FILES TO LOCAL TEMPORARY DIRECTORY
% Six file groups: *.alk, *.phy, *.eco, *.no3, *.dox, *.ocr
% ************************************************************************
if isempty(FileList{2}.list)
    disp(['No float .phy files found to process for float ',MBARI_ID_str]);
    disp(['Last processed cast: ',sprintf('%03.0f',last_cast)])
    return
end

disp(['Copying all relevant BGC-SOLO files to ', dirs.temp, '  .......'])

% COPY *.MSG files
% if ~isempty(ls([dirs.temp,'*'])) % JP 06/15/22 - This tries to delete folders to
%     delete([dirs.temp,'*']) % Clear any files in temp dir (need to specify extensions??)
% end

for i6 = 1:length(FileList)
    if ~isempty(FileList{i6}.list)

        % Check for old files in temp dir & delete if need be JP 06/15/22
        % add file dletion check in loop
        fn_tmp = FileList{i6}.list{1,1}; % get a file name to get extension
        ext    = regexp(fn_tmp,'\w+$','match','once');
        if ~isempty(ls([dirs.temp,'*.', ext]))
            delete([dirs.temp,'*.', ext])
        end

        for i = 1:size(FileList{i6}.list,1)
            fp = fullfile(FileList{i6}.list{i,2}, FileList{i6}.list{i,1});
            copyfile(fp, dirs.temp);
        end
    end
end
clear fn_tmp ext
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

% CHECK FOR EXISTING WMO DIR, CREATE IF NOT THERE
if exist([dirs.mat,WMO,filesep],'dir') ~= 7
    status = mkdir([dirs.mat,WMO,filesep]);
    if status == 0
        disp(['Directory could not be created at: ', ...
            [dirs.mat,WMO,filesep]]);
        tf_float.status = 0;
        return
    end
end

% IF UPDATE STR = ALL, CLEAR MAT FILES IN WMO DIR
if strcmp(update_str,'all')
    file_chk = ls([dirs.mat,WMO,'/*.mat']);
    if ~isempty(file_chk)
        disp(['Updating all files, clearing existing *.mat files from ',...
            dirs.mat,WMO,'\ first!'])
        delete([dirs.mat,WMO,'/*.mat']);
    end
end

% ************************************************************************
% PROCESS B-SOLO FLAT FILES FOR GIVEN FLOAT
% ************************************************************************

% ***********************************************************************
% GET DIR LIST AS STRUCTURE - THIS WILL BE USED TO SEND AN EMAIL OF
% ARRIVING NEW MSG FILES
mdir = dir([dirs.temp,'*.phy']); %use ctd ".phy" flat files as the foundational file - need ctd for everything!
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
    end
    tf_float.new_messages = new_msgs(1:ct); % list of new msgs from float
    %clear mdir rr new_msgs str
end

% ***********************************************************************

msg_list   = ls([dirs.temp,'*.phy']); % get list of file names to process

phy_ind_chk = 0; % indice toggle - find indices once per float
hr_ind_chk = 0; % toggle

master_FLBB = 0; % some floats (i.e.7558) change mode, start 1 , always 1

% GET BAD SENSOR LIST FOR FLOAT IF ANY
BSL = dirs.BSL;
disp(['Processing SOLO float ' cal.info.name, '.........'])
%------------------------------------Tanya temp
if str2num(WMO) == 1902370
    M = cellstr(msg_list);
    mm = regexp(M,'\d+(?=\.phy)','once','match');
    x2 = str2num(cell2mat(mm));
    s_ind = find(x2>=6 & x2<=12);
    msg_list(s_ind,:) = []; %remove from processing loop
    tf_float.new_messages(s_ind)=[]; %remove from email notificaiton list
end
%------------------------------------END Tanya temp

for msg_ct = 1:size(msg_list,1)
    clear LR HR INFO
    msg_file = strtrim(msg_list(msg_ct,:));
    disp(['PROCESSING CTD .phy FILE ', msg_file])
    %sensor file cell structure: {file expected; string tag; cal variable}
    sensor_files(1,:) = {regexprep(msg_file,'phy','no3'),'no3','N'}; %NO3 file;
    sensor_files(2,:)  = {regexprep(msg_file,'phy','alk'),'alk','pH'}; %pH file
    sensor_files(3,:)  = {regexprep(msg_file,'phy','dox'),'dox','O'}; %O2 file
    sensor_files(4,:)  = {regexprep(msg_file,'phy','eco'),'eco','CHL'}; %ECO file; use CHL for cal file chk
    sensor_files(5,:)  = {regexprep(msg_file,'phy','ocr'),'ocr','OCR'}; %OCR file

    % find block of numbers then look ahead to see if '.phy' follows
    cast_num = regexp(msg_file,'\d+(?=\.phy)','once','match');
    if str2num(cast_num) == 0
        continue
    end

    %TM 3/2/22; issue with cycle 44 pressures returned; In touch with John
    %Gilson about this.  Breaking our code, but waiting to hear more
    %details before addressing a proper fix.  Skip over this cycle for now.
    %     if contains(cast_num,'44')
    %         disp('TEMPORARILY SKIPPING CYCLE 44 DUE TO PRESSURE PROBLEM THAT IS CURRENTLY BEING INVESTIGATED BY IDG...')
    %         continue
    %     end

    fileinfos = dir([dirs.temp,msg_file]);
    timestamps = fileinfos.date;
    timediff = now-datenum(timestamps); % use file-stamp
    %%%cycdate = datenum(timestamps); %use date when file came in, not internal cycle date (used for DATA-MODE assignments)

    %----------------------------------------------------------------------
    %- Make sure ALL B-SOLO flat files are present if float includes all 6 sensors,
    %- otherwise, end processing for this cycle (but if msg file is > 20 days old, and still not complete, then process as usual)
    % timediff = 50; % FOR MANUAL PROCESSING OVERRIDE, UNCOMMENT THIS LINE

    %%%%TM NOTE: if the B-SOLO deployments require any 'exceptions' down the road
    %(ie floats that never return any NO3 files, ie something unique happened
    %to the SUNA), then add those exceptions here (refer to
    %Process_APEX_float.m)
    Loop_CHK = 0;
    for isf = 1:length(sensor_files)

        if strcmp(sensor_files{isf,2},'alk') == 1 && ~isempty(regexp(MBARI_ID_str, bad_ALK_filter,'once'))
            continue
        end
        if (exist([dirs.temp, sensor_files{isf,1}],'file')==0 && isfield(cal,sensor_files{isf,3}) && ...
                timediff<=20)
            disp('*******************************************************')
            disp(['WARNING: .',sensor_files{isf,2},' FILE IS MISSING: ',sensor_files{isf,1}])
            disp('PROCESSING HALTED UNTIL ALL FILES ARE PRESENT.')
            disp('*******************************************************')
            Indext1 = strfind(tf_float.new_messages,msg_file);
            t1 = find(not(cellfun('isempty',Indext1)));
            if ~isempty(t1)
                tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' ',sensor_files{isf,2},' lag']}; %mark on list as in limbo
            end
            Loop_CHK = 1;
            continue
        end
    end
    if Loop_CHK == 1 % some file type is missing!  skip this cycle, so need to break out of second loop
        continue
    end

        INFO.ice_flag = 0;

    %----------------------------------------------------------------------
    %%%%% NOW PARSE DATA, CHECK FOR EMPTY RETURN, AND BUILD REVERSE INDICES.
    %%%%% ARGO WANTS SHALLOW TO DEEP
    % A)------------
    % Start with phy (ctd) file
    Dphy = parse_BGCsio([dirs.temp,msg_file], 1, 'phy');
    if (isempty(Dphy.data))% || size(Dphy.data,1)<=3) % CHECK FOR DATA.  %%%TM NOTE: may want to add additional check on length of ctd data returned?
        disp(['No CTD data in .phy file for ', ...
            strtrim(msg_list(msg_ct,:))])
        LAGstr = 'phy';
        Loop_CHK = 1;
        %         continue
    else
        phy_rows = size(Dphy.data,1);
        IX      = (phy_rows: -1 : 1)'; % build reverse indices - ARGO WANTS SHALLOW TO DEEP
        phy_d    = Dphy.data(IX,:); % ARGO WANTS SHALLOW TO DEEP
        clear B IX phy_rows
    end
    % B)------------
    % Now deal with bgc files.
    for isf2 = 1:length(sensor_files)
        if strcmp(sensor_files{isf2,2},'no3')==1
            Dno3 = parse_NO3sio([dirs.temp,sensor_files{isf2,1}], 1);
            if (isempty(Dno3.data))% || size(Dphy.data,1)<=3) % CHECK FOR DATA.  %%%TM NOTE: may want to add additional check on length of ctd data returned?
                disp(['No NO3 data in .no3 file for ', ...
                    [dirs.temp,sensor_files{isf2,1}]])
                LAGstr = sensor_files{isf2,2};
                Loop_CHK = 1;
                %                 continue
            else
                % If missing spectra...wait.  TM. 11/6/23 solo 002.103;
                % 002.111
                %                 keyboard
                tmpN = Dno3.data;
                xtmpN = find(isnan(tmpN));
                if ~isempty(xtmpN)
                    tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' ','WAIT: Missing NO3 spectra.']}; %mark on list as in limbo
                    continue
                end
                no3_rows = size(Dno3.data,1);
                IX      = (no3_rows: -1 : 1)'; % build reverse indices - ARGO WANTS SHALLOW TO DEEP
                no3_d    = Dno3.data(IX,:); % ARGO WANTS SHALLOW TO DEEP
                clear B IX no3_rows
            end
        else
            d = parse_BGCsio([dirs.temp,sensor_files{isf2,1}],1,sensor_files{isf2,2});
            if (isempty(d.data))% || size(Dphy.data,1)<=3) % CHECK FOR DATA.  %%%TM NOTE: may want to add additional check on length of ctd data returned?
                disp(['WARNING: No ',sensor_files{isf2,2},'data in file for ', ...
                    [dirs.temp,sensor_files{isf2,1}]])
                LAGstr = sensor_files{isf2,2};
                Loop_CHK = 1;
                %                 continue
                %             end
            else % if data also header variable
                bgc_rows = size(d.data,1);
                IX      = (bgc_rows: -1 : 1)'; % build reverse indices - ARGO WANTS SHALLOW TO DEEP
                tmp_d    = d.data(IX,:); % ARGO WANTS SHALLOW TO DEEP
                clear B IX bgc_rows
                %             end
                if strcmp(sensor_files{isf2,2},'dox')==1
                    dox_InAir = d.InAir.data;
                    iPhaseAir = find(strcmp('PhaseDelay',   d.InAir.hdr) == 1); % Phase optode
                    iToAir = find(strcmp('TempOpt',   d.InAir.hdr) == 1); % Temp optode
                    iPumpAir = find(strcmp('pumpIndicator',   d.InAir.hdr) == 1); % Air pump on?
                    iESecAir = find(strcmp('ElapsedSec',   d.InAir.hdr) == 1); % Time (sec) elapsed since pump on
                end
                eval([sensor_files{isf2,2},'_d = tmp_d;']); % simplified data array, with proper indexing
                eval(['D',sensor_files{isf2,2},'=d;']);
                clear tmp_d d
            end
        end
    end
    if Loop_CHK == 1 &&  isempty(regexp(MBARI_ID_str, bad_ALK_filter,'once'))% some file type is empty!  skip this cycle, so need to break out of msg file loop
        disp(['SKIPPING CYCLE PROCESSING FOR ',msg_file,' DUE TO EMPTY DATA ON ONE OR MORE PERSSURE AXES!'])
        Indext1 = strfind(tf_float.new_messages,msg_file);
        t1 = find(not(cellfun('isempty',Indext1)));
        if ~isempty(t1)
            tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' ',LAGstr,' lag']}; %mark on list as in limbo
        end
        continue
    end

    % Add ice detection true/false (as of May, 2024, there is currently no ice algorithm on solo floats...so if will default to INFO.ice_flag = 0.)
%     INFO.ice_flag = 0;
% %     if isfield(d,'ice_flag')
% %         if d.ice_flag
% %             INFO.ice_flag = 1;
% %         end
% %     end

    %----------------------------------------------------------------------
    %%%%TM NOTE:
    %%%%% SAVE THIS BLOCK FOR NOW... DIRECT USE DOESN'T MAKE SENSE SINCE
    %%%%% PHY FILE PRESSURE AXES DON'T TRANSLATE TO BGC.  PERFORM CHECKS ON
    %%%%% DOXY ETC?
    %     d = parse_APEXmsg4ARGO([dirs.temp,msg_file]);
    %     % IF MSG FILE EXIST BUT NO LR DATA REMOVE FROM NEW LIST AND
    %     % PUT IN BAD LIST. IF NOTHING GOOD LEFT, EXIT
    %         if isempty(d.lr_d) || size(d.lr_d,1)<=3
    %             t1 = ~cellfun(@isempty,regexp(tf_float.new_messages, msg_file,'once'));
    %             tf_float.bad_messages = [tf_float.bad_messages; ...
    %                 tf_float.new_messages(t1)];
    %             tf_float.new_messages(t1) =[];
    %             if isempty(tf_float.new_messages)
    %                 tf_float.status = 0; % NO VALID MESSAGES LEFT TO PROCESS
    %                 disp(['No complete float message files found to process', ...
    %                     'for float ',MBARI_ID_str]);
    %                 disp(['Last processed cast: ',sprintf('%03.0f',last_cast)])
    %                 return
    %             end
    %         end

    % If made it here -- then cycle msg_ct has all necessary B-SOLO flat
    % files.  Now, time to parse all files and check for data!

    %----------------------------------------------------------------------
    %%%%%
    % GATHER INFO FOR BUILDING ODV FILE LATER AND JUST TO MAKE LIFE
    % EASIER
    % % % THESE COMMENTED OUT INFO FIELDS ARE NOT PRESENT IN THE BSOLO FILES, SO DO NOT FILL THESE FOR NOW.  I DO NOT THINK THEY WERE ACTIVELY BEING USED FOR APEX ANYWAY...
    % % %     INFO.CpActivationP = d.CpActivationP;
    % % %     INFO.FwRev    = d.FwRev;
    % % %     INFO.CTDtype  = d.CTDtype;
    % % %     INFO.CTDsn    = d.CTDsn;
    INFO.FlbbMode = 1; %%%%TM NOTE: set this explicitly for now...keep this code here in case can be used later.
    if isfinite(INFO.FlbbMode)  && INFO.FlbbMode >0
        master_FLBB = 1;
    end

    INFO.sdn      = Dphy.EOP_sdn;  %use end-of-profile time; within the ctd phy file.
    cycdate = INFO.sdn; %for DATA mode assignment.  try using cycle date.  For apex, we use date cycle came in but since John is always reprocessing this doesn't work.
    % % % %     % CHECK FOR BAD GPS TIME 10/30/19 -jp %%% TM NOTE: keep this in here
    % % % %     % for now.  Doesn't hurt?
    % % % %     % WE ARE GETTING TIME FROM PROFILE TERMINATION TIME NOT GPS FIX TIME!!!
    % % % %     if ~isempty(first_sdn) && ~isempty(INFO.sdn)
    % % % %         dvec = datevec(INFO.sdn);
    % % % %         if INFO.sdn > first_sdn + 365*20 && dvec(1) == 2099 % 20 yrs from start?
    % % % %             disp(['GPS time for this profile is > 20 years past start ', ...
    % % % %                 '- gps week day number bug?!!'])
    % % % %             %dvec = datevec(INFO.sdn); % I think this should be moved above 2nd IF? jp 06/20/21
    % % % %             dvec(1)  = 1999; % per aoml suggestion don't quite understand jump to 2099
    % % % %             INFO.sdn = datenum(dvec) + 1024*7;
    % % % %         elseif INFO.sdn < first_sdn % bad gps time fix 10/30/19
    % % % %             disp('GPS time for this profile is unreasonable - gps week day number bug!!')
    % % % %             disp(['days since first profile = ',num2str((INFO.sdn - ...
    % % % %                 first_sdn),'%0.0f')]);
    % % % %             INFO.sdn = INFO.sdn + 1024*7;
    % % % %         end
    % % % %     end
    INFO.cast     = Dphy.cast;
    INFO.gps      = Dphy.gps;
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

    % % %     % NEED TO CORRECT GPS SDN FOR WEEK DAY ROLLOVER BUG TOO (ua9634 & ua8482)
    % % %     if ~isempty(first_sdn) && ~isempty(INFO.gps)
    % % %         dvec = datevec(INFO.gps(:,1));
    % % %         if INFO.gps(1,1) > first_sdn + 365*20 && dvec(1,1) == 2099 % 20 yrs from start?
    % % %             disp(['GPS time for this profile is > 20 years past start ', ...
    % % %                 '- gps week day number bug?!!'])
    % % %             dvec(:,1)  = 1999; % per aoml suggestion don't quite understand jump to 2099
    % % %             INFO.gps(:,1)= datenum(dvec) + 1024*7;
    % % %         elseif INFO.gps(1,1) < first_sdn % bad gps time fix 10/30/19
    % % %             disp('GPS time for this profile is unreasonable - gps week day number bug!!')
    % % %             disp(['days since first profile = ',num2str((INFO.gps(1,1) - ...
    % % %                 first_sdn),'%0.0f')]);
    % % %             INFO.gps(:,1) = INFO.gps(:,1) + 1024*7;
    % % %         end
    % % %     end

    INFO.INST_ID  = cal.info.INST_ID;
    INFO.name   = cal.info.name;
    INFO.WMO_ID = cal.info.WMO_ID;
    INFO.float_type = float_type;
    INFO.file_name = msg_file; %phy file name
    INFO.BOP_sdn = Dphy.BOP_sdn;
    INFO.EOP_sdn = Dphy.EOP_sdn;
    INFO.EOP_sdn_str = Dphy.EOP_sdn_str;
    INFO.sensors = Dphy.sensors;
    INFO.pres_axes_ct = Dphy.pres_axes_ct;

    %%%%INFO.EOT    = d.EOT; %save for now.

    %----------------------------------------------------------------------
    % DO A BUNCH OF ONE TIME TASKS - ONLY NEED TO DO THIS ONCE PER FLOAT
    % GET NUMBER OF PRESSURE AXES AND SORT OUT INDEXING!!
    % GET VARIABLE INDICES

    %if phy_ind_chk <500 % some indexing issues due to sensor file format mods...?  Just check every time while SOLO output under dev...
    phy_ind_chk = phy_ind_chk+1;  %TM try check ind twice, in case faulty msg ie 000 skips index chk

    %ctd inds: -------------------------------------------------------------
    iPhr   = find(strcmp('PRESSURE (dbar)',      Dphy.hdr) == 1); % CTD HR P
    iThr   = find(strcmp('TEMPERATURE-90 (degC)',      Dphy.hdr) == 1); % CTD HR T
    iShr   = find(strcmp('SALINITY (psu)',      Dphy.hdr) == 1); % CTD HR S
    iPlr   = find(strcmp('PRESSURE2 (dbar)',      Dphy.hdr) == 1); % CTD LR P
    % At the start...check that low-res pressure exists!  If not, spit out
    % a warning... ie ss0004.067, see John email from July 10 2023.
    %     if isempty(iPlr)
    %         Indext12 = strfind(tf_float.new_messages,msg_file);
    %         t12 = find(not(cellfun('isempty',Indext12)));
    %         tf_float.new_messages(t12) ={[char(tf_float.new_messages(t12)),' ','SKIP: No valid LR PRES']}; %mark on list as in limbo
    %         continue
    %     end
    iTlr   = find(strcmp('TEMPERATURE2-90 (degC)',      Dphy.hdr) == 1); % CTD HR T
    iSlr   = find(strcmp('SALINITY2 (psu)',      Dphy.hdr) == 1); % CTD HR S
    %pres inds: -------------------------------------------------------------
    %(these will all always be in the phy file!  TM NOTE:
    %will they also be duplicated in the bgc flat files??  If so do we
    %want to use both to perform sanity checks per Josh's suggestion?
    %If so, where should this be done?
    %
    % Use the sensors cell variable to loop through and identify each
    % BGC axis.  (maybe this little loop section isn't necessary and the pressure axis
    % indexing embedded down in the T/S interpolation loop?)
    FLTsensors = Dphy.sensors;
    ctdN = strfind(INFO.sensors,'SBE41');
    ctdN2 = find((cellfun('isempty',ctdN)));
    nBGCpres = length(ctdN2);
    CTD_Nax = length(FLTsensors)-nBGCpres;
    %         nBGCpres = length(FLTsensors)-CTD_Nax; %always 2 ctd axes?
    iPbgc = NaN(length(nBGCpres),1);
    for inP = 1:nBGCpres
        if strcmp(FLTsensors{inP+CTD_Nax},'NO3') == 1
            NO3AX = inP; %need this later on in interp block
        end
        tmpind = find(not(cellfun('isempty',strfind(Dphy.hdr,FLTsensors{inP+CTD_Nax}))));
        iPbgc(inP) = tmpind;
        clear tmpind
    end
    %sbe83 inds: -------------------------------------------------------------
    %         iPdox  = find(strcmp('PRESSURE (dbar)', Ddox.hdr) ==1); % SBE83 pres (if returned!)
    if exist('Ddox','var')
        iPhase  = find(strcmp('PHASE DELAY DOXY (uSec)',   Ddox.hdr) == 1); % SBE83 O2 phase delay
        iTo = find(strcmp('TEMP DOXY (degC)', Ddox.hdr) == 1); % SBE83 Temp (degC)
    end
    %pH inds: -------------------------------------------------------------
    %         iPpH  = find(strcmp('PRESSURE (dbar)',  Dalk.hdr) == 1); % pH pres (if returned!)
    if exist('Dalk','var')
        ipH  = find(strcmp('Vrs (V)',  Dalk.hdr) == 1); % pH pres (if returned!)
        iVk  = find(strcmp('Vk (V)',  Dalk.hdr) == 1); % pH Vk
        iIk  = find(strcmp('Ik (nA)',  Dalk.hdr) == 1); % pH Ik
        iIb  = find(strcmp('Ib (nA)',  Dalk.hdr) == 1); % pH Ib
    end
    %eco: -------------------------------------------------------------
    %         iPchl = find(strcmp('PRESSURE (dbar)', Deco.hdr) == 1); % eco Pres (if returned!)
    if exist('Deco','var')
        iChl = find(strcmp('CHLA FLUOR. (count)',   Deco.hdr) == 1); % CHL fluor
        iBb  = find(strcmp('BETA BACKSCAT700 (count)',  Deco.hdr) == 1); % Backscatter
        iCdm = find(strcmp('CDOM FLUOR. (count)',   Deco.hdr) == 1); % CDOM
    end
    %no3: -------------------------------------------------------------
    %         iPNO3 = find(strcmp('PRESSURE (dbar)',    Dno3.hdr) == 1); % Nitrate pres (if returned!)
    %         iNO3 = find(strcmp('NITRATE (mmole/kg)',    Dno3.hdr) == 1); % Onboard nitrate no longer returned in ss0002!  so this index was coming empty
    if exist('Dno3','var')
        iNO3 = find(strcmp('UV_INT_DARK (count)',    Dno3.hdr) == 1); % Use Dark header - always needed
    end
    %ocr: -------------------------------------------------------------
    %         iPocr = find(strcmp('PRESSURE (dbar)',  Docr.hdr) == 1); % OCR pres (if returned!)
    % 3 ocr irr channels.  Check wavelengths listed in cal file and prep
    % accordingly.
    if exist('Docr','var')
        for iocrch = 1:4
            eval(['twv = cal.OCR.CH0',num2str(iocrch),'.WL;']);
            if strcmp(twv,'PAR')
                iocrpar = find(strcmp('RAW_DOWNWELL_PAR (count)',  Docr.hdr) == 1);
                indocr(iocrch) = iocrpar;
            else
                iocr = find(strcmp(['RAW_DOWN_IRR',twv,' (count)'],  Docr.hdr) == 1); % OCR pres (if returned!)
                indocr(iocrch) = iocr; %eval(['iocr',num2str(iocrch),' = iocr;'])
            end
        end
    end

    % % % % % %         % GET FLOATVIZ DATA - REG AND QC & GET HR DATA TOO
    % % % % % %         % WILL BE USED TO EXTRACT QF DATA FLAGS LATER
    % % % % % %         % DATA IS BEING PULLED FROM SIROCCO
    % % % % % %         FVQC_flag   = 1;
    % % % % % %         %         mymockSirocco = 'C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING_siroccoMock\DATA\FLOATVIZ\';
    % % % % % %         %         FV_data     = get_FloatViz_data([mymockSirocco,INFO.WMO_ID,'.TXT']);
    % % % % % %         %         FV_QCdata   = get_FloatViz_data([mymockSirocco,'QC\',INFO.WMO_ID,'QC.TXT']);
    % % % % % %         %         FV_HRdata   = get_FloatViz_data([mymockSirocco,'HR\',INFO.WMO_ID,'_HR.TXT']);
    % % % % % %         %         FV_HRQCdata = get_FloatViz_data([mymockSirocco,'HRQC\',INFO.WMO_ID,'_HRQC.TXT']);
    % % % % % %         FV_data     = get_FloatViz_data(INFO.WMO_ID);
    % % % % % %         FV_QCdata   = get_FloatViz_data([INFO.WMO_ID,'QC']);
    % % % % % %         FV_HRdata   = get_FloatViz_data([INFO.WMO_ID,'_HR']);
    % % % % % %         FV_HRQCdata = get_FloatViz_data([INFO.WMO_ID,'_HRQC']);
    % % % % % %
    % % % % % %         if isempty(FV_QCdata)
    % % % % % %             FV_QCdata = FV_data;
    % % % % % %             FV_HRQCdata = FV_HRdata;
    % % % % % %             FVQC_flag = 0;
    % % % % % %         end
    % end
    INFO.bgc_pres_axes_ct = nBGCpres;

    % *************************************************************************
    % ******************************CTD DATA***********************************
    % *************************************************************************

    % *********************************************************************
    % START WITH HR CTD DATA;  DEFINE ARGO PARAMETERS----------------------
    % *********************************************************************

    if isempty(phy_d) % CHECK FOR DATA
        disp(['No ctd data in message file for ', ...
            strtrim(msg_list(msg_ct,:))])
        HR.PRES = [];
        HR.PRES_ADJUSTED = [];
        HR.PSAL = [];
        HR.PSAL_ADJUSTED = [];
        HR.TEMP = [];
        HR.TEMP_ADJUSTED = [];
        HR.NBIN_CTD = [];  %%% TM NOTE: where will we get this info for solo?
    else % if data also header variable
        HR.PRES = phy_d(:,iPhr);
        HR.PRES_ADJUSTED = HR.PRES;
        HR.PSAL = phy_d(:,iShr);
        HR.PSAL_ADJUSTED = HR.PSAL; % NO CORRECTIONS DONE FOR S OR T
        HR.TEMP = phy_d(:,iThr);
        HR.TEMP_ADJUSTED = HR.TEMP;
        %         HR.NBIN_CTD = hr_d(:,iBIN);

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
                    for i2 = 1:length(singleBADs)
                        xxtmp = find(HR.PRES == singleBADs{1});
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
                for i3 = 1:length(rangeBADs)
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

        % % % % % %         % CHECK FOR ANY HIGH RESOLUTION RAW QUALITY FLAGS TO ADD
        % % % % % %         % TM NOTE:  LOOK BACK AT PROCESS_APEX_FLOAT WHEN READY TO
        % % % % % %         % REINSERT THIS CODE BLOCK FOR GRABBING SIROCCO
        % % % % % %         % HR QUALITY FLAGS...(trying to reduce clutter for now!)
        % % % % % %         if ~isempty(FV_HRdata)
        % % % % % %         ... ; end

        clear HRQF_P HRQF_S HRQF_T

    end


    % *********************************************************************
    % NOW DEFINE LR CTD DATA;  DEFINE ARGO PARAMETERS----------------------
    % *********************************************************************
    % % %     LR.PRES      = phy_d(:,iPlr);
    % % %     tpNan = isnan(LR.PRES);
    % % %     LR.PRES = LR.PRES(~tpNan);
    % % %     % MAKE AN ARRAY OF ZEROS FOR FILLING ARRAYS LATER
    % % %     fill0  = ones(size(LR.PRES))* 0;
    % % %
    % % %     LR.PRES_QC          = fill0 + fv.QC;
    % % %     LR.PRES_ADJUSTED    = LR.PRES;
    % % %     LR.PRES_ADJUSTED_QC = fill0 + fv.QC;
    % % %
    % % %     LR.PSAL             = phy_d(~tpNan,iSlr);
    % % %     LR.PSAL_QC          = fill0 + fv.QC;
    % % %     LR.PSAL_ADJUSTED    = LR.PSAL;
    % % %     LR.PSAL_ADJUSTED_QC = fill0 + fv.QC;
    % % %
    % % %     LR.TEMP             = phy_d(~tpNan,iTlr);
    % % %     LR.TEMP_QC          = fill0 + fv.QC;
    % % %     LR.TEMP_ADJUSTED    = LR.TEMP;
    % % %     LR.TEMP_ADJUSTED_QC = fill0 + fv.QC;

    LR.PRES      = phy_d(:,iPlr);
    tpNan = isnan(LR.PRES);
    LR.PRES(tpNan) = fv.bio;
    %     LR.PRES = LR.PRES(~tpNan);
    % MAKE AN ARRAY OF ZEROS FOR FILLING ARRAYS LATER
    fill0  = ones(size(LR.PRES))* 0;
    LR.PRES_QC          = fill0 + fv.QC;
    LR.PRES_ADJUSTED    = LR.PRES;
    LR.PRES_ADJUSTED_QC = fill0 + fv.QC;

    LR.PSAL      = phy_d(:,iSlr);
    tsNan = isnan(LR.PSAL);
    LR.PSAL(tsNan) = fv.bio;
    LR.PSAL_QC          = fill0 + fv.QC;
    LR.PSAL_ADJUSTED    = LR.PSAL;
    LR.PSAL_ADJUSTED_QC = fill0 + fv.QC;

    LR.TEMP      = phy_d(:,iTlr);
    ttNan = isnan(LR.TEMP);
    LR.TEMP(ttNan) = fv.bio;
    LR.TEMP_QC          = fill0 + fv.QC;
    LR.TEMP_ADJUSTED    = LR.TEMP;
    LR.TEMP_ADJUSTED_QC = fill0 + fv.QC;

    % if INFO.cast  == 110
    %     keyboard
    % end
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
                    xxtmp = find(LR.PRES == singleBADs{1});
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

    clear LRQF_P LRQF_S LRQF_T
    % *************************************************************************
    % ***************************END CTD DATA**********************************
    % *************************************************************************


    % *************************************************************************
    % NOW DO NECESSARY T & S INTERPOLATIONS TO BGC PRESSURE AXES!
    % *************************************************************************


    % If NO pressure levels exist in HR, can't proceed!!  (ie
    % 4903026.181...)
    Anan = isnan(HR.PRES);
    if sum(Anan)==length(HR.PRES)
        disp(['Pressure is missing at all levels for HR Pressure axis!  CANNOT PROCEED for cycle ', strtrim(msg_list(msg_ct,:))])
        %now remove from the final processed list
        Indext11 = strfind(tf_float.new_messages,msg_file);
        t1 = find(not(cellfun('isempty',Indext11)));
        if ~isempty(t1)
            tf_float.new_messages(t1) ={[char(tf_float.new_messages(t1)),' ','SKIP: No valid PRES']}; %mark on list as in limbo
        end
        continue
    end
    for iB = 1:nBGCpres
        tmpKpres = phy_d(:,iPbgc(iB));
        [maxKP,indKP] = nanmax(tmpKpres);
        if isnan(maxKP) %ALL entries are NaN!!  Sensor not turned off, but no data was returned (so nan-filled file.  ie ss0001.047 no3)

            Aint = int16.empty(0,1);
            Bdub = double.empty(size(Aint));
            eval(['BGC0',num2str(iB),' = struct("PRES",Bdub,"PRES_QC",Bdub,"PRES_ADJUSTED",Bdub,"PRES_ADJUSTED_QC",Bdub,"TEMPi",Bdub,"TEMPi_QC",Bdub,"PSALi",Bdub,"PSALi_QC",Bdub);'])
            fill0=[];
        else
            kidpres = tmpKpres(1:indKP); %remove nans beyond max pres
            inanp = isnan(kidpres);
            tpkNan = kidpres(~inanp); %remove any other nans but only for the interpolation
            %         kidpres = kidpres(~tpkNan);
            fill0  = ones(size(kidpres))* 0;

            kidpres_QC          = fill0 + fv.QC;
            kidpres_ADJUSTED    = kidpres;
            kidpres_ADJUSTED_QC = fill0 + fv.QC;

            % CHECK FOR BAD PRESS VALUES
            kidQF_P = kidpres < RCR.P(1) | kidpres > RCR.P(2);
            kidpres_QC(kidQF_P)  = 4;  % BAD
            kidpres_QC(~kidQF_P & kidpres ~= fv.bio) = 1; % GOOD

            kidpres_ADJUSTED_QC(kidQF_P)  = 4;  % BAD
            kidpres_ADJUSTED_QC(~kidQF_P & kidpres_ADJUSTED ~= fv.bio) = 1; % GOOD
            if iB == NO3AX %for nitrate
                offset = 1.1; %positive offset.  ADD pressure to NO3 axis (SUNA sensor is deeper)
                disp('*********************************************')
                disp(['**NOTE: Adding 1.1m sensor offset to ',FLTsensors{iB+CTD_Nax},' pressure axis during T/S interpolation.'])
                disp('*********************************************')
            else
                offset = 0;
            end
            tol = 2;
            PST = [HR.PRES HR.PRES_QC HR.TEMP HR.TEMP_QC HR.PSAL HR.PSAL_QC];
            %             keyboard
            newPTS = interp_SOLO_ST(PST,tpkNan,offset,tol);
            % TM aug22, now take care to match up pressures (want to preserve
            % missing pressures below max pres.)
            foxSal = NaN(length(kidpres),1);
            foxTmp = NaN(length(kidpres),1);
            foxSalQ = NaN(length(kidpres),1);
            foxTmpQ = NaN(length(kidpres),1);
            [~, fox, ~] = intersect(kidpres,tpkNan);
            %             keyboard
            foxSal(fox) = newPTS.PSALi;
            foxSalQ(fox) = newPTS.PSALi_QC;
            foxTmp(fox) = newPTS.TEMPi;
            foxTmpQ(fox) = newPTS.TEMPi_QC;
            %now fill each BGC structure (one per pressure axis) with Argo
            %variables (and intermediate MBARI interpolated T,S variables)
            %         eval(['BGC0',num2str(iB),' = struct("PRES",kidpres,"PRES_QC",kidpres_QC,"PRES_ADJUSTED",kidpres_ADJUSTED,"PRES_ADJUSTED_QC",kidpres_ADJUSTED_QC,"TEMPi",newPTS.TEMPi,"TEMPi_QC",newPTS.TEMPi_QC,"PSALi",newPTS.PSALi,"PSALi_QC",newPTS.PSALi_QC);'])
            eval(['BGC0',num2str(iB),' = struct("PRES",kidpres,"PRES_QC",kidpres_QC,"PRES_ADJUSTED",kidpres_ADJUSTED,"PRES_ADJUSTED_QC",kidpres_ADJUSTED_QC,"TEMPi",foxTmp,"TEMPi_QC",foxTmpQ,"PSALi",foxSal,"PSALi_QC",foxSalQ);'])
        end
        clear tpkNan fill0
    end
    %     save('temp.mat')

    % *************************************************************************
    % NOW PROCESS EACH BGC PARAMETER IN TURN.  AT THE BEGINNING, PERFORM A
    % CHECK ON WHICH BGC AXIS THE PARAMTER BELONGS TO, AND PROCEED ACCORDINGLY.
    % *************************************************************************

    if isfield(cal,'O')
        BGCIND = find(not(cellfun('isempty',strfind(FLTsensors,'DOX'))))-CTD_Nax; %#ok<STRCL1>
        eval(['DOXarray = BGC0',num2str(BGCIND),';']); %to reduce the number of eval calls, rename to temporary structure until all desired fields are populated.
        fill0  = ones(size(DOXarray.PRES))* 0;
        % CALCULATE POTENTIAL DENSITY (STILL NEED TO UPDATE TO TEOS 10 !)
        potT   = theta(DOXarray.PRES, DOXarray.TEMPi, DOXarray.PSALi,0); % potential temp, p=0
        % USING POTENTIAL DENSITY ( per: Processing Argo OXYGEN data at the
        % DAC level Version 2.0 October 1st 2015, Chap 5)
        DOX_den = density(DOXarray.PSALi, potT); % potential density, kg /m^3

        % ****************************************************************
        % CALCULATE OXYGEN CONCENTRATION
        % ****************************************************************
        if ~isempty(iPhase) % Phase data exists so O2 data should exist.  TM NOTE: should add a check on ~isempty(PSALi)?
            t_nan = isnan(dox_d(:,iPhase));
            O2_chk = dox_d(:,iPhase) < RCR.OP(1) | ...
                dox_d(:,iPhase) > RCR.OP(2) | ...
                dox_d(:,iTo) < RCR.OT(1) | ...
                dox_d(:,iTo) > RCR.OT(2); % Flag phase, Temp out of range
            if any(O2_chk)
                disp(['Out of range phase or optode T detected for ', ...
                    msg_file])
            end

            if strcmp(cal.O.type,'SBE83')
                DOXarray.PHASE_DELAY_DOXY            = fill0 + fv.bio;
                DOXarray.PHASE_DELAY_DOXY_QC         = fill0 + fv.QC;
                DOXarray.PHASE_DELAY_DOXY(~t_nan)    = dox_d(~t_nan,iPhase);
                DOXarray.PHASE_DELAY_DOXY_QC(~t_nan) =  fv.QC;
                %Initialize In-Air data variables:
                zero_fill = dox_InAir(:,iPhaseAir)* 0; % make array of zeros
                TRAJ.InAir.RAW = dox_InAir;
                TRAJ.InAir.TEMP_DOXY = zero_fill+fv.bio;
                TRAJ.InAir.PHASE_DELAY_DOXY = zero_fill+fv.bio;
                TRAJ.InAir.PPOX_DOXY = zero_fill+fv.bio;
                TRAJ.InAir.TEMP_DOXY_QC = zero_fill+fv.QC;
                TRAJ.InAir.PHASE_DELAY_DOXY_QC = zero_fill+fv.QC;
                TRAJ.InAir.PPOX_DOXY_QC = zero_fill+fv.QC;
                TRAJ.InAir.IN_AIR_LOGICAL = zero_fill+1;
                TRAJ.InAir.IN_AIR_LOGICAL = dox_InAir(:,1); %TM 3/28/22 add in-air logical per John's code modifications
                if ~isempty(dox_InAir) & sum(dox_InAir(:)) ~= 0 % 0 means ice detection on
                    myadata = [zero_fill dox_InAir(:,iToAir) zero_fill dox_InAir(:,iPhaseAir) dox_InAir(:,iToAir)]; %use optode T here for CTD T input?
                    [ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myadata, cal.O, 0); % jp 02/0924 looks like all SOLO's report Optode T  in deg C
                    %[ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myadata, cal.O,cal.info.float_type);
                    TRAJ.InAir.TEMP_DOXY = dox_InAir(:,iToAir);
                    TRAJ.InAir.PHASE_DELAY_DOXY = dox_InAir(:,iPhaseAir);
                    TRAJ.InAir.PPOX_DOXY = ppoxdoxy;
                    TRAJ.InAir.PH2O = pH2O;
                    po2BAD = TRAJ.InAir.PPOX_DOXY < RCR.PO2(1) | TRAJ.InAir.PPOX_DOXY > RCR.PO2(2);
                    tempBAD = TRAJ.InAir.TEMP_DOXY < RCR.OT(1) | TRAJ.InAir.TEMP_DOXY > RCR.OT(2);
                    tphaseBAD = TRAJ.InAir.PHASE_DELAY_DOXY < RCR.OP(1) | TRAJ.InAir.PHASE_DELAY_DOXY > RCR.OP(2);
                    TRAJ.InAir.PPOX_DOXY_QC(po2BAD)=4;
                    TRAJ.InAir.PPOX_DOXY_QC(~po2BAD)=1;
                    TRAJ.InAir.TEMP_DOXY_QC(tempBAD)=4;
                    TRAJ.InAir.TEMP_DOXY_QC(~tempBAD)=1;
                    TRAJ.InAir.PHASE_DELAY_DOXY_QC(tphaseBAD)=4;
                    TRAJ.InAir.PHASE_DELAY_DOXY_QC(~tphaseBAD)=1;
                    clear O2 po2BAD tempBAD tphaseBAD
                end

            else
                disp('SBE83 optode not detected!!! Error in SOLO oxygen processing...')
            end
            DOXarray.DOXY                = fill0 + fv.bio;
            DOXarray.DOXY_QC             = fill0 + fv.QC;
            DOXarray.TEMP_DOXY           = fill0 + fv.bio;
            DOXarray.TEMP_DOXY_QC        = fill0 + fv.QC;
            DOXarray.DOXY_ADJUSTED       = fill0 + fv.bio;
            DOXarray.DOXY_ADJUSTED_QC    = fill0 + fv.QC;
            DOXarray.DOXY_ADJUSTED_ERROR = fill0 + fv.bio;
            INFO.DOXY_SCI_CAL_EQU  = 'not applicable';
            INFO.DOXY_SCI_CAL_COEF = 'not applicable';
            INFO.DOXY_SCI_CAL_COM  = 'not applicable';
            INFO.DOXY_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
            if strcmp(cal.O.type,'SBE83')
                myOdata = [DOXarray.PRES(~t_nan) DOXarray.TEMPi(~t_nan) DOXarray.PSALi(~t_nan) dox_d((~t_nan),iPhase),dox_d((~t_nan),iTo)];
                [ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myOdata, cal.O, 0); % jp 02/0924 looks like all SOLO's report Optode T  in deg C
                %[ppoxdoxy, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(myOdata, cal.O,cal.info.float_type);
                DOXarray.DOXY(~t_nan) = O2_uM ./ DOX_den(~t_nan) .* 1000; % umol/kg
            else
                disp('SBE83 not detected!!! Error in SOLO oxygen processing!')
            end
            tDOXY = abs(DOXarray.DOXY) > crazy_val & ~t_nan; % Unrealistic bad value
            DOXarray.DOXY(tDOXY) = crazy_val; % SET TO crazy bad value
            % find any final NaN in computed doxy (ie if phase value good
            % but computed O2 missing due to NaN in psal or temp used?)
            t_nan_DOX = isnan(DOXarray.DOXY);
            DOXarray.DOXY(t_nan_DOX) = fv.bio;
            DOXarray.DOXY_QC(~t_nan & ~t_nan_DOX)      = 3;
            % JP: QC assignment below not needed - crazy value greater than range limits
            %LR.DOXY_QC(tDOXY) = 4; % set crazy bad value QF to 4
            DOXarray.TEMP_DOXY(~t_nan)    = dox_d(~t_nan,iTo);
            DOXarray.TEMP_DOXY_QC(~t_nan) = fv.QC;
            % Save O2sol, pO2?
            clear O2


            if isfield(QC,'O')
                % !! ONE TIME GAIN CORRECTION ONLY !!
                %             LR.DOXY_ADJUSTED(~t_nan)  = LR.DOXY(~t_nan) .* QC.O.steps(1,3);
                QCD = [DOXarray.PRES(~t_nan & ~t_nan_DOX), DOXarray.TEMPi(~t_nan & ~t_nan_DOX), DOXarray.PSALi(~t_nan & ~t_nan_DOX), DOXarray.DOXY(~t_nan & ~t_nan_DOX)];
                DOXarray.DOXY_ADJUSTED(~t_nan & ~t_nan_DOX) = ...
                    apply_QC_corr(QCD, INFO.sdn, QC.O);
                tDOXY_ADJ = abs(DOXarray.DOXY_ADJUSTED) > crazy_val & ~t_nan & ~t_nan_DOX; % Unrealistic bad value
                DOXarray.DOXY_ADJUSTED(tDOXY_ADJ) = crazy_val; % SET TO crazy bad value
                DOXarray.DOXY_ADJUSTED_QC(~t_nan & ~t_nan_DOX) = 1; % set=1 9/27/16 vs 2 = probably good
                % JP: QC assignment below not needed - crazy value greater than range limits
                %LR.DOXY_ADJUSTED_QC(tDOXY_ADJ) = 4; %set crazy val QF to bad
                %LR.DOXY_ADJUSTED_ERROR(~t_nan) = LR.DOXY_ADJUSTED(~t_nan) * 0.01;
                INFO.DOXY_SCI_CAL_EQU  = 'DOXY_ADJUSTED=DOXY*G; G = G_INIT + G_DRIFT*(JULD_PROF - JULD_INIT)/365';
                % QC matrix entry relevant to current cycle.
                steptmp = find(QC.O.steps(:,2)<=INFO.cast,1,'last');
                juld_prof = INFO.sdn-datenum(1950,01,01); %convert to JULD
                juld_init = QC.O.steps(steptmp,1)-datenum(1950,01,01); %convert to JULD
                juld_end = QC.date; %date at last DMQC
                if juld_prof<juld_end
                    if ~isempty(dox_InAir)
                        DOXarray.DOXY_ADJUSTED_ERROR(~t_nan & ~t_nan_DOX) = convert_O2mb_error_to_conc(DOXarray.TEMPi(~t_nan & ~t_nan_DOX),DOXarray.PSALi(~t_nan & ~t_nan_DOX),2); %2mbar error spec for air-calibrated
                        %                         LAST_DOXY_ERROR = DOXarray.DOXY_ADJUSTED_ERROR;
                    else
                        DOXarray.DOXY_ADJUSTED_ERROR(~t_nan & ~t_nan_DOX) = convert_O2mb_error_to_conc(DOXarray.TEMPi(~t_nan & ~t_nan_DOX),DOXarray.PSALi(~t_nan & ~t_nan_DOX),5); %4-6mbar error spec for air-calibrated
                        %                         LAST_DOXY_ERROR = DOXarray.DOXY_ADJUSTED_ERROR;
                    end
                else
                    extra_error_ppox = 1.*(juld_prof-juld_end)./365; %1 mb per year error inflation per Argo rec
                    O2error = convert_O2mb_error_to_conc(DOXarray.TEMPi(~t_nan & ~t_nan_DOX),DOXarray.PSALi(~t_nan & ~t_nan_DOX),extra_error_ppox);
                    DOXarray.DOXY_ADJUSTED_ERROR(~t_nan & ~t_nan_DOX) = DOXarray.DOXY_ADJUSTED_ERROR(~t_nan & ~t_nan_DOX)+O2error;
                end

                INFO.DOXY_SCI_CAL_COEF = ['G_INIT = ', ...
                    num2str(QC.O.steps(steptmp,3),'%0.4f'),...
                    '; G_DRIFT = ',num2str(QC.O.steps(steptmp,5),'%0.4f'),...
                    '; JULD_PROF = ',num2str(juld_prof,'%9.4f'),...
                    '; JULD_INIT = ',num2str(juld_init,'%9.4f')];

                if ~isempty(Ddox.InAir)
                    INFO.DOXY_SCI_CAL_COM  = ['G determined from ',...
                        'float  measurements in air. See Johnson et al.,2015,', ...
                        'doi:10.1175/JTECH-D-15-0101.1'];
                else
                    INFO.DOXY_SCI_CAL_COM  = ['G determined by surface' ...
                        ' measurement comparison to World Ocean Atlas 2009.', ...
                        'See Takeshita et al.2013,doi:10.1002/jgrc.20399'];
                end
            end


            % DO SOME FINAL CHECK ON VALUES, IF BAD SET QF = 4
            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'O');
            t_bio = DOXarray.DOXY ~= fv.bio; % Non fill value samples
            tST   = DOXarray.PSALi_QC == 4 | DOXarray.TEMPi_QC == 4 | DOXarray.PRES_QC == 4; % Bad S or T will affect O2
            t_chk = DOXarray.DOXY < RCR.O(1)| DOXarray.DOXY > RCR.O(2); % range check
            t_chk = t_chk | O2_chk | tST; % BAD O2 T or phase

            t_chk = t_chk & t_bio; % & not MVI either

            if isfield(DOXarray,'TPHASE_DOXY')
                DOXarray.TPHASE_DOXY_QC(t_bio) = DOXarray.TPHASE_DOXY_QC(t_bio) *  ...
                    ~BSLflag + BSLflag*theflag;
                DOXarray.TPHASE_DOXY_QC(t_chk) = 4;
            end

            % VERY BAD PSAL & TEMP = BAD O2 TOO
            DOXarray.DOXY_QC(t_bio) = DOXarray.DOXY_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            DOXarray.DOXY_QC(t_chk) = 4;

            if isfield(QC,'O')
                t_bio = DOXarray.DOXY_ADJUSTED ~= fv.bio;
                %             tST   = DOXarray.PSALi_ADJUSTED_QC == 4 | DOXarray.TEMPi_ADJUSTED_QC == 4 ...
                %                 | DOXarray.PRESi_ADJUSTED_QC == 4;
                t_chk = DOXarray.DOXY_ADJUSTED < RC.O(1)| ...
                    DOXarray.DOXY_ADJUSTED > RC.O(2) | O2_chk;
                t_chk = (t_chk | tST) & t_bio;

                DOXarray.DOXY_ADJUSTED_QC(t_bio) = DOXarray.DOXY_ADJUSTED_QC(t_bio) * ...
                    ~BSLflag + BSLflag*theflag;
                DOXarray.DOXY_ADJUSTED_QC(t_chk) = 4;
            end

            if yesBSAML == 1 && yesBSAMLcyc==1
                SbsIND = find(strcmp('O',tmpBSAML.list(:,ibsSENS)));
                if ~isempty(SbsIND)
                    TMPsbs = tmpBSAML.list(SbsIND,:);
                    singleBADs = TMPsbs(:,ibsD);
                    singleBADSflags = TMPsbs(:,ibsFL);
                    rangeBADs = TMPsbs(:,ibsDB);
                    rangeBADsflags = TMPsbs(:,ibsFL);
                    if ~isempty(singleBADs{1})
                        for i2=1:length(singleBADs)
                            xxtmp = find(DOXarray.PRES == singleBADs{1});
                            if isempty(xxtmp)
                                disp('WARNING!! DOXY PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                            else
                                if DOXarray.DOXY~=fv.bio
                                    DOXarray.DOXY_QC(xxtmp) = str2double(singleBADSflags{i2});
                                end
                                if DOXarray.DOXY_ADJUSTED~=fv.bio
                                    DOXarray.DOXY_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                                end
                            end
                        end
                        clear i2
                    end
                    for i3=1:length(rangeBADs)
                        DOXarray.DOXY_QC(DOXarray.PRES>=rangeBADs{i3}(1) & DOXarray.PRES<=rangeBADs{i3}(2) & DOXarray.DOXY~=fv.bio) = str2double(rangeBADsflags{i3});
                        DOXarray.DOXY_ADJUSTED_QC(DOXarray.PRES>=rangeBADs{i3}(1) & DOXarray.PRES<=rangeBADs{i3}(2)& DOXarray.DOXY_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                    end
                    clear i3
                end
            end


            % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            % DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
            % (BAD).  BGC_spiketest screens for nans and fill values and pre-identified bad values.
            %
            % RUN TEST ON RAW DOXY
            % however....if Arctic float, do not perform O2 spiketest
            %(very large gradient near surface...does not work well in this region!)
            %             if regexp(MBARI_ID_str ,no_o2_spike,'once')
            %                 disp([MBARI_ID_str,' is on the no DOXY spike test list'])
            %             else
            QCscreen_O = DOXarray.DOXY_QC == 4; % screen for BAD data already assigned.
            [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str, ...
                str2double(cast_num),[DOXarray.PRES DOXarray.DOXY],'O2',dirs.cal,fv.bio,QCscreen_O);
            if ~isempty(spike_inds)
                DOXarray.DOXY_QC(spike_inds) = quality_flags;
                disp(['DOXY QUALITY FLAGS ADJUSTED FOR IDENTIFIED ', ...
                    'SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
            end
            %
            % RUN TEST ON QC DOXY_ADJUSTED
            QCscreen_Oadj = DOXarray.DOXY_ADJUSTED_QC == 4; % screen for BAD data already assigned.
            [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str, ...
                str2double(cast_num),[DOXarray.PRES DOXarray.DOXY_ADJUSTED],'O2',...
                dirs.cal,fv.bio,QCscreen_Oadj);
            if ~isempty(spike_inds)
                DOXarray.DOXY_ADJUSTED_QC(spike_inds) = quality_flags;
                disp(['DOXY_ADJUSTED QUALITY FLAGS ADJUSTED FOR ', ...
                    'IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
            end
            % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            clear QCscreen_O QCscreenOadj spike_inds quality_flags
            %             end

            % OXYGEN DATA-MODE
            XEMPT_O = find(DOXarray.DOXY ~= 99999,1); % if empty, then cycle has no data
            %             if isempty(QC) || isempty(XEMPT_O) % no adjustments have been made yet, or all data is 99999
            %                 INFO.DOXY_DATA_MODE = 'R';
            if ~isempty(QC) && isfield(QC,'O') && ~isempty(XEMPT_O)
                %                 if cycdate > QC.date

                if str2num(cast_num) > QC.N.steps(end,2) %us last cycle assessed in Nitrate for O2 as well (O2 cycle listed in qc matrix is typically '1' if avg gain used)
                    INFO.DOXY_DATA_MODE = 'A';
                else
                    INFO.DOXY_DATA_MODE = 'D';
                end
            else
                INFO.DOXY_DATA_MODE = 'R';
                INFO.DOXY_SCI_CAL_EQU  = 'not applicable';
                INFO.DOXY_SCI_CAL_COEF = 'not applicable';
                INFO.DOXY_SCI_CAL_COM  = 'not applicable';
            end
            eval(['BGC0',num2str(BGCIND),' = DOXarray;']); %reassign array structure name.
            clear BGCIND fill0
        end
    end

    % ****************************************************************
    % CALCULATE CHLOROPHYLL CONCENTRATION (ug/L or mg/m^3)
    % ****************************************************************

    if (~isempty(iChl) && master_FLBB ~= 0) || ...
            (~isempty(iChl) && ~isempty(regexp(MBARI_ID_str, bad_chl_filter,'once')))
        BGCIND = find(not(cellfun('isempty',strfind(FLTsensors,'ECO'))))-CTD_Nax; %#ok<STRCL1>
        eval(['ECOarray = BGC0',num2str(BGCIND),';']); %to reduce the number of eval calls, rename to temporary structure until all desired fields are populated.
        fill0  = ones(size(ECOarray.PRES))* 0;
        % Not sure why master_FLBB variable part of the statement, seems having
        % just ~istempty(iChl) would suffice?
        % JP 03/10/20201 - There are a handful of earlier floats where the msg
        % file header indicates that an FLBB exists but there isn't one on the
        % float. For these floats the FLBB_mode = 0 from the start maybe there
        % is a better way...

        % Predim variables with fill values then adjust as needed
        ECOarray.FLUORESCENCE_CHLA    = fill0 + fv.bio;
        ECOarray.FLUORESCENCE_CHLA_QC = fill0 + fv.QC;
        ECOarray.CHLA                 = fill0 + fv.bio;
        ECOarray.CHLA_QC              = fill0 + fv.QC;
        ECOarray.CHLA_ADJUSTED        = fill0 + fv.bio;
        ECOarray.CHLA_ADJUSTED_QC     = fill0 + fv.QC;
        ECOarray.CHLA_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CHLA_SCI_CAL_EQU   = 'not applicable';
        INFO.CHLA_SCI_CAL_COEF  = 'not applicable';
        INFO.CHLA_SCI_CAL_COM   = 'not applicable';
        INFO.CHLA_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)
		ECOarray.CHLA_FLUORESCENCE                 = fill0 + fv.bio;
        ECOarray.CHLA_FLUORESCENCE_QC              = fill0 + fv.QC;
        ECOarray.CHLA_FLUORESCENCE_ADJUSTED        = fill0 + fv.bio;
        ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC     = fill0 + fv.QC;
        ECOarray.CHLA_FLUORESCENCE_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CHLA_FLUORESCENCE_SCI_CAL_EQU   = 'not applicable';
        INFO.CHLA_FLUORESCENCE_SCI_CAL_COEF  = 'not applicable';
        INFO.CHLA_FLUORESCENCE_SCI_CAL_COM   = 'not applicable';
        INFO.CHLA_FLUORESCENCE_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

        t_nan = isnan(eco_d(:,iChl)); % NaN's in data if any

        ECOarray.FLUORESCENCE_CHLA(~t_nan)    = eco_d(~t_nan,iChl);
        ECOarray.FLUORESCENCE_CHLA_QC(~t_nan) = fv.QC;

        if isfield(cal,'CHL') % Sensor could be bad so maybe no cal info
            ECOarray.CHLA_FLUORESCENCE(~t_nan) = (eco_d(~t_nan,iChl) - cal.CHL.ChlDC) .* ...
                cal.CHL.ChlScale;
            ECOarray.CHLA_FLUORESCENCE_QC(~t_nan) =  3; % 3 do not use w/o adjusting

            ECOarray.CHLA(~t_nan) = (eco_d(~t_nan,iChl) - cal.CHL.ChlDC) .* ...
                cal.CHL.ChlScale;
            ECOarray.CHLA_QC(~t_nan) =  3; % 3 do not use w/o adjusting
			
            % ADJUSTED DATA BASED ON ADMT18 CONCENSUS -WILL BE UPDATED
            % WITHIN THE YEAR - jp 12/13/2017

            % 1st CHECK FOR IN SITU DARKS & TRY AND GET THEM IF NOT THERE
            if ~isfield(cal.CHL, 'SWDC')
                %SWDC = get_CHLdark(MBARI_ID_str, dirs, 5); % 1st 5 good profiles
                SWDC = get_CHLdark(cal, dirs, 5); % 1st 5 good profiles
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
            ECOarray.CHLA_FLUORESCENCE_ADJUSTED(~t_nan) = (eco_d(~t_nan,iChl) - ...
                CHL_DC) .* cal.CHL.ChlScale; %no factor of 2 scale correction for CHLA_FLUORESCENCE!
            ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC(~t_nan) =  2;
			
            ECOarray.CHLA_ADJUSTED(~t_nan) = (eco_d(~t_nan,iChl) - ...
                CHL_DC) .* cal.CHL.ChlScale ./ 2;
            ECOarray.CHLA_ADJUSTED_QC(~t_nan) =  2;
            % NPQ NEXT
            NPQ_CHL = ECOarray.CHLA_ADJUSTED;
            NPQ_CHL(t_nan) = NaN; % fill back to NaN

            NPQ = get_NPQcorr(nanmean(INFO.gps,1), ...
                [[ECOarray.PRES,ECOarray.TEMPi,ECOarray.PSALi],NPQ_CHL], dirs);
            NPQ.data(t_nan,2:end) = fv.bio; % nan back to fill

            if ~isempty(NPQ.data)
                iXing   = find(strcmp('Xing_MLD',NPQ.hdr) == 1);
                iSPIKE  = find(strcmp('CHLspike',NPQ.hdr) == 1);
                tNPQ = ECOarray.PRES <= NPQ.XMLDZ;
                ECOarray.CHLA_ADJUSTED(tNPQ & ~t_nan) = ...
                    NPQ.data(tNPQ & ~t_nan,iXing) + ...
                    NPQ.data(tNPQ & ~t_nan,iSPIKE);
                ECOarray.CHLA_ADJUSTED_QC(tNPQ & ~t_nan) =  5;

            end
			
            ECOarray.CHLA_ADJUSTED_ERROR(~t_nan) = ...
                abs(ECOarray.CHLA_ADJUSTED(~t_nan) * 2);

            INFO.CHLA_SCI_CAL_EQU  = ['CHLA_ADJUSTED=CHLA/A, '...
                'NPQ corrected (Xing et al., 2012), spike profile ', ...
                'added back in'];
            INFO.CHLA_SCI_CAL_COEF = 'A=2';
            INFO.CHLA_SCI_CAL_COM  =['A is best estimate ', ...
                'from Roesler et al., 2017, doi: 10.1002/lom3.10185'];
				
			INFO.CHLA_FLUORESCENCE_SCI_CAL_EQU  = ['CHLA_FLUORESCENCE_ADJUSTED = ((FLUORESCENCE_CHLA-DARK"_CHLA)*SCALE_CHLA)'];
            INFO.CHLA_FLUORESCENCE_SCI_CAL_COEF = ['DARK"_CHLA = ',num2str(CHL_DC),', SCALE_CHLA = ',num2str(cal.CHL.ChlScale)];
            INFO.CHLA_FLUORESCENCE_SCI_CAL_COM  =['CHLA_FLUORESCENCE RT adj specified in http://dx.doi.org/10.13155/35385'];



            % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
            [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CHL');
            t_bio   = ECOarray.CHLA ~= fv.bio;
            ECOarray.FLUORESCENCE_CHLA_QC(t_bio) = ECOarray.FLUORESCENCE_CHLA_QC(t_bio) ...
                * ~BSLflag + BSLflag*theflag;
            ECOarray.CHLA_QC(t_bio) = ECOarray.CHLA_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            ECOarray.CHLA_FLUORESCENCE_QC(t_bio) = ECOarray.CHLA_FLUORESCENCE_QC(t_bio) * ~BSLflag + BSLflag*theflag;
            t_chk = t_bio & (ECOarray.CHLA < RCR.CHL(1)|ECOarray.CHLA > RCR.CHL(2));
            ECOarray.CHLA_QC(t_chk) = 4;
            ECOarray.FLUORESCENCE_CHLA_QC(t_chk) = 4;
            ECOarray.CHLA_FLUORESCENCE_QC(t_chk) = 4;


            if isfield(cal.CHL, 'SWDC')
                t_bio   = ECOarray.CHLA_ADJUSTED ~= fv.bio;
                ECOarray.CHLA_ADJUSTED_QC(t_bio) = ECOarray.CHLA_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
				ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC(t_bio) = ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC(t_bio) ...
                    * ~BSLflag + BSLflag*theflag;
                t_chk = t_bio & ...
                    (ECOarray.CHLA_ADJUSTED < RC.CHL(1)|ECOarray.CHLA_ADJUSTED > RC.CHL(2));
                ECOarray.CHLA_ADJUSTED_QC(t_chk) = 4;
				ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC(t_chk) = 4;
            end

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
                            xxtmp = find(ECOarray.PRES == singleBADs{1});
                            if isempty(xxtmp)
                                disp('WARNING!! CHLA PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                            else
                                if ECOarray.CHLA~=fv.bio
                                    ECOarray.CHLA_QC(xxtmp) = str2double(singleBADSflags{i2});
                                    ECOarray.CHLA_FLUORESCENCE_QC(xxtmp) = str2double(singleBADSflags{i2});									
                                end
                                if ECOarray.CHLA_ADJUSTED~=fv.bio
                                    ECOarray.CHLA_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                                    ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});									
                                end
                            end
                        end
                    end
                    for i3 = 1:length(rangeBADs)
                        ECOarray.CHLA_FLUORESCENCE_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2) & ECOarray.CHLA_FLUORESCENCE~=fv.bio) = str2double(rangeBADsflags{i3});
                        ECOarray.CHLA_FLUORESCENCE_ADJUSTED_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2)& ECOarray.CHLA_FLUORESCENCE_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});					
                        ECOarray.CHLA_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2) & ECOarray.CHLA~=fv.bio) = str2double(rangeBADsflags{i3});
                        ECOarray.CHLA_ADJUSTED_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2)& ECOarray.CHLA_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                    end
                    clear i3
                end
            end
        end

        % CHL DATA_MODE
        if isfield(cal,'CHL')
            %         if isfield(cal.CHL,'SWDC') && isfield(cal.CHL.SWDC,'DC') && sum(LR.CHLA_ADJUSTED<99999)>0  % median dark count has been quantified (and there is data for that profile) --> adjustment has been made
            if sum(ECOarray.CHLA_ADJUSTED<99999)>0  % median dark count has been quantified (and there is data for that profile) --> adjustment has been made
                INFO.CHLA_DATA_MODE = 'A';
				INFO.CHLA_FLUORESCENCE_DATA_MODE = 'A';
            else
                INFO.CHLA_DATA_MODE = 'R';
                INFO.CHLA_FLUORESCENCE_DATA_MODE = 'R';				
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
        ECOarray.BETA_BACKSCATTERING700    = fill0 + fv.bio;
        ECOarray.BETA_BACKSCATTERING700_QC = fill0 + fv.QC;
        ECOarray.BBP700                    = fill0 + fv.bio;
        ECOarray.BBP700_QC                 = fill0 + fv.QC;
        ECOarray.BBP700_ADJUSTED           = fill0 + fv.bio;
        ECOarray.BBP700_ADJUSTED_QC        = fill0 + fv.QC;
        ECOarray.BBP700_ADJUSTED_ERROR     = fill0 + fv.bio;
        INFO.BBP700_SCI_CAL_EQU   = 'not applicable';
        INFO.BBP700_SCI_CAL_COEF  = 'not applicable';
        INFO.BBP700_SCI_CAL_COM   = 'not applicable';
        INFO.BBP700_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

        t_nan = isnan(eco_d(:,iBb)); % NaN's in data if any

        if isfield(cal,'BB') % Sensor could be bad so maybe no cal info
            % VOLUME SCATTERING FUNCTION (VSF) (m^-1 sr^-1)
            VSF(~t_nan) = (eco_d(~t_nan, iBb) - cal.BB.BetabDC) ...
                .* cal.BB.BetabScale;
            % NOW CALCULATE PARTICLE BACKSCATTER COEFFICIENT
            %X       = 1.13*2*pi; % (Barnes and Antoine, 2014)
            %X       = 1.17*2*pi; % (email from E. Boss 24 May 2016)

            % (Proc. Bio-Argo particle backscattering at the DAC level
            % Version 1.2, July 21th 2016
            X      = 1.076*2*pi; % FLBBCD (three channel sensor) chi is 1.076 (for FLBBAP2 chi is 1.097)
            LAMBDA = 700;
            DELTA  = 0.039;      % depolarization ratio
            THETA  = 124;        % FLBBCD 124 (FLBBAP2 is 142)
            BETA_SW_ind = find(t_nan == 0);
            if ~isempty(BETA_SW_ind)
                for ct = 1:size(BETA_SW_ind,1)
                    [BETA_SW(BETA_SW_ind(ct)), b90sw, bsw] = ...
                        betasw_ZHH2009(LAMBDA,ECOarray.TEMPi(BETA_SW_ind(ct)), ...
                        THETA, ECOarray.PSALi(BETA_SW_ind(ct)), DELTA);
                end
            end
            %LR.BETA_BACKSCATTERING700 = VSF; % (DOUBLE CHECK SPEC say counts???)
            ECOarray.BETA_BACKSCATTERING700(~t_nan) = eco_d(~t_nan,iBb); % counts
            ECOarray.BETA_BACKSCATTERING700_QC(~t_nan) = fv.QC; % counts
            ECOarray.BBP700(~t_nan)    = (VSF(~t_nan) - BETA_SW(~t_nan)) * X; %b_bp m^-1
            t_nan_BBP = isnan(ECOarray.BBP700);
            ECOarray.BBP700(t_nan_BBP) = fv.bio;
            ECOarray.BBP700_QC(~t_nan & ~t_nan_BBP) = 2; % 3 do not use w/o adjusting ... 6/10/21 modify qcraw flag from 3 to 2.

            % CALCULATE ADJUSTED DATA _ NO ADJUSTMENTS AT THIS TIME
            % BUT...6/10/21 START POPULATING BBP700_ADJUSTED WITH BBP
            % DIRECTLY!!!  KEEP ORIGINAL BLOCK IN CASE AN ACTUAL ADJUSTMENT
            % GETS IMPLEMENTED/STORED IN THE QC MATRIX.
            ECOarray.BBP700_ADJUSTED(~t_nan& ~t_nan_BBP) = ECOarray.BBP700(~t_nan& ~t_nan_BBP);
            ECOarray.BBP700_ADJUSTED_QC(~t_nan& ~t_nan_BBP) = 1;
            ECOarray.BBP700_ADJUSTED_ERROR(~t_nan& ~t_nan_BBP) = fv.bio; % PLACE HOLDER FOR NOW
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
        t_bio   = ECOarray.BBP700 ~= fv.bio;
        ECOarray.BETA_BACKSCATTERING700_QC(t_bio) = ...
            ECOarray.BETA_BACKSCATTERING700_QC(t_bio) * ~BSLflag + BSLflag*theflag;
        ECOarray.BBP700_QC(t_bio) = ECOarray.BBP700_QC(t_bio) * ~BSLflag + BSLflag*theflag;


        t_chk = t_bio & ...
            (ECOarray.BBP700 < RCR.BB700(1)| ECOarray.BBP700 > RCR.BB700(2));
        ECOarray.BBP700_QC(t_chk) = 4;
        ECOarray.BETA_BACKSCATTERING700_QC(t_chk) = 4;

        t_bio = ECOarray.BBP700_ADJUSTED ~= fv.bio;
        ECOarray.BBP700_ADJUSTED_QC(t_bio) = ECOarray.BBP700_ADJUSTED_QC(t_bio) ...
            * ~BSLflag + BSLflag*theflag;
        t_chk = t_bio & (ECOarray.BBP700_ADJUSTED < RC.BB700(1)| ...
            ECOarray.BBP700_ADJUSTED > RC.BB700(2));
        ECOarray.BBP700_ADJUSTED_QC(t_chk) = 4;

        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('BBP',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(ECOarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! BBP700 PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if ECOarray.BBP700~=fv.bio
                                ECOarray.BBP700_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if ECOarray.BBP700_ADJUSTED~=fv.bio
                                ECOarray.BBP700_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    ECOarray.BBP700_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2) & ECOarray.BBP700~=fv.bio) = str2double(rangeBADsflags{i3});
                    ECOarray.BBP700_ADJUSTED_QC(ECOarray.PRES>=rangeBADs{1}(1) & ECOarray.PRES<=rangeBADs{1}(2)& ECOarray.BBP700_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end

        % ASSIGN ECO PARAMETER DATA MODES:
        % BBP700
        if isfield(cal,'BB')
            if sum(ECOarray.BBP700_ADJUSTED<99999)>0  % there is data for that profile --> adjustment has been made
                INFO.BBP700_DATA_MODE = 'A';
            else
                INFO.BBP700_DATA_MODE = 'R';
            end
        end

        clear BETA_SW X VSF ct b90sw bsw
    end

    % ****************************************************************
    % CALCULATE CDOM (ppb)
    % ****************************************************************
    if ~isempty(iCdm)
        ECOarray.FLUORESCENCE_CDOM    = fill0 + fv.bio;
        ECOarray.FLUORESCENCE_CDOM_QC = fill0 + fv.QC;
        ECOarray.CDOM                 = fill0 + fv.bio;
        ECOarray.CDOM_QC              = fill0 + fv.QC;
        ECOarray.CDOM_ADJUSTED        = fill0 + fv.bio;
        ECOarray.CDOM_ADJUSTED_QC     = fill0 + fv.QC;
        ECOarray.CDOM_ADJUSTED_ERROR  = fill0 + fv.bio;
        INFO.CDOM_SCI_CAL_EQU   = 'not applicable';
        INFO.CDOM_SCI_CAL_COEF  = 'not applicable';
        INFO.CDOM_SCI_CAL_COM   = 'not applicable';
        INFO.CDOM_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

        t_nan = isnan(eco_d(:,iCdm)); % NaN's in data if any
        if isfield(cal,'CDOM') % Sensor could be bad so maybe no cal info

            ECOarray.FLUORESCENCE_CDOM(~t_nan)    = eco_d(~t_nan,iCdm);
            ECOarray.FLUORESCENCE_CDOM_QC(~t_nan) = fv.QC;
            ECOarray.CDOM(~t_nan)  = (eco_d(~t_nan,iCdm) - cal.CDOM.CDOMDC) ...
                .* cal.CDOM.CDOMScale;
            ECOarray.CDOM_QC(~t_nan) =  3; % 3 do not use w/o adjusting

            % NO ADJUSTED DATA YET !!!! THESE ARE PLACE HOLDERS FOR NOW !!!!
            if isfield(QC,'CDOM')
                ECOarray.CDOM_ADJUSTED(~t_nan)        = fv.bio;
                ECOarray.CDOM_ADJUSTED_QC(~t_nan)     = fv.QC;
                ECOarray.CDOM_ADJUSTED_ERROR(~t_nan) =  fv.bio; % PLACE HOLDER FOR NOW
            end
        end

        % DO A FINAL RANGE CHECK ON VALUES, IF BAD SET QF = 4
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'CDOM');
        t_bio = ECOarray.CDOM ~= fv.bio;
        ECOarray.FLUORESCENCE_CDOM_QC(t_bio) = ECOarray.FLUORESCENCE_CDOM_QC(t_bio) ...
            * ~BSLflag + BSLflag*theflag;
        ECOarray.CDOM_QC(t_bio) = ECOarray.CDOM_QC(t_bio) * ~BSLflag + BSLflag*theflag;

        t_chk = t_bio & (ECOarray.CDOM < RCR.CDOM(1)|ECOarray.CDOM > RCR.CDOM(2));
        ECOarray.CDOM_QC(t_chk) = 4;
        ECOarray.FLUORESCENCE_CDOM_QC(t_chk) = 4;

        t_bio = ECOarray.CDOM_ADJUSTED ~= fv.bio;
        ECOarray.CDOM_ADJUSTED_QC(t_bio) = ECOarray.CDOM_ADJUSTED_QC(t_bio)  ...
            * ~BSLflag + BSLflag*theflag;
        t_chk = t_bio & ...
            (ECOarray.CDOM_ADJUSTED < RC.CDOM(1)|ECOarray.CDOM_ADJUSTED > RC.CDOM(2));
        ECOarray.CDOM_ADJUSTED_QC(t_chk) = 4;

        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('CDOM',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(ECOarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! CDOM PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if ECOarray.CDOM~=fv.bio
                                ECOarray.CDOM_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if ECOarray.CDOM_ADJUSTED~=fv.bio
                                ECOarray.CDOM_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    ECOarray.CDOM_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2) & ECOarray.CDOM~=fv.bio) = str2double(rangeBADsflags{i3});
                    ECOarray.CDOM_ADJUSTED_QC(ECOarray.PRES>=rangeBADs{i3}(1) & ECOarray.PRES<=rangeBADs{i3}(2)& ECOarray.CDOM_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end

        % CDOM PARAMETER_DATA_MODE
        if isfield(cal,'CDOM')
            INFO.CDOM_DATA_MODE = 'R';
        end
    end


    if exist('ECOarray','var')
        eval(['BGC0',num2str(BGCIND),' = ECOarray;']); %reassign array structure name.
        clear BGCIND fill0
    end


    % ****************************************************************
    % CALCULATE pH (umol / kg scale)
    % ****************************************************************
    alkN = strfind(INFO.sensors,'ALK');
    alkN2 = find(~(cellfun('isempty',alkN)));
    if exist('ipH','var') && ~isempty(ipH) &&  ~isempty(alkN2)
        BGCIND = find(not(cellfun('isempty',strfind(FLTsensors,'ALK'))))-CTD_Nax; %#ok<STRCL1>
        eval(['ALKarray = BGC0',num2str(BGCIND),';']); %to reduce the number of eval calls, rename to temporary structure until all desired fields are populated.
        fill0  = ones(size(ALKarray.PRES))* 0;

        ALKarray.VRS_PH                          = fill0 + fv.bio;
        ALKarray.VRS_PH_QC                       = fill0 + fv.QC;
        ALKarray.PH_IN_SITU_FREE                 = fill0 + fv.bio;
        ALKarray.PH_IN_SITU_FREE_QC              = fill0 + fv.QC;
        ALKarray.PH_IN_SITU_TOTAL                = fill0 + fv.bio;
        ALKarray.PH_IN_SITU_TOTAL_QC             = fill0 + fv.QC;
        ALKarray.IB_PH                           = fill0 + fv.bio;
        ALKarray.IB_PH_QC                        = fill0 + fv.QC;
        ALKarray.IK_PH                           = fill0 + fv.bio;
        ALKarray.IK_PH_QC                        = fill0 + fv.QC;
        ALKarray.VK_PH                           = fill0 + fv.bio;
        ALKarray.VK_PH_QC                        = fill0 + fv.QC;
        ALKarray.PH_IN_SITU_TOTAL_ADJUSTED       = fill0 + fv.bio;
        ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC    = fill0 + fv.QC;
        ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR = fill0 + fv.bio;
        INFO.PH_SCI_CAL_EQU  = 'not applicable';
        INFO.PH_SCI_CAL_COEF = 'not applicable';
        INFO.PH_SCI_CAL_COM  = 'not applicable';
        INFO.PH_DATA_MODE  = 'R'; %"not applicable" is not acceptable for this field, should be 'R' (per Coriolis)

        t_nan = isnan(alk_d(:,ipH)); %NaN's?
        t_nan_diag = isnan(alk_d(:,iVk)); %need to check separately for NaNs in the diagnostics (sample targets can vary...ie decimated diagnostic sampling on cycles 1-12 for solo1)... % TM, 10/20/22

        ALKarray.VRS_PH(~t_nan)      = alk_d(~t_nan,ipH); % I param
        ALKarray.VRS_PH_QC(~t_nan)   = fv.QC;
        if isfield(cal,'pH')
            [phfree,phtot] = phcalc(ALKarray.VRS_PH(~t_nan), ...
                ALKarray.PRES(~t_nan), ALKarray.TEMPi(~t_nan), ALKarray.PSALi(~t_nan), ...
                cal.pH.k0, cal.pH.k2, cal.pH.pcoefs);
            ALKarray.PH_IN_SITU_FREE(~t_nan)    = phfree; % I param
            ALKarray.PH_IN_SITU_FREE_QC(~t_nan) = fv.QC;

            ALKarray.PH_IN_SITU_TOTAL(~t_nan)    = phtot;
            t_nan_pH = isnan(ALKarray.PH_IN_SITU_TOTAL);
            ALKarray.PH_IN_SITU_TOTAL(t_nan_pH) = fv.bio;
            ALKarray.PH_IN_SITU_TOTAL_QC(~t_nan & ~t_nan_pH) = 3;

            LR_inf = isinf(ALKarray.PH_IN_SITU_FREE); % happens if S = 0
            ALKarray.PH_IN_SITU_FREE(LR_inf)    = 20.1; %UNREAL #
            ALKarray.PH_IN_SITU_FREE_QC(LR_inf) = 4;
            ALKarray.PH_IN_SITU_TOTAL(LR_inf)    = 20.1; %UNREAL #
            ALKarray.PH_IN_SITU_TOTAL_QC(LR_inf) = 4;
        else % data returned  but no cal
            disp(['pH data exists but no cal file for ',INST_ID_str])
        end

        % Diagnostics
        ALKarray.IB_PH(~t_nan_diag & ~t_nan_pH)    = alk_d(~t_nan_diag & ~t_nan_pH,iIb);% * 1e9; % nano amps
        ALKarray.IB_PH_QC(~t_nan_diag & ~t_nan_pH) = fv.QC;
        ALKarray.IK_PH(~t_nan_diag & ~t_nan_pH) = alk_d(~t_nan_diag & ~t_nan_pH, iIk);% * 1e9; % nano amps
        ALKarray.IK_PH_QC(~t_nan_diag & ~t_nan_pH) = fv.QC;
        ALKarray.VK_PH(~t_nan_diag & ~t_nan_pH) = alk_d(~t_nan_diag & ~t_nan_pH, iVk); % nano amps
        ALKarray.VK_PH_QC(~t_nan_diag & ~t_nan_pH) = fv.QC;

        clear ipH_p ipH_t ipH_t ipH_Ib pH_data pH_hdr
        %     end

        if isfield(cal,'pH') && isfield(QC,'pH')
            %             QCD = [ALKarray.PRES(~t_nan & ~t_nan_pH), ...
            %                 ALKarray.TEMPi(~t_nan & ~t_nan_pH), ...
            %                 ALKarray.PSALi(~t_nan & ~t_nan_pH), phtot(~t_nan_pH)];
            QCD = [ALKarray.PRES(~t_nan & ~t_nan_pH), ...
                ALKarray.TEMPi(~t_nan & ~t_nan_pH), ...
                ALKarray.PSALi(~t_nan & ~t_nan_pH), ...
                ALKarray.PH_IN_SITU_TOTAL(~t_nan & ~t_nan_pH)]; % jp 06/16/2022

            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED(~t_nan & ~t_nan_pH) = ...
                apply_QC_corr(QCD, INFO.sdn, QC.pH);
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_nan & ~t_nan_pH) = 1; % set=1 9/27/16 vs
            %ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR(~t_nan) = 0.01;
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED(LR_inf)    = 20.1; %UNREAL #
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(LR_inf) = 4;
            step_tmpPH = find(QC.pH.steps(:,2)<=INFO.cast,1,'last');
            juld_prof = INFO.sdn-datenum(1950,01,01); %convert to JULD
            juld_init = QC.pH.steps(step_tmpPH,1)-datenum(1950,01,01); %convert to JULD
            juld_end = QC.date-datenum(1950,01,01); %date at last DMQC
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR = zeros(size(ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR,1),1)+0.01 + DOXarray.DOXY_ADJUSTED_ERROR(1).*0.0016;
            if juld_prof<juld_end
                ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR = ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR + 0.03.*(juld_prof-juld_end)./365;
            end
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR(t_nan & t_nan_pH) = fv.bio;
            INFO.PH_SCI_CAL_EQU  = ['PH_IN_SITU_TOTAL_ADJUSTED=', ...
                '[PH_IN_SITU_TOTAL+[PUMP_OFFSET - [OFFSET + DRIFT(JULD-JULD_PIVOT)/365]*TCOR]]/GAIN;',...
                'TCOR=(TREF+273.15)./(TEMP+273.15);  TREF = TEMP at 1500m.'];
            INFO.PH_SCI_CAL_COEF = ['PUMP_OFFSET = ',num2str(QC.pH.pHpumpoffset),...
                '; OFFSET = ',num2str(QC.pH.steps(step_tmpPH,4),'%6.4f'),...
                '; DRIFT = ',num2str(QC.pH.steps(step_tmpPH,5),'%6.4f'),...
                '; GAIN = ',num2str(QC.pH.steps(step_tmpPH,3),'%6.4f'),...
                '; JULD = ',num2str(juld_prof,'%9.4f'),...
                '; JULD_PIVOT = ',num2str(juld_init,'%9.4f')];
            INFO.PH_SCI_CAL_COM  =['PUMP_OFFSET derived manually, applied to data above 980m.  OFFSET and DRIFT derived following '...
                'Maurer et al., 2021 (https://doi.org/10.3389/fmars.2021.683207).'...
                'Contact: Tanya Maurer (tmaurer@mbari.org.'];
            if regexp(MBARI_ID_str, bad_O2_filter, 'once')
                [ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR, INFO.PH_SCI_CAL_COM,~,~] = Reassign_Argos_LIReqn8(MBARI_ID_str,INFO.cast,ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_ERROR,INFO.PH_SCI_CAL_COM,0,0);
            end

        end


        % DO A FINAL QC CHECK
        % 04/07/2020 JP - Now that the BSL flag can be 3 or 4 it has the
        % potential to overwrite a 4 with a 3 so don't let this happen
        t_bio = ALKarray.PH_IN_SITU_TOTAL ~= fv.bio;

        %6/11/20 TM need this incase no diagnostic data, 0 if no data.
        % added the ".*1e9" to this line.  It was lacking in Josh's update
        % from 4/7/20 and causing erroneous flagging for certain cases, ie
        % float 9634 surface samples of cycles 67 and 97.
        %t_diag = LR.IB_PH ~= fv.bio.*1e9;

        % JP 12/17/20 slight modification to TM's fix. If the dura file is
        % too small no Ib & Ik tests happen & LR.IB_PH, LR.IK_PH is never
        % reassigned so you need to check for fv.bio*1e9 & fv.bio as fill
        % values
        %
        % TM NOTE: The diagnostics for the BSOLO are only returned every
        % other sample.  Leaving the diag-checks for flagging for now, but
        % we may want to add some smarter checks for full profile flagging
        % if many samples are flagged?  OR remove diag checks?  Or?
        t_diag = ~(ALKarray.IB_PH == fv.bio.*1e9 | ALKarray.IB_PH == fv.bio); % JP 12/17/20

        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'PH');
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (ALKarray.PSALi_QC == 4 | ALKarray.TEMPi_QC == 4 | ALKarray.PRES_QC == 4)*4; % Bad STP will affect pH
        tRC   = (ALKarray.PH_IN_SITU_TOTAL < RCR.PH(1)|ALKarray.PH_IN_SITU_TOTAL > RCR.PH(2))*4; % Range check
        tVRS  = (ALKarray.VRS_PH < RCR.PHV(1) | ALKarray.VRS_PH > RCR.PHV(2))*4;
        tIB3  = (ALKarray.IB_PH < RCR.IB(1) | ALKarray.IB_PH > RCR.IB(2)).*t_diag*3; % DIAGNOSTIC
        tIK3  = (ALKarray.IK_PH < RCR.IK(1) | ALKarray.IK_PH > RCR.IK(2)).*t_diag*3; % DIAGNOSTIC
        tIB4  = (ALKarray.IB_PH < RCR.IB(1)*25 | ALKarray.IB_PH > RCR.IB(2)*25).*t_diag*4; % DIAGNOSTIC
        tIK4  = (ALKarray.IK_PH < RCR.IK(1)*25 | ALKarray.IK_PH > RCR.IK(2)*25).*t_diag*4; % DIAGNOSTIC


        tALL  = max([ALKarray.PH_IN_SITU_TOTAL_QC, tBSL, tSTP, tRC, ...
            tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        ALKarray.PH_IN_SITU_TOTAL_QC(t_bio) = tALL(t_bio);
        ALKarray.PH_IN_SITU_FREE_QC(t_bio)  = tALL(t_bio);

        tALL  = max([ALKarray.VRS_PH_QC, tBSL, tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        ALKarray.VRS_PH_QC(t_bio) = tALL(t_bio);

        tALL  = max([ALKarray.IB_PH_QC, tBSL, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        ALKarray.IB_PH_QC(t_bio)  = tALL(t_bio);
        ALKarray.IK_PH_QC(t_bio)  = tALL(t_bio);
        ALKarray.VK_PH_QC(t_bio)  = tALL(t_bio);

        % 04/07/2020 JP - Now that the BSL flag can be 3 or 4 it has the
        % potential to overwrite a 4 with a 3 so don't let this happen
        t_bio = ALKarray.PH_IN_SITU_TOTAL_ADJUSTED ~= fv.bio;
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (ALKarray.PSALi_QC == 4 | ALKarray.TEMPi_QC == 4 | ....
            ALKarray.PRES_ADJUSTED_QC == 4)*4; % Bad STP will affect pH
        tRC   = (ALKarray.PH_IN_SITU_TOTAL_ADJUSTED < RC.PH(1) | ...
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED > RC.PH(2))*4; % Range check
        tVRS  = (ALKarray.VRS_PH < RC.PHV(1) | ALKarray.VRS_PH > RC.PHV(2))*4;

        % double check fill value QC flags - 9634 returns fv's from qc_adjustmet
        % due to bad gps date bug - jp 10/04/19
        ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(~t_bio)  = fv.QC;

        tALL  = max([ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC, tBSL, tSTP, tRC, ...
            tVRS, tIB3, tIK3, tIB4, tIK4],[],2); % get highest flag
        ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(t_bio) = tALL(t_bio);


        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        % FINALLY DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
        % (BAD).  BGC_spiketest screens for nans and fill values.
        %
        % RUN TEST ON RAW PH_IN_SITU_TOTAL
        QCscreen_phT = ALKarray.PH_IN_SITU_TOTAL_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[ALKarray.PRES ALKarray.PH_IN_SITU_TOTAL],'PH',dirs.cal,fv.bio,QCscreen_phT);
        if ~isempty(spike_inds)
            % Could add optional code to compare qf already assigned (ie in range
            % checks above), and qf assigned during spiketest.  Keep
            % whichever is greater.  (At this point, I'm pre-screening the vector using
            % current quality flags.  Screening occurs within Function.  Currently assigning 4 (bad) to spikes, although may
            % be modified in the future).  5/22/18. TM.
            %spikeQF = LR.PH_INSITU_FREE_QC(spike_inds) < quality_flags;
            %LR.PH_INSITU_FREE_QC(spike_inds(spikeQF)) = quality_flags(spikeQF);
            ALKarray.PH_IN_SITU_TOTAL_QC(spike_inds) = quality_flags;
            disp(['ALKarray.PH_IN_SITU_TOTAL QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
            %         else
            %             disp(['NO LR.PH_IN_SITU_TOTAL SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        %
        % RUN TEST ON QC PH_IN_SITU_TOTAL_ADJUSTED
        QCscreen_phTadj = ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[ALKarray.PRES ALKarray.PH_IN_SITU_TOTAL_ADJUSTED],'PH',dirs.cal,fv.bio,QCscreen_phTadj);
        if ~isempty(spike_inds)
            ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(spike_inds) = quality_flags;
            disp(['PH_IN_SITU_TOTAL_ADJUSTED QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
            %         else
            %             disp(['NO LR.PH_IN_SITU_TOTAL_ADJUSTED SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        % if INFO.cast == 51; keyboard;end
        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('PH',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(ALKarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! PH PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if ALKarray.PH_IN_SITU_TOTAL~=fv.bio
                                ALKarray.PH_IN_SITU_TOTAL_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if ALKarray.PH_IN_SITU_TOTAL_ADJUSTED~=fv.bio
                                ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    ALKarray.PH_IN_SITU_TOTAL_QC(ALKarray.PRES>=rangeBADs{i3}(1) & ALKarray.PRES<=rangeBADs{i3}(2) & ALKarray.PH_IN_SITU_TOTAL~=fv.bio) = str2double(rangeBADsflags{i3});
                    ALKarray.PH_IN_SITU_TOTAL_ADJUSTED_QC(ALKarray.PRES>=rangeBADs{1}(1) & ALKarray.PRES<=rangeBADs{1}(2)& ALKarray.PH_IN_SITU_TOTAL_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end

        % PH DATA-MODE
        XEMPT_PH = find(ALKarray.PH_IN_SITU_TOTAL ~= 99999,1); % if empty, then cycle has no data
        %         if isempty(QC) || isempty(XEMPT_PH) % no adjustments have been made yet, or all data is 99999
        %             INFO.PH_DATA_MODE = 'R';
        if ~isempty(QC) && isfield(QC,'pH') && ~isempty(XEMPT_PH)
            %             if cycdate > QC.date
            if str2num(cast_num) > QC.pH.steps(end,2)
                INFO.PH_DATA_MODE = 'A';
            else
                INFO.PH_DATA_MODE = 'D';
            end
        else
            INFO.PH_DATA_MODE = 'R';
            INFO.PH_SCI_CAL_EQU  = 'not applicable';
            INFO.PH_SCI_CAL_COEF = 'not applicable';
            INFO.PH_SCI_CAL_COM  = 'not applicable';
        end

        eval(['BGC0',num2str(BGCIND),' = ALKarray;']); %reassign array structure name.
        clear BGCIND fill0 phfree phtot QCD dura t_chk1 t_chk2 t_chk3 QCscreen_phF QCscreen_phT QCscreen_phTadj QCscreen_O QCscreenOadj spike_inds quality_flags LR_inf

        %if INFO.cast == 30,pause,end % TESTING

    end

    % ****************************************************************
    % CALCULATE NITRATE (umol / kg scale)
    % DO DEPTH CORRECTION 1st
    % CONVERT TO umol/kg
    % ****************************************************************

    if ~isempty(iNO3)
        BGCIND = find(not(cellfun('isempty',strfind(FLTsensors,'NO3'))))-CTD_Nax; %#ok<STRCL1>
        eval(['NO3array = BGC0',num2str(BGCIND),';']); %to reduce the number of eval calls, rename to temporary structure until all desired fields are populated.
        fill0  = ones(size(NO3array.PRES))* 0;

        NO3array.NITRATE                      = fill0 + fv.bio;
        NO3array.NITRATE_QC                   = fill0 + fv.QC;
        NO3array.UV_INTENSITY_DARK_NITRATE    = fill0 + fv.bio;
        NO3array.UV_INTENSITY_DARK_NITRATE_QC = fill0 + fv.QC;
        NO3array.UV_INTENSITY_NITRATE         = fill0 + fv.bio;
        NO3array.UV_INTENSITY_NITRATE_QC      = fill0 + fv.QC;

        NO3array.NITRATE_ADJUSTED             = fill0 + fv.bio;
        NO3array.NITRATE_ADJUSTED_QC          = fill0 + fv.QC;
        NO3array.NITRATE_ADJUSTED_ERROR       = fill0 + fv.bio;
        INFO.NITRATE_SCI_CAL_EQU  = 'not applicable';
        INFO.NITRATE_SCI_CAL_COEF = 'not applicable';
        INFO.NITRATE_SCI_CAL_COM  = 'not applicable';
        INFO.NITRATE_DATA_MODE = 'R';

        % IF ISUS FILE & CAL INFO EXIST, TRY AND PROCESS IT
        if isfield(cal,'N') && Dno3.spectra_pix_range(1)~=999
            spec = Dno3; %everything is included in the original parsed .no3 file aside from a few structure variables.  Add them here.
            % The raw spec data read in from the no3 file is still reversed.
            % Flip the Argo interp arrays such that the NO3 processing is
            % operating on properly aligned arrays.
            spec.P = flipud(NO3array.PRES);
            spec.T = flipud(NO3array.TEMPi);
            spec.S = flipud(NO3array.PSALi);

            UV_INTEN = spec.UV_INTEN;

            spec.SDN = spec.DC*NaN;
            if ~isempty(UV_INTEN) && ~isempty(NO3array.PRES)
                % [SDN, DarkCur, Pres, Temp, Sal, NO3, BL_int,BL_slope,
                %  RMS_ER, Wl~240, ABS~240] !!! NITRATE STILL umol/L !!!

                NO3  = calc_FLOAT_NO3(spec, cal.N, 1); % ESW P corr
                %NO3  = calc_APEX_NO3_JP(spec, cal.N, 0); % NO P corr

                IX = (size(NO3,1):-1:1)'; % FLIP SHALLOW TO DEEP FOR ARGO
                %[B,IX]   = sort(NO3(:,3)); % SORT SHALLOW TO DEEP FOR ARGO
                NO3      = NO3(IX,:); % umol/L
                UV_INTEN = UV_INTEN(IX,:);
                clear B IX
            end

        else
            disp('No calibration info to process nitrate with')
        end

        if exist('NO3', 'var') == 1
            % TM NOTE: I don't think this is needed as the pres t/s interp
            % makes this check obsolete.
            %--
            % MATCH NITRATE DEPTHS TO CTD DEPTHS IF LENGTHS DIFFERENT
            % MISSING SPECTRA  AND NO3 CONC. SET TO NaN
            % could also  replace missing NO3 from *.isus with *.msg

            % % %             if size(LR.PRES,1) ~= size(NO3,1) % Mismatch in sample lengths
            % % %                 cS       = size(spec.UV_INTEN,2); % cols of WL's
            % % %                 rP       = size(LR.PRES,1); % # of CTD samples (*.msg)
            % % %                 [rN, cN] = size(NO3); % RN = # of isus samples (*.isus)
            % % %                 disp(['*.isus(',num2str(rN),') & *.msg(',num2str(rP), ...
            % % %                     ') sample counts not equal for ',msg_file])
            % % %                 NO3tmp   = ones(rP,cN)* NaN;
            % % %                 spectmp  = ones(rP,cS)* fv.bio;
            % % %                 % fill msg file data then replace with NO3 data
            % % %                 %NO3tmp(:,3:6) = [LR.PRES, LR.TEMP, LR.PSAL lr_d(:,iNO3)];
            % % %                 NO3tmp(:,3:5) = [LR.PRES, LR.TEMP, LR.PSAL];
            % % %
            % % %                 for i = 1: rP
            % % %                     ind  = find(NO3(:,3) == LR.PRES(i),1,'first');
            % % %                     if ~isempty(ind)
            % % %                         NO3tmp(i,:)  = NO3(ind,:);
            % % %                         spectmp(i,:) = UV_INTEN(ind,:);
            % % %                     end
            % % %                 end
            % % %                 NO3      = NO3tmp;
            % % %
            % % %                 UV_INTEN = spectmp; % n x m matrix
            % % %                 clear rP rN cN NO3tmp i
            % % %             end

            t_nan  = isnan(NO3(:,6));
            % CHECK FIT ERROR AND BASELINE ABSORBANCE
            tABS08 = NO3(:,11) > 0.8; %ABS240 > 0.8 (QF=3)
            tABS11 = NO3(:,9) > 0.003 | NO3(:,11) > 1.1; % RMS & ABS240 (QF=4)

            NO3array.UV_INTENSITY_DARK_NITRATE    = NO3(:,2);
            NO3array.UV_INTENSITY_DARK_NITRATE(t_nan) = fv.bio;
            NO3array.UV_INTENSITY_DARK_NITRATE_QC = fill0 + fv.QC;

            NO3array.UV_INTENSITY_NITRATE    =  UV_INTEN;
            NO3array.UV_INTENSITY_NITRATE_QC =  UV_INTEN * 0 + fv.QC;

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

            NO3_p0_kg         = NO3_p0 ./ N_den * 1000; % umol/kg
            NO3array.NITRATE        = NO3_p0_kg;
            NO3array.NITRATE(t_nan) = fv.bio;
            % %             t_nan_NIT = isnan(DOXarray.DOXY);  %Not Needed for NO3??

            % QC RAW
            NO3array.NITRATE_QC = fill0 + 3; % questionable to start
            NO3array.NITRATE_QC(t_nan) = fv.QC;
            NO3array.NITRATE_QC(~t_nan & tABS11) = 4;

            % ********************************************************
            % APPLY QC CORRECTIONS
            % CORRECTIONS DETERMINED ON umol/L scale so adjust on that
            % scale and then convert
            if isfield(cal,'N') && isfield(QC,'N')
                QCD = [NO3array.PRES, NO3array.TEMPi, NO3array.PSALi,NO3_p0_kg];
                NO3array.NITRATE_ADJUSTED = fill0 + fv.bio;
                NO3array.NITRATE_ADJUSTED(~t_nan) = ...
                    apply_QC_corr(QCD(~t_nan,:), INFO.sdn, QC.N);

                NO3array.NITRATE_ADJUSTED_QC = fill0 + 1; % set=1 9/27/16
                NO3array.NITRATE_ADJUSTED_QC(t_nan) = fv.QC; % nan's to fill value

                NO3array.NITRATE_ADJUSTED_QC(~t_nan & tABS08) = 3; %Bad ABS240
                if sum(tABS11) > 0
                    disp(['Nitrate ABS@240 > 1.1 or RMS > 0.003 ',...
                        'detected. Setting QF = 4 for flagged values'])
                    NO3array.NITRATE_ADJUSTED_QC(~t_nan & tABS11) = 4; % Really bad abs or RMS
                end

                %LR.NITRATE_ADJUSTED_ERROR = (abs(LR.NITRATE - ...
                %    LR.NITRATE_ADJUSTED)) * 0.1 + 0.5;
                %LR.NITRATE_ADJUSTED_ERROR(t_nan) = fv.bio;
                step_tmpN = find(QC.N.steps(:,2)<=INFO.cast,1,'last');
                juld_prof = INFO.sdn-datenum(1950,01,01); %convert to JULD
                juld_init = QC.N.steps(step_tmpN,1)-datenum(1950,01,01); %convert to JULD
                juld_end = QC.date-datenum(1950,01,01); %date at last DMQC
                NO3array.NITRATE_ADJUSTED_ERROR = ones(size(NO3array.NITRATE_ADJUSTED_ERROR,1),1) + DOXarray.DOXY_ADJUSTED_ERROR(1)./10;
                if juld_prof<juld_end
                    NO3array.NITRATE_ADJUSTED_ERROR = NO3array.NITRATE_ADJUSTED_ERROR + 1.*(juld_prof-juld_end)./365;
                end
                NO3array.NITRATE_ADJUSTED_ERROR(t_nan) = fv.bio;
                INFO.NITRATE_SCI_CAL_EQU  = ['NITRATE_ADJUSTED=', ...
                    '[NITRATE-[OFFSET + DRIFT(JULD-JULD_PIVOT)/365]]/GAIN'];
                INFO.NITRATE_SCI_CAL_COEF = ['OFFSET = ', ...
                    num2str(QC.N.steps(step_tmpN,4),'%6.4f'),...
                    '; DRIFT = ',num2str(QC.N.steps(step_tmpN,5),'%6.4f'),...
                    '; GAIN = ',num2str(QC.N.steps(step_tmpN,3),'%6.4f'),...
                    '; JULD = ',num2str(juld_prof,'%9.4f'),...
                    '; JULD_PIVOT = ',num2str(juld_init,'%9.4f')];
                INFO.NITRATE_SCI_CAL_COM  =['Adjustments derived following '...
                    'Maurer et al., 2021 (https://doi.org/10.3389/fmars.2021.683207).'...
                    'Contact Tanya Maurer (tmaurer@mbari.org) ',...
                    'for more information.'];
                if regexp(MBARI_ID_str, bad_O2_filter, 'once')
                    [~, ~, NO3array.NITRATE_ADJUSTED_ERROR, INFO.NITRATE_SCI_CAL_COM] = Reassign_ArgoSpecs_LIReqn8(MBARI_ID_str,INFO.cast,0,0,NO3array.NITRATE_ADJUSTED_ERROR,INFO.NITRATE_SCI_CAL_COM);
                end
            end
            clear QCD UV_INTEN
        end

        % DO A FINAL RANGE CHECK BASED ON BSL & OTHER DEPENDANT PARAMETERS
        % MODIFIED 04/07/2020 BY JP BECAUSE BSL FLAG OF 3 COULD OVERWRITE A 4
        t_bio = NO3array.NITRATE ~= fv.bio; % non fill values
        [BSLflag, theflag] = isbadsensor(BSL, MBARI_ID_str, INFO.cast, 'N');
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (NO3array.PSALi_QC == 4 | NO3array.TEMPi_QC == 4 | NO3array.PRES_QC == 4)*4; % Bad STP will affect nitrate
        tRC   = (NO3array.NITRATE < RCR.NO3(1)| NO3array.NITRATE > RCR.NO3(2))*4; % RAnge check

        tALL  = max([NO3array.NITRATE_QC, tBSL, tSTP, tRC],[],2); % get highest flag
        NO3array.NITRATE_QC(t_bio) = tALL(t_bio);

        tALL  = max([NO3array.UV_INTENSITY_DARK_NITRATE_QC, tBSL, tRC],[],2); % get highest flag
        NO3array.UV_INTENSITY_DARK_NITRATE_QC(t_bio) = tALL(t_bio);


        t_bio = NO3array.NITRATE_ADJUSTED ~= fv.bio;
        tBSL  = ones(size(t_bio)) * ~BSLflag + BSLflag*theflag; % flag from BSL
        tSTP  = (NO3array.PSALi_QC == 4 | NO3array.TEMPi_QC == 4 | ...
            NO3array.PRES_ADJUSTED_QC == 4)*4;
        tRC   = (NO3array.NITRATE_ADJUSTED < RC.NO3(1)| NO3array.NITRATE_ADJUSTED > RC.NO3(2))*4;

        tALL  = max([NO3array.NITRATE_ADJUSTED_QC, tBSL, tSTP, tRC],[],2); % get highest flag
        NO3array.NITRATE_ADJUSTED_QC(t_bio) = tALL(t_bio);

        % double check fill value QC flags - 9634 returns fv's from qc_adjustmet
        % due to bad gps date bug - jp 10/04/19
        NO3array.NITRATE_ADJUSTED_QC(~t_bio)  = fv.QC;

        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        % DO A SPIKE TEST ON VALUES, IF SPIKES ARE IDENTIFIED, SET QF TO 4
        % (BAD).  BGC_spiketest screens for nans and fill values.
        %
        % RUN TEST ON RAW NITRATE
        QCscreen_N = NO3array.NITRATE_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[NO3array.PRES NO3array.NITRATE],'NO3',dirs.cal,fv.bio,QCscreen_N);
        if ~isempty(spike_inds)
            NO3array.NITRATE_QC(spike_inds) = quality_flags;
            disp(['NITRATE QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
            %         else
            %             disp(['NO LR.NITRATE SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        %
        % RUN TEST ON QC NITRATE_ADJUSTED
        QCscreen_Nadj = NO3array.NITRATE_ADJUSTED_QC == 4; % screen for BAD data already assigned.
        [spike_inds, quality_flags] = BGC_spiketest(MBARI_ID_str,str2double(cast_num),[NO3array.PRES NO3array.NITRATE_ADJUSTED],'NO3',dirs.cal,fv.bio,QCscreen_Nadj);
        if ~isempty(spike_inds)
            NO3array.NITRATE_ADJUSTED_QC(spike_inds) = quality_flags;
            disp(['LR.NITRATE_ADJUSTED QUALITY FLAGS ADJUSTED FOR IDENTIFIED SPIKES, PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
            %         else
            %             disp(['NO LR.NITRATE_ADJUSTED SPIKES IDENTIFIED FOR PROFILE ',cast_num,' FLOAT ',MBARI_ID_str,'.'])
        end
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('N',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(NO3array.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! PH PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if NO3array.NITRATE~=fv.bio
                                NO3array.NITRATE_QC(xxtmp) = str2double(singleBADSflags{i3});
                            end
                            if NO3array.NITRATE_ADJUSTED~=fv.bio
                                NO3array.NITRATE_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i3});
                            end
                        end
                    end
                    clear i2
                end
                for i3 = 1:length(rangeBADs)
                    NO3array.NITRATE_QC(NO3array.PRES>=rangeBADs{i3}(1) & NO3array.PRES<=rangeBADs{i3}(2) & NO3array.NITRATE~=fv.bio) = str2double(rangeBADsflags{i3});
                    NO3array.NITRATE_ADJUSTED_QC(NO3array.PRES>=rangeBADs{1}(1) & NO3array.PRES<=rangeBADs{1}(2)& NO3array.NITRATE_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end

        % NITRATE DATA MODE
        XEMPT_N = find(NO3array.NITRATE ~= 99999,1); % if empty, then cycle has no data %% TM 2/14/23; I don't think this is right!  Can still be D-mode if no data!
        %if isempty(QC) || isempty(XEMPT_N) % no adjustments have been made yet, or all data is 99999
        %    INFO.NITRATE_DATA_MODE = 'R';
        if ~isempty(QC) && isfield(QC,'N') && ~isempty(XEMPT_N)
            %             if cycdate > QC.date
            if str2num(cast_num) > QC.N.steps(end,2)
                INFO.NITRATE_DATA_MODE = 'A';
            else
                INFO.NITRATE_DATA_MODE = 'D';
            end
        else
            INFO.NITRATE_DATA_MODE = 'R';
            INFO.NITRATE_SCI_CAL_EQU  = 'not applicable';
            INFO.NITRATE_SCI_CAL_COEF = 'not applicable';
            INFO.NITRATE_SCI_CAL_COM  = 'not applicable';
        end
        eval(['BGC0',num2str(BGCIND),' = NO3array;']); %reassign array structure name.
        clear BGCIND tchk tABS08 tABS11 NO3 QCscreen_N QCscreen_Nadj spike_inds quality_flags
    end


    % ****************************************************************
    % NOW DEAL WITH OCR DATA EXCEPTION
    % ****************************************************************
    if isfield(cal,'OCR')
        BGCIND = find(not(cellfun('isempty',strfind(FLTsensors,'OCR'))))-CTD_Nax; %#ok<STRCL1>
        eval(['OCRarray = BGC0',num2str(BGCIND),';']); %to reduce the number of eval calls, rename to temporary structure until all desired fields are populated.
        fill0  = ones(size(OCRarray.PRES))* 0;
        for ocr_i = 1:4
            eval(['tmpwv = cal.OCR.CH0',num2str(ocr_i),'.WL;']);
            eval(['a0 = cal.OCR.CH0',num2str(ocr_i),'.a0;']);
            eval(['a1 = cal.OCR.CH0',num2str(ocr_i),'.a1;']);
            eval(['im = cal.OCR.CH0',num2str(ocr_i),'.im;']);
            if strcmp(tmpwv,'PAR')
                INFO = setfield(INFO,['DOWNWELLING_',tmpwv,'_SCI_CAL_EQU'],'not applicable');
                INFO = setfield(INFO,['DOWNWELLING_',tmpwv,'_SCI_CAL_COEF'],'not applicable');
                INFO = setfield(INFO,['DOWNWELLING_',tmpwv,'_SCI_CAL_COM'],'not applicable');
                INFO = setfield(INFO,['DOWNWELLING_',tmpwv,'_DATA_MODE'],'R');
                OCRarray = setfield(OCRarray,['RAW_DOWNWELLING_',tmpwv],fill0 + fv.bio);
                OCRarray = setfield(OCRarray,['RAW_DOWNWELLING_',tmpwv,'_QC'],fill0 + fv.QC);
                OCRarray = setfield(OCRarray,['DOWNWELLING_',tmpwv],fill0 + fv.bio);
                OCRarray = setfield(OCRarray,['DOWNWELLING_',tmpwv,'_QC'],fill0 + fv.QC);
                OCRarray = setfield(OCRarray,['DOWNWELLING_',tmpwv,'_ADJUSTED'],fill0 + fv.bio);
                OCRarray = setfield(OCRarray,['DOWNWELLING_',tmpwv,'_ADJUSTED_QC'],fill0 + fv.QC);
                iIND = indocr(ocr_i);
                Tnan = isnan(ocr_d(:,iIND));
                eval(['OCRarray.RAW_DOWNWELLING_',tmpwv,'(~Tnan) = ocr_d(~Tnan,iIND);']);
                eval(['OCRarray.RAW_DOWNWELLING_',tmpwv,'_QC(~Tnan) = fv.QC;']);
                eval(['OCRarray.DOWNWELLING_',tmpwv,'(~Tnan) = (OCRarray.RAW_DOWNWELLING_',tmpwv,'(~Tnan)-a0).*a1.*im;']);
                eval(['OCRarray.DOWNWELLING_',tmpwv,'_QC(~Tnan) = 2;']);
                eval(['tmpvarn = OCRarray.DOWNWELLING_',tmpwv,';']);
                eval(['RngChktmp =RCR.OCR',tmpwv,';']);
                t_bio = tmpvarn ~= fv.bio; % Non fill value samples
                t_chk = tmpvarn < RngChktmp(1)| tmpvarn > RngChktmp(2); % range check
                t_chk = t_chk & t_bio; % & not MVI either
                eval(['OCRarray.DOWNWELLING_',tmpwv,'_QC(t_chk) = 4;']);
                clear RngChktmp t_bio t_chk tmpwv
            else
                INFO = setfield(INFO,['DOWN_IRRADIANCE',tmpwv,'_SCI_CAL_EQU'],'not applicable');
                INFO = setfield(INFO,['DOWN_IRRADIANCE',tmpwv,'_SCI_CAL_COEF'],'not applicable');
                INFO = setfield(INFO,['DOWN_IRRADIANCE',tmpwv,'_SCI_CAL_COM'],'not applicable');
                INFO = setfield(INFO,['DOWN_IRRADIANCE',tmpwv,'_DATA_MODE'],'R');
                OCRarray = setfield(OCRarray,['RAW_DOWNWELLING_IRRADIANCE',tmpwv],fill0 + fv.bio);
                OCRarray = setfield(OCRarray,['RAW_DOWNWELLING_IRRADIANCE',tmpwv,'_QC'],fill0 + fv.QC);
                OCRarray = setfield(OCRarray,['DOWN_IRRADIANCE',tmpwv],fill0 + fv.bio);
                OCRarray = setfield(OCRarray,['DOWN_IRRADIANCE',tmpwv,'_QC'],fill0 + fv.QC);
                OCRarray = setfield(OCRarray,['DOWN_IRRADIANCE',tmpwv,'_ADJUSTED'],fill0 + fv.bio);
                OCRarray = setfield(OCRarray,['DOWN_IRRADIANCE',tmpwv,'_ADJUSTED_QC'],fill0 + fv.QC);
                iIND = indocr(ocr_i);
                Tnan = isnan(ocr_d(:,iIND));
                eval(['OCRarray.RAW_DOWNWELLING_IRRADIANCE',tmpwv,'(~Tnan) = ocr_d(~Tnan,iIND);']);
                eval(['OCRarray.RAW_DOWNWELLING_IRRADIANCE',tmpwv,'_QC(~Tnan) = fv.QC;']);
                eval(['OCRarray.DOWN_IRRADIANCE',tmpwv,'(~Tnan) = (OCRarray.RAW_DOWNWELLING_IRRADIANCE',tmpwv,'(~Tnan)-a0).*a1.*im.*0.01;']); %units W/m2/nm
                eval(['OCRarray.DOWN_IRRADIANCE',tmpwv,'_QC(~Tnan) = 2;']);
                eval(['tmpvarn = OCRarray.DOWN_IRRADIANCE',tmpwv,';']);
                eval(['RngChktmp =RCR.OCR',tmpwv,';']);
                t_bio = tmpvarn ~= fv.bio; % Non fill value samples
                t_chk = tmpvarn < RngChktmp(1)| tmpvarn > RngChktmp(2); % range check
                t_chk = t_chk & t_bio; % & not MVI either
                eval(['OCRarray.DOWN_IRRADIANCE',tmpwv,'_QC(t_chk) = 4;']);
                clear RngChktmp t_bio t_chk tmpwv
            end
        end
        %% Brute force add bad sample list for OCR!-------------------------
        % TM 11/21/22 Add OCR to bad sample list processing.  MAKE THIS SMARTER!!
        if yesBSAML == 1 && yesBSAMLcyc==1
            SbsIND = find(strcmp('OCR380',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2 = 1:length(singleBADs)
                        xxtmp = find(OCRarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! OCR PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if OCRarray.DOWN_IRRADIANCE380~=fv.bio
                                OCRarray.DOWN_IRRADIANCE380_QC(xxtmp) = str2double(singleBADSflags{i3});
                            end
                            if OCRarray.DOWN_IRRADIANCE380_ADJUSTED~=fv.bio
                                OCRarray.DOWN_IRRADIANCE380_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i3});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    OCRarray.DOWN_IRRADIANCE380_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2) & OCRarray.DOWN_IRRADIANCE380~=fv.bio) = str2double(rangeBADsflags{i3});
                    OCRarray.DOWN_IRRADIANCE380_ADJUSTED_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2)& OCRarray.DOWN_IRRADIANCE380_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end


            SbsIND = find(strcmp('OCR412',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(OCRarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! OCR PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if OCRarray.DOWN_IRRADIANCE412~=fv.bio
                                OCRarray.DOWN_IRRADIANCE412_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if OCRarray.DOWN_IRRADIANCE412_ADJUSTED~=fv.bio
                                OCRarray.DOWN_IRRADIANCE412_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    OCRarray.DOWN_IRRADIANCE412_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2) & OCRarray.DOWN_IRRADIANCE412~=fv.bio) = str2double(rangeBADsflags{i3});
                    OCRarray.DOWN_IRRADIANCE412_ADJUSTED_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2)& OCRarray.DOWN_IRRADIANCE412_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end


            SbsIND = find(strcmp('OCR490',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(OCRarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! OCR PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if OCRarray.DOWN_IRRADIANCE490~=fv.bio
                                OCRarray.DOWN_IRRADIANCE490_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if OCRarray.DOWN_IRRADIANCE490_ADJUSTED~=fv.bio
                                OCRarray.DOWN_IRRADIANCE490_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3=1:length(rangeBADs)
                    OCRarray.DOWN_IRRADIANCE490_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2) & OCRarray.DOWN_IRRADIANCE490~=fv.bio) = str2double(rangeBADsflags{i3});
                    OCRarray.DOWN_IRRADIANCE490_ADJUSTED_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2)& OCRarray.DOWN_IRRADIANCE490_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end


            SbsIND = find(strcmp('OCRPAR',tmpBSAML.list(:,ibsSENS)));
            if ~isempty(SbsIND)
                TMPsbs = tmpBSAML.list(SbsIND,:);
                singleBADs = TMPsbs(:,ibsD);
                singleBADSflags = TMPsbs(:,ibsFL);
                rangeBADs = TMPsbs(:,ibsDB);
                rangeBADsflags = TMPsbs(:,ibsFL);
                if ~isempty(singleBADs{1})
                    for i2=1:length(singleBADs)
                        xxtmp = find(OCRarray.PRES == singleBADs{1});
                        if isempty(xxtmp)
                            disp('WARNING!! OCR PRESSURES INDICATED AS BAD ON BAD SAMPLE LIST DO NOT EXIST IN FILE.  NO FLAGGING CHANGED.')
                        else
                            if OCRarray.DOWNWELLING_PAR~=fv.bio
                                OCRarray.DOWNWELLING_PAR_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                            if OCRarray.DOWNWELLING_PAR_ADJUSTED~=fv.bio
                                OCRarray.DOWNWELLING_PAR_ADJUSTED_QC(xxtmp) = str2double(singleBADSflags{i2});
                            end
                        end
                    end
                    clear i2
                end
                for i3 = 1:length(rangeBADs)
                    OCRarray.DOWNWELLING_PAR_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2) & OCRarray.DOWNWELLING_PAR~=fv.bio) = str2double(rangeBADsflags{i3});
                    OCRarray.DOWNWELLING_PAR_ADJUSTED_QC(OCRarray.PRES>=rangeBADs{i3}(1) & OCRarray.PRES<=rangeBADs{i3}(2)& OCRarray.DOWNWELLING_PAR_ADJUSTED~=fv.bio) = str2double(rangeBADsflags{i3});
                end
                clear i3
            end
        end
        %------------------------------------------------------------------------


        eval(['BGC0',num2str(BGCIND),' = OCRarray;']); %reassign array structure name.
        clear fill0
    end
    % end
    % % %
    % % % 	% TM 2/13/23; Now perform final check to ensure no bad data was propagated to adjusted fields.
    % % % 	% Note 4903026 cycle 81 had no pH data returned
    % % % 	if exist('BGC01')
    % % %     [BGC01] = screen_bad_adjusted(BGC01);
    % % % 	end
    % % % 	if exist('BGC02')
    % % % 	[BGC02] = screen_bad_adjusted(BGC02);
    % % % 	end
    % % % 	if exist('BGC03')
    % % % 	[BGC03] = screen_bad_adjusted(BGC03);
    % % % 	end
    % % % 	if exist('BGC04')
    % % % 	[BGC04] = screen_bad_adjusted(BGC04);
    % % % 	end
    % % % 	if exist('BGC05')
    % % % 	[BGC05] = screen_bad_adjusted(BGC05);
    % % % 	end


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

    cast_NUM = cast_num(end-2:end); %TM NOTE: Keep with MBARI standardized formatting for now.  SOLOs are using 4char cast number, our mat files use 3char.
    save_str = [dirs.mat, WMO,filesep, WMO,'.', cast_NUM,'.mat'];

    % clear any remaining extraneous variables that would slip in with regexp
    clear BGCIND
    if exist('TRAJ','var')
        save(save_str,'-regexp','^BGC0*','LR','HR','INFO','TRAJ');
    else
        save(save_str,'-regexp','^BGC0*','LR','HR','INFO');
    end
    %     if msg_ct == 1
    %         copyfile(fp_cal, [dirs.mat, WMO,filesep]); % copy over cal file
    %     end
    clear BGC*
end
copyfile(fp_cal, [dirs.mat, WMO,filesep]) % copy over cal file
%fprintf('\r\n')

% *********************************************************************
% CLEAN UP
% Check for msg files in temp dir & delete if there JP 06/15/22
for i6 = 1:length(FileList)
    if ~isempty(FileList{i6}.list)
        fn_tmp = FileList{i6}.list{1,1}; % get a file name to get extension
        ext    = regexp(fn_tmp,'\w+$','match','once');
        if ~isempty(ls([dirs.temp,'*.', ext]))
            delete([dirs.temp,'*.', ext])
        end
    end
end

% if ~isempty(ls([dirs.temp,'*.???']))
%     delete([dirs.temp,'*.???']);
% end
%     tf_float.status = 1;
%end












