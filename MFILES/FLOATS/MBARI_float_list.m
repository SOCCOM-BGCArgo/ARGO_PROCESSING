function d = MBARI_float_list(dirs)
%
% OVERVIEW:
%    This function builds a master list structure for MBARI profiling
%    float processing:
%
% INPUT:
%   dirs - a structure of directories & paths or []. If empty defaults
%          are used
%
% OUTPUT:
%   d - s tructure with 3 fields (hdr, list, and log)
%       hdr  - 1 x n cell array describing the n colums in "list"
%       list - master list for MBARI float processing
%       log  - cell array of strings noting any problems
%
%       A *.mat & *.txt file will be saved to dirs.cal
%
% EXAMPLES:
%       d = MBARI_float_list(dirs);
%       d = MBARI_float_list([]);
%
% CHANGE HISTORY:
%    01/09/2017 - JP,if UW site down function won't break now
%    01/09/2020 - JP, Complete over haul, removed obsolete tasks, condensed
%                 & cleaned up
%    06/09/2020 - JP, fixed bugs in message string if no connection was made
%                 to UW
%    11/04/2020 - JP, changed source URL and handling for UW WMO listing &
%                 added new exception blocks for 1117SOOCN, 0948STNP &
%                 0948STNP2
%    01/12/2021 - JP, Complete overhual & merged with Float type list in
%                 prep for GO_GBC. GETTING WMO & INSTITUTE ID'S FROM UW,
%                 SIO, WHOI
%    01/28/2021 - JP, switched to using FloatConfig dir to define master list
%    02/25/2021 - JP, Version 2, reconfigured for new mbari id naming
%                 convention & GO-BGC
%    03/01/2021 - JP, un1173 fixed incorrect WMO assignment to 3xO2 float
%                 due to older duplicate UW ID with assigned WMO
%    03/18/2021 - JP added code to flag dead SOCCOM floats from Sharon /
%                 Tanya's float status table & identify if ph, no3 or chl
%                 on board
%    03/23/2021 - JP - changes due to changes in Sharon's stats table & set
%                 tf_bfile for temporary WMO = -1

% *************************************************************************
% TESTING
% dirs = [];
% *************************************************************************
%
% *************************************************************************
% *************************************************************************
% SET DIRS, FILENAMES, & FILTERS
% *************************************************************************
% *************************************************************************
% Initialize output
%out_name = 'new_MBARI_float_list'; % save name
tic
out_name = 'MBARI_float_list'; % save name
d.list = {};
d.hdr  = {};
d.log  = {};

log = cell(100,1);
log_ct = 0;

if isempty(dirs)
    user_dir   = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir   = [user_dir, '\Documents\MATLAB\'];
    dirs.mat   = [user_dir,'ARGO_PROCESSING\DATA\FLOATS\'];
    dirs.cal   = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
    dirs.temp  = 'C:\temp\'; % path to temporary working dir
    dirs.msg   = '\\seaecho.shore.mbari.org\floats\';
    dirs.bfile = [user_dir,'ARGO_MBARI2AOML\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

% PROGRAM LIST FILE PATHS
dirs.config = [dirs.cal,'FLOAT_CONFIG\']; % MBARI master list here
dirs.ftp    = '\\atlas\ftp\pub\SOCCOM\RawFloatData\Cal_files\';

% ADD URL's TO GET WMO #'s FROM INSTITUTIONS
% dirs.UW_html = ['http://runt.ocean.washington.edu/argo/', ...
%     'heterographs/rollcall.html'];
dirs.UW_URL   = 'http://runt.ocean.washington.edu/swift/WmoIdMap';
dirs.SIO_URL  = 'http://sio-argo.ucsd.edu/historical.html';
%dirs.WHOI_URL = 'http://argo.whoi.edu/gdac_index.html';
%dirs.WHOI_DOM = 'http://db.whoifloatgroup.org/';
dirs.WHOI_DOM = 'https://db.whoifloatgroup.org/';
api_call_str = '%sapi/%s?FLOAT_SERIAL_NO=%s&PLATFORM_TYPE=%s';



uw_wmo_fn   = [dirs.cal,'UW_WMO_list.mat'];
sio_wmo_fn  = [dirs.cal,'SIO_WMO_list.mat'];
whoi_wmo_fn = [dirs.cal,'WHOI_WMO_list.mat'];

% URL TO SHARON & TANYA FLOAT STATUS TABLES
%dirs.sio_stats    = 'http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html';
dirs.sio_stats    = 'http://soccom.ucsd.edu/SOCCOM_float_performance.html';
dirs.mbari_stats  = ['https://www3.mbari.org/chemsensor/', ...
    'MBARI_float_table/mbari_float_table.html'];
dead_float_exp    = '^D|^dead|^MPD|^recover|^R'; % change as needed

% *************************************************************************
%                  REGULAR EXPRESSION FILTERS
% *************************************************************************
% MBARI ID starts with w (WHOI), s (SIO), or u (UW) followed by floatype
% char [ans] then at least 3 numbers
mbari_id_exp = '^[uws][ans]\d{3}\d*';  % MBARI ID

% MBARI FLOAT NAMES TO EXCLUDE ON FLoatVizconfig.txt - Ignore these on master list
exclude_expr = 'MRY|^old|^Cert|ua6962|ua6963|';

% SPECIAL CASE FLOATS [mbari_id, WMO, msg_dir, tf_bfile
special_case(1,:) = {'un0412' , 'NO_WMO_' , 'duplicate\UW\n0412\', 0};
special_case(2,:) = {'ua8501C', '5904680', 'duplicate\UW\f8501\', 1};
special_case(3,:) = {'ua8514H', '5904172', 'duplicate\UW\f8514_2\', 1};
special_case(4,:) = {'ua6966' , 'NO_WMO_' , 'duplicate\UW\f6966\' ,0};
special_case(5,:) = {'un0948' , 'NO_WMO_' , 'duplicate\UW\n0948\', 0};
special_case(6,:) = {'un0948B', 'NO_WMO_' , 'UW\n0948\', 0};
special_case(7,:) = {'un1173' , 'NO_WMO_' , '\\atlas\chem\tripleO2float\n1173\', 0};


% *************************************************************************
% *************************************************************************
% CHECK FOR CONFIG DIR, MSG DIR, INSTITUTION WMO LISTINGS
% IF NOT FOUND (i.e networks down) EXIT, OTHERWISE BUILD MASTER LIST
% *************************************************************************
% *************************************************************************
% TRY AND GET MBARI MASTER LIST
if ~isfolder(dirs.config)
    disp('Float Config directory could not be accessed - network may be down!')
    fprintf('%s\n',dirs.config);
    return
else
    tmp = dir([dirs.config,'*FloatConfig.txt']);
    config_names  = {tmp.name}';
    float_names   = regexprep(config_names,'_FloatConfig.txt','');
    
    % FILTER MASTER LIST WITH LOGICAL TESTS & PULL INST_ID FROM MBARI NAME
    t1 = ~cellfun(@isempty,regexp(float_names, mbari_id_exp,'once')); % desired floats
    t2 =  cellfun(@isempty,regexpi(float_names,exclude_expr,'once')); %exluded floats
    float_names  = float_names(t1&t2,:);
    config_names = config_names(t1&t2,:);
end

if ~isfolder(dirs.msg)
    disp(['Could not access root float *.msg dirs to determine float', ...
        'type at: ', dirs.msg]);
    return
end

% ************************************************************************
% **** GET INSTITUTION WMO LOOKUP TABLES ****
% ************************************************************************
disp('Getting institution ID & WMOs lookup tables....')

% TRY & GET UW WMO LISTING FROM PLAIN TEXT FILE: WMO...spaces...uw_id
% FOR DUPLICATE UW ID FLOAT DEPLOYMENTS OLDER UW ID'S ARE COMMENTED OUT IN
% DANA's LIST SO ONLY MOST RECENT WMO WILL BE RETRIEVED. OLDER ONES ARE
% HARDWIRED
try
    str = webread(dirs.UW_URL);
    txt_rows = regexp(str,'\n','split')'; % text rows to strings in cell array
    % Find comment lines, old duplicate float lines & remove first
    t1 = cellfun(@isempty, regexp(txt_rows,'^/|^\s+/|^\#|^\$','once'));
    txt_rows = txt_rows(t1);
    uw_wmo   = [regexp(txt_rows,'(?<=\d{7}\s+)\d+','match','once'), ...
        regexp(txt_rows,'\d{7}(?=\s+\d+)','match','once')]; % inst_id, wmo
    uw_wmo = uw_wmo(~cellfun(@isempty, uw_wmo(:,2)),:); %clear empty rows
    uw_wmo = uw_wmo(~cellfun(@isempty, uw_wmo(:,1)),:); %clear empty rows
    save(uw_wmo_fn,'uw_wmo')
    clear str txt_rows
catch
    if isfile(uw_wmo_fn)
        str = ['Could not access current UW WMO listing at: ', dirs.UW_URL];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        disp(['USING CACHED FILE: ',uw_wmo_fn])
        load(uw_wmo_fn)
    else
        disp('Could not find any UW WMO listing - EXITING!')
        return
    end
end

% TRY & GET SIO WMO LISTING HTML FILE (WMO in HREF string block)
try
    str       = webread(dirs.SIO_URL); %one huge string!
    tbl_entry = regexp(str,'</td>','split')'; % html table data strings to cell array
    %tbl_rows = regexp(str,'<tr>','split')';
    sio_wmo   = [regexp(tbl_entry,'\d+(?=</a>)','match','once'), ...
        regexp(tbl_entry,'\d{7}(?=a\.html)','match','once')]; % inst_id, wmo
    sio_wmo = sio_wmo(~cellfun(@isempty, sio_wmo(:,2)),:); %clear empty wmo rows
    sio_wmo = sio_wmo(~cellfun(@isempty, sio_wmo(:,1)),:); %clear empty sio id rows
    save(sio_wmo_fn,'sio_wmo')
    clear str tbl_entry
catch
    if isfile(sio_wmo_fn)
        str = ['Could not access current UW WMO listing at: ', dirs.SIO_URL];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        disp(['USING CACHED FILE: ',sio_wmo_fn])
        load(sio_wmo_fn)
    else
        disp('Could not find any SIO WMO listing - EXITING!')
        return
    end
end

% TRY & GET WHOI WMO LISTING (WMO in HREF string block)
% only want WHOI float ID block from html page which is the last table block
% try
%     str = webread(dirs.WHOI_URL);
%     ind = regexp(str,'Float Serial Number','once'); % start WHOI ID block
%     str = str(ind:end);
%     tbl_entry = regexp(str,'</td>','split')'; % html table row strings to cell array
%     whoi_wmo  = [regexp(tbl_entry,'\d+(?=</a>)','match','once'), ...
%         regexp(tbl_entry,'(?<=wmo/)\d{7}','match','once')]; % inst_id, wmo
%     whoi_wmo   = whoi_wmo(~cellfun(@isempty, whoi_wmo(:,2)),:); %clear empty wmo rows
%     whoi_wmo   = whoi_wmo(~cellfun(@isempty, whoi_wmo(:,1)),:); %clear empty wmo rows
%     save(whoi_wmo_fn,'whoi_wmo')
%     clear str ind tbl_entry
% catch
%     if isfile(whoi_wmo_fn)
%         str = ['Could not access current WHOI WMO listing at: ', dirs.WHOI_URL];
%         log_ct = log_ct+1;
%         log{log_ct} = str;
%         disp(str)
%         disp(['USING CACHED FILE: ',whoi_wmo_fn])
%         load(whoi_wmo_fn)
%     else
%         disp('Could not find any WHOI WMO listing - EXITING!')
%         return
%     end
% end

% ************************************************************************
%  LOAD FLOAT STATUS TABLE FOR DEAD FLOAT ASSIGNMENT SO FAR ONLY SOCCOM
%  & MBARI STATUS TABLES BUT OTHER TABLES CAN BE PARSED & ADDED WHEN
%  AVAILABLE / NEEDED
% ************************************************************************

% SHARONS'S HTML TABLE AT SIO
disp('Parsing & loading SOCCOM float stats table')
d = parse_html_table(dirs.sio_stats);
stats_hdr  = d.hdr;
stats_data = d.data;
clear d
sWMO   = strcmp(stats_hdr,'MBARI ODV_file WMO');
sSTAT  = strcmp(stats_hdr,'Status'); % using this for dead, op flag OK too?
sOP    = strcmp(stats_hdr,'Op Flag');
t_dead = ~cellfun(@isempty, regexp(stats_data(:,sSTAT), dead_float_exp,'once'));
dead_WMO = stats_data(t_dead, sWMO); % dead SOCCOM float WMO #'s


% TANYA'S HTML TABLE AT MBARI
disp('Parsing & loading MBARI float stats table')
d = parse_html_table(dirs.mbari_stats);
stats_hdr  = d.hdr;
stats_data = d.data;
clear d
sWMO   = strcmp(stats_hdr,'Argo Data WMO Number');
sSTAT  = strcmp(stats_hdr,'Status'); % using this for dead, op flag OK too?
t_dead = ~cellfun(@isempty, regexp(stats_data(:,sSTAT), dead_float_exp,'once'));
dead_WMO = [dead_WMO;stats_data(t_dead, sWMO)];

% ************************************************************************
% PULL INST_ID FROM MBARI NAME, SORT BY FLOAT NAMES BY INSITUTE & THEN
% BY INSTITUTE ID NUMBER. PREDIM VARIABLES & GET INDICES, LOAD SOCCOM &
% BGC LIST FILES
% ************************************************************************
INST_ID  = regexp(float_names,'\d+','match', 'once');
[~,IX]   = sort(str2double(INST_ID));
tmp_flt  = float_names(IX);
tmp_cfg  = config_names(IX);
tmp_inst = INST_ID(IX);
tUW      = ~cellfun(@isempty,regexp(tmp_flt,'^u','once'));
tSIO     = ~cellfun(@isempty,regexp(tmp_flt,'^s','once'));
tWH      = ~cellfun(@isempty,regexp(tmp_flt,'^w','once'));
mbari_id = [tmp_flt(tSIO); tmp_flt(tUW); tmp_flt(tWH)];
INST_ID  = [tmp_inst(tSIO); tmp_inst(tUW); tmp_inst(tWH)];
cfg_fnames = [tmp_cfg(tSIO); tmp_cfg(tUW); tmp_cfg(tWH)];
r           = size(mbari_id,1);

clear float_names IX tmp_flt tmp_inst tmp_cfg tSIO tWN t1 t2

% PREDIMNSION ARRAYS FOR NEXT STEPS
WMO_ID       = cell(r,1); % Predim for next step
float_type   = cell(r,1); % Predim, i.e. APEX or NAVIS or SOLO
message_dirs = cell(r,1); % Predim, i.e. APEX or NAVIS or SOLO
NC_template  = ones(r,1)*-1; % value code for netcdf skeleton to use
program      = cell(r,1); % program name
region       = program; % ocean region descriptor
tf_Bfile     = ones(r,1); % 1 = send to aoml
GPS_fix      = ones(r,3) * NaN; % lon lat SDN
tf_dead      = ones(r,1)*0; % 1 = dead float
tf_NO3       = ones(r,1)*0; % 1 = send to aoml
tf_pH        = ones(r,1)*0; % 1 = send to aoml
tf_Chl       = ones(r,1)*0; % 1 = send to aoml

MaxMsgF      = ones(r,1)* NaN; 
MaxMsgFD     = ones(r,1)* NaN; 
LastMsgF     = ones(r,1)* NaN; 
LastMsgFD    = ones(r,1)* NaN; 
MaxProcCyc   = ones(r,1)* NaN; 
MaxProcCycD  = ones(r,1)* NaN;


hdr = {'MBARI ID' 'INST ID' 'WMO' 'float type' '1st lon' '1st lat' ...
    '1st date' 'msg dir' 'NC template' 'Program' 'Region', ...
    'tf Bfile' 'tf Dead' 'tf NO3' 'tf pH' 'tf Chl'};
hdr = [hdr, 'max cycle msg file' 'max cycle file date' 'latest cycle msg file', ...
    'latest cycle file date' 'max cycle proc' 'max cycle proc date'];

I.MBA = strcmp('MBARI ID', hdr);  % DEFINE INDICES
I.INS = strcmp('INST ID', hdr);
I.WMO = strcmp('WMO', hdr);
I.TYP = strcmp('float type', hdr);
I.LON = strcmp('1st lon', hdr);
I.LAT = strcmp('1st lat', hdr);
I.SDN = strcmp('1st date', hdr);
I.MSG = strcmp('msg dir', hdr);
I.NC  = strcmp('NC template', hdr);
I.PRG = strcmp('Program', hdr);
I.REG = strcmp('Region', hdr);
I.BFL = strcmp('tf Bfile', hdr);
I.D   = strcmp('tf Dead', hdr);
I.NO3 = strcmp('tf NO3', hdr);
I.PH  = strcmp('tf pH', hdr);
I.CHL = strcmp('tf Chl', hdr);

I.MMF  = strcmp('max cycle msg file', hdr);
I.MMFD = strcmp('max cycle file date', hdr);
I.LMF  = strcmp('latest cycle msg file', hdr);
I.LMFD = strcmp('latest cycle file date', hdr);
I.MPC  = strcmp('max cycle proc', hdr);
I.MPCD = strcmp('max cycle proc date', hdr);

out   = cell(r, size(hdr,2));

clear ans fid d t1 t2 float_exp exclude_exp FloatViz_names IX
clear inst_num tSIO tWH float_names

% ************************************************************************

% *************************************************************************
% *************************************************************************
% STEP THROUGH INSTITUTE ID's & FILL IN WMO #'s & PROGRAM & REGION
disp('Assigning WMOs to MBARI IDs .....')
for i = 1:r
    float_name = mbari_id{i};
    config_name = cfg_fnames{i};
    inst_num   = regexp(float_name,'\d+','once','match');
    wmo        = '';
    
    tSpecial = strcmp(special_case(:,1), float_name);
    
    if sum(tSpecial) == 1
        if strcmp(special_case{tSpecial,2},'NO_WMO_')
            wmo = {[special_case{tSpecial,2}, float_name]};
            str = ['NO WMO_ID for ',float_name,'. WMO was never assigned'];
            log_ct = log_ct+1;
            log{log_ct} = str;
            disp(str)
        else
            wmo = special_case(tSpecial,2);
        end
    elseif regexp(float_name,'^u','once') % Univ of Wash WMO lookup
        t1  = strcmp(inst_num, uw_wmo(:,1));
        wmo = uw_wmo(t1,2);
%     elseif regexp(float_name,'^w','once') % WHOI WMO lookup
%         t1  = strcmp(inst_num, whoi_wmo(:,1));
%         wmo = whoi_wmo(t1,2);
    elseif regexp(float_name,'^w','once') % GET WMO WITH API CALL
        if float_name(2) == 'n'
            flt_type = 'NAVIS_EBR';
        elseif float_name(2) == 'a'
            flt_type = 'APEX'; % CHECK IF MORE TO STRING?
        end
        
        % TRY AND GET WMO USING INSTITUTE ID & FLOAT TYPE
        api_call = sprintf(api_call_str, dirs.WHOI_DOM,'wmo', inst_num, flt_type);
        try
            info = webread(api_call);
            wmo = {info.WMO};
        catch
            wmo ={['NO_WMO_', mbari_id]};
        end
    elseif regexp(float_name,'^s','once') % SIO WMO lookup
        t1  = strcmp(inst_num, sio_wmo(:,1));
        wmo = sio_wmo(t1,2);
    end
    
    if ~isempty(wmo) % remove duplicates but keep order (WHOI issue)
        wmo = unique(wmo,'stable');
        WMO_ID{i} = wmo{1,1};
    else
        str = ['WMO_ID for ',float_name,' is not assigned yet'];
        wmo = {['NO_WMO_', float_name]}; % assign temporary
        WMO_ID{i} = wmo{1,1};
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
    end
    
    % CHECK TO SEE IF FLOAT IS DEAD (ONLY TESTING SOCCOM FOR NOW)
    if sum(strcmp(dead_WMO, wmo{1,1}),1) == 1
        tf_dead(i) = 1;
    end
    
    % CHECK TO SEE IF NO3 CAL EXISTS
    if isfile([dirs.cal,'NO3_CAL\',float_name,'.cal'])
        tf_NO3(i) = 1;
    end
    
    % LOOK FOR PROGRAM & REGION STRINGS & SENSOR INFO IN FLOAT CONFIG FILES
    fid = fopen([dirs.cal,'FLOAT_CONFIG\',config_name]);
    tline = '';
    while ischar(tline)
        if regexp(tline,'^Program','once')
            program{i} = regexp(tline,'(?<=^Program\:\s+)[\w+/-\s]+', ...
                'match','once');
        elseif regexp(tline,'^Region','once')
            region{i} = regexp(tline,'(?<=^Region\:\s*)\w+','match','once');
        elseif regexp(tline,'^pH sensor|number of cal','once')
            tf_pH(i) = 1;
        elseif regexp(tline,'^MCOMS|FLBB','once')
            tf_Chl(i) = 1;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end
clear wmo tline fid t1

out(:,I.MBA) = mbari_id;
out(:,I.INS) = INST_ID;
out(:,I.WMO) = WMO_ID;
out(:,I.PRG) = program;
out(:,I.REG) = region;
out(:,I.D)   = num2cell(tf_dead);
out(:,I.NO3) = num2cell(tf_NO3);
out(:,I.PH)  = num2cell(tf_pH);
out(:,I.CHL) = num2cell(tf_Chl);

% *************************************************************************
% DETERMINE FLOAT TYPE AND ADD TO LIST BY CHECKING MSG DIR NAME
% LOOP THROUGH THE FLOAT NAMES AGAIN
% *************************************************************************
disp('Assigning float types (i.e APEX, NAVIS or SOLO .....')

% GET LIST OF MSG FILE DIRS FOR EACH INSTITUTION
d       = dir([dirs.msg,'UW\']);
uw_dirs = regexpi({d.name}','^[fns]\d.+','once','match');
uw_dirs = uw_dirs(~cellfun(@isempty, uw_dirs)); % remove empty cells

d         = dir([dirs.msg,'WHOI\']);
whoi_dirs = regexpi({d.name}','^w[ans]\d+','once','match');
whoi_dirs = whoi_dirs(~cellfun(@isempty, whoi_dirs)); % remove empty cells

d        = dir([dirs.msg,'SIO\']);
sio_dirs = regexpi({d.name}','^[ans]\d+','once','match');
sio_dirs = sio_dirs(~cellfun(@isempty, sio_dirs)); % remove empty cells
clear d

% STEP THROUGH MASTER LIST & DO SOME MATCH UPS
for i = 1:r
    float_name = mbari_id{i};
    inst_name  = INST_ID{i};
    clear msg_dirs dp mdir
    
    % IDENTIFY INSTITUTION PATH
    if regexp(float_name,'^u','once') % Univ of Wash msg dir
        msg_dirs  = uw_dirs;
        dp = [dirs.msg,'UW\'];
        inst_str = '';
    elseif regexp(float_name,'^w','once') % WHOI msg dir
        msg_dirs  = whoi_dirs;
        dp = [dirs.msg,'WHOI\'];
        inst_str = 'w';
    elseif regexp(float_name,'^s','once') % SIO msg dir
        msg_dirs  = sio_dirs;
        dp = [dirs.msg,'SIO\'];
        inst_str = '';
    end
    
    tf_dir   = ~cellfun(@isempty, regexp(msg_dirs,[inst_str,'[fnsa]', ...
        inst_name], 'once'));
    tSpecial = strcmp(special_case(:,1), float_name);
    
    if sum(tSpecial) == 1 % special case floats take priority
        special_dir = special_case{tSpecial,3};
        if strncmp(special_dir,'\\',2) % Entire path
            message_dirs{i} = special_dir;
        else
%            message_dirs{i} = [dp, special_dir];
            message_dirs{i} = [dirs.msg, special_dir]; %TM 4/12/21; duplicates are in .../duplicate/UW/... (not .../UW/duplicate/...)
        end
        float_type{i} = regexp(special_dir, '[fans](?=\d{3}\d+)', ...
            'once','match');
    elseif sum(tf_dir) == 1
        mdir = msg_dirs{tf_dir};
        message_dirs{i} = [dp, mdir,'\'];
        float_type{i} = regexp(mdir, '[fans](?=\d{2}\d+)','once','match');
    elseif sum(tf_dir) == 0 % No directory match
        str = ['could not find float directory on seaecho for : ',float_name];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        message_dirs{i} = 'UNKNOWN';
        float_type{i}   = 'UNKNOWN';
    end
end

float_type = regexprep(float_type,'wa|f|a','APEX');
float_type = regexprep(float_type,'wn|n','NAVIS');
float_type = regexprep(float_type,'s','SOLO');

out(:,I.TYP) = float_type;
out(:,I.MSG) = message_dirs;

% *************************************************************************
% *************************************************************************
% DETERMINE WHICH NETCDF BFILE SKELETON TO USE. TYPES ARE BASED ON SENSOR
% & FLOAT DESGIN
% *************************************************************************
% *************************************************************************
hdr_start = '$ p t s '; % ALL APEX & NAVIS HEADERS START WITH THESE VARS

% ************************************************************************
% APEX FLOAT DEFINITIONS: annie type, hdr vars, FlbbMode, PalMode, flt filter
% -99 is a fill value
% ************************************************************************
% APEX TYPE 1:
% Aanderaa 4330, ISUS, DURA, FLBB (Example: 9313)
APEX_types(1,:) = {1,'Topt TPhase RPhase no3 pH(V) FSig BbSig TSig', 1, -99, ''};

% APEX TYPE 2:
% Aanderaa 4330, ISUS, DURA FLBB data rows but FlbbMode(0)(Example: 9662)
APEX_types(2,:) = {2,'Topt TPhase RPhase no3 pH(V) FSig BbSig TSig', 0, -99, ''};

% APEX TYPE 3:
% Aanderaa 4330, ISUS, FLBB (NO pH, Example: 7567)
APEX_types(3,:) = {3,'Topt TPhase RPhase no3 FSig BbSig TSig',1, -99, ''};

% APEX TYPE 4:
% Aanderaa 4330, DURA (NO NO3 or FLBB, Example: 9018)
APEX_types(4,:) = {4,'Topt TPhase RPhase pH(V)', -99, -99, ''};

% APEX TYPE 5:
% Aanderaa 3830, ISUS , FLBB (older O2 optode & NO pH, Example: 7613)
APEX_types(5,:) = {5,'BPhase Topt no3 FSig BbSig TSig',1, -99, ''};

% APEX TYPE 6:
% Aanderaa 4330 & ISUS (NO pH or FLBB, Example: 7550)
APEX_types(6,:) = {6,'TPhase Topt no3',-99, -99, ''};

% APEX TYPE 7:
% Aanderaa 3830 & ISUS (older O2 optode, NO pH or FLBB, Example: 6972)
APEX_types(7,:) = {7,'BPhase Topt no3',-99, -99, ''};

% APEX TYPE 8:
% Aanderaa 4330, ISUS, FLBB (NO pH, Example: 7601, uses type 3 skeleton)
APEX_types(8,:) = {8,'TPhase Topt no3 FSig BbSig TSig',1, -99, ''};

% APEX TYPE 9:
% Aanderaa 4330 & ISUS (NO pH or FLBB, Example: 7663, uses type 6 skeleton)
APEX_types(9,:) = {9,'Topt TPhase RPhase no3',-99, -99, ''};

% APEX TYPE 10:
% Aanderaa 4330 & DURA (NO NO3 or FLBB, Example: 8374, uses type 4 skeleton)
APEX_types(10,:) = {10,'TPhase Topt pH(V)',-99, -99, ''};

% APEX TYPE 11: NO SKELETON YET - NOT GOING TO ARGO
% Aanderaa 4330 ONLY (NO NO3,FLBB or pH, Example: 12472)
APEX_types(11,:) = {11,'Topt TPhase RPhase',-99, -99, ''};

% APEX TYPE 12: NO SKELETON YET - NOT GOING TO ARGO
% Aanderaa 4330, DURA, FLBB (NO NO3, Example: 12792, PAL sensor equiped)
APEX_types(12,:) = {12,'Topt TPhase RPhase pH(V) FSig BbSig TSig',1, 1, ''};

% ************************************************************************
% NAVIS FLOAT DEFINITIONS: annie type, hdr vars, FlbbMode, PalMode, flt filter
% ************************************************************************
% NAVIS TYPE 1:
% NAVIS SBE63 O2, SUNA, DURA, MCOMS (CHL & BBP & FDOM), Example: 0508
NAVIS_types(1,:) = {1,'no3 O2ph O2tV Fl Bb Cdm phV phT', -99, -99, ''};

% NAVIS TYPE 2:
% NAVIS SBE63 O2, SUNA, MCOMS (CHL & BBP700 & FDOM), NO pH, Example: 0037
NAVIS_types(2,:) = {2,'no3 O2ph O2tV Fl Bb Cd', -99, -99, ''};

% NAVIS TYPE 3: BBP532 instead of FDOM
% NAVIS SBE63, SUNA, DURA, MCOMS (CHL & BBP700 & BBP532), Example: 0565
NAVIS_types(3,:) = {3,'no3 O2ph O2tV Fl Bb Cdm phV phT', -99, -99, '0565'};

% NAVIS TYPE 4: Andrea's OSP floats. can use type 1 skeleton
% NAVIS SBE63, SUNA, DURA, MCOMS (CHL & BBP700 & BBP532), No pH,Example:0949
NAVIS_types(4,:) = {4,'no3 O2ph O2tV Mch1 Mch2 Mch3 phV phT', -99, -99, ''};

% NAVIS TYPE 5:
% NAVIS SBE63, SUNA, DURA, MCOMS (CHL & BBP700 & FDOM), No NO3,Example:0412
NAVIS_types(5,:) = {5,'O2ph O2tV Fl Bb Cdm phV phT', -99, -99, ''};

% NAVIS TYPE 6: Newer NAVIS deployments, beginning 111X series (Oct2020).
% Can use type I skeleton. (May include I-params for pH in the future).
% NAVIS SBE63, SUNA, DURA, MCOMS (CHL & BBP700 & BBP532), includes pH
% diagnostics on spot samples.  pHT column also removed from msg file.
NAVIS_types(6,:) = {6,'no3 O2ph O2tV Mch1 Mch2 Mch3 phVrs phVk phIb pHIk', -99, -99, ''};


% ************************************************************************
% SOLO FLOAT DEFINITIONS: annie type, hdr vars, FlbbMode, PalMode, flt filter
% ************************************************************************
% THIS IS JUST A PLACE HOLDER FOR NOW
SOLO_types(1,:) = {6,'no3 O2ph O2tV Mch1 Mch2 Mch3 phVrs phVk phIb pHIk', -99, -99, ''};


% ************************************************************************
% STEP THROUGH MASTER LIST & ASSIGN ANNIE FLOAT TYPES FOR BR FILES
% ************************************************************************
disp(' ')
disp('Assigning Bfile skeleton types.....')
annie_types  = ones(size(out(:,1)))*0; % predim


for fct = 1:size(out,1)
    
    flbb_mode = -99;
    pal_mode  = -99;
    tf_type   = 0;
    tf_hdr    = 0;
    tf_gps    = 0;
    
    float_name = out{fct,I.MBA};
    inst_name  = out{fct,I.INS};
    msg_folder = out{fct,I.MSG};
    wmo        = out{fct,I.WMO};
    disp(float_name)
    
    tSpecial = strcmp(special_case(:,1), float_name);
    if sum(tSpecial) == 1 %special case
        tf_Bfile(fct) = special_case{tSpecial,4};
    elseif strncmp(wmo, 'NO_WMO',6) % no wmo assigned yet
        tf_Bfile(fct) = -1;
    end
    
    if strcmp(msg_folder, 'UNKNOWN')
        str = ['msg file directory unknown for: ',float_name, ...
            ' skipping type ID'];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        continue
    elseif ~isfolder(msg_folder)
        str = ['No msg directory found for: ',float_name, ...
            ' skipping type ID'];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        continue
    end
    
    if strcmp(out{fct,I.TYP},'APEX')
        FLOAT_types = APEX_types;
    elseif strcmp(out{fct,I.TYP},'NAVIS')
        FLOAT_types = NAVIS_types;
    elseif strcmp(out{fct,I.TYP},'SOLO')
        FLOAT_types = SOLO_types;
    else
        str = ['Unknown float type - skipping type ID for: ', float_name];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        continue
    end
    
    hdr_lines   = regexprep(FLOAT_types(:,2),' ',''); % remove spaces from hdr patterns
    
    % CHECK FOR SPECIAL CASE MCOMS (2 BBP SENSORS)
    % Check if float is an MCOMS special case if not isempty = 1
    tf_flt_filt = strcmp(FLOAT_types(:,5), inst_name);
    if sum(tf_flt_filt) == 0
        tf_flt_filt = cellfun(@isempty,FLOAT_types(:,5));
    end
    
    
    % GET MSG LIST WHICH WILL BE USED TO GET SOME FLOAT CYCLE STATS / LISTS
    % AND ANNIE FLOAt TYPE
    
    cal_info.msg_dir = msg_folder;
    cal_info.name    = float_name;
    d        = get_msg_list(cal_info, 'msg');
    dlist    = d.list;
    if isempty(dlist)
        str = ['msg directory exists for: ',float_name,' but no msg ',...
            'files inside - skiping type ID'];
        log_ct = log_ct+1;
        log{log_ct} = str;
        disp(str)
        continue
    end
    
    cycles   = str2double(regexp(dlist(:,1),'\d{3}(?=\.msg)','match','once'));
    tmax     = cycles == max(cycles); % max msg file cycle in respository
    fmod_sdn = cell2mat(dlist(:,3));
    tlast    = fmod_sdn == max(fmod_sdn); % % last msg file cycle in respository
    if sum(tlast) > 1 % multiple files with same time stamp ie ua8501C
        disp([float_name,' file listing has multiple identical last file ',...
            'time stamps'])
        tlast = tlast & cycles == max(cycles(tlast));
    end
    
    tAnnie   = cycles > 0 & cycles < 6; % look for 1st 6 cycles
    
    MaxMsgF(fct)   = cycles(tmax);
    MaxMsgFD(fct)  = fmod_sdn(tmax);
    LastMsgF(fct)  = cycles(tlast);
    LastMsgFD(fct) = fmod_sdn(tlast);
    clear cal_info d cycles tmax fmod_sdn tlast
    
    % GET LAST PROCESSED CYCLE NUMBER & DATE
    fp_mat = [dirs.mat, wmo,'\'];
    if isfolder(fp_mat)
        tmp = ls([fp_mat, wmo,'*.mat']);
        if ~isempty(tmp)
                    ProcCycList = cellstr(ls([fp_mat, wmo,'*.mat']));
            mat_cycles  = str2double(regexp(ProcCycList,'\d{3}(?=\.mat)', ...
                'match','once'));
            tmax        =   mat_cycles == max(mat_cycles);
            load([fp_mat, ProcCycList{tmax}],'INFO');
            
            MaxProcCyc(fct) = mat_cycles(tmax);
            MaxProcCycD(fct) = INFO.sdn;
            clear fp_mat ProcCycList mat_cycles tmax INFO
        else
            clear tmp
        end
    end
    
    % LOOK IN MSG FILES FOR INFO
    msg_files = dlist(tAnnie,1);
    %     if isempty(msg_files)
    %         str = ['msg directory exists for: ',float_name,' but no msg ',...
    %             'files inside - skiping type ID'];
    %         log_ct = log_ct+1;
    %         log{log_ct} = str;
    %         disp(str)
    %         continue
    %     end
    
    
    for ct = 1:size(msg_files,1) % Try 1st 3 msg files for header line
        msg_fn = msg_files{ct};
        fp = [msg_folder, msg_fn];
        if isfile(fp)
            temp_fp = [dirs.temp, msg_fn];
            copyfile(fp, temp_fp);
            fid = fopen(temp_fp);
            tline = ' ';
            while ischar(tline)
                % TEST FOR FLBBMODE LINE
                tflbb = regexp(tline,'(?<=^\$\s+FlbbMode\()\d+','match','once');
                if ~isempty(tflbb)
                    flbb_mode = str2double(tflbb);
                    if flbb_mode > 1 %17350, FlbbMode(255)
                        flbb_mode = 1;
                    end
                end
                % TEST FOR PALMODE LINE
                tpal = regexp(tline,'(?<=^\$\s+PalMode\()\d+','match','once');
                if ~isempty(tpal)
                    pal_mode = str2double(tpal);
                end
                
                % TEST FOR DISCRETE SAMPLE HEADER LINE
                thdr = regexp(tline,'^\$\s+p\s+t\s+s\s+', 'once','end');
                if ~isempty(thdr)
                    hdr_line = tline;
                    hdr_start = thdr;
                    tf_hdr = 1;
                    %break % header line found so stop looking for it
                end
                
                % TEST FOR GPS FIX LINE with valid position fix
                if regexp(tline, '(?<=^Fix:\s+)[-\d]','once')
                    gps_line = tline;
                    tf_gps = 1;
                    gps = regexp(gps_line,'\s+','split');
                end
                
                if tf_hdr == 1 && tf_gps == 1
                    break
                end
                
                tline = fgetl(fid);
            end
            fclose(fid);
            delete(temp_fp)
            clear tflbb tpal fid
            
            if tf_hdr == 0 && tline == -1 % no header & GPS line found try another msg file
                continue
            end
            
            % IF YOU GET HERE ALL IS GOOD - NOW TRY AND MATCH TO ANNIE TYPES
            % AND ADD PROFILE GPS INFO
            if size(gps,2) >= 3
                GPS_fix(fct,1:2) = str2double(gps(2:3)); % lon and lat recovered
            end
            
            if size(gps,2) >=5
                dstr = [gps{4},' ',gps{5}]; % date & time to one string
                GPS_fix(fct,3) = datenum(dstr,'mm/dd/yyyy HHMMSS');
            end
            
            hdr_str     = regexprep(hdr_line(hdr_start+1:end),' ',''); %remove spaces
            tf_flt_type = strcmp(hdr_lines, hdr_str);  % match header params
            tflbb       = cell2mat(FLOAT_types(:,3)) == flbb_mode;
            tpal        = cell2mat(FLOAT_types(:,4)) == pal_mode;
            ttype       = tf_flt_type & tflbb & tpal & tf_flt_filt;
            
            if sum(ttype) == 1 % A MATCH!
                annie_types(fct) = cell2mat(FLOAT_types(ttype,1));
                tf_type   = 1; % SUCCESS!
                break % get out of msg file loop & go to next float
            end
        else
            str = ['Msg file ', fp, ' not found'];
            log_ct = log_ct+1;
            log{log_ct} = str;
            disp(str)
        end
    end
    if tf_type == 0
        if tf_hdr == 1
            str = ['Header line found but no match for ',float_name, ...
                '. Float type = 0'];
            log_ct = log_ct+1;
            log{log_ct} = str;
            disp(str)
        else
            str = ['No header line found for ', float_name, ...
                ' (possibly no msg files yet). Float type = 0'];
            log_ct = log_ct+1;
            log{log_ct} = str;
            disp(str)
        end
    end
end
out(:,I.NC)  = num2cell(annie_types);
out(:,I.BFL) = num2cell(tf_Bfile);
out(:,I.LON) = num2cell(GPS_fix(:,1));
out(:,I.LAT) = num2cell(GPS_fix(:,2));
out(:,I.SDN) = num2cell(GPS_fix(:,3));

out(:,I.MMF)  = num2cell(MaxMsgF);
out(:,I.MMFD) = num2cell(MaxMsgFD);
out(:,I.LMF)  = num2cell(LastMsgF);
out(:,I.LMFD) = num2cell(LastMsgFD);
out(:,I.MPC)  = num2cell(MaxProcCyc);
out(:,I.MPCD) = num2cell(MaxProcCycD);



% *************************************************************************
% BUILD FLOAT LIST CELL ARRAY
% *************************************************************************
d.list = out;
d.hdr  = hdr;
d.log  = log(1:log_ct);

%clearvars -except d dirs out_name hdr

% *************************************************************************
% SAVE AS *.MAT & *.TXT
% *************************************************************************
save([dirs.cal, out_name,'.mat'],'d');

% FOR TEXT FILE CONVERT SDN TO DATE STR
% 1st profile SDN
plist = d.list;
tmp   = cell2mat(plist(:, I.SDN));
tnan  = isnan(tmp);
dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
plist(~tnan, I.SDN) = cellstr(dstr);
plist(tnan, I.SDN)  = cellstr('');

% Max msg cycle file mod SDN
tmp   = cell2mat(plist(:, I.MMFD));
tnan  = isnan(tmp);
dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
plist(~tnan, I.MMFD) = cellstr(dstr);
plist(tnan, I.MMFD)  = cellstr('');

% Last msg cycle file mod SDN
tmp   = cell2mat(plist(:, I.LMFD));
tnan  = isnan(tmp);
dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
plist(~tnan, I.LMFD) = cellstr(dstr);
plist(tnan, I.LMFD)  = cellstr('');

% Last processed cycle date
tmp   = cell2mat(plist(:, I.MPCD));
tnan  = isnan(tmp);
dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
plist(~tnan, I.MPCD) = cellstr(dstr);
plist(tnan, I.MPCD)  = cellstr('');

% FOR TEXT FILE PRINTING, CONVERT NaN's in LAT & LON & TO EMPTY
mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.LON));
plist(mask,I.LON) = {[]};
plist(mask,I.LAT) = {[]};

mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.MMF));
plist(mask,I.MMF) = {[]};

mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.LMF));
plist(mask,I.LMF) = {[]};

mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.MPC));
plist(mask,I.MPC) = {[]};

% hdr = [hdr, 'max cycle msg file' 'max cycle file date' 'last cycle msg file', ...
%     'last cycle file date' 'max cycle proc' 'max cycle proc date'];

r = size(d.list,1);
fid = fopen([dirs.cal,out_name,'.txt'], 'w');
fmt = '';
for i = 1:size(d.hdr,2)-1
    fprintf(fid,'%s\t', d.hdr{i});
    if regexp(d.hdr{i},'^NC|^tf|cycle msg|cycle proc(?!\s)','once')
        fmt = [fmt,'%0.0f\t'];
    elseif regexp(d.hdr{i},'^1st l','once')
        fmt = [fmt,'%0.3f\t'];
    else
        fmt = [fmt,'%s\t'];
    end
end
fmt = [fmt,'%s\r\n'];
fprintf(fid,'%s\r\n',d.hdr{i+1});

D = plist';
fprintf(fid, fmt, D{:}); %
fclose(fid);

stop_time = toc;
fprintf('MBARI master float list built in %2.1f minutes\n',stop_time/60);
% COPY LIST FILES TO MBARI BFIL & FTP DIRS
copyfile([dirs.cal,out_name,'*'], dirs.bfile);
copyfile([dirs.cal,out_name,'*'], dirs.ftp); % ACTIVATE ONCE LIVE




