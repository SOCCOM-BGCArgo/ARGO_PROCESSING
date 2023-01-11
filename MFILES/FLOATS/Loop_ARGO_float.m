function Loop_ARGO_float(update_str,filter_str,exclude_floats)

% ************************************************************************
% ------------------------------------------------------------------------
% ************************************************************************

% FUNCTION TO PROCESSES BGC ARGO FLOATS TO CREATE *.MAT FILES FOR EVERY PROFILE,
% CREATES ODV COMPATIBLE *TXT FILES FOR EACH FLOAT AND THEN COPIES THESE
% FILES TO CHEM AND SIRROCO FOR BOTH APEX AND NAVIS FLOATS

% INPUTS: update_str: 'update', or 'all' (use 'all' to reprocess all
%                     floats!)
%         filter_str: float subset when update mode is set to all (only
%                      specified floats are processed)
%         exclude_floats: floats to exclude from processing

% USE FILTER STRING TO SUBSET MAIN FLOAT LIST USING REGULAR EXPRESSIONS
% ONLY FLOATS MATCHING EXPRESSION WILL BE PROCESSED
%filter_str = '^9\d{3}|^12\d{3}|8514|8501';
%filter_str = '^0\d{3}'; % THIS WILL GET NAVIS FLOATS


% BUILD EXCLUSION LIST
% DUPLICATE OR ODD BALL FLOATS THAT ARE RUN ON THEIR OWN
% THIS COULD ALSO BE DEAD FLOATS THAT DON'T NEED TO BE RUN ANY MORE
% exlude_floats ={'ua8501C' 'ua8514H' 'ua8374'};


% AUTHOR: Josh Plant, MBARI
% UPDATES: 5/25/17, TLM, converted code from script to function
%   	08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
%       06/02/2020 - added code section at the end to process .dura files
%       into diagnostic text files for use in Eviz http://www3.mbari.org/chemsensor/eviz.htm
%       03/08/2021 TM Modifications to bring in line with the new MBARI master
%        float list and switch to WMO for processed file names.  Also, modified code structure in 'update' mode
%		 such that missing files that come in out of sequence actually get processed in real time (previously a full
%		 float refresh done manually was required for such cases).  Also enhanced the notification system to better
% 		 differentiate between such cases and msg files received in sequence.
%        5/16/2022: TM Remove MLR "flavor" of data files for upcoming snapshot archive!

% ************************************************************************
% ------------------------------------------------------------------------
% ************************************************************************

if ~isempty(filter_str)
    disp(['A filter string exists - only a subset of floats will be ',...
        'processed']);
    disp(['filter str = "',filter_str,'"']);
    disp('*********************  WARNING  **************************')
    disp(' DID YOU SHUT DOWN AUTOMATIC UPDATE-MODE PROCESSING JOB?? ')
    disp('*********************  WARNING  **************************')
    disp('')
    disp('PRESS ANY KEY TO CONTINUE / CTRL C to EXIT');
    
    pause
end

%initialize floats processed list:
floats_processed = [];
float_types = [];
kk=1;

% if ~exist('update_str', 'var') % if both lines commented, default to update
%     disp('Could not find update_str var: setting update_str = update')
%     update_str = 'update';
% end


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

% PATHS & NAMES - THESE ARE FOR MBARI PROCESSING PC!!!
% WHEN DEALING WITH ODD-BALL FLOATS SET THE msg DIR TO C:\temp\(UW_ID)\ AND
% THEN DUMP THE MESSAGE FILES IN THE DIR TO PROCESS

user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir = [user_dir, '\Documents\MATLAB\'];

dirs.mat       = [user_dir,'ARGO_PROCESSING\DATA\FLOATS\'];
dirs.cal       = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
%dirs.NO3config = [user_dir,'ARGO_PROCESSINGDATA\CAL\'];
dirs.FVlocal   = [user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ\'];
dirs.FV        = [user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ\'];
dirs.CANY = [user_dir,'ARGO_PROCESSING\DATA\CANYON\'];
%dirs.pk        = [dirs.FV,'PARK\'];
dirs.pk       = [user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ_ParkData\'];



%dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';

dirs.QCadj     = [user_dir,'ARGO_PROCESSING\DATA\CAL\QC_LISTS\'];
dirs.temp      = 'C:\temp\';
dirs.msg       = '\\seaecho.shore.mbari.org\floats\';
% dirs.msg       = '\\atlas\ChemWebData\floats\';
%dirs.msg       = 'C:\temp\';
%dirs.msg       = '\\atlas\chemwebdata\floats\duplicate\';
%dirs.config    = '\\atlas\Chem\ISUS\Argo\';

dirs.log = [user_dir,'ARGO_PROCESSING\DATA\Processing_logs\'];
dirs.bat = [user_dir,'batchfiles\']; %for STERNA & ARGOSY.  Is different on Josh's system.
dirs.float_stats = 'http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html';

% ************************************************************************
% ************************************************************************
% DO SOME PREP WORK
% ************************************************************************
% ************************************************************************

% ************************************************************************
% START LOG FILE
diary off % just to be sure it closed properly
log_file = ['ARGO_log_',datestr(now,'yyyymmddHHMM.txt')];
diary([dirs.log,log_file]); % start a processing log
disp(['Processing mode = ', update_str])
start_time = now;

% ************************************************************************
% CHECK FOR STUCK JOBS AND KILL IF NEEDED
d = get_system_jobs('MATLAB.exe');

iS = strcmp('Status',d.hdr);
iT = strcmp('Window Title',d.hdr);
iP  = strcmp('PID',d.hdr);

tS1 = strcmp('Running',d.list{1,iS}); % can fail in either state
tS0 = strcmp('Not Responding',d.list{1,iS});
tT = strcmp('MATLAB Command Window',d.list{1,iT}); % *.bat Matlab launch

PID = d.list{1,iP}((tS1|tS0)&tT); % Get process ID's

if size(PID,1) > 1 % at least one stale instance
    rPID = size(PID,1);
    sprintf('%f stuck Matlab processes found - trying to kill ...',rPID-1)
    kill_str = '';
    for i = 1:rPID-1
        kill_str = [kill_str,' /pid ',PID{i,1}];
    end
    kill_str =['taskkill',kill_str];
    system(kill_str);
else
    disp('NO OTHER MATLAB INSTANCES ARE RUNNING!')
end
clear d iS iT iP tS1 tS0 tT PID rPID kill_str

% ************************************************************************
% CONFIGURE EMAIL
setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
setpref('Internet','E_mail','tmaurer@mbari.org'); % define sender
%email_list ={'johnson@mbari.org';'tmaurer@mbari.org';'jplant@mbari.org'}; %5/16/19, Ken asked to be removed from processing email list.
email_list ={'tmaurer@mbari.org';'jplant@mbari.org';'eclark@mbari.org';'johnson@mbari.org'}; %TM 8/16/21 added Emily to recipients.
% email_list ={'tmaurer@mbari.org'};
new_msg_fn = ['new_argo_msgs',datestr(now,'yyyymmddHHMM.txt')];
bad_msg_fn = ['bad_argo_msgs',datestr(now,'yyyymmddHHMM.txt')];

% ************************************************************************
% DISPLAY SOME STARTING INFO
disp(' ')
disp(['ARGO PROCESSING STARTED ON ', ...
    datestr(start_time,'mm/dd/yyyy HH:MM:SS')]);

if strcmp(update_str, 'all')
    disp(['Processing mode is: "', update_str, ...
        '" All profiles will be processed.'])
elseif strcmp(update_str, 'update')
    disp(['Processing mode is: "', update_str, ...
        '" Only new profiles will be processed.'])
else
    disp('Processing mode is: unknown');
    return
end

% ************************************************************************
% GENERATE LIST OF FLOATS TO PROCESS & GENERATE PH VERSION LIST FOR DURA
% DIAGNOSTIC EVIZ CODE.
try
    d = MBARI_float_list([]);
    FLOAT_LIST = d.list;
catch
    disp('****MBARI FLOAT LIST BUILDER FAILED.****')
    disp('Trying to load pre-existing file...')
    load([dirs.cal,'MBARI_float_list.mat']);
    %     load('C:\temp_tanya\MBARI_float_list_NAVIS.mat') %TEMP for testing Navis
    FLOAT_LIST = d.list;
end
if ~isempty(d) % network down - couldn't make list?
    % get some indices
    iMSG =find(strcmp('msg dir', d.hdr)  == 1);
    iWMO = find(strcmp('WMO',d.hdr) == 1);
    iMB  = find(strcmp('MBARI ID',d.hdr) == 1);
    iINST = find(strcmp('INST ID',d.hdr) == 1);
    iFLT  = find(strcmp('float type',d.hdr) == 1);
    iTFPH   = find(strcmp('tf pH',d.hdr) == 1); % master list test
    iTFNO3  = find(strcmp('tf NO3',d.hdr) == 1);
    clear d
else
    disp('****UNABLE TO BUILD (OR LOAD PRE-EXISTING) MBARI MASTER FLOAT LIST.  EXITING!!!!****')
    return
end

% ************************************************************************
% TEMORARY CODE BLOCK TO EXCLUDE NAVIS FROM REPROCESSING LIST
% tf_NAVIS   = strcmp(FLOAT_LIST(:,4),'APEX'); % TESTING - EXCLUDE NAVIS
% FLOAT_LIST = FLOAT_LIST(~tf_NAVIS,:); % TESTING - EXLUDE NAVIS
% ************************************************************************

%dura pH ver code:
try
    disp('TRYING TO GET PH VERSION INFO FOR MBARI-CALIBRATED SENSORS.')
    dirs.ph_ver = get_ph_version;
catch
    disp('***get_ph_version FAILED!!** LOADING STATIC MAT FILE FROM ....\Documents\MATLAB\ARGO_PROCESSING\DATA\PH_DIAGNOSTICS')
    dirs.ph_ver = load([user_dir,'\ARGO_PROCESSING\DATA\PH_DIAGNOSTICS\MBARI_pHsensor_versions.mat']);
end

% % SELECT A SPECIFIC FLOAT OR SUBSET OF FLOATS
if ~isempty(filter_str)
    tm = regexpi(FLOAT_LIST(:,iWMO),filter_str);
    tm2 = cellfun(@isempty,tm);
    FLOAT_LIST_wmo = FLOAT_LIST(~tm2,:);
    if ~isempty(FLOAT_LIST_wmo)
        disp('WMOs DETECTED AS INPUT TO FILTER STR.')
    end
    jp = regexpi(FLOAT_LIST(:,iMB),filter_str);
    t1 = cellfun(@isempty,jp);
    FLOAT_LIST_instid = FLOAT_LIST(~t1,:);
    if ~isempty(FLOAT_LIST_instid)
        disp('INST_IDs DETECTED AS INPUT TO FILTER STR.')
    end
    FLOAT_LIST = [FLOAT_LIST_wmo;FLOAT_LIST_instid];
    clear jp t1 tm tm2 FLOAT_LIST_wmo FLOAT_LIST_istid
end

% GET BAD SENSOR LIST
bad_sensor_list = parse_bad_sensor_list([dirs.cal,'bad_sensor_list.txt']);
iM   = find(strcmp('MBARI ID STR',bad_sensor_list.hdr) == 1);

% ************************************************************************
% ************************************************************************
% PROCESS MSG FILES FOR FLOATS IN FLOAT LIST
% ************************************************************************
% ************************************************************************
new_msgs = {};
bad_msgs = {};
start_num = 1; % Starting row of FLOATLIST
stop_num  = size(FLOAT_LIST,1); % ENDING ROW OF FLOATLIST
for loop_ctr = start_num: stop_num
    try
        INST_ID_str    = FLOAT_LIST{loop_ctr,iINST}; % Institution ID str
        MBARI_ID_str = FLOAT_LIST{loop_ctr,iMB}; % MBARI ID str
        FLOAT_TYPE   = FLOAT_LIST{loop_ctr,iFLT}; % APEX or NAVIS or SOLO?
        WMO_ID_str   = FLOAT_LIST{loop_ctr,iWMO}; % WMO ID str
        
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
        
        % Exclude by MBARI name OR WMO - Skip if in exlude_floats
        if ~isempty(exclude_floats) && ...
                any(strcmp(FLOAT_LIST{loop_ctr,iMB},exclude_floats)) || any(strcmp(FLOAT_LIST{loop_ctr,iWMO},exclude_floats))
            continue
        end
        
        if strcmp('APEX',FLOAT_LIST{loop_ctr,iFLT})
            tf_float = Process_APEX_float(MBARI_ID_str, dirs, update_str);
        elseif strcmp('NAVIS',FLOAT_LIST{loop_ctr,iFLT})
            tf_float = Process_NAVIS_float(MBARI_ID_str, dirs, update_str);
        elseif strcmp('SOLO',FLOAT_LIST{loop_ctr,iFLT})
            tf_float = Process_SOLO_float(MBARI_ID_str, dirs, update_str);
        else
            disp(['Unknown float type for ',MBARI_ID_str,', WMO: ',WMO_ID_str,...
                '! processing next float'])
            continue
        end
        
        if tf_float.status == 1
            disp(['NEW FLOAT PROFILES WERE PROCESSED FOR ',MBARI_ID_str,', WMO: ',WMO_ID_str])
            floats_processed{kk} = MBARI_ID_str;
            float_types{kk} = FLOAT_TYPE;
            kk = kk+1;
            if ~isempty(tf_float.new_messages)
                new_msgs =[new_msgs;tf_float.new_messages];
            end
            
            if ~isempty(tf_float.bad_messages)
                bad_msgs =[bad_msgs;tf_float.bad_messages];
            end
        end
        
    catch ME
        disp(['ARGO PROCESSING FAILED ON FLOAT ',MBARI_ID_str, ...
            ' at ', datestr(now,'mm/dd/yyyy HH:MM:SS')]);
        new_msgs = [new_msgs; ['ARGO PROCESSING ERROR: ',MBARI_ID_str]];
        new_msgs = [new_msgs; ME.message];
        
        for k=1:length(ME.stack)
            ME.stack(k)
            new_msgs = [new_msgs; [ME.stack(k).file, ' LINE: ', ...
                num2str(ME.stack(k).line)]];
        end
        fclose('all');
    end
end
end_mat_time = now;



% ************************************************************************
% ************************************************************************
% CREATE ODV COMPATIBLE TXT FILES
% ************************************************************************
% ************************************************************************
for loop_ctr = start_num : stop_num % TM 3/3/21 - why is this a separate loop?
    try
        INST_ID_str    = FLOAT_LIST{loop_ctr,iINST};
        MBARI_ID_str = FLOAT_LIST{loop_ctr,iMB}; % MBARI ID str
        WMO_ID       = FLOAT_LIST{loop_ctr,iWMO};
        
        %TM 6/14/21 temporary(?) exclusion for OCR-test float.  Don't yet
        %write ODV files.
        % TM 2/7/22, added the first BSOLO to this exclusion...will likely
        % remove the exclusion soon!
       % if strcmp(WMO_ID,'5906446')==1 || strcmp(WMO_ID,'5906320')==1 %|| strcmp(WMO_ID,'4903026')==1
        %    continue
       % end
        
        if ~isempty(exclude_floats) && ...
                any(strcmp(FLOAT_LIST{loop_ctr,iMB},exclude_floats)) || any(strcmp(FLOAT_LIST{loop_ctr,iWMO},exclude_floats))
            continue
        end
        
        ODV_tf = argo2odv_LIAR(WMO_ID, dirs, update_str,1);
        if ODV_tf == 1
            disp(['A NEW ODV FILE WAS CREATED FOR ',MBARI_ID_str,', WMO: ',WMO_ID,'.']);
            
            %** IF THE ODV PROFILE FILE WAS UPDATED, UPDATE THE PARK DATA FILE TOO **
            try
                out = merge_PARK_mat(WMO_ID, dirs, 1);
                disp(['A NEW PARK DATA ODV FILE WAS CREATED FOR ',MBARI_ID_str,', WMO: ',WMO_ID,'.']);
                clear out
            catch ME
                disp(['ARGO PARK ODV FILE CREATION FAILED ON FLOAT ',MBARI_ID_str,', WMO: ',WMO_ID, ...
                    ' at ', datestr(now,'mm/dd/yyyy HH:MM:SS')]);
                new_msgs = [new_msgs; ['PARK ODV FILE CREATION ERROR: ',MBARI_ID_str,', WMO: ',WMO_ID]];
                new_msgs = [new_msgs; ME.message];
                
                for k=1:length(ME.stack)
                    ME.stack(k)
                    new_msgs = [new_msgs; [ME.stack(k).file, ' LINE: ', ...
                        num2str(ME.stack(k).line)]];
                end
            end
            
        end
        
    catch ME
        disp(['ARGO ODV(LIAR) FILE CREATION FAILED ON FLOAT ',MBARI_ID_str,', WMO: ',WMO_ID, ...
            ' at ', datestr(now,'mm/dd/yyyy HH:MM:SS')]);
        new_msgs = [new_msgs; ['ODV(LIAR) FILE CREATION ERROR: ',MBARI_ID_str,', WMO: ',WMO_ID]];
        new_msgs = [new_msgs; ME.message];
        
        for k=1:length(ME.stack)
            ME.stack(k)
            new_msgs = [new_msgs; [ME.stack(k).file, ' LINE: ', ...
                num2str(ME.stack(k).line)]];
        end
        fclose('all');
    end
    
    try
        ODV_tf = argo2odv_CANY(WMO_ID, dirs, update_str,1);
        if ODV_tf == 1
            disp(['A NEW ODV FILE WAS CREATED USING ',...
                'CANYON ALKALINITY']);
        end
        
    catch ME
        disp(['ARGO ODV(CANYON) FILE CREATION FAILED ON FLOAT ',MBARI_ID_str,', WMO: ',WMO_ID, ...
            ' at ', datestr(now,'mm/dd/yyyy HH:MM:SS')]);
        new_msgs = [new_msgs; ['ODV(CANYON) FILE CREATION ERROR: ',MBARI_ID_str,', WMO: ',WMO_ID]];
        new_msgs = [new_msgs; ME.message];
        
        for k=1:length(ME.stack)
            ME.stack(k)
            new_msgs = [new_msgs; [ME.stack(k).file, 'LINE: ', ...
                ME.stack(k).line]];
        end
        fclose('all');
    end
%% TM 5/16/22: Remove MLR "flavor" of data files for upcoming snapshot archive!!!    
% % % %     try
% % % %         ODV_tf = argo2odv_MLR(WMO_ID, dirs, update_str,1);
% % % %         if ODV_tf == 0
% % % %             %continue
% % % %         else
% % % %             disp(['A NEW ODV FILE WAS CREATED USING ',...
% % % %                 'Williams MLR ALKALINITY']);
% % % %         end
% % % %         
% % % %     catch ME
% % % %         disp(['ARGO ODV FILE CREATION FAILED ON FLOAT ',MBARI_ID_str,', WMO: ',WMO_ID, ...
% % % %             ' at ', datestr(now,'mm/dd/yyyy HH:MM:SS')]);
% % % %         new_msgs = [new_msgs; ['ODV(MLR) FILE CREATION ERROR: ',MBARI_ID_str,', WMO: ',WMO_ID]];
% % % %         new_msgs = [new_msgs; ME.message];
% % % %         
% % % %         for k=1:length(ME.stack)
% % % %             ME.stack(k)
% % % %             new_msgs = [new_msgs; [ME.stack(k).file, 'LINE: ', ...
% % % %                 ME.stack(k).line]];
% % % %         end
% % % %         fclose('all');
% % % %     end
end
end_txt_time = now;

% ************************************************************************
% ************************************************************************
% CREATE DURA TXT FILES FOR E-VIZ (code added June 2, 2020) TM
% CREATE ISUS TXT FILES FOR N-VIZ (code added Aug 17, 2020) TM
% ************************************************************************
% ************************************************************************
%
% If running in "all" mode, will want to reprocess all.
if strcmp(update_str, 'all')
    floats_processed = FLOAT_LIST(:,iMB); %MBARI ID str
    %     float_types = FLOAT_LIST(:,iFLT);
end

% If running in "update" mode, only reprocess the floats for which new msg
% files come in.
 if ~isempty(floats_processed)
     for ifp = 1:length(floats_processed)
         try
             MBARI_id_str = floats_processed{ifp}; % MBARI ID str
             tg = strcmp(FLOAT_LIST(:,iMB), MBARI_id_str) & cell2mat(FLOAT_LIST(:,iTFPH)) == 1 & ~strcmp( FLOAT_LIST(:,iFLT),'NAVIS');
             if sum(tg) == 1
                 disp(['NOW PROCESSING .DURA FILES FOR FLOAT ',...
                     MBARI_id_str,' AND GENERATING PH DIAGNOSTIC TXT FILE FOR USE IN E-VIZ.'])
                 dirs.save = [user_dir,'ARGO_PROCESSING\', ...
                     'DATA\PH_DIAGNOSTICS\'];
                 Merge_dura_msgs(MBARI_id_str, dirs)
             end
             
             %             thefloattype = float_types{ifp};
             disp(['DONE PROCESSING .DURA FILES FOR FLOAT ',...
                 MBARI_id_str,' AND GENERATING PH DIAGNOSTIC TXT FILE FOR USE IN E-VIZ.'])
             %             Loop_ARGO_dura(MBARI_id_str,thefloattype)
         catch ME
             disp(['Merge_dura_msgs failed on float ',MBARI_id_str,'.']);
             for k=1:length(ME.stack)
                 ME.stack(k)
                 ME.stack(k).file
                 ME.stack(k).line
             end
         end
         
         try
             MBARI_id_str = floats_processed{ifp}; % MBARI ID str
             %             thefloattype = float_types{ifp};
             tg = strcmp(FLOAT_LIST(:,iMB), MBARI_id_str) & cell2mat(FLOAT_LIST(:,iTFNO3)) == 1;
             if sum(tg) == 1
                 disp(['NOW PROCESSING .ISUS FILES FOR FLOAT ',...
                     MBARI_id_str,' AND GENERATING NO3 DIAGNOSTIC TXT FILE FOR USE IN N-VIZ.'])
                 dirs.save = [user_dir,'ARGO_PROCESSING\', ...
                     'DATA\NO3_DIAGNOSTICS\'];
                 Merge_isus_msgs(MBARI_id_str, dirs)
             end
             disp(['DONE PROCESSING .ISUS FILES FOR FLOAT ',...
                 MBARI_id_str,' AND GENERATING NO3 DIAGNOSTIC TXT FILE FOR USE IN N-VIZ.'])
             %             Loop_ARGO_isus(MBARI_id_str)
         catch ME
             disp(['Loop_ARGO_isus failed on float ',MBARI_id_str,'.']);
             for k=1:length(ME.stack)
                 ME.stack(k)
                 ME.stack(k).file
                 ME.stack(k).line
             end
         end
     end
 end

% ************************************************************************
% NOW REGENERATE ALL THE FLOATVIZ HTML AND FLOATVIZ LIST FILES.
% Im setting this up to get rebuilt with each processing, such that any new
% floats are automatically ingested into the web-based FloatVIZ system.
% However, we could easily modify this block to only run once per day, or
% perhaps only if there's a cycle '001' that gets processed?  But, it
% really doesn't take long at all to run, so just keep it simple for now...
% (TM, 03.24.21)
% ************************************************************************
try
    status = rewrite_XViz_html([]);
    stat_ind = find(cell2mat(status(:,2))==0);
    if ~isempty(stat_ind)
        for kstat = 1:length(stat_ind)
            theerror = status{stat_ind(kstat),1};
            disp('**********')
            disp(['** ERROR IN RE-GENERATING ',theerror,'. **'])
            disp('** TRY RUNNING rewrite_XViz_html LOCALLY. **')
            disp('**********')
        end
    else
        disp('** CONGRATS, ALL FLOATVIZ HTML FILES WERE RE-GENERATED SUCCESSFULLY. **')
    end
catch
    disp('*** ERROR RUNNING rewrite_XViz_html ***')
    disp('*** CHECK FLOATVIZ HTML FILES FOR ACCURACY! ***')
end
try
    status = rewrite_XViz_txt([]);
    stat_ind = find(cell2mat(status(:,2))==0);
    if ~isempty(stat_ind)
        for kstat = 1:length(stat_ind)
            theerror = status{stat_ind(kstat),1};
            disp('**********')
            disp(['** ERROR IN RE-GENERATING ',theerror,'. **'])
            disp('** TRY RUNNING rewrite_XViz_html LOCALLY. **')
            disp('**********')
        end
    else
        disp('** CONGRATS, ALL FLOATVIZ LIST FILES WERE RE-GENERATED SUCCESSFULLY. **')
    end
catch
    disp('*** ERROR RUNNING rewrite_XViz_txt ***')
    disp('*** CHECK FLOATVIZ LIST TXT FILES FOR ACCURACY! ***')
end

% ************************************************************************
% IF IN UPDATE MODE EMAIL LIST OF NEW FLOATS & CREATE TXT FILES
% ************************************************************************
% disp(['THESE ARE THE NEW MSGS, TESTING!!'])
if strcmp(update_str, 'update') && ~isempty(new_msgs)
    fid = fopen([dirs.log,new_msg_fn], 'w');
    if fid ~= -1
        fprintf(fid,'Name\tDate\tSize(kb)\t\r\n'); % header line
        for i = 1:size(new_msgs,1)
            fprintf(fid,'%s\r\n',new_msgs{i});
        end
    else
        disp(['Could not print to: ',[dirs.log,new_msg_fn]])
        fprintf('fid = %f and rows in new msgs variable = %f', ...
            fid,size(new_msgs,1));
    end
    disp('sending mail now')
    sendmail(email_list,'ARGOSY: NEW BGC ARGO MESSAGES PROCESSED', new_msgs)
end

% bad messages
if strcmp(update_str, 'update') && ~isempty(bad_msgs)
    fid = fopen([dirs.log,bad_msg_fn], 'w');
    if fid ~= -1
        fprintf(fid,'Name\tDate\tSize(kb)\t\r\n'); % header line
        for i = 1:size(bad_msgs,1)
            fprintf(fid,'%s\r\n',bad_msgs{i});
        end
    else
        disp(['Could not print to: ',[dirs.log,bad_msg_fn]])
        fprintf('fid = %f and rows in new msgs variable = %f', ...
            fid,size(bad_msgs,1));
    end
end

end_email_time = now;

% ************************************************************************
% LAST STEP - COPY FILES TO THE NETWORK
% ************************************************************************
disp(' ');
disp('COPYING FILES TO THE NETWORK........');
str = [ dirs.bat,'copyARGO2network_xcopy.bat'];
disp(str)
status = system(str); % Comment this line out if you don't want to copy to Chem!

end_copy_time = now;

disp(' ')
disp(['Message files Processed in ', ...
    num2str((end_mat_time-start_time)*24*60,'%0.2f'), ' minutes'])
disp(['FloatViz *.txt files created in ', ...
    num2str((end_txt_time-end_mat_time)*24*60,'%0.2f'), ' minutes'])
disp(['Email & new/bad message file lists created in ', ...
    num2str((end_email_time-end_txt_time)*24*60,'%0.2f'), ' minutes'])
disp(['Network copy jobs finished in  ', ...
    num2str((end_copy_time-end_email_time)*24*60, '%0.2f'), ' minutes'])
disp(['Overall update Processes completed in ', ...
    num2str((end_copy_time-start_time)*24*60,'%0.2f'), ' minutes'])
diary off








