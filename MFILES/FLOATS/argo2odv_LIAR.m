function tf_odv = argo2odv_LIAR(WMO_ID, dirs, update_str, HR_flag)

% ************************************************************************
% PURPOSE:
%    This function creates ODV compaitble text files used in FloatViz
%    and SOCCOMViz using the *.mat profile files. If pH data exists the
%    LIAR approach will be used to estimate alkalinity
%
%
% USAGE:
%	tf_odv = argo2odv10(WMO_ID, dirs, update_str)
%
% INPUTS:
%   WMO_ID  = WMO ID, as a string
%   update_str = "all" or "update"
%                 all    - to process all available msg files
%                 update - only process new msg files
%
%   HR_flag = 1 or 0, 1 to create HR txt files
%
%   dirs       = Either an empty variable ora structure with directory
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
%   dirs.QCadj      = path to QC adjustment list for all floats
%   dirs.FVlocal   = path to Floatviz file made by matlab go here
%
% CHANGE HISTORY
% 04/21/2017 - added code to include S & T ADJUSTED & QC variables for LR
%   and HR data.
% 04/25/2107 - add updated LIAR_V2 calc line - one more output variable
%   (unused) and function name change.
% 6/29/2017 - added code to change O2 values > 99990 to 99990, since such values were getting
%	reassigned to "fill-value" with erroneous QC flags, and thus getting rejected by the GDAC.
%   (See associated email correspondence with Annie Wong on 6/27)
% 07/19/17 added code to print comments from cal file if the exist
% 08/28/17 added code to re-assign dirs.msg directory for special case
%   floats which include duplicate UW ID floats and floats with NO WMO
% 10/09/2017, %add check for bad ALK values, bogging CO2SYS down
% 11/06/2017 made change to bias pH calculation - set Pres for any bad pH
%   equal to NaN during calculation so bias depth would not coincide with
%   a bad pH value - jp
% 12/18/2017 - Added code to deal with ARGO QC flag == 5 & to process chl
%   for in situ DC, NPQ corr, and chl_a corr with diff slope factore S of
%   30S. Also calculating POC & turning off ftp code to UMaine - jp
% 01/04/2017 - Fixed BBP532 bug: for backscatter there is no QC data so
%   fill adjusted data with raw data. This was done for BBP700 but not for
%   BBP532 (Presently this only affects float 0565). - JP
% 11/14/18 - Added "HR_flag" input to function to controll wether or not
%   HR txt files are created, -jp
% 01/10/19 - TM, Added code to deal with floats with all missing position
%   fixes (ie 12768, deployed too close to ice).  If no position info, try to
%   grab the position from the 000.msg file as estimate for first cycle.
%   This only goes in the ODV files.  QF will be = 4 (questionable).
% 11/05/19 - JP added code to fix minor BIAS pH calc bug. I previously assumed
%   if pressure value existed so would pH. Now check for valid pH value
%   too. This is a fix for spotty pH data in floats like 12888.
% 02/26/20 - JP Changed index searches (iP & look up tables)  argo2ODV*
%   for adjusted data from 'PRES' to 'PRES_ADJUSTED to account for changes
%   in Process_APEX_float & Process_NAVIS_float
% 05/11/2020 - TM modification to fix implemented on 01/10/19 (deals with floats with missing position fixes for cycle 1.  Now 2 floats in this category, 12768soocn and 12783eqpac). Bug fix.
% 08/25/2020 - TM; enhancement to header description (flagging for NPQ data is addressed).
% 12/13/20 - TM & JP - Forced all fopen writes to UTF-8, because that is the
%    new default for Matlab 2020 and better cross platform sharing
% 12/21/20 - JP, Added header line to txt files to alert ODV that format is UTF-8, //<Encoding>UTF-8</Encoding>
% 3/5/21 - TM, Code modifications to switch to WMO-based filename.
% 02/23/22 - JP added code to include BGC-SOLO, OCR APEX & 2X02 floats
% 05/17/22 - JP modified code to eliminate Chl_a_corr[mg/m^3], b_bp_corr[1/m]
%            & their associated QC columns from the generated text file
% 07/21/2022 - JP Updated QC Flagging code for ALK, DIC, pH25 for "sprof" type
%     floats (i.e SOLO). Bug related to lots of fill values + QF = 8 
%     not properly assigning ODV QF flags
% 07/25/2022 - TM Officially switched to CO2SYSv3 from "CO2SYSSOCCOM", although all constants used are the same (the latter inclued hard-wired inputs from N. Williams).
% 02/15/2023 - JP Added code block to look in master list for 1st cycle
%              position if missing & not found in 000.msg file
% 05/31/2023 - TM Added additional parameters for OCR443 and CHLA435
% 06/21/2023 - TM, incorporated ability to use psal-proxy (ARMOR3D) for LIAR calculation (for use in pco2)
% 02/01/2024 TM, modifications in support of ss4003 (different combo of OCR channels.
% 02/03/2024 JP, added code so raw OCR gets propagated to ADJ ODV file
%               similar to CHL, BBP & some general code cleanup/houskeeping
%5/2024 	TM - added ice evasion record column to ODV processed file.

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% TESTING
% WMO_ID     = '4903026'; % s0001
% WMO_ID     = '5906446'; %'ua19314' % OCR 5906446	APEX
% WMO_ID     = '5905131'; % ua12733	12733	5905131	APEX
% %WMO_ID = '5906481'; % ua19842 2XO2 float
% %WMO_ID = '5906482'; % ua19298 2XO2 float
% WMO_ID = '5906767'; 
% WMO_ID = '5906044';
% WMO_ID = '5906522';
% WMO_ID = '5906495';

% WMO_ID = 'NO_WMO_un0948';
% dirs       = [];
% update_str = 'all';
% HR_flag    = 1;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



% ************************************************************************
% SET UP DIRECTORIES AND PATHS
% ************************************************************************
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    dirs.NO3config = [user_dir,'CAL\'];
    dirs.FVlocal   = [user_dir,'FLOATVIZ\'];
    dirs.FV        = [user_dir,'FLOATVIZ\'];
    %dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    
    dirs.QCadj     = [user_dir,'CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    dirs.log = [user_dir,'Processing_logs\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end


% ************************************************************************
% LOAD MBARI FLOAT LIST
list_path  = [dirs.cal, 'MBARI_float_list.mat']; % file path
if exist(list_path, 'file') % load list variables
    load(list_path);
    FLOAT_LIST = d.list;
else
    disp('BUILDING FLOAT ID LIST ...')
    %[float_name, UW_ID, WMO_ID, float type]
    d = MBARI_float_list(dirs);
    FLOAT_LIST = d.list;
end

%Get Some indices:
iMSG = find(strcmp('msg dir', d.hdr)  == 1);
iWMO = find(strcmp('WMO',d.hdr) == 1);
iMB  = find(strcmp('MBARI ID',d.hdr) == 1);
iFLT = find(strcmp('float type',d.hdr) == 1);
iPRJ = find(strcmp('Program',d.hdr) == 1);
iREG = find(strcmp('Region',d.hdr) == 1);

iLON1 = find(strcmp('1st lon',d.hdr) == 1);
iLAT1 = find(strcmp('1st lat',d.hdr) == 1);
iSDN1 = find(strcmp('1st date',d.hdr) == 1);
clear d

% get index of float in master list
% i_tmp    = strfind(FLOAT_LIST(:,3),WMO_ID); %02/15/23 JP this is getting 2 matches for 0949
% iTHISFLT = find(not(cellfun('isempty',i_tmp)));
iTHISFLT = strcmp(FLOAT_LIST(:,3),WMO_ID);

% Set some definitions
dirs.msg     = FLOAT_LIST{iTHISFLT,iMSG};
MBARI_ID_str = FLOAT_LIST{iTHISFLT,iMB};
FLT_type     = FLOAT_LIST{iTHISFLT,iFLT};
PROJ_Name    = FLOAT_LIST{iTHISFLT,iPRJ};
REGION       = FLOAT_LIST{iTHISFLT,iREG};

% only gets used if missing 1st position not found by other means - jp
% 02/15/23
first_pos    = [1, cell2mat(FLOAT_LIST(iTHISFLT,[iSDN1,iLON1,iLAT1]))];

% ************************************************************************

tf_odv = 0;
% ************************************************************************
% CHECK IF FILE NEEDS TO BE UPDATED
% OR IF FLOAT IS PROBABLY DEAD
% ************************************************************************
if strcmp(update_str, 'update')
    % GET WMO# TO BUILD PATH TO *.MAT FILES FOR FLOAT
    
    data_dir   = [dirs.mat, WMO_ID,'\'];
    mat_files  = dir([data_dir,WMO_ID,'*.mat']);
    mat_cycles = regexp({mat_files.name},'(?<=\.)\d+(?=\.)','match');
    N_mat      = length(mat_cycles);
    
    clear list float_names ind ind data_dir
    clear mat_files mat_cycles
    
    % NOW CHECK FLOATVIZ FILE
    fv_name = [WMO_ID, '.TXT'];
    if exist([dirs.FVlocal, fv_name],'file')
        old_d = get_FloatViz_data([dirs.FVlocal, fv_name]);
        all_FV = unique(old_d.data(:,2));
        max_FV = max(old_d.data(:,2));
        
        %         if max_mat <= max_FV
        if N_mat == length(all_FV)
            fprintf(['LIAR FloatViz file: %s is current ', ...
                '(last cycle #: %03.0f)\r\n'], fv_name, max_FV);
            return
        end
    end
end

% ************************************************************************
% NEED TO DO RANGE CHECK FOR T & S NOT SAVED in *.mat files - THIS QC'ing
% HAS BEEN DONE ON Annie's end -- VALID BIO ARGO RANGE CHECK [MIN MAX]
RC.S    = [2 41]; % from argo parameter list
RC.T    = [-2.5 40]; % from argo parameter list
RC.P    = [0 12000]; % from argo parameter list
% ************************************************************************
% MERGE INDIVIDUAL MAT FILES FOR GIVEN BIO-ARGO FLOAT
% returns a structure: INFO, hdr, data
% ************************************************************************
tf_APEX_OCR = strcmp(WMO_ID,'5906320') | strcmp(WMO_ID,'5906446'); % ua19191, ua19314

if  ~strcmp(FLT_type,'SOLO') & ~tf_APEX_OCR % APEX, NAVIS but not OCR ONLY APEX
    d = merge_ARGO_mat(WMO_ID, dirs);
else % SOLO
    d = merge_SOLO_mat(WMO_ID, dirs);
end

if isempty(d)
    disp(['Could not merge data for ',WMO_ID,'. Probably no data '...
        'to merge - Process float data first.'])
    return
end
info    = d.INFO;
rhdr    = d.rhdr; % raw data
% keyboard
rdata   = d.rdata;

ahdr    = d.ahdr; % adjusted data
adata   = d.adata;

[dr,dc] = size(rdata); % raw and adjusted are the same size

hrrdata     = d.hrrdata; % if not an APEX float this will be empty
[hrdr,hrdc] = size(hrrdata);

% ************************************************************************
% ADD LAT QF COL AND INTERPOLATE MISSING POSITIIONS IF POSSIBLE
% QC FLAGS STILL ARGO AT THIS POINT
% INTERP QF FLAG = 3, POSITIONS STILL MISSING = NaN, QF = 99, GOOD = 1
% LOOK FOR MISSING PROFILES TOO
% ************************************************************************

% First, generate explanation of the IceEvRec column.  This goes in every file, regardless of Ice Algorithm presence or status!
iceEvRec_hdr_str = '"IceEvRec" = 1 indicates that the float was unable to surface as a result of the ice-avoidance algorithm. Note that this result does not necessarily signify the presence of sea-ice.';
iceEvRec_hdr_str2 = '"IceEvRec" = 0 indicates that the float was able to surface as a result of the ice-avoidance algorithm (and/or the absence of ice-avoidance algorithm software, ie for low-latitude floats)';
iceEvRec_hdr_str3 = 'For more information on the technical aspects of the ice-avoidance algorithm onboard profiling floats, please see Riser et al (2018)  https://doi.org/10.1002/2017JC013419';
fill_0     = ones(dr,1)*0; % for adding QF arrays
hr_fill_0  = ones(hrdr,1)*0; % for adding QF arrays

rhdr    = [rhdr(1:4),'Lat_QC', rhdr(5:dc)]; % raw data
rdata   = [rdata(:,1:4),fill_0+1 , rdata(:,5:dc)]; % ARGO QF's, START GOOD

ahdr    = [ahdr(1:4),'Lat_QC', ahdr(5:dc)]; % adjusted data
adata   = [adata(:,1:4),fill_0+1 , adata(:,5:dc)]; % ARGO QF's

if strcmp(info.float_type, 'APEX') & ~tf_APEX_OCR
    hrrdata   = [hrrdata(:,1:4),hr_fill_0+1 , hrrdata(:,5:hrdc)];
end

[~,ia,~] = unique(rdata(:,1)); % unique casts
pos_fix  = rdata(ia,1:4);   % cast, SDN, LON, LAT (position subset)

% ************************************************************************
% SCAN FOR MISSING FLOAT PROFILES
missing_profile_str = 'No missing float profiles';
missing_profiles =[];
for i = min(pos_fix(:,1)) : max(pos_fix(:,1))
    t1 = sum(pos_fix(:,1) == i,1);
    if t1 == 0
        missing_profiles = [missing_profiles, i];
    end
end
if ~isempty(missing_profiles)
    missing_profile_str = ['Missing Float profile(s) for station(s): ', ...
        sprintf('%0.0f ',missing_profiles)];
    disp(missing_profile_str);
end
clear i t1 missing_profiles

% ************************************************************************
% SCAN FOR MISSING POSITION FIXES

missing_pos_str = '';
interp_pos_str  = '';
%t_nan = isnan(pos_fix(:,4)); % any nan's in LAT?
cy1 = pos_fix(pos_fix(:,1)==1,:);
t_nan = isnan(cy1(:,4)); % any nans in first cycle LAT?
if ~isempty(cy1)
    if sum(t_nan,1)==length(cy1(:,4))  %first cycle missing position fix
        disp('POSITION INFO MISSING AT CYCLE 1...ATTEMPTING TO GRAB LAT/LON FROM 000.msg file...')
        %Try to retrieve the position info from cycle 000 if available.  Use
        %this as an estimate of the first cycle.
        %load cal info:
        load([dirs.cal,'cal',MBARI_ID_str,'.mat']);
        msglisting = get_msg_list(cal.info, 'msg');
        ilistMSG = find(strcmp('msg file', msglisting.hdr)  == 1);
        i_000 = strfind(msglisting.list(:,ilistMSG),'000.msg');
        itmp2 = find(not(cellfun('isempty',i_000)));
        if ~isempty(itmp2)
            msg000 = msglisting.list{itmp2,ilistMSG};
            %         if ~isempty(msg000)
            [LON,LAT] = get000msg_position([cal.info.msg_dir,msg000]);
            if ~isempty(LON) && ~isempty(LAT) %successful position grab
                pos_fix(1,3) = LON;
                pos_fix(1,4) = LAT;
                rdata(rdata(:,1)==1,3) = LON;
                rdata(rdata(:,1)==1,4) = LAT;
                rdata(rdata(:,1)==1,5) = 3; %if using 000.msg position info, mark lat/lon QF as "questionable"
                disp('Lat/Lon successfully extracted from 000.msg file for cycle 1');
            end
        elseif ~isnan(first_pos(:,3)) && ~isnan(first_pos(:,4)) % 1st fix in master list?
            pos_fix(1,3) = first_pos(:,3);
            pos_fix(1,4) = first_pos(:,4);
            rdata(rdata(:,1)==1,3) = first_pos(:,3);
            rdata(rdata(:,1)==1,4) = first_pos(:,4);
            rdata(rdata(:,1)==1,5) = 3; %if using 000.msg position info, mark lat/lon QF as "questionable"
            disp('Lat/Lon successfully extracted from master processing list')
        end
    end
end

t_nan = isnan(pos_fix(:,4)); % any nan's in LAT? Try again now that checked for 000 position info.
if sum(t_nan,1) > 0 && sum(t_nan,1)~=length(pos_fix(:,4)) %if not all nans, try to interpolate (ie if never a pos fix then no interpolation, move on)
    diff_nan   = [0;diff(t_nan)]; % 1 = start of NaN, -1 = start of good
    nan_start  = find(diff_nan == 1);
    nan_end    = find(diff_nan == -1) - 1;
    
    % some interp may be possible - resize to complete bounds only
    % # of start indices should equal # of end indices
    if size(nan_start,1) > size(nan_end,1) % more starts than ends, no last profile
        nan_start = nan_start(1:size(nan_end,1));
    elseif size(nan_start,1) < size(nan_end,1) % more ends than starts
        nan_end = nan_end(2:end); % 1st profile doesn't have a position
    elseif nan_end(1)<nan_start(1) && nan_start(end)>nan_end(end) % same number of starts, ends, but don't line up cuz missing positions at both beginning and end
        nan_end = nan_end(2:end);
        nan_start = nan_start(1:end-1);
    end
    
    for i = 1 : size(nan_start,1)
        % no upper bound, can't interpolate, move on
        if isempty(nan_end)
            break
        end
        
        % all good now, nan's bounded at this point, procede to interp
        if size(nan_start,1) == size(nan_end,1)
            bnds = pos_fix([nan_start(i)-1,nan_end(i)+1],:); % bounds
            bait = pos_fix(nan_start(i):nan_end(i),:); % data to interp
            % Meridian crossing - float unlikely to move more than 180 deg
            if abs(bnds(2,3) - bnds(1,3)) > 180
                t1 = bnds(:,3) < 180;
                bnds(t1,3) = bnds(t1,3) + 360; % temp add 360 to small side
            end
            missing_pos = [bait(:,1:2), ...
                interp1(bnds(:,2),bnds(:,3:4),bait(:,2))];
            % bring meridian crossing back to reality if need be
            t1 = missing_pos(:,3) > 360;
            missing_pos(t1,3) = missing_pos(t1,3) - 360;
            pos_fix(nan_start(i):nan_end(i),3:4) = missing_pos(:,3:4);
            for j = 1:size(missing_pos,1) % step through casts
                t1 = rdata(:,1) == missing_pos(j,1);
                rdata(t1,2:4) = ones(sum(t1),1) * missing_pos(j,2:4); % matrix
                rdata(t1,5)  = 3; %set QF to 3, ODV coversion takes it to 4
                if strcmp(info.float_type, 'APEX') & ~tf_APEX_OCR
                    t1 = hrrdata(:,1) == missing_pos(j,1);
                    hrrdata(t1,2:4) = ones(sum(t1),1) * missing_pos(j,2:4);
                    hrrdata(t1,5)  = 3; %set QF to 3, ODV coversion takes it to 4
                end
                
            end
        end
    end
    clear diff_nan nan_start nan_end bnds bait i j t1
    
    % NOW ASSES WHAT HAS BEEN INTERPOLATED AND WHAT HAS NOT
    t_nan2 = isnan(pos_fix(:,4)); % Still any nan's in LAT? NO INTERP YET
    t_nan3 = t_nan & ~t_nan2; % used to be nan's but now interpolated
    if sum(t_nan3 > 0) % Interpolated values
        interp_pos_str = ['Missing Float position interpolated for ', ...
            'station(s): ', sprintf('%0.0f ',pos_fix(t_nan3,1))];
        interp_pos_str = [interp_pos_str,'\r\n//Latitude quality flag = 4',...
            ' for interpolated float positions'];
        disp(interp_pos_str)
    end
    if sum(t_nan2 > 0) % Could not interpolate these casts
        missing_pos_str = ['No position for station(s): ', ...
            sprintf('%0.0f ',pos_fix(t_nan2,1))];
        disp(missing_pos_str)
        rdata(isnan(rdata(:,4)),5) = 99;  % BIO ARGO MISSING DATA QF VALUE
        if strcmp(info.float_type, 'APEX') & ~tf_APEX_OCR
            hrrdata(isnan(hrrdata(:,4)),5) = 99;
        end
        
    end
    adata(:,3:5) = rdata(:,3:5); % ADJUSTED POSITIONS EQUAL RAW
    %clear t_nan t_nan2 t_nan3
end
[dr,dc]     = size(rdata); % raw and adjusted are the same size
[hrdr,hrdc] = size(hrrdata);

% ************************************************************************
% LOOK FOR BIO ARGO MISSING VALUES & REPLACE WITH NaN
% ************************************************************************
t_99999 = rdata == 99999; % BIO ARGO MISSING DATA VALUE
rdata(t_99999) = NaN;

t_99999 = adata == 99999;
adata(t_99999) = NaN;

t_99999 = hrrdata == 99999;
hrrdata(t_99999) = NaN;

clear t_99999 t_99999_sum tmp;

% ************************************************************************
% ************************************************************************
% NOW ADD SOME DATA COLUMNS AND HEADERS NEEDED FOR THE ODV FILE
% PRES QC, TEMP QC, PSAL QC, density,  O2sat, pH25
% ************************************************************************
% ************************************************************************
% GET SOME INDICES (will be the same for raw and adjusted data)
iL  = find(strncmp('Lat [',rhdr,5) == 1); % set order =lat, lat QF, p, t, s
iP  = find(strcmp('PRES',rhdr)     == 1);
iT  = find(strcmp('TEMP',rhdr)     == 1);
iS  = find(strcmp('PSAL',rhdr)     == 1);

% 02/11/20 CHECK FOR PRESS_QC FIELD - IT MAY NOT EXIST IN OLDER MAT FILES
% THIS IS  MESSY WORK AROUND - mfiles need an overhaul
if sum(strcmp(rhdr,'PRES_QC'),2) == 0 % NO PRES QC YET - old file
    PRES_QF = fill_0 + 1;             % ARGO QF = GOOD
    PRES_QF(isnan(rdata(:,iP))) = 99; % BIO ARGO QF MISSING VALUE
else
    PRES_QF = max([rdata(:,iP+1) adata(:,iP+1)],[],2);
end

% ADJUSTED & RAW ARE THE SAME VALUES BUT FLAGS MAY ONLY BE SET IN ADJUSTED
% DATA SO LOOK AT BOTH AND TAKE THE GREATEST VALUE
TEMP_QF = max([rdata(:,iT+1) adata(:,iT+1)],[],2);
PSAL_QF = max([rdata(:,iS+1) adata(:,iS+1)],[],2);


potT   = theta(rdata(:,iP),rdata(:,iT),rdata(:,iS),0);
den    = density(rdata(:,iS),potT) -1000; % pot den anom (sigma-theta)
den_QF = fill_0 + 1;
den_QF(isnan(den)) = 99;
den_QF(TEMP_QF == 4 | PSAL_QF == 4) = 4;% out of range T or S so den bad too

% CALC DEPTH
float_z  = ones(size(rdata(:,iP)))*NaN;
nan_lat  = isnan(rdata(:,iL));
mean_lat = mean(rdata(:,iL),'omitnan');
% If pos fix available
float_z(~nan_lat) = sw_dpth(rdata(~nan_lat,iP),rdata(~nan_lat,iL));
% otherwise use average lat
float_z(nan_lat) = sw_dpth(rdata(nan_lat,iP),rdata(nan_lat,iP)*0+mean_lat);

%float_z_QF = fill_0 + 1;
float_z_QF = PRES_QF; %TM, ie cycle 7593.032.msg has bad pres values and they weren't getting propagated to depth(z)-qc.  Just use PRES_QC?? 8/1/22 
float_z_QF(isnan(float_z)) = 99;
clear nan_lat mean_lat

% ADD P,T,S QF's DENSITY AND DEPTH TO THE DATA SETS
raw_hdr = [rhdr(1:iL+2),rhdr(iP),'PRES_QC',rhdr(iT:iT+3), ...
    'SIGMA_THETA','SIGMA_THETA_QC','DEPTH','DEPTH_QC',rhdr(iS+2:dc)];

raw_data = [rdata(:,1:iL+2), rdata(:,iP), PRES_QF, rdata(:,iT), ...
    TEMP_QF, rdata(:,iS), PSAL_QF, den, den_QF, float_z, float_z_QF, ...
    rdata(:,iS+2:dc)];

adj_hdr = [ahdr(1:iL+2),ahdr(iP),'PRES_ADJUSTED_QC',ahdr(iT:iT+3), ...
    'SIGMA_THETA','SIGMA_THETA_QC','DEPTH','DEPTH_QC',ahdr(iS+2:dc)];

adj_data = [adata(:,1:iL+2), adata(:,iP), PRES_QF, adata(:,iT), ...
    TEMP_QF, adata(:,iS), PSAL_QF, den, den_QF, float_z, float_z_QF, ...
    adata(:,iS+2:dc)];

%%%%%%%%% START PSAL PROXY BLOCK  %%%%%%%%%%
% Now modify adjusted psal data, depending on whether a "psal proxy" float.  These columns will get reverted to original psal at the very end!

if isfield(info,'PSAL_PROXY_USED')
    if info.PSAL_PROXY_USED
	disp(['%%%%% NOTE: USING ARMOR3D PSAL PROXY FOR COMPUTING DERIVED PARAMS, FLOAT ',WMO_ID,'. %%%%'])
	cast_inds = ~isnan(d.psalprox); %these are indices of entries where psal-proxy was assigned in the merger.
	iaS  = find(strcmp('PSAL_ADJUSTED',adj_hdr) == 1); %let's be sure we grab the right index
% 	iacyc  = find(strcmp('Station',adj_hdr) == 1); %let's be sure we grab the right index
% 	xcycind = find(adj_data(:,iacyc)>cast_num_start);
% 	mycycs = unique(adj_data(xcycind,iacyc));
	psal_orig = adj_data(:,iaS); %for replacing later, after all the calls to LIAR, CO2SYS have been done.
	psalqc_orig = adj_data(:,iaS+1);
	adj_data(cast_inds,iaS) = d.psalprox(cast_inds);
    adj_data(cast_inds,iaS+1) = 1; %Argo Quality flag still at this point... make all = 1 (good) if using ARMOR3D
    end
end
%%%%%%%%% END PSAL PROXY BLOCK  %%%%%%%%%%

if strcmp(info.float_type, 'APEX') & ~tf_APEX_OCR%
    hrPRES_QF = hr_fill_0 + 1;             % ARGO QF = GOOD
    hrPRES_QF(isnan(hrrdata(:,iP))) = 99; % BIO ARGO QF MISSING VALUE
    hrPRES_QF(hrrdata(:,iP) < RC.P(1) | hrrdata(:,iP) > RC.P(2)) = 4;
    
    hrTEMP_QF = hrrdata(:,iT+1);
    hrTEMP_QF(isnan(hrrdata(:,iT))) = 99;
    hrTEMP_QF(hrrdata(:,iT) < RC.T(1) | hrrdata(:,iT) > RC.T(2)) = 4;
    hrPSAL_QF = hrrdata(:,iS+1);
    hrPSAL_QF(isnan(hrrdata(:,iS))) = 99;
    hrPSAL_QF(hrrdata(:,iS) < RC.S(1) | hrrdata(:,iS) > RC.S(2)) = 4;
    
    hrpotT   = theta(hrrdata(:,iP),hrrdata(:,iT),hrrdata(:,iS),0);
    hrden    = density(hrrdata(:,iS),hrpotT) -1000; % pot den anom (sigma-theta)
    hrden_QF = hr_fill_0 + 1;
    hrden_QF(isnan(hrden)) = 99;
    hrden_QF(hrTEMP_QF == 4 | hrPSAL_QF == 4) = 4;
    
    % CALC DEPTH
    hrfloat_z = sw_dpth(hrrdata(:,iP),hrrdata(:,iL));
    hrfloat_z_QF = hr_fill_0 + 1;
    hrfloat_z_QF(isnan(hrfloat_z)) = 99;
    
    hrraw_data = [hrrdata(:,1:iL+2), hrrdata(:,iP), hrPRES_QF, hrrdata(:,iT), ...
        hrTEMP_QF, hrrdata(:,iS), hrPSAL_QF, hrden, hrden_QF, hrfloat_z, hrfloat_z_QF, ...
        hrrdata(:,iS+2:hrdc)];
    
    [hrdr,hrdc] = size(hrraw_data); % get new size
    clear hrrdata hrPRES_QF hrTEMP_QF hrPSAL_QF hrden_QF hrden
    clear hrfloat_z_QF
end


[dr,dc] = size(raw_data); % get new size

clear t_MVI rhdr rdata ahdr adata den float_z float_z_QF
clear PRES_QF TEMP_QF PSAL_QF den_QF

% REDO INDICES
iP  = find(strcmp('PRES',raw_hdr) == 1);
iT  = find(strcmp('TEMP',raw_hdr) == 1);
iS  = find(strcmp('PSAL',raw_hdr) == 1);
iO  = find(strcmp('DOXY',raw_hdr) == 1);

% ************************************************************************
% ADD O2 % SAT IF O2 exists
if ~isempty(iO)

    
    %(umol/kg) / (umol/kg) * 100
    crazyO = raw_data(:,iO) == 99990; % If "crazy O2" conc value has been set (in Process_APEX/NAVIS_float), make O2 sat 99990 as well
    O2sat(crazyO,:) = 99990;
    O2sat(~crazyO,:)     = raw_data(~crazyO,iO) ./ oxy_sol(raw_data(~crazyO,iT), ...
        raw_data(~crazyO,iS),0) *100;
    O2sat_QF  = raw_data(:,iO+1);  %use oxygen conc QF
    
    crazyOadj = adj_data(:,iO) == 99990;
    O2sat_adj(crazyOadj,:) = 99990;
    O2sat_adj(~crazyOadj,:)    = adj_data(~crazyOadj,iO) ./ oxy_sol(adj_data(~crazyOadj,iT), ...
        adj_data(~crazyOadj,iS),0) *100;
    O2sat_adj_QF = adj_data(:,iO+1);
    
    raw_hdr = [raw_hdr(1:iO+1),'DOXY_%SAT', 'DOXY_%SAT_QC', ...
        raw_hdr(iO+2:dc)];
    raw_data = [raw_data(:,1:iO+1), O2sat, O2sat_QF, raw_data(:,iO+2:dc)];
    
    adj_hdr = [adj_hdr(1:iO+1),'DOXY_%SAT_ADJUSTED', ...
        'DOXY_%SAT_ADJUSTED_QC', adj_hdr(iO+2:dc)];
    adj_data = [adj_data(:,1:iO+1), O2sat_adj, O2sat_adj_QF, ...
        adj_data(:,iO+2:dc)];
    
    [dr,dc] = size(raw_data); % get new size
    clear O2sat O2sat_QF O2sat_adj O2sat_adj_QF
end

% CHECK FOR 2XO2 APEX FLOATS
iO2  = find(strcmp('DOXY2',raw_hdr) == 1);
if ~isempty(iO2)
    %(umol/kg) / (umol/kg) * 100
    crazyO = raw_data(:,iO2) == 99990; % If "crazy O2" conc value has been set (in Process_APEX/NAVIS_float), make O2 sat 99990 as well
    O2sat(crazyO,:) = 99990;
    O2sat(~crazyO,:)     = raw_data(~crazyO,iO2) ./ oxy_sol(raw_data(~crazyO,iT), ...
        raw_data(~crazyO,iS),0) *100;
    O2sat_QF  = raw_data(:,iO2+1);  %use oxygen conc QF
    
    crazyOadj = adj_data(:,iO2) == 99990;
    O2sat_adj(crazyOadj,:) = 99990;
    O2sat_adj(~crazyOadj,:)    = adj_data(~crazyOadj,iO2) ./ oxy_sol(adj_data(~crazyOadj,iT), ...
        adj_data(~crazyOadj,iS),0) *100;
    O2sat_adj_QF = adj_data(:,iO2+1);
    
    raw_hdr = [raw_hdr(1:iO2+1),'DOXY2_%SAT', 'DOXY2_%SAT_QC', ...
        raw_hdr(iO2+2:dc)];
    raw_data = [raw_data(:,1:iO2+1), O2sat, O2sat_QF, raw_data(:,iO2+2:dc)];
    
    adj_hdr = [adj_hdr(1:iO2+1),'DOXY2_%SAT_ADJUSTED', ...
        'DOXY2_%SAT_ADJUSTED_QC', adj_hdr(iO2+2:dc)];
    adj_data = [adj_data(:,1:iO2+1), O2sat_adj, O2sat_adj_QF, ...
        adj_data(:,iO2+2:dc)];
    
    [dr,dc] = size(raw_data); % get new size
    clear O2sat O2sat_QF O2sat_adj O2sat_adj_QF
end


[raw_r,raw_c] = size(raw_data); % get new size
[adj_r,adj_c] = size(adj_data); % get new size

% ************************************************************************
% ************************************************************************
% CACULATE VARIOUS FLAVORS OF ALKALINITY
% ************************************************************************
% ************************************************************************

% BUILD SOME INDICES AGAIN
iP   = find(strcmp('PRES_ADJUSTED',adj_hdr)   == 1); % order always lat, p, t, s
iT   = find(strcmp('TEMP_ADJUSTED',adj_hdr)   == 1);
iS   = find(strcmp('PSAL_ADJUSTED',adj_hdr)   == 1);
iZ   = find(strcmp('DEPTH',adj_hdr)  == 1);
iCHL = find(strcmp('CHLA_ADJUSTED',adj_hdr)    == 1); 
iCHL435  = find(strcmp('CHLA435_ADJUSTED',adj_hdr)    == 1);

iBBP    = find(strcmp('BBP700_ADJUSTED',adj_hdr) == 1);
iBBP532 = find(strcmp('BBP532_ADJUSTED',adj_hdr) == 1);
iCDOM   = find(strcmp('CDOM_ADJUSTED',adj_hdr)   == 1);

iO   = find(strcmp('DOXY_ADJUSTED',adj_hdr)    == 1);
iO2  = find(strcmp('DOXY2_ADJUSTED',adj_hdr)    == 1);
iN   = find(strcmp('NITRATE_ADJUSTED',adj_hdr) == 1);
iST  = find(strcmp('SIGMA_THETA',adj_hdr) == 1);
iL   = find(strncmp('Lat [',adj_hdr,5)    == 1);
iPH  = find(strcmp('PH_IN_SITU_TOTAL_ADJUSTED',adj_hdr)   == 1); % pH
iOCR = find(~cellfun(@isempty,regexp(adj_hdr,'^DOWN.+[ED]$','once')) == 1); %1 x 4, all channels, logical

LIAR_alk_str   = '';
CANYON_alk_str = '';
MLR_alk_str    = '';

% ************************************************************************
% CALCULATE LIAR ALKALINITY FOR THE ADJUSTED DATA SET
% ************************************************************************
if ~isempty(iPH) % empty will be false
    LIAR_pos = [adj_data(:,3), adj_data(:,4), adj_data(:,iZ)];
    ptemp   = theta(adj_data(:,iP),adj_data(:,iT),adj_data(:,iS),0);
    
    % TEST FOR pH+O2 & pH+NO O2 jp 07/21/2022, predim ALK array
    % ONLY GENERATE ALK where pH data exist
    tPH_O    = ~isnan(adj_data(:,iPH)) & ~isnan(adj_data(:,iO)); % ph & o2 exist
    tPH_NO_O = ~isnan(adj_data(:,iPH)) & ... 
        (isnan(adj_data(:,iO)) | adj_data(:,iO+1)>2);  % ph but no valid O2
    Alk        = adj_data(:,iPH) * NaN; % pre dim
    Alk_QF     = ones(size(adj_data,1),1) * 99; % pre dim
    Alk_uncert = Alk;
   
    % 07/21/2022 JP
    % With SOLO floats QC "maxing" is tricker because of QF 8 & lots of
    % fill values. Value added only going to ODV files so convert QC 8's to
    % nearest & largest QC non interp flag first before "maxing"
    % PREP QC FLAGS FOR SOLO TYPE FLOATS TO GET FLAGGING RIGHT
    % NO O2 ALK EST. WILL ONLY USE 1st COLUMN
    QF_tmp = [adj_data(:,iS+1) adj_data(:,iO+1)];
    chk8   = QF_tmp == 8; % only triggered for sprof type files
    if sum(chk8(:)) > 1
        ct     = (1:size(QF_tmp,1))';
        for i = 1:size(QF_tmp,2)
            tmp    = QF_tmp(:,i);
            t8     = tmp == 8;
            tnan   = isnan(tmp); % just in case
            i_next = interp1(ct(~t8 & ~tnan), tmp(~t8 & ~tnan), ct(t8),...
                'next');
            i_prev = interp1(ct(~t8 & ~tnan), tmp(~t8 & ~tnan), ct(t8),...
                'previous');
            QF_tmp(t8,i) = max([i_next, i_prev],[],2);
        end
    end

    % �mol/kg unless MolalityTF is set to 0
    
    % MeasIDVec Parameter Key:                  Equation Key:
    % 1. Salinity                               5.  S, Theta, N, AOU
    % 2. Potential temperature                  6.  S, Theta, N
    % 3. Nitrate                                7.  S, Theta, AOU
    % 4. AOU                                    8.  S, Theta
    % 5. Silicate                               13. S, N, AOU
    % 6. O2                                     14. S, N
    % 7. Temperature                            15. S, AOU
    %                                           16. S
    
    %     MeasIDVec    = [1 2 6]; % PSAL, Pot_TEMP, DOXY_ADJ,
    %     Measurements = [adj_data(:,iS), ptemp, adj_data(:,iO)];
    
    MeasIDVec    = [1 6 7]; % PSAL, DOXY_ADJ, TEMP
    Measurements = [adj_data(:,iS), adj_data(:,iO), adj_data(:,iT)];
    
    if any(tPH_O) % pH & O normal alkalinity estimate
        Equations     = 7; % S, Theta, AOU
        Alk_QF(tPH_O) = max(QF_tmp(tPH_O,:),[],2);

        LIAR_alk_str = ['TALK_LIAR = Alkalinity estimated with the LIAR V2.2 ', ...
            'algorithim Eqn7, MeasIDVec = [1 6 7] (PSAL, DOXY_ADJ, TEMP), ', ...
            '(Carter et al., 2016,doi:10.1002/lom3.10087; Carter et al., 2017, https://doi.org/10.1002/lom3.10232).', ...
            'LIAR code used was downloaded on 10/31/2018 from https://github.com/BRCScienceProducts/LIRs.'];
        LIAR_alk_str2 = '';
        disp('Calculating Alkalinity using LIAR_v2.2 ........')
        
        [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos(tPH_O,:), ...
            Measurements(tPH_O,:), MeasIDVec,'Equations', Equations); % update Jp 07/21/22
        Alk(tPH_O)        = Alk_Est;
        Alk_uncert(tPH_O) = Uncert_Est;
    end
    % IF O2 DATA MISSING OR BAD DO ALK ESTIMATE WITH S AND potT ONLY
    if any(tPH_NO_O) % MISSING O2 DATA, REDO ALK AND QC FLAGS = 4
        MeasIDVec     = [1 7]; % PSAL, TEMP % updated 8/11/17 jp
        Equations     = 8; % S, Theta
        LIAR_alk_str2  = ['SOME O2 DATA IS MISSING! Using only salinity ',...
            'and potential temperature in LIAR estimate for these samples'];
        NO_O_QC             = max(QF_tmp(tPH_NO_O,1),[],2);
        NO_O_QC(NO_O_QC <3) = 3; % All No O2 estimates will have a QC of 3 at best
        Alk_QF(tPH_NO_O)    = NO_O_QC;

        [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos(tPH_NO_O,:), ...
            Measurements(tPH_NO_O,[1,3]), MeasIDVec,'Equations', Equations); % update JP 07/21/22
        Alk(tPH_NO_O)        = Alk_Est;
        Alk_uncert(tPH_NO_O) = Uncert_Est;
    end    
% % % % %     % IF O2 DATA MISSING OR BAD DO ALK ESTIMATE WITH S AND potT ONLY
% % % % %     if any(tPH_NO_O) % MISSING O2 DATA, REDO ALK AND QC FLAGS = 4
% % % % %         MeasIDVec     = [1 7]; % PSAL, TEMP % updated 8/11/17 jp
% % % % %         Equations     = 8; % S, Theta
% % % % %         LIAR_alk_str2  = ['SOME O2 DATA IS MISSING! Using only salinity ',...
% % % % %             'and potential temperature in LIAR estimate for these samples'];
% % % % %         Alk_QF(tPH_NO_O) = max(QF_tmp(tPH_NO_O,1),[],2);
% % % % %         Alk_QF(Alk_QF <3) = 3;
% % % % % 
% % % % %         [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos(tPH_NO_O,:), ...
% % % % %             Measurements(tPH_NO_O,[1,3]), MeasIDVec,'Equations', Equations); % update JP 07/21/22
% % % % %         Alk(tPH_NO_O)        = Alk_Est;
% % % % %         Alk_uncert(tPH_NO_O) = Uncert_Est;
% % % % %     end
    
    t1 = Alk_uncert > 25 & Alk_QF <3; % CAN SET TO BETTER LIMITS AT SOME POINT
    t2 = Alk_uncert > 50 & Alk_QF <4;
    
    Alk_QF(t1) = 3; % ARGO QF's STILL
    Alk_QF(t2) = 4;
    
    t_nan = isnan(Alk);
    Alk_QF(t_nan) = 99; % ARGO QF missing values
    
    % NOW ADD TO HEADER AND DATA VARIABLES
    adj_hdr = [adj_hdr, 'LIAR_ALK' 'LIAR_ALK_QC'];
    adj_data = [adj_data, Alk, Alk_QF];
    
    [raw_r,raw_c] = size(raw_data); % get new size
    [adj_r,adj_c] = size(adj_data); % get new size
    
    clear Uncert_Est MinUncert_Equ Measurements Uncert_Est ptemp
    clear  pos_fix  t1 t2 Alk Alk_uncert Alk_est Alk_QF
end


% ************************************************************************
% NOW GET PH TOTAL @ 25C&0P & DIC ESTIMATES USING CO2SYSv3
% ONLY USE CO2SYS WITH IN RANGE pH VALUES - OTHERWISE INF LOOP IN CALC
% FOR ODV FILES WITH REALLY BAD pH, THE pH WILL BE SEt TO MVI & QF == 1
% THIS IS BECAUSE CO2SYS WILL BOG DOWN ON UNREALISTIC PH VALUES
% ************************************************************************
% NEED A FEW MORE INDICES

iAL  = find(strcmp('LIAR_ALK',adj_hdr)   == 1); % LIAR ESTIMATE
% iAC  = find(strcmp('CANYON_ALK',adj_hdr) == 1); % CANYON estimate
% iAM  = find(strcmp('MLR_ALK',adj_hdr)    == 1); % MLR estimate

if ~isempty(iPH) % pH data exists
    tQC_PH        = adj_data(:,iPH) >= 6.5 & adj_data(:,iPH) <= 9;
    pHSCALEIN     = 1;  % TOTAL SCALE
    K1K2CONSTANTS = 10; % Lueker et al, 2000
    KSO4CONSTANTS = 3;  % KSO4 of Dickson & TB of Lee 2010 (USE THIS ONE !!!)
    %KSO4CONSTANTS = 1;  % KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED)
    
    % IN SITU CONDITIONS - USING pH & TALK TO ESTIMATE DIC
    SAL      = adj_data(:,iS);
    TEMPIN   = adj_data(:,iT);
    PRESIN   = adj_data(:,iP);
    
    % NEED TO ESTIMATE SILICATE AND PHOSPHATE
    % USE REDFILED TO APROXIMATE IF GOOD NITRATE EXISTS, OTHERWISE 0
    %     SI       = ones(size(SAL)) + 0; % set to ) since no data
    %     PO4      = ones(size(SAL)) + 0;
    SI  = ones(size(adj_data(:,iP))) * 0; % start with zeros
    PO4 = SI;
    if ~isempty(iN) % check for nitrate data
        tgood       = adj_data(:,iN+1) == 1;   % good nitrate data = 1, REST =0
        SI          = adj_data(:,iN) * 2.5;  % APROX WITH REDFIELD
        PO4         = adj_data(:,iN) / 16;
        SI(~tgood)  = 0;
        PO4(~tgood) = 0;
    end
    
    TEMPOUT  = TEMPIN;
    PRESOUT  = PRESIN;
    PAR1     = adj_data(:,iPH);
    PAR1(~tQC_PH) = NaN;
    PAR1TYPE = 3; % pH
    PAR2TYPE = 1; % Alkalinity
    
    if ~isempty(iAL) % LIAR ALKALINITY DIC
        PAR2     = adj_data(:,iAL);
        tQC_PAR2 = adj_data(:,iAL+1) ==4; %add check for bad ALK values, bogging CO2SYS down 10/9/2017
        PAR2(tQC_PAR2) = nan;
        
% % % %         [OUT] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
% % % %             TEMPIN, 25, PRESIN, 0, SI, PO4, ...
% % % %             pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);
                [OUT] = CO2SYSv3(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
            TEMPIN, 25, PRESIN, 0, SI, PO4, 0, 0, ...
            pHSCALEIN, K1K2CONSTANTS, 1,2,2);
        
        LIAR_PH25C =  OUT(:,20); % Col 18 in OUT is pH 25C; TM 20 in Sharp v3
        LIAR_PH25C(LIAR_PH25C==-999)=NaN; % CO2SYSv3 sets all NaN to -999...set them back!
        
        % MAKE A BIASED PH DATA ARRAY FOR CALCULATING pCO2 IN ACCORDANCE
        % WITH WILLIAMS et al 2017, SECT 3.4 EQUATION 3
        % CALCULATE BIAS AT CALIBRATION DEPTH (1500m or 1000m) FOR EACH
        % PROFILE AND ADD THIS TO THE PROFILE PH TO MAKE THE BIASED ARRAY
        cycles   = unique(adj_data(:,1));
        BIAS_PH  = ones(size(adj_data(:,iPH)))* NaN;
        pres_tol = 150; % max lookup offset in meters
        ind      = [];
        disp('Calculating biased pH for CO2 calculation using CO2SYS')
        
        for i = 1:size(cycles,1)
            ind = [];
            tcycle  = adj_data(:,1) == cycles(i); % flag profile
            
            % get P, pH, pH25C for profile
            tmp     = [adj_data(tcycle,[iP,iPH,iPH+1]), LIAR_PH25C(tcycle)];
            tbad    = tmp(:,3) == 4 | isnan(tmp(:,2));
            tmp(tbad,:) = [];
            
            if isempty(tmp)
                disp(['Could not calculate biased pH for cycle ', ...
                    num2str(cycles(i)), ' no usable deep pH found']);
                continue
            end
            
            p1500   = abs(tmp(:,1) - 1500);
            min1500 = min(p1500);
            if min1500 < pres_tol % look for 1500m sample first
                ind = find(p1500 == min1500,1);
            else
                p1000   = abs(tmp(:,1) - 1000); % look for 1000m sample
                min1000 = min(p1000);
                if min1000 < pres_tol % look for 1000m sample next
                    ind = find(p1000 == min1000,1);
                end
            end
            
            if ~isempty(ind)
                pH_bias = 0.034529 * tmp(ind,4) - 0.26709;
                BIAS_PH(tcycle) = adj_data(tcycle,iPH) - pH_bias;
            else
                disp(['Could not calculate biased pH for cycle ', ...
                    num2str(cycles(i))]);
            end
        end
        
        % SET OUT OF RANGE pH = NAN SO CO2SYS DOESN't BOG DOWN
        BIAS_PH(~tQC_PH) = NaN;
        
        % REDO FOR pCO2 using biased pH at INSITU P & T
% % % %         [OUT] = CO2SYSSOCCOM(BIAS_PH, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
% % % %             TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
% % % %             pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);

        [OUT] = CO2SYSv3(BIAS_PH, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
            TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, 0, 0, ...
            pHSCALEIN, K1K2CONSTANTS, 1, 2, 2);
        
        LIAR_pCO2 = OUT(:,21); % Col 19 in OUT is pCO2, (uatm); TM, 21 in Sharp V3
        LIAR_pCO2(LIAR_pCO2==-999)=NaN; % CO2SYSv3 sets all NaN to -999...set them back!
        
        
        
        clear cycles pres_tol i tcycle tmp p1500 min1500 ind
        clear  p1000 min1000 BIAS_PH
        
        
        % REDO FOR DIC at INSITU P & T
% % % %         [OUT] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
% % % %             TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
% % % %             pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);

        [OUT] = CO2SYSv3(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
            TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, 0, 0, ...
            pHSCALEIN, K1K2CONSTANTS, 1 ,2 ,2);
        
        LIAR_DIC  = OUT(:,2); % Col 2 in OUT is TCO2, (umol/kgSW); TM same for Sharp v3.
        LIAR_DIC(LIAR_DIC==-999)=NaN; % CO2SYSv3 sets all NaN to -999...set them back!
        %         DIC_QF = max([adj_data(:,iS+1) adj_data(:,iO+1) ...
        %                       adj_data(:,iAL+1) adj_data(:,iPH+1)],[],2);
%         DIC_QF = max([adj_data(:,iS+1) ... %remove Oxygen flag as part of max calc.  It's incorporated into the ALK flag (above) already, and including here again causes problems if value is missing (would trump a 'quest' value cuz missing in argo is 99 (larger than a 3)).
%             adj_data(:,iAL+1) adj_data(:,iPH+1)],[],2);
        
        % 07/21/2022 JP
        % With SOLO floats QC "maxing" is tricker because of QF 8 & lots of
        % fill values, & ALK goes full resolution but pH does not. Value
        % added only going to ODV files so convert QC 8's to nearest &
        % largest QC non interp flag first before "maxing"
        
        QF_tmp = [adj_data(:,iS+1), adj_data(:,iAL+1), adj_data(:,iPH+1)];
        chk8   = QF_tmp == 8; % only triggered for sprof type files
        if sum(chk8(:)) > 1
            ct     = (1:size(QF_tmp,1))';
            for i = 1:size(QF_tmp,2)
                tmp    = QF_tmp(:,i);
                t8     = tmp == 8;
                tnan   = isnan(tmp); % just in case
                i_next = interp1(ct(~t8 & ~tnan), tmp(~t8 & ~tnan), ct(t8),...
                    'next');
                i_prev = interp1(ct(~t8 & ~tnan), tmp(~t8 & ~tnan), ct(t8),...
                    'previous');
                QF_tmp(t8,i) = max([i_next, i_prev],[],2);
            end
        end
        DIC_QF = max(QF_tmp,[],2);
        clear QF_tmp chk8 ct i tmp t8 tnan i_next i_prev

       
        adj_hdr  = [adj_hdr, 'LIAR PHTOT25C' 'LIAR PHTOT25C_QC', ...
            'LIAR_DIC' 'LIAR_DIC_QC' 'LIAR_pCO2' 'LIAR_pCO2_QC'];
        adj_data = [adj_data, LIAR_PH25C, DIC_QF, LIAR_DIC, DIC_QF, ...
            LIAR_pCO2, DIC_QF];
        [adj_r,adj_c] = size(adj_data); % get new size
        
        % DO A LAST CHECK FOR MISSING VALUES
        iPH25   = find(strcmp('LIAR PHTOT25C',adj_hdr)   == 1);
        iDIC    = find(strcmp('LIAR_DIC',adj_hdr)   == 1);
        iPCO2   = find(strcmp('LIAR_pCO2',adj_hdr)   == 1);
        
        adj_data(isnan(adj_data(:,iPH25)), iPH25+1) = 99; % MVI QF VALUE
        adj_data(isnan(adj_data(:,iDIC)), iDIC+1)   = 99;
        adj_data(isnan(adj_data(:,iPCO2)), iPCO2+1) = 99;
        

    end
end

clear iP iT iS %iO iN


% ************************************************************************
% CONVERT ARGO DATA QUALITY FLAGS TO ODV DATA QUALITY FLAGS
% DO THIS BEFORE CHECKING FOR BOSS DATA AT U MAINE
% U MAINE QUALITY FLAGS ALREADY ODV STYLE
% ************************************************************************

for i = 1:raw_c
    % SETTING ALL QF's TO 1, EXCEPT FOR PTSZ & OBVIOUSLY BAD
    if regexp(raw_hdr{i},'\_QC', 'once')
        tmp1 = raw_data(:,i);
        
        % solo sprof merge can have 8's, replace with highest bounding QC
        t8 = tmp1 == 8; 
        if sum(t8,1) > 0 % sprof can have interpolated QC values  (8)
            ct     = (1:size(tmp1,1))';
            tnan   = isnan(tmp1);
            i_next = interp1(ct(~t8 & ~tnan), tmp1(~t8 & ~tnan), ct(t8),...
                'next');
            i_prev = interp1(ct(~t8 & ~tnan), tmp1(~t8 & ~tnan), ct(t8),...
                'previous');
            tmp1(t8) = max([i_next, i_prev],[],2);
            clear i_next i_prev t8 ct tnan
        end

        tmp1(tmp1 == 5) = 3; % This is for NPQ corrected data
        tmp1(tmp1 == 4) = 8;
        tmp1(tmp1 == 3) = 4; % This will catch interp lat 2/9/17 jp
        tmp1(tmp1 == 0) = 10; % temporary
        
        if strncmp(raw_hdr{i},'BBP',3)
            tmp1(tmp1 == 1) = 0;
            tmp1(tmp1 == 2) = 1; % BGC Argo flags raw bbp qc =2
        else
            tmp1(tmp1 == 1 |tmp1 == 2) = 0; % P T S Z 2/17/22 jp and ==2 piece here
        end
        tmp1(tmp1 == 10 | tmp1 == 99) = 1; % NO QC or Missing value
        raw_data(:,i) = tmp1;
    end
end


for i = 1:adj_c
    if regexp(adj_hdr{i},'\_QC', 'once')
        tmp2 = adj_data(:,i);
        
        % solo sprof merge can have 8's, replace with highest bounding QC
        t8 = tmp2 == 8;
        if sum(t8,1) > 0 % sprof has interpolated values
            ct     = (1:size(tmp2,1))';
            tnan   = isnan(tmp2);
            i_next = interp1(ct(~t8 & ~tnan), tmp2(~t8 & ~tnan), ct(t8),...
                'next');
            i_prev = interp1(ct(~t8 & ~tnan), tmp2(~t8 & ~tnan), ct(t8),...
                'previous');
            tmp2(t8) = max([i_next, i_prev],[],2);
            clear i_next i_prev t8 ct tnan
        end
        
        tmp2(tmp2 == 5) = 3; % This is for NPQ corrected data
        tmp2(tmp2 == 4) = 8;
        tmp2(tmp2 == 3) = 4;
        tmp2(tmp2 == 0) = 10; % temporary
        tmp2(tmp2 == 2 | tmp2 == 1) = 0;
        tmp2(tmp2 == 10 | tmp2 == 99) = 1; % NO QC or Missing value
        adj_data(:,i) = tmp2;
    end
end



if strcmp(info.float_type, 'APEX') & ~tf_APEX_OCR
    for i = 1:hrdc
        if regexp(raw_hdr{i},'\_QC', 'once')
            tmp2 = hrraw_data(:,i);
            tmp2(tmp2 == 5) = 3; % This is for NPQ corrected data
            tmp2(tmp2 == 4) = 8;
            tmp2(tmp2 == 3) = 4;
            tmp2(tmp2 == 0) = 10; % temporary
            tmp2(tmp2 == 2 | tmp2 == 1) = 0;
            tmp2(tmp2 == 10 | tmp2 == 99) = 1; % NO QC or Missing value
            hrraw_data(:,i) = tmp2;
        end
    end
end

clear i tmp1 tmp2


% ************************************************************************
% ************************************************************************
% NOW DEAL WITH Bio-optical data - IF NO ADJUSTED DATA EXISTS
% FILL WITH RAW FOR TXT FILES
% iChl    = find(strcmp('CHLA_ADJUSTED',adj_hdr)   == 1); %Are these needed if laready defined? jp 02/04/24
% iBbp    = find(strcmp('BBP700_ADJUSTED',adj_hdr) == 1);
% iCdm    = find(strcmp('CDOM_ADJUSTED',adj_hdr)   == 1);
% iBbp532 = find(strcmp('BBP532_ADJUSTED',adj_hdr) == 1);
% OCR_cols = ~cellfun(@isempty, regexp(adj_hdr,'E\d{3}(?!.+QC)|PAR(?!.+QC)','once')); % This will get IRR & PAR & NOT QC
% iOCR = find(OCR_cols == 1);

info.Chl_data_str    = '';
info.Chl435_data_str = '';
info.Bbp_data_str    = '';
info.Cdm_data_str    = '';
info.Bbp532_data_str = '';
info.OCR_data_str    = '';

% These work but the param index really should be checked for non empty first. Early JP
% work!! JP 02/04/24
if ~isempty(iCHL) && sum(~isnan(adj_data(:,iCHL)),1) == 0 % No adjusted data
    adj_data(:,iCHL) = raw_data(:,iCHL);
    adj_data(:,iCHL+1) = raw_data(:,iCHL+1);
    adj_data(adj_data(:,iCHL+1)~= 8,iCHL+1) = 1; % uninspected
    info.Chl_data_str = 'Adjusted Chlorophyll data = raw Chlorophyll data!';
    disp(info.Chl_data_str)
end

if ~isempty(iCHL435) && sum(~isnan(adj_data(:,iCHL435)),1) == 0 % No adjusted data
    adj_data(:,iCHL435) = raw_data(:,iCHL435);
    adj_data(:,iCHL435+1) = raw_data(:,iCHL435+1);
    adj_data(adj_data(:,iCHL435L+1)~= 8, iCHL435+1) = 1; % uninspected
    info.Chl435_data_str = 'Adjusted Chlorophyll4535 data = raw Chlorophyll435 data!';
    disp(info.Chl435_data_str)
end

if ~isempty(iBBP) && sum(~isnan(adj_data(:,iBBP)),1) == 0 % No adjusted data
    adj_data(:,iBBP)   = raw_data(:,iBBP);
    adj_data(:,iBBP+1) = raw_data(:,iBBP+1);
    adj_data(adj_data(:,iBBP+1)~= 8,iBBP+1) = 1; % uninspected
    info.Bbp_data_str = 'Adjusted backscatter data (700) = raw backscatter data (700)!';
    disp(info.Bbp_data_str)
end

if ~isempty(iCDOM) && sum(~isnan(adj_data(:,iCDOM)),1) == 0 % No adjusted data
    adj_data(:,iCDOM) = raw_data(:,iCDOM);
    adj_data(:,iCDOM+1) = raw_data(:,iCDOM+1);
    adj_data(adj_data(:,iCDOM+1)~= 8,iCDOM+1) = 1; % uninspected
    info.Cdm_data_str = 'Adjusted FDOM data = raw FDOM data!';
    disp(info.Cdm_data_str)
end

if ~isempty(iBBP532) && sum(~isnan(adj_data(:,iBBP532)),1) == 0 % No adjusted data
    adj_data(:,iBBP532)   = raw_data(:,iBBP532);
    adj_data(:,iBBP532+1) = raw_data(:,iBBP532+1);
    adj_data(adj_data(:,iBBP532+1)~= 8,iBBP532+1) = 1; % uninspected
    info.Bbp532_data_str = 'Adjusted backscatter (532) data = raw backscatter data (532)!';
    disp(info.Bbp532_data_str)
end

if ~isempty(iOCR) && all(isnan(adj_data(:,iOCR)),'all') % No adjusted data
    adj_data(:,iOCR)   = raw_data(:,iOCR);
    t8 = raw_data(:,iOCR+1) == 8;
    t8(t8 ==1) = 8; % set flags
    t8(t8 == 0) = 1; % set flags, good data are uninspected = 1
    adj_data(:,iOCR+1) = t8;
    info.OCR_data_str = 'Adjusted OCR radiometer values = raw OCR radiometer values!';
    disp(info.OCR_data_str)
end


% ************************************************************************
% % ADD "VALUE ADDED" BIO-OPTICS PRODUCTS

if ~isempty(iBBP)
    % using Johnson et al.,2017, doi:10.1002/2017JC012838
    POC = (3.12e4 * adj_data(:,iBBP) + 3.04) / 12; % mmol/m^3
    adj_hdr    = [adj_hdr, 'POC[mmol/m^3]' 'QF']; % ADD TO HEADER
    adj_data = [adj_data, POC, adj_data(:,iBBP+1)];
end

[raw_r,raw_c] = size(raw_data); % get new size
[adj_r,adj_c] = size(adj_data); % get new size

clear NPQ iNPQ tNPQ tLAT CHL_EXP CHL_new i t1 tmp cycle_ct

% ************************************************************************
% CHECK FOR MISSING VALUES - SHOULD ALL BE SET TO NaN's NOW
% ************************************************************************
raw_missing_data_str = 'No missing data found';
nan_sum =(sum(isnan(raw_data), 2))>0; % cols so I can get casts later

if strcmp(FLT_type,'APEX') && ~tf_APEX_OCR && sum(nan_sum,1) > 0
    tmp = unique(raw_data(nan_sum,1));
    raw_missing_data_str = ['Missing Float data detected for raw data', ...
        ' station(s): ', sprintf('%0.0f ',tmp)];
    disp(raw_missing_data_str);
end
clear nan_sum tmp

adj_missing_data_str = 'No missing data found';
nan_sum = sum(isnan(adj_data), 2)>0; % cols so I can get casts later
if strcmp(FLT_type,'APEX') && ~tf_APEX_OCR && sum(nan_sum,1) > 0
    tmp = unique(adj_data(nan_sum,1));
    adj_missing_data_str = ['Missing Float data detected for adjusted ',...
        'station(s): ', sprintf('%0.0f ',tmp)];
    disp(adj_missing_data_str);
end
clear nan_sum tmp

% ************************************************************************
%          BUILD COMBINED LR + HR DATA SET FOR "HR" FLOATVIZ FILE
%                 !!! APEX FLOATS ONLY, HR DATA = P T S !!!
% ************************************************************************
if strcmp(info.float_type, 'APEX') && ~tf_APEX_OCR
    cycles = unique(raw_data(:,1));
    allraw_data = ones(hrdr + raw_r, raw_c) * NaN;
    alladj_data = ones(hrdr + adj_r, adj_c) * NaN;
    allraw_ct = 1;
    alladj_ct = 1;
    iP     = find(strcmp('PRES',raw_hdr) == 1);
    
    %COMBINE HR & LR AND SORT
    for i = 1 : size(cycles,1)
        t_hr   = hrraw_data(:,1) == cycles(i);
        hr_ct  = sum(t_hr);
        t_raw  = raw_data(:,1)   == cycles(i);
        raw_ct = sum(t_raw);
        t_adj  = adj_data(:,1)   == cycles(i);
        adj_ct = sum(t_adj);
        
        tmp_raw = ones(raw_ct + hr_ct, raw_c)* NaN;
        tmp_raw(1:raw_ct, 1:raw_c) = raw_data(t_raw,:);
        tmp_raw(raw_ct+1: raw_ct+hr_ct, 1:hrdc) = hrraw_data(t_hr,:);
        [~,IX]  = sort(tmp_raw(:,iP),'descend');
        tmp_raw = tmp_raw(IX,:);
        allraw_data(allraw_ct:allraw_ct + raw_ct + hr_ct-1,:) = tmp_raw;
        allraw_ct = allraw_ct + raw_ct + hr_ct;
        
        tmp_adj = ones(adj_ct + hr_ct, adj_c)* NaN;
        tmp_adj(1:adj_ct, 1:adj_c) = adj_data(t_adj,:);
        tmp_adj(adj_ct+1: adj_ct+hr_ct, 1:hrdc) = hrraw_data(t_hr,:);
        [~,IX]  = sort(tmp_adj(:,iP),'descend');
        tmp_adj = tmp_adj(IX,:);
        alladj_data(alladj_ct:alladj_ct + adj_ct + hr_ct-1,:) = tmp_adj;
        alladj_ct = alladj_ct + adj_ct + hr_ct;
        
        clear t_hr hr_ct t_raw raw_ct t_adj adj_ct tmp_raw tmp_adj
    end
    
    % NOW SET ALL QF NaN's to 1
    for i = 1:raw_c % SET QF NaN's to 1
        if regexp(raw_hdr{i},'\_QC', 'once')
            tmp2 = allraw_data(:,i);
            tmp2(isnan(tmp2)) = 1; % NO QC or Missing value
            allraw_data(:,i) = tmp2;
        end
    end
    for i = 1:adj_c % SET QF NaN's to 1
        if regexp(adj_hdr{i},'\_QC|^QF', 'once')
            tmp2 = alladj_data(:,i);
            tmp2(isnan(tmp2)) = 1; % NO QC or Missing value
            alladj_data(:,i) = tmp2;
        end
    end
end


% ************************************************************************
% ************************************************************************
% PRINT DATA TO FILE
% http://blogs.mathworks.com/loren/2006/04/19/high-performance-file-io/
% Change 'w' to 'W' for the file feopens to save time
% ************************************************************************
% ************************************************************************
load([dirs.cal,'cal',MBARI_ID_str,'.mat']);
notes_flag = 0;
if isfield(cal.info,'notes')
    notes_flag = 1;
    notes = cal.info.notes;
    %     clear cal %TM 3/3/21; DON'T CLEAR YET; IT'S USED JUST A BIT FURTHER
    %     DOWN.
end
% keyboard
MVI_str = '-1e10'; % MISSING VALUE INDICATOR FOR ODV

Bio_optics_str = ['See: Boss, E.B. and N. Ha�ntjens, 2016. ', ...
    'http://soccom.princeton.edu/sites/default/files/files/', ...
    'SOCCOM_2016-1_Bio-optics-primer.pdf'];

CO2SYS_str = ['//NOTE ON CO2SYS CARBONATE SYSTEM CALCULATIONS:\r\n',...
    '//All carbonate system variables calculated with CO2SYS for Matlab\r\n',...
    '//(Sharp et al., 2020; see also Lewis and Wallace 1998)\r\n',...
    '//used the following conditions: pH was reported on the total scale.\r\n',...
    '//K1 and K2 dissociation constants were from Lueker et al., 2000, doi:\r\n', ...
    '//10.1016/S0304-4203(00)00022-0. The KSO4 dissociation constant was\r\n',...
    '//from Dickson, 1990, doi: 10.1016/0021-9614(90)90074-Z. The KF dissociation\r\n',...
    '//constant was from Perez and Fraga 1987, doi: 10.1016/0304-4203(87)90036-3.\r\n'....
    '//The borate to salinity ratio was from Lee et al., 2010,\r\n', ...
    '//doi:10.1016/j.gca.2009.12.027. Silicate and Phosphate were not\r\n', ...
    '//measured by the float, but estimates based on Redfieldian ratios improved\r\n',...
    '//the carbonate system estimates. If a nitrate value was considered to\r\n',...
    '//be of good quality silicate = nitrate*2.5 and phosphate = nitrate/16,\r\n',...
    '//otherwise the best estimate for both was considered to be 0. When pCO2\r\n',...
    '//was estimated from TALK_LIAR and pHinsitu, a bias was first added to pHinsitu\r\n',...
    '//following Williams et al., 2017, doi: https://doi.org/10.1002/2016GB005541 , section 3.4, equation 3.\r\n',...
	'//This correction is not necessary for DIC and DIC is computed with the reported pH and the TALK_LIAR value.\r\n'];

CHL_str = ['//NOTE ON Chl_a [mg/m^3] CONCENTRATION:\r\n',...
    '//There is community-established calibration bias of 2 for the WET Labs 413',...
    ' ECO-series fluorometers\r\n//(Roesler et al, 2017, doi: 10.1002/lom3.10185).',...
    '\r\n//Chl_a has been recalculated using in situ measured dark counts. ', ...
    'Chl_a is then divided by the Roesler factor of 2.\r\n//Lastly, profiles ', ...
    'with sun elevaltion > 0 are corrected for NPQ (Xing et al, 2012, doi: ', ...
    '10.4319/lom.2012.10.483)\r\n//with uncorrected spikes (raw-filtered data)', ...
    ' added back to the corrected profile.\r\n'];
    
%     //Chl_a_corr uses a modified Xing ',...
%     'approach to correct for NPQ: the reference depth is the\r\n//',...
%     'shallower of the mixed layer depth or the 1 percent light depth (Kd based on ', ...
%     '\r\n//Kim et al., 2015, doi:10.5194/bg-12-5119-2015).No spike profile added.',...
%     '\r\n//Note that all NPQ-corrected data receives an Argo quality flag of 5 ("value changed")',...
%     '\r\n//and an ODV-style quality flag of 4 ("questionable").',...
%     '\r\n//South of 30S a slope correction of 6 was used(E. Boss unpublished data).',...
%     '\r\n//This correction scheme was decided upon at the 18th Argo Data Management Team ',...
%     '\r\n//meeting in Hamburg, Germany (Nov, 2017), and is subject to change as research on optimal ',...
%     '\r\n//correction methods for float data from FLBB sensors continues.\r\n//'];

SOLO_info_str = ['//\r\n//THIS IS a BGC SOLO FLOAT!\r\n',...
    '//Each BGC sensor may have a unique pressure axis. If so all data have been aligned\r\n',...
    '//on to a common pressure axis using a modified version of code from Henry Bittig\r\n',...
    '//which was originally used to build the Sprof profile files. More information at:\r\n',...
    '//https://archimer.ifremer.fr/doc/00445/55637/  and  ',...
    'https://www.seanoe.org/data/00345/45589/\r\n',...
    '//(\\decArgo_YYYYMMDD_XXXX\\decArgo_soft\\soft\\util\\sub_foreign', ...
    '\\ARGO_simplified_profile.m\r\n', ...
    '//Data here may not exactly match that of the Sprof profile file ',...
    'due to slight QC flagging differences\r\n//\r\n'];

TEST_OCR_info_str = ['//\r\n//THIS IS AN APEX TEST FLOAT\r\n//It has one of the ', ...
    'first installations of an OCR 504 Radiometer on an APEX float\r\n//The OCR data have ',...
    'their own pressure axis so data have been aligned\r\n',...
    '//on to a common pressure axis using a modified version of code from Henry Bittig\r\n',...
    '//which was originally used to build the Sprof profile files. Prior to alignment \r\n',...
    '//non unique data shallower than 4 dbar were replaced by the median of these data points.\r\n', ...
    '//More information at:\r\n',...
    '//https://archimer.ifremer.fr/doc/00445/55637/  and  ',...
    'https://www.seanoe.org/data/00345/45589/\r\n',...
    '//(\\decArgo_YYYYMMDD_XXXX\\decArgo_soft\\soft\\util\\sub_foreign', ...
    '\\ARGO_simplified_profile.m\r\n', ...
    '//Data here may not exactly match that of the Sprof profile file ',...
    'due to slight QC flagging differences\r\n//\r\n'];

% LOAD CAL FILE - will use info in meta data headers
% % load([dirs.cal,'cal',info.FloatViz_name,'.mat']);  % TM 3/3/21; THE CAL
% FILE IS LOADED ON LINE ~1020, DOES THIS NEED TO BE LOADED AGAIN??
% ************************************************************************
% BUILD LOOK UP CELL ARRAY to match variables and set format string
% ************************************************************************
%RAW ODV FILE
ODV_raw(1,:)  = {'Pressure[dbar]'        '%0.2f' 'PRES' '' '' ''}; % ?
ODV_raw(2,:)  = {'Temperature[�C]'       '%0.4f' 'TEMP' '' '' ''};
ODV_raw(3,:)  = {'Salinity[pss]'         '%0.4f' 'PSAL' '' '' ''};
ODV_raw(4,:)  = {'Sigma_theta[kg/m^3]'   '%0.3f' 'SIGMA_THETA' '' '' ''};
ODV_raw(5,:)  = {'Depth[m]'              '%0.3f' 'DEPTH' '' '' ''};
ODV_raw(6,:)  = {'Oxygen[�mol/kg]'       '%0.2f' 'DOXY' '' '' ''};
ODV_raw(7,:)  = {'OxygenSat[%]'          '%0.1f' 'DOXY_%SAT' '' '' ''};
ODV_raw(8,:)  = {'Nitrate[�mol/kg]'      '%0.2f' 'NITRATE' '' '' ''};
ODV_raw(9,:)  = {'Chl_a[mg/m^3]'         '%0.4f' 'CHLA' '' '' ''};
ODV_raw(10,:)  = {'Chl_a435[mg/m^3]'         '%0.4f' 'CHLA435' '' '' ''};
ODV_raw(11,:) = {'b_bp700[1/m]'          '%0.6f' 'BBP700' '' '' ''};
ODV_raw(12,:) = {'pHinsitu[Total]'       '%0.4f' 'PH_IN_SITU_TOTAL' '' '' ''};
ODV_raw(13,:) = {'b_bp532[1/m]'          '%0.6f' 'BBP532' '' '' ''};
ODV_raw(14,:) = {'CDOM[ppb]'             '%0.2f' 'CDOM' '' '' ''};

% OCR VARIABLES
ODV_raw(15,:) = {'DOWN_IRRAD380[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE380' '' '' ''};
ODV_raw(16,:) = {'DOWN_IRRAD412[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE412' '' '' ''};
ODV_raw(17,:) = {'DOWN_IRRAD443[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE443' '' '' ''};
ODV_raw(18,:) = {'DOWN_IRRAD490[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE490' '' '' ''};
ODV_raw(19,:) = {'DOWN_IRRAD555[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE555' '' '' ''};
ODV_raw(20,:) = {'DOWNWELL_PAR[�mol Quanta/m^2/sec]'  '%0.6f' 'DOWNWELLING_PAR' '' '' ''};

% FOR DUAL 4330 O2 APEX
ODV_raw(21,:)  = {'Oxygen2[�mol/kg]'       '%0.2f' 'DOXY2' '' '' ''};
ODV_raw(22,:)  = {'Oxygen2Sat[%]'          '%0.1f' 'DOXY2_%SAT' '' '' ''};

% ************************************************************************
% ADD VARIABLE DESCRIPTORS TO RAW CELL LOOKUP TABLE - SENSOR TYPE SN COMMENT
if strcmp(info.float_type,'APEX') && ~tf_APEX_OCR
    O2sensor = 'Aanderaa ';
else
    O2sensor = '';
end

ODV_raw(1,4:6) = {info.CTDtype info.CTDsn ''}; %P
ODV_raw(2,4:6) = {info.CTDtype info.CTDsn ''}; %T
ODV_raw(3,4:6) = {info.CTDtype info.CTDsn ''}; %S
ODV_raw(4,4:6) = {'' '' 'Potential density at the sea surface'}; %Sigma-T
ODV_raw(5,4:6) = {'' '' 'Depth calculated from pressure and latitude'}; %Z
if sum(strcmp(ODV_raw{6,3},raw_hdr))
    ODV_raw(6,4:6) = {[O2sensor,cal.O.type] cal.O.SN ''}; %O2
    ODV_raw(7,4:6) = {'' '' ['Calculation assumes atmospheric pressure',...
        '= 1013.25 mbar']}; % O2 % sat
end
if sum(strcmp(ODV_raw{8,3},raw_hdr))
    if isfield(cal,'N')
        ODV_raw(8,4:6) = {cal.N.type cal.N.SN ''}; %Nitrate
    else
        ODV_raw(8,4:6) = {'!!!NO NITRATE CALIBRATION DATA !!!' '' ''}; %Nitrate
    end
end
if sum(strcmp(ODV_raw{9,3},raw_hdr))
    ODV_raw(9,4:6) = {cal.CHL.type cal.CHL.SN ''}; %CHL
    ODV_raw(10,4:6) = {cal.CHL.type cal.CHL.SN ''}; %CHL
end
if sum(strcmp(ODV_raw{11,3},raw_hdr))
    ODV_raw(11,4:6) = {cal.BB.type cal.BB.SN ''}; %Back scatter 700
end
if sum(strcmp(ODV_raw{12,3},raw_hdr))
    ODV_raw(12,4:6) = {cal.pH.type cal.pH.SN ''}; %pH
end
if sum(strcmp(ODV_raw{13,3},raw_hdr))
    ODV_raw(13,4:6) = {cal.CDOM.type cal.CDOM.SN ''}; %bbp532
end
if sum(strcmp(ODV_raw{14,3},raw_hdr))
    ODV_raw(14,4:6) = {cal.CDOM.type cal.CDOM.SN ''}; %CDOM
end
if sum(strcmp(ODV_raw{15,3},raw_hdr))
    ODV_raw(15:20,4) = {cal.OCR.type};
    ODV_raw(15:20,5) = {cal.OCR.SN};
    ODV_raw(15:20,6) = {''};

    % ODV_raw(15,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_raw(16,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_raw(17,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_raw(18,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_raw(19,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_raw(20,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
end
if sum(strcmp(ODV_raw{21,3},raw_hdr))
    ODV_raw(21,4:6) = {[O2sensor,cal.O2.type] cal.O2.SN ''}; %O2
    ODV_raw(22,4:6) = {'' '' ['Calculation assumes atmospheric pressure',...
        '= 1013.25 mbar']}; % O2 % sat
end

% ************************************************************************
% JP 02/04/24 I think order could be improved so raw & adj match up to
% value added parameters?
% ************************************************************************
% QC ODV FILE
ODV_adj(1,:)  = {'Pressure[dbar]'        '%0.2f' 'PRES_ADJUSTED' '' '' ''};
ODV_adj(2,:)  = {'Temperature[�C]'       '%0.4f' 'TEMP_ADJUSTED' '' '' ''};
ODV_adj(3,:)  = {'Salinity[pss]'         '%0.4f' 'PSAL_ADJUSTED' '' '' ''};
ODV_adj(4,:)  = {'Sigma_theta[kg/m^3]'   '%0.3f' 'SIGMA_THETA' '' '' ''};
ODV_adj(5,:)  = {'Depth[m]'              '%0.3f' 'DEPTH' '' '' ''};
ODV_adj(6,:)  = {'Oxygen[�mol/kg]'       '%0.2f' 'DOXY_ADJUSTED' '' '' ''};
ODV_adj(7,:)  = {'OxygenSat[%]'          '%0.1f' 'DOXY_%SAT_ADJUSTED' '' '' ''};
ODV_adj(8,:)  = {'Nitrate[�mol/kg]'      '%0.2f' 'NITRATE_ADJUSTED' '' '' ''};
ODV_adj(9,:)  = {'Chl_a[mg/m^3]'         '%0.4f' 'CHLA_ADJUSTED' '' '' ''};
ODV_adj(10,:) = {'Chl_a435[mg/m^3]'         '%0.4f' 'CHLA435_ADJUSTED' '' '' ''};
ODV_adj(11,:) = {'b_bp700[1/m]'          '%0.6f' 'BBP700_ADJUSTED' '' '' ''};
ODV_adj(12,:) = {'POC[mmol/m^3]'         '%0.2f' 'POC[mmol/m^3]' '' '' ''};
ODV_adj(13,:) = {'pHinsitu[Total]'       '%0.4f' 'PH_IN_SITU_TOTAL_ADJUSTED' '' '' ''};
ODV_adj(14,:) = {'pH25C[Total]'          '%0.4f' 'LIAR PHTOT25C' '' '' ''};
ODV_adj(15,:) = {'TALK_LIAR[�mol/kg]'    '%4.0f' 'LIAR_ALK' '' '' ''};
ODV_adj(16,:) = {'DIC_LIAR[�mol/kg]'     '%4.0f' 'LIAR_DIC' '' '' ''};
ODV_adj(17,:) = {'pCO2_LIAR[�atm]'       '%4.1f' 'LIAR_pCO2' '' '' ''};
ODV_adj(18,:) = {'b_bp532[1/m]'          '%0.6f' 'BBP532_ADJUSTED' '' '' ''};
ODV_adj(19,:) = {'CDOM[ppb]'             '%0.2f' 'CDOM_ADJUSTED' '' '' ''};

% OCR VARIABLES - JP 02/04/24 added adjusted to merged param header names
% (col 3)
ODV_adj(20,:) = {'DOWN_IRRAD380[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE380_ADJUSTED' '' '' ''};
ODV_adj(21,:) = {'DOWN_IRRAD412[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE412_ADJUSTED' '' '' ''};
ODV_adj(22,:) = {'DOWN_IRRAD443[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE443_ADJUSTED' '' '' ''};
ODV_adj(23,:) = {'DOWN_IRRAD490[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE490_ADJUSTED' '' '' ''};
ODV_adj(24,:) = {'DOWN_IRRAD555[W/m^2/nm]'            '%0.6f' 'DOWN_IRRADIANCE555_ADJUSTED' '' '' ''};
ODV_adj(25,:) = {'DOWNWELL_PAR[�mol Quanta/m^2/sec]'  '%0.6f' 'DOWNWELLING_PAR_ADJUSTED' '' '' ''};

% FOR DUAL 4330 O2 APEX
ODV_adj(26,:)  = {'Oxygen2[�mol/kg]'       '%0.2f' 'DOXY2' '' '' ''};
ODV_adj(27,:)  = {'Oxygen2Sat[%]'          '%0.1f' 'DOXY2_%SAT' '' '' ''};

% ************************************************************************
% ADD VARIABLE DESCRIPTORS TO ADJ CELL LOOKUP TABLE - SENSOR TYPE SN COMMENT

ODV_adj(1:9,4:6) = ODV_raw(1:9,4:6);
if sum(strcmp(ODV_adj{9,3},adj_hdr))
    ODV_adj(9,4:6) = {cal.CHL.type cal.CHL.SN Bio_optics_str}; %chl
    ODV_adj(10,4:6) = {cal.CHL.type cal.CHL.SN Bio_optics_str}; %chl
end
if sum(strcmp(ODV_adj{11,3},adj_hdr))
    ODV_adj(11,4:6) = {cal.BB.type cal.BB.SN Bio_optics_str}; %bbp700
end

if sum(strcmp(ODV_adj{12,3},adj_hdr))
    ODV_adj(12,4:6) = {cal.BB.type cal.BB.SN Bio_optics_str}; %POC
end

if sum(strcmp(ODV_adj{13,3},adj_hdr))
    ODV_adj(13,4:6) = {cal.pH.type cal.pH.SN ''}; %pH in situ
end
if sum(strcmp(ODV_adj{14,3},adj_hdr))
    str =['estimated with CO2SYS(TALK_LIAR,pHinsitu) for Matlab ', ...
        ' see note below'];
    ODV_adj(14,4:6) = {cal.pH.type cal.pH.SN str}; %pH 25C total
end
if sum(strcmp(ODV_adj{15,3},adj_hdr))
    ODV_adj(15,4:6) = {'' '' LIAR_alk_str}; %TALK LIAR
end
if sum(strcmp(ODV_adj{16,3},adj_hdr))
    str =['estimated with CO2SYS(TALK_LIAR,pHinsitu) for Matlab ', ...
        '(see note below)'];
    ODV_adj(16,4:6) = {'' '' str}; %LIAR DIC
end
if sum(strcmp(ODV_adj{17,3},adj_hdr))
    str =['estimated with CO2SYS(TALK_LIAR,biased pHinsitu) for Matlab ', ...
        '(see note below)'];
    ODV_adj(17,4:6) = {'' '' str}; %LIAR pCO2 DIC
end

if sum(strcmp(ODV_adj{18,3},adj_hdr))
    ODV_adj(18,4:6) = {cal.CDOM.type cal.CDOM.SN ''}; % BBP532
end
if sum(strcmp(ODV_adj{19,3},adj_hdr))
    ODV_adj(19,4:6) = {cal.CDOM.type cal.CDOM.SN ''}; %CDOM
end
if sum(strcmp(ODV_adj{20,3},adj_hdr))
    ODV_adj(20:25,4) = {cal.OCR.type}; %OCR
    ODV_adj(20:25,5) = {cal.OCR.SN}; %OCR
    ODV_adj(20:25,6) = {''}; %OCR


    % ODV_adj(20,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_adj(21,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_adj(22,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_adj(23,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_adj(24,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
    % ODV_adj(25,4:6) = {cal.OCR.type cal.OCR.SN ''}; %OCR
end

if sum(strcmp(ODV_adj{26,3},adj_hdr))
    ODV_adj(26,4:6) = {[O2sensor,cal.O2.type] cal.O2.SN ''}; %O2
    ODV_adj(27,4:6) = {'' '' ['Calculation assumes atmospheric pressure',...
        '= 1013.25 mbar']}; % O2 % sat
end

% ************************************************************************
% ODV_RAW & ADV_ADJ DEFINE THE ORDER OF PARAMETERS PRINTED, SUBSET THIS
% MASTER LIST TO ONLY THE BGC PARAMETERS IN THE DATA FILE, NO LONGER A
% STATIC ORDER WITH FILL VALUES IF PARAMETERS ARE MISSING
%REMOVE SPECIFC VARIABLES FROM THE LOOKUP CELL ARRAY IF NOT PRESENT IN
% DATA

Lia_raw = ismember(ODV_raw(:,3),raw_hdr); %find data params existing in ODV_raw
ODV_raw = ODV_raw(Lia_raw,:);

Lia_adj = ismember(ODV_adj(:,3),adj_hdr); %find data params existing in ODV_adj
ODV_adj = ODV_adj(Lia_adj,:);

% 
% 
% if isempty(iPH)
%     t1 = strcmp('pHinsitu[Total]',ODV_raw(:,1));
%     ODV_raw(t1,:) = [];
% 
%     t1 = strcmp('pHinsitu[Total]',ODV_adj(:,1));
%     t2 = strcmp('pCO2_LIAR[�atm]',ODV_adj(:,1));
%     ODV_adj(t1|t2,:) = [];
% end
% 
% if isempty(iCDOM)
%     t1 = strncmp('CDOM[ppb]',ODV_raw(:,1),4);
%     ODV_raw(t1,:) = [];
% 
%     t1 = strncmp('CDOM[ppb]',ODV_adj(:,1),4);
%     ODV_adj(t1,:) = [];
% end
% 
% % JP 02/04/24 I believe un0565 is the only float with Bbp532?
% if ~strcmp(info.MBARI_ID_str, 'un0565') 
% %if ~strcmp(info.float_type, 'NAVIS') %TMNEW
%     t1 = strcmp('b_bp532[1/m]',ODV_raw(:,1));
%     ODV_raw(t1,:) = [];
% 
%     t1 = strcmp('b_bp532[1/m]',ODV_adj(:,1));
%     ODV_adj(t1,:) = [];
% end
% 
% if isempty(iCHL435)
%     t1 = strncmp('Chl_a435',ODV_raw(:,1),8);
%     ODV_raw(t1,:) = [];
% 
%     t1 = strncmp('Chl_a435',ODV_adj(:,1),8);
%     ODV_adj(t1,:) = [];
% end
% 
% % OCR: get logicals, OCR in print hdr list, all OCR channels start with "DOWN"
% t1_raw = strncmp('DOWN',ODV_raw(:,3),4); 
% t1_adj = strncmp('DOWN',ODV_adj(:,3),4);
% if isempty(iOCR)
%     ODV_raw(t1_raw,:) = [];
%     ODV_adj(t1_adj,:) = [];
% else % figure out which channels to keep. JP 02/04/24
%     OCR_parms = raw_hdr(iOCR); %iOCR defined earlier (1x4 logical)
%     Lia       = ismember(ODV_raw(:,3), OCR_parms); % keep these (float OCR PARAMS)
%     ODV_raw(t1_raw & ~Lia,:) = [];
% 
%     OCR_parms = adj_hdr(iOCR);
%     Lia       = ismember(ODV_adj(:,3), OCR_parms); % keep these (float OCR PARAMS)
%     ODV_adj(t1_adj & ~Lia,:) = [];
% end
% 
% % if isempty(iOCR412)
% %     ind1 = find(strncmp('DOWN_IRRAD412',ODV_raw(:,1),13)   == 1);
% %     ODV_raw(ind1,:) = [];
% %     
% %     ind1 = find(strncmp('DOWN_IRRAD412',ODV_adj(:,1),13)   == 1);
% %     ODV_adj(ind1,:) = [];
% % end
% % 
% % if isempty(iOCR443)
% %     ind1 = find(strncmp('DOWN_IRRAD443',ODV_raw(:,1),13)   == 1);
% %     ODV_raw(ind1,:) = [];
% %     
% %     ind1 = find(strncmp('DOWN_IRRAD443',ODV_adj(:,1),13)   == 1);
% %     ODV_adj(ind1,:) = [];
% % end
% % 
% % if isempty(iOCR555)
% %     ind1 = find(strncmp('DOWN_IRRAD555',ODV_raw(:,1),13)   == 1);
% %     ODV_raw(ind1,:) = [];
% %     
% %     ind1 = find(strncmp('DOWN_IRRAD555',ODV_adj(:,1),13)   == 1);
% %     ODV_adj(ind1,:) = [];
% % end
% % 
% % if isempty(iOCRPAR)
% %     ind1 = find(strncmp('DOWNWELL_PAR',ODV_raw(:,1),12)   == 1);
% %     ODV_raw(ind1,:) = [];
% %     
% %     ind1 = find(strncmp('DOWNWELL_PAR',ODV_adj(:,1),12)   == 1);
% %     ODV_adj(ind1,:) = [];
% % end
% 
% 
% if tf_APEX_OCR % APEX OCR TEST FLOATS
%     tg = cellfun(@isempty, regexp(ODV_raw(:,1),'^Oxy|^Nit|^Chl|^b_bp','once'));
%     ODV_raw = ODV_raw(tg,:);
% 
%     tg = cellfun(@isempty, regexp(ODV_adj(:,1),'^Oxy|^Nit|^Chl|^b_bp|^POC','once'));
%     ODV_adj = ODV_adj(tg,:);
% end
% 
% if isempty(iO2)
%     t1 = strncmp('Oxygen2',ODV_raw(:,1),7);
%     ODV_raw(t1,:) = [];
% 
%     t1 = strncmp('Oxygen2',ODV_adj(:,1),7);
%     ODV_adj(t1,:) = [];
% else
%     tg = cellfun(@isempty, regexp(ODV_raw(:,1),'^Nit|^Chl|^b_bp','once'));
%     ODV_raw = ODV_raw(tg,:);
% 
%     tg = cellfun(@isempty, regexp(ODV_adj(:,1),'^Nit|^Chl|^b_bp|^POC','once'));
%     ODV_adj = ODV_adj(tg,:);
% end

raw_var_ct = size(ODV_raw,1);
adj_var_ct = size(ODV_adj,1);


% ************************************************************************
% ************************************************************************
% CREATE RAW ACII FILE - LOW RES
% ************************************************************************
% ************************************************************************

% PRINT META DATA HEADER LINES FIRST
disp(['Printing raw data to: ',dirs.FVlocal, info.FloatViz_name, '.txt']);
fid_raw  = fopen([dirs.FVlocal, info.FloatViz_name,'.TXT'],'W','n','UTF-8');

fprintf(fid_raw,'//0\r\n');
fprintf(fid_raw,'//<Encoding>UTF-8</Encoding>\r\n');

fprintf(fid_raw,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
    '\r\n']);
fprintf(fid_raw,['//WMO ID: ',info.WMO,'\r\n']);
fprintf(fid_raw,['//Institution ID: ',info.INST_ID_num,'\r\n']);
fprintf(fid_raw,['//MBARI ID: ',info.MBARI_ID_str,'\r\n']);
fprintf(fid_raw,['//Project Name: ',PROJ_Name,'\r\n']);
fprintf(fid_raw,['//Region: ',REGION,'\r\n']);

if strcmp(FLT_type,'SOLO')
    fprintf(fid_raw, SOLO_info_str);
end

if tf_APEX_OCR
    fprintf(fid_raw, TEST_OCR_info_str);
end

fprintf(fid_raw,['//', missing_profile_str,'\r\n']);
if ~isempty(interp_pos_str)
    fprintf(fid_raw,['//',interp_pos_str,'\r\n']);
end
if ~isempty(missing_pos_str)
    fprintf(fid_raw,['//',missing_pos_str,'\r\n']);
end
fprintf(fid_raw,['//',raw_missing_data_str,'\r\n']);
	fprintf(fid_raw,['//',iceEvRec_hdr_str,'\r\n']);
	fprintf(fid_raw,['//',iceEvRec_hdr_str2,'\r\n']);
	fprintf(fid_raw,['//',iceEvRec_hdr_str3,'\r\n']);
    fprintf(fid_raw,['//<MetaVariable>label="IceEvRec" value_type="INTEGER" significant_digits="0"</MetaVariable>','\r\n']);

% PRINT OUT SPECIAL NOTES IF THEY EXIST
if notes_flag == 1
    fprintf(fid_raw,'//\r\n//SPECIAL NOTICE:\r\n');
    for i = 1: size(notes,1)
        fprintf(fid_raw,'//%s\r\n',notes{i,1});
    end
end

% PRINT OUT FLOAT VARIABLE INFO
fprintf(fid_raw,'//\r\n//FLOAT VARIABLES:\r\n');
fprintf(fid_raw,'//Variable\tSensor\tSerial number\tComment\r\n');
for i = 1:raw_var_ct
    fprintf(fid_raw,'//%s\t%s\t%s\t%s\r\n',ODV_raw{i,1},ODV_raw{i,4}, ...
        ODV_raw{i,5},ODV_raw{i,6});
end
fprintf(fid_raw,'//\r\n');

% if strncmp('Miss',raw_missing_data_str,4)
%     fprintf(fid_raw,['//Missing data value = ',MVI_str,'\r\n']);
% end
fprintf(fid_raw,['//Missing data value = ',MVI_str,'\r\n']);

fprintf(fid_raw,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
    '1=Missing or not inspected \r\n']);
fprintf(fid_raw,'//Note: all timestamps are in GMT. \r\n');

% NOW PRINT THE RAW DATA HEADER
std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
    'Lon [�E]' 'Lat [�N]' 'QF' 'IceEvRec'}; % SIZE = 9
std_size = size(std_ODV_vars,2);

for i = 1:std_size % PRINT STANDARD HEADER VARS
    fprintf(fid_raw,'%s\t',std_ODV_vars{1,i}); % std vars
end

for i = 1:raw_var_ct % PRINT FLOAT SPECIFIC HEADER VARS
    if i < raw_var_ct
        fprintf(fid_raw,'%s\t%s\t',ODV_raw{i,1},'QF'); % std vars
    else
        fprintf(fid_raw,'%s\t%s\r\n',ODV_raw{i,1},'QF'); % std vars
    end
end

% BUILD DATA FORMAT STRING & CREATE DUMMY MATRIX TO DUMP DATA TO
% THIS CATCHES FLOATS WHICH DON'T HAVE SENSORS FOR GIVEN VARIABLES
% THESE COLUMNS WILL GET FILLED WITH -1e10 and QC FLAG = 1
dummy_out  = ones(raw_r, raw_var_ct) * NaN;
fill_MVI = ones(raw_r, 1) * -1e10;
fill_QC  = ones(raw_r, 1);
ODV_std_f = '%s\t%0.0f\t%s\t%s\t%s\t%0.3f\t%0.3f\t%0.0f\t%0.0f\t'; %std_vars
ODV_raw_f =''; ODV_space_f = '';
for i = 0:raw_var_ct-1
    c_ct = i*2+1; % Need to add data and QC col
    ind = find(strcmp(ODV_raw{i+1,3}, raw_hdr) == 1);
    
    if ~isempty(ind)
        if i < raw_var_ct-1
            ODV_raw_f   = [ODV_raw_f, ODV_raw{i+1,2},'\t%0.0f\t'];
            ODV_space_f = [ODV_space_f,'\t\t'];
        else
            ODV_raw_f   = [ODV_raw_f, ODV_raw{i+1,2},'\t%0.0f\r\n'];
            ODV_space_f = [ODV_space_f,'\t\r\n'];
        end
        dummy_out(:,c_ct:c_ct+1)  = raw_data(:,ind:ind+1);
    else
        if i < raw_var_ct-1
            ODV_raw_f   = [ODV_raw_f, '%1.0E\t%0.0f\t'];
            ODV_space_f = [ODV_space_f,'\t\t'];
        else
            ODV_raw_f   = [ODV_raw_f, '%1.0E\t%0.0f\r\n'];
            ODV_space_f = [ODV_space_f,'\t\r\n'];
        end
        dummy_out(:,c_ct:c_ct+1)  = [fill_MVI,fill_QC]; % No data, but need cols
    end
end

% NOW PRINT DATA LINES TO FILE
cast_num   = 0; %initalize
line_ct    = 0;
for sample_ct = 1 : raw_r
    if raw_data(sample_ct,1) - cast_num > 0 % build standard cast part
        if sample_ct > 1
            fprintf(fid_raw,[std_str,ODV_space_f]); % add profile spacer line
            line_ct = line_ct+1;
        end
        cast_num = raw_data(sample_ct,1);
        date_str = datestr(raw_data(sample_ct,2),'mm/dd/yyyy');
        time_str = datestr(raw_data(sample_ct,2),'HH:MM');
        std_str  = sprintf(ODV_std_f, info.WMO, cast_num, 'C', ...
            date_str, time_str, raw_data(sample_ct,3), ...
            raw_data(sample_ct,4), raw_data(sample_ct,5), raw_data(sample_ct,6));
        std_str = regexprep(std_str,'NaN',MVI_str);
    end
    data_str = sprintf(ODV_raw_f, dummy_out(sample_ct,:));
    out_str = [std_str,data_str];
    % replace NaN' w/ missing value indicator
    out_str = regexprep(out_str,'NaN',MVI_str);  % replace NaN' w/ missing value indicator
    fprintf(fid_raw, '%s', out_str);
    line_ct = line_ct+1;
    clear out_str data_str date_str time_str
    if sample_ct == raw_r
        fprintf(fid_raw,[std_str,ODV_space_f]); % add profile spacer line
        line_ct = line_ct+1;
    end
end

fclose(fid_raw);
clear fid_raw cast_num sample_ct

% MAKE CONFIG FILE
%fid_raw  = fopen([dirs.FVlocal, info.FloatViz_name, '.CFG'],'w');
fid_raw  = fopen([dirs.FVlocal, info.FloatViz_name, '.CFG'],'W','n','UTF-8');

fprintf(fid_raw,'//%0.0f\r\n',line_ct);
fclose(fid_raw);
clear fid_raw line_ct dummy_out

% ************************************************************************
% ************************************************************************
% CREATE ADJUSTED ACII FILE - LOW RES
% ************************************************************************
% ************************************************************************

% START REPLACEMENT OF PSAL -----------------------
% Ok, first thing's first -- if a psal-proxy float, re-insert original psal!  (we don't want the proxy in the final data files)
if isfield(info,'PSAL_PROXY_USED')
    if info.PSAL_PROXY_USED %we never deleted this index...
	disp(['%%%%% NOTE: RE-INSERTIN ORIGINAL PSAL FOR FINAL WRITE-TO-FILE, FLOAT ',WMO_ID,'. %%%%'])
	adj_data(:,iaS) = psal_orig;
    % OK need to convert to ODV flag...this was already done earlier for
    % the other arrays. Keep it simple and clean for this block.
    psalqc_orig(psalqc_orig == 5) = 3; 
    psalqc_orig(psalqc_orig == 4) = 8;
    psalqc_orig(psalqc_orig == 3) = 4; 
    psalqc_orig(psalqc_orig == 1) = 0; 
	adj_data(:,iaS+1) = psalqc_orig;
    end
end
% END REPLACEMENT OF PSAL -----------------------

% TEST FOR QC ADJUSTMENT CONSTANTS, GET IF THEY EXIST
QC = get_QC_adjustments(info.FloatViz_name,dirs);
QC_check = 1; % Adjustments exist
if isempty(QC)
    QC_check = 0; % NO Adjustments
end

if QC_check == 1
    % PRINT HEADER FIRST
    disp(['Printing adjusted data to: ',dirs.FVlocal, 'QC\', ...
        info.FloatViz_name,'QC', '.txt']);
    fid_adj  = fopen([dirs.FVlocal,'QC\', info.FloatViz_name,'QC','.TXT'],'W','n','UTF-8');
    fprintf(fid_adj,'//0\r\n');
    fprintf(fid_adj,'//<Encoding>UTF-8</Encoding>\r\n');
    fprintf(fid_adj,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
        '\r\n']);
    fprintf(fid_adj,'//!! ADJUSTED DATA FILE !!\r\n');
    
    fprintf(fid_adj,['//WMO ID: ',info.WMO,'\r\n']);
    fprintf(fid_adj,['//Institution ID: ',info.INST_ID_num,'\r\n']);
    fprintf(fid_adj,['//MBARI ID: ',info.MBARI_ID_str,'\r\n']);
    fprintf(fid_adj,['//Project Name: ',PROJ_Name,'\r\n']);
    fprintf(fid_adj,['//Region: ',REGION,'\r\n']);
    
    if strcmp(FLT_type,'SOLO')
        fprintf(fid_adj, SOLO_info_str);
    end
    
    if tf_APEX_OCR
        fprintf(fid_raw, TEST_OCR_info_str);
    end
    
    fprintf(fid_adj,['//', missing_profile_str,'\r\n']);
    if ~isempty(interp_pos_str)
        fprintf(fid_adj,['//',interp_pos_str,'\r\n']);
    end
    if ~isempty(missing_pos_str)
        fprintf(fid_adj,['//',missing_pos_str,'\r\n']);
    end
    fprintf(fid_adj,['//',adj_missing_data_str,'\r\n']);
    fprintf(fid_adj,['//',iceEvRec_hdr_str,'\r\n']);
	fprintf(fid_adj,['//',iceEvRec_hdr_str2,'\r\n']);
	fprintf(fid_adj,['//',iceEvRec_hdr_str3,'\r\n']);
    fprintf(fid_adj,['//<MetaVariable>label="IceEvRec" value_type="INTEGER" significant_digits="0"</MetaVariable>','\r\n']);

    % PRINT OUT SPECIAL NOTES IF THEY EXIST
    if notes_flag == 1
        fprintf(fid_adj,'//\r\n//SPECIAL NOTICE:\r\n');
        for i = 1: size(notes,1)
            fprintf(fid_adj,'//%s\r\n',notes{i,1});
        end
    end
    
    % PRINT OUT FLOAT VARIABLE INFO
    fprintf(fid_adj,'//\r\n//FLOAT VARIABLES:\r\n');
    fprintf(fid_adj,'//Variable\tSensor\tSerial number\tComment\r\n');
    for i = 1:adj_var_ct
        fprintf(fid_adj,'//%s\t%s\t%s\t%s\r\n',ODV_adj{i,1},ODV_adj{i,4}, ...
            ODV_adj{i,5},ODV_adj{i,6});
    end
    
    if ~isempty(iPH)
        fprintf(fid_adj,'//\r\n');
        fprintf(fid_adj,CO2SYS_str); % Print out CO2SYS statement
        %fprintf(fid_adj,'//\r\n');
    end
    
    if ~isempty(iCHL)
        fprintf(fid_adj,'//\r\n');
        fprintf(fid_adj,CHL_str); % Print out CHL statement
    end
    
    
    % PRINT OUT FLOAT VARIABLE QC CORRECTION INFO
    fprintf(fid_adj,'//\r\n//QUALITY CONTROLLED DATA CORRECTIONS:\r\n');
    fprintf(fid_adj,'//Measurement\tStation\tGain\tOffset\tDrift\r\n');
    possible_fields = {'O' 'N' 'pH' 'CHL' 'BB' 'CDOM'};
    possible_hdr_names = {'Oxygen' 'Nitrate' 'pH' 'Chl' 'BB' 'CDOM'};
    for i = 1 : size(possible_fields,2)
        if isfield(QC, possible_fields{i}) %
            tmp = QC.(possible_fields{i});
            if isfield(tmp,'steps')
                steps = tmp.steps; % QC constants in here
            else
                continue
            end
        else
            continue
        end
        if ~isempty(steps) 
            [QCr,QCc] = size(steps);
            for j = 1: QCr
                fprintf(fid_adj, ['//',possible_hdr_names{i},'\t', ...
                    '%0.0f\t%0.4f\t%0.4f\t%0.4f\r\n'], steps(j,2:QCc));
            end
        end
    end
    fprintf(fid_adj,'//\r\n');

    % Add note if raw data are used to populated adjusted data
    tmp_fields ={'Chl_data_str' 'Chl435_data_str' 'Bbp_data_str', ...
        'Cdm_data_str' 'Bbp532_data_str' 'OCR_data_str'};
    for i = 1: size(tmp_fields,2)
        if ~isempty(info.(tmp_fields{i}))
            fprintf(fid_adj, '//%s\r\n', info.(tmp_fields{i}));
        end
    end



    info.Chl_data_str    = '';
    info.Chl435_data_str = '';
    info.Bbp_data_str    = '';
    info.Cdm_data_str    = '';
    info.Bbp532_data_str = '';
    info.OCR_data_str    = '';

    
    if strncmp('Miss',adj_missing_data_str,4)
        fprintf(fid_adj,['//Missing data value = ',MVI_str,'\r\n']);
    end
    fprintf(fid_adj,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
        '1=Missing or not inspected \r\n']);
    fprintf(fid_adj,['//Note: all timestamps are in GMT. \r\n']);
    
std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
    'Lon [�E]' 'Lat [�N]' 'QF' 'IceEvRec'}; % SIZE = 9
    std_size = size(std_ODV_vars,2);
    
    % ************************************************************************
    % PRINT THE ADJUSTED DATA HEADER
    std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
        'Lon [�E]' 'Lat [�N]' 'QF' 'IceEvRec'}; % SIZE = 9
    std_size = size(std_ODV_vars,2);
    
    for i = 1:std_size % PRINT STANDARD HEADER VARS
        fprintf(fid_adj,'%s\t',std_ODV_vars{1,i}); % std vars
    end
    
    for i = 1:adj_var_ct % PRINT FLOAT SPECIFIC HEADER VARS
        if i < adj_var_ct
            fprintf(fid_adj,'%s\t%s\t',ODV_adj{i,1},'QF'); % std vars
        else
            fprintf(fid_adj,'%s\t%s\r\n',ODV_adj{i,1},'QF'); % std vars
        end
    end
    
    % BUILD DATA FORMAT STRING & CREATE DUMMY MATRIX TO DUMP DATA TO
    % THIS CATCHES FLOATS WHICH DON'T HAVE SENSORS FOR GIVEN VARIABLES
    % THESE COLUMNS WILL GET FILLED WITH -1e10 and QC FLAG = 1
    dummy_out  = ones(adj_r, adj_var_ct) * NaN;
    fill_MVI = ones(adj_r, 1) * -1e10;
    fill_QC  = ones(adj_r, 1);
    ODV_std_f = '%s\t%0.0f\t%s\t%s\t%s\t%0.3f\t%0.3f\t%0.0f\t%0.0f\t'; %std_vars
    ODV_adj_f =''; ODV_space_f = '';
    
    for i = 0:adj_var_ct-1
        c_ct = i*2+1; % Need to add data and QC col
        ind = find(strcmp(ODV_adj{i+1,3}, adj_hdr) == 1);
        
        if ~isempty(ind)
            if i < adj_var_ct-1
                ODV_adj_f   = [ODV_adj_f, ODV_adj{i+1,2},'\t%0.0f\t'];
                ODV_space_f = [ODV_space_f,'\t\t'];
            else
                ODV_adj_f   = [ODV_adj_f, ODV_adj{i+1,2},'\t%0.0f\r\n'];
                ODV_space_f = [ODV_space_f,'\t\r\n'];
            end
            dummy_out(:,c_ct:c_ct+1)  = adj_data(:,ind:ind+1);
        else
            if i < adj_var_ct-1
                ODV_adj_f   = [ODV_adj_f, '%1.0E\t%0.0f\t'];
                ODV_space_f = [ODV_space_f,'\t\t'];
            else
                ODV_adj_f   = [ODV_adj_f, '%1.0E\t%0.0f\r\n'];
                ODV_space_f = [ODV_space_f,'\t\r\n'];
            end
            dummy_out(:,c_ct:c_ct+1)  = [fill_MVI,fill_QC]; % No data, but need cols
        end
    end
    
    % NOW PRINT DATA LINES TO FILE
    cast_num   = 0; %initalize
    line_ct    = 0;
    for sample_ct = 1 : adj_r
        if adj_data(sample_ct,1) - cast_num > 0 % build standard cast part
            if sample_ct > 1
                fprintf(fid_adj,[std_str,ODV_space_f]); % add profile spacer line
                line_ct = line_ct+1;
            end
            cast_num = adj_data(sample_ct,1);
            date_str = datestr(adj_data(sample_ct,2),'mm/dd/yyyy');
            time_str = datestr(adj_data(sample_ct,2),'HH:MM');
            std_str  = sprintf(ODV_std_f, info.WMO, cast_num, 'C', ...
                date_str, time_str, adj_data(sample_ct,3), ...
                adj_data(sample_ct,4), adj_data(sample_ct,5), adj_data(sample_ct,6));
            std_str = regexprep(std_str,'NaN',MVI_str);
        end
        data_str = sprintf(ODV_adj_f, dummy_out(sample_ct,:));
        out_str = [std_str,data_str];
        % replace NaN' w/ missing value indicator
        out_str = regexprep(out_str,'NaN',MVI_str);  % replace NaN' w/ missing value indicator
        fprintf(fid_adj, '%s', out_str);
        line_ct = line_ct+1;
        clear out_str data_str date_str time_str
        if sample_ct == adj_r
            fprintf(fid_adj,[std_str,ODV_space_f]); % add profile spacer line
            line_ct = line_ct+1;
        end
    end
    
    fclose(fid_adj);
    clear fid_adj cast_num sample_ct dummy_out
    
    % MAKE CONFIG FILE
    %fid_adj  = fopen([dirs.FVlocal,'QC\', info.FloatViz_name,'QC', '.CFG'],'w');
    fid_adj  = fopen([dirs.FVlocal,'QC\', info.FloatViz_name,'QC', '.CFG'],'W','n','UTF-8');
    
    fprintf(fid_adj,'//%0.0f\r\n',line_ct);
    fclose(fid_adj);
    clear fid_adj line_ct adj_out
end

% ************************************************************************
% ************************************************************************
% CREATE COMBINED HR+LR RAW ACII FILE
% ************************************************************************
% ************************************************************************
% PRINT HEADER FIRST
if strcmp(info.float_type, 'APEX') && ~tf_APEX_OCR && HR_flag == 1
    [allraw_r,allraw_c] = size(allraw_data);
    disp(['Printing HR+LR raw data to: ',dirs.FVlocal, 'HR\' ...
        info.FloatViz_name, '_HR.txt']);
    %     fid_raw  = fopen([dirs.FVlocal,'HR\' info.FloatViz_name,'_HR.TXT'],'W');
    fid_raw  = fopen([dirs.FVlocal,'HR\' info.FloatViz_name,'_HR.TXT'],'W','n','UTF-8');
    fprintf(fid_raw,'//0\r\n');
    fprintf(fid_raw,'//<Encoding>UTF-8</Encoding>\r\n');
    fprintf(fid_raw,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
        '\r\n']);
    fprintf(fid_raw,['//WMO ID: ',info.WMO,'\r\n']);
    fprintf(fid_raw,['//Institution ID: ',info.INST_ID_num,'\r\n']);
    fprintf(fid_raw,['//MBARI ID: ',info.MBARI_ID_str,'\r\n']);
    fprintf(fid_raw,['//Project Name: ',PROJ_Name,'\r\n']);
    fprintf(fid_raw,['//Region: ',REGION,'\r\n']);
    fprintf(fid_raw,['//', missing_profile_str,'\r\n']);
    if ~isempty(interp_pos_str)
        fprintf(fid_raw,['//',interp_pos_str,'\r\n']);
    end
    if ~isempty(missing_pos_str)
        fprintf(fid_raw,['//',missing_pos_str,'\r\n']);
    end
    fprintf(fid_raw,['//',raw_missing_data_str,'\r\n']);
	fprintf(fid_raw,['//',iceEvRec_hdr_str,'\r\n']);
	fprintf(fid_raw,['//',iceEvRec_hdr_str2,'\r\n']);
	fprintf(fid_raw,['//',iceEvRec_hdr_str3,'\r\n']);
    fprintf(fid_raw,['//<MetaVariable>label="IceEvRec" value_type="INTEGER" significant_digits="0"</MetaVariable>','\r\n']);
    % PRINT OUT SPECIAL NOTES IF THEY EXIST
    if notes_flag == 1
        fprintf(fid_raw,'//\r\n//SPECIAL NOTICE:\r\n');
        for i = 1: size(notes,1)
            fprintf(fid_raw,'//%s\r\n',notes{i,1});
        end
    end
    
    % PRINT OUT FLOAT VARIABLE INFO
    fprintf(fid_raw,'//\r\n//FLOAT VARIABLES:\r\n');
    fprintf(fid_raw,'//Variable\tSensor\tSerial number\tComment\r\n');
    for i = 1:raw_var_ct
        fprintf(fid_raw,'//%s\t%s\t%s\t%s\r\n',ODV_raw{i,1},ODV_raw{i,4}, ...
            ODV_raw{i,5},ODV_raw{i,6});
    end
    fprintf(fid_raw,'//\r\n');
    
    if strncmp('Miss',raw_missing_data_str,4)
        fprintf(fid_raw,['//Missing data value = ',MVI_str,'\r\n']);
    end
    fprintf(fid_raw,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
        '1=Missing or not inspected \r\n']);
    fprintf(fid_raw,['//Note: all timestamps are in GMT. \r\n']);
    
    % NOW PRINT THE RAW DATA HEADER
    std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
        'Lon [�E]' 'Lat [�N]' 'QF' 'IceEvRec'}; % SIZE = 9
    std_size = size(std_ODV_vars,2);
    
    for i = 1:std_size % PRINT STANDARD HEADER VARS
        fprintf(fid_raw,'%s\t',std_ODV_vars{1,i}); % std vars
    end
    
    for i = 1:raw_var_ct % PRINT FLOAT SPECIFIC HEADER VARS
        if i < raw_var_ct
            fprintf(fid_raw,'%s\t%s\t',ODV_raw{i,1},'QF'); % std vars
        else
            fprintf(fid_raw,'%s\t%s\r\n',ODV_raw{i,1},'QF'); % std vars
        end
    end
    
    % BUILD DATA FORMAT STRING & CREATE DUMMY MATRIX TO DUMP DATA TO
    % THIS CATCHES FLOATS WHICH DON'T HAVE SENSORS FOR GIVEN VARIABLES
    % THESE COLUMNS WILL GET FILLED WITH -1e10 and QC FLAG = 1
    dummy_out = ones(allraw_r, raw_var_ct) * NaN;
    fill_MVI  = ones(allraw_r, 1) * -1e10;
    fill_QC   = ones(allraw_r, 1);
    ODV_std_f = '%s\t%0.0f\t%s\t%s\t%s\t%0.3f\t%0.3f\t%0.0f\t%0.0f\t'; %std_vars
    ODV_raw_f =''; ODV_space_f = '';
    for i = 0:raw_var_ct-1
        c_ct = i*2+1; % Need to add data and QC col
        ind = find(strcmp(ODV_raw{i+1,3}, raw_hdr) == 1);
        
        if ~isempty(ind)
            if i < raw_var_ct-1
                ODV_raw_f   = [ODV_raw_f, ODV_raw{i+1,2},'\t%0.0f\t'];
                ODV_space_f = [ODV_space_f,'\t\t'];
            else
                ODV_raw_f   = [ODV_raw_f, ODV_raw{i+1,2},'\t%0.0f\r\n'];
                ODV_space_f = [ODV_space_f,'\t\r\n'];
            end
            dummy_out(:,c_ct:c_ct+1)  = allraw_data(:,ind:ind+1);
        else
            if i < raw_var_ct-1
                ODV_raw_f   = [ODV_raw_f, '%1.0E\t%0.0f\t'];
                ODV_space_f = [ODV_space_f,'\t\t'];
            else
                ODV_raw_f   = [ODV_raw_f, '%1.0E\t%0.0f\r\n'];
                ODV_space_f = [ODV_space_f,'\t\r\n'];
            end
            dummy_out(:,c_ct:c_ct+1)  = [fill_MVI,fill_QC]; % No data, but need cols
        end
    end
    
    % NOW PRINT DATA LINES TO FILE
    cast_num   = 0; %initalize
    line_ct    = 0;
    for sample_ct = 1 : allraw_r
        if allraw_data(sample_ct,1) - cast_num > 0 % build standard cast part
            if sample_ct > 1
                fprintf(fid_raw,[std_str,ODV_space_f]); % add profile spacer line
                line_ct = line_ct+1;
            end
            cast_num = allraw_data(sample_ct,1);
            date_str = datestr(allraw_data(sample_ct,2),'mm/dd/yyyy');
            time_str = datestr(allraw_data(sample_ct,2),'HH:MM');
            std_str  = sprintf(ODV_std_f, info.WMO, cast_num, 'C', ...
                date_str, time_str, allraw_data(sample_ct,3), ...
                allraw_data(sample_ct,4), allraw_data(sample_ct,5), allraw_data(sample_ct,6));
            std_str = regexprep(std_str,'NaN',MVI_str);
        end
        data_str = sprintf(ODV_raw_f, dummy_out(sample_ct,:));
        out_str = [std_str,data_str];
        % replace NaN' w/ missing value indicator
        out_str = regexprep(out_str,'NaN',MVI_str);  % replace NaN' w/ missing value indicator
        fprintf(fid_raw, '%s', out_str);
        line_ct = line_ct+1;
        clear out_str data_str date_str time_str
        if sample_ct == allraw_r
            fprintf(fid_raw,[std_str,ODV_space_f]); % add profile spacer line
            line_ct = line_ct+1;
        end
    end
    
    fclose(fid_raw);
    clear fid_raw cast_num sample_ct
    
    % MAKE CONFIG FILE
    %fid_raw  = fopen([dirs.FVlocal,'HR\', info.FloatViz_name, '_HR.CFG'],'w');
    fid_raw  = fopen([dirs.FVlocal,'HR\', info.FloatViz_name, '_HR.CFG'],'W','n','UTF-8');
    
    fprintf(fid_raw,'//%0.0f\r\n',line_ct);
    fclose(fid_raw);
    clear fid_raw line_ct dummy_out
end

% ************************************************************************
% ************************************************************************
% CREATE COMBINED HR+LR ADJUSTED ACII FILE
% ************************************************************************
% ************************************************************************
if QC_check == 1 && HR_flag == 1
    % PRINT HEADER FIRST
    if strcmp(info.float_type, 'APEX') && ~tf_APEX_OCR
        [alladj_r,alladj_c] = size(allraw_data);
        disp(['Printing HR+LR adjusted data to: ',dirs.FVlocal, ...
            'HRQC\', info.FloatViz_name, '_HRQC', '.txt']);
        %         fid_adj  = fopen([dirs.FVlocal,'HRQC\', info.FloatViz_name, ...
        %             '_HRQC','.TXT'],'W');
        fid_adj  = fopen([dirs.FVlocal,'HRQC\', info.FloatViz_name, ...
            '_HRQC','.TXT'],'W','n','UTF-8');
        
        fprintf(fid_adj,'//0\r\n');
        fprintf(fid_adj,'//<Encoding>UTF-8</Encoding>\r\n');
        fprintf(fid_adj,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
            '\r\n']);
        fprintf(fid_adj,'//!! ADJUSTED DATA FILE !!\r\n');
        
        fprintf(fid_adj,['//WMO ID: ',info.WMO,'\r\n']);
        fprintf(fid_adj,['//Institution ID: ',info.INST_ID_num,'\r\n']);
        fprintf(fid_adj,['//MBARI ID: ',info.MBARI_ID_str,'\r\n']);
        fprintf(fid_adj,['//Project Name: ',PROJ_Name,'\r\n']);
        fprintf(fid_adj,['//Region: ',REGION,'\r\n']);
        
        fprintf(fid_adj,['//', missing_profile_str,'\r\n']);
        if ~isempty(interp_pos_str)
            fprintf(fid_adj,['//',interp_pos_str,'\r\n']);
        end
        if ~isempty(missing_pos_str)
            fprintf(fid_adj,['//',missing_pos_str,'\r\n']);
        end
        fprintf(fid_adj,['//',adj_missing_data_str,'\r\n']);
        fprintf(fid_adj,['//',iceEvRec_hdr_str,'\r\n']);
		fprintf(fid_adj,['//',iceEvRec_hdr_str2,'\r\n']);
		fprintf(fid_adj,['//',iceEvRec_hdr_str3,'\r\n']);
        fprintf(fid_adj,['//<MetaVariable>label="IceEvRec" value_type="INTEGER" significant_digits="0"</MetaVariable>','\r\n']);
        % PRINT OUT SPECIAL NOTES IF THEY EXIST
        if notes_flag == 1
            fprintf(fid_adj,'//\r\n//SPECIAL NOTICE:\r\n');
            for i = 1: size(notes,1)
                fprintf(fid_adj,'//%s\r\n',notes{i,1});
            end
        end
        
        % PRINT OUT FLOAT VARIABLE INFO
        fprintf(fid_adj,'//\r\n//FLOAT VARIABLES:\r\n');
        fprintf(fid_adj,'//Variable\tSensor\tSerial number\tComment\r\n');
        for i = 1:adj_var_ct
            fprintf(fid_adj,'//%s\t%s\t%s\t%s\r\n',ODV_adj{i,1},ODV_adj{i,4}, ...
                ODV_adj{i,5},ODV_adj{i,6});
        end
        
        if ~isempty(iPH)
            fprintf(fid_adj,'//\r\n');
            fprintf(fid_adj,CO2SYS_str); % Print out CO2SYS statement
            %fprintf(fid_adj,'//\r\n');
        end
        
        if ~isempty(iCHL)
            fprintf(fid_adj,'//\r\n');
            fprintf(fid_adj,CHL_str); % Print out CHL statement
        end
        
        % PRINT OUT FLOAT VARIABLE QC CORRECTION INFO
        fprintf(fid_adj,'//\r\n//QUALITY CONTROLLED DATA CORRECTIONS:\r\n');
        fprintf(fid_adj,'//Measurement\tStation\tGain\tOffset\tDrift\r\n');
        possible_fields = {'O' 'N' 'pH' 'CHL' 'BB' 'CDOM'};
        possible_hdr_names = {'Oxygen' 'Nitrate' 'pH' 'Chl' 'BB' 'CDOM'};
        for i = 1 : size(possible_fields,2)
            if isfield(QC, possible_fields{i}) %
                tmp = QC.(possible_fields{i});
                if isfield(tmp,'steps')
                    steps = tmp.steps; % QC constants in here
                else
                    continue
                end
            else
                continue
            end
            if ~isempty(steps) %TM 3/3/21 again...is this block necessary?
                [QCr,QCc] = size(steps);
                if strcmp(possible_fields{i}, 'CHL') % ALL CHL: gain = 2, offest 0
                    steps(:,3) = 1/2;
                    steps(:,4) = 0;
                end
                for j = 1: QCr
                    fprintf(fid_adj, ['//',possible_hdr_names{i},'\t', ...
                        '%0.0f\t%0.4f\t%0.4f\t%0.4f\r\n'], steps(j,2:QCc));
                end
            end
        end
        fprintf(fid_adj,'//\r\n');
        
        
        
        if strncmp('Miss',adj_missing_data_str,4)
            fprintf(fid_adj,['//Missing data value = ',MVI_str,'\r\n']);
        end
        fprintf(fid_adj,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
            '1=Missing or not inspected \r\n']);
        fprintf(fid_adj,['//Note: all timestamps are in GMT. \r\n']);
        
        % NOW PRINT THE ADJUSTED DATA HEADER
std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
    'Lon [�E]' 'Lat [�N]' 'QF' 'IceEvRec'}; % SIZE = 9
        std_size = size(std_ODV_vars,2);
        
        for i = 1:std_size % PRINT STANDARD HEADER VARS
            fprintf(fid_adj,'%s\t',std_ODV_vars{1,i}); % std vars
        end
        
        for i = 1:adj_var_ct % PRINT FLOAT SPECIFIC HEADER VARS
            if i < adj_var_ct
                fprintf(fid_adj,'%s\t%s\t',ODV_adj{i,1},'QF'); % std vars
            else
                fprintf(fid_adj,'%s\t%s\r\n',ODV_adj{i,1},'QF'); % std vars
            end
        end
        
        % BUILD DATA FORMAT STRING & CREATE DUMMY MATRIX TO DUMP DATA TO
        % THIS CATCHES FLOATS WHICH DON'T HAVE SENSORS FOR GIVEN VARIABLES
        % THESE COLUMNS WILL GET FILLED WITH -1e10 and QC FLAG = 1
        dummy_out = ones(alladj_r, adj_var_ct) * NaN;
        fill_MVI  = ones(alladj_r, 1) * -1e10;
        fill_QC   = ones(alladj_r, 1);
        ODV_std_f = '%s\t%0.0f\t%s\t%s\t%s\t%0.3f\t%0.3f\t%0.0f\t%0.0f\t'; %std_vars
        ODV_adj_f =''; ODV_space_f = '';
        
        for i = 0:adj_var_ct-1
            c_ct = i*2+1; % Need to add data and QC col
            ind = find(strcmp(ODV_adj{i+1,3}, adj_hdr) == 1);
            
            if ~isempty(ind)
                if i < adj_var_ct-1
                    ODV_adj_f   = [ODV_adj_f, ODV_adj{i+1,2},'\t%0.0f\t'];
                    ODV_space_f = [ODV_space_f,'\t\t'];
                else
                    ODV_adj_f   = [ODV_adj_f, ODV_adj{i+1,2},'\t%0.0f\r\n'];
                    ODV_space_f = [ODV_space_f,'\t\r\n'];
                end
                dummy_out(:,c_ct:c_ct+1)  = alladj_data(:,ind:ind+1);
            else
                if i < adj_var_ct-1
                    ODV_adj_f   = [ODV_adj_f, '%1.0E\t%0.0f\t'];
                    ODV_space_f = [ODV_space_f,'\t\t'];
                else
                    ODV_adj_f   = [ODV_adj_f, '%1.0E\t%0.0f\r\n'];
                    ODV_space_f = [ODV_space_f,'\t\r\n'];
                end
                dummy_out(:,c_ct:c_ct+1)  = [fill_MVI,fill_QC]; % No data, but need cols
            end
        end
        
        
        % NOW PRINT DATA LINES TO FILE
        cast_num   = 0; %initalize
        line_ct    = 0;
        for sample_ct = 1 : alladj_r
            if alladj_data(sample_ct,1) - cast_num > 0 % build standard cast part
                if sample_ct > 1
                    fprintf(fid_adj,[std_str,ODV_space_f]); % add profile spacer line
                    line_ct = line_ct+1;
                end
                cast_num = alladj_data(sample_ct,1);
                date_str = datestr(alladj_data(sample_ct,2),'mm/dd/yyyy');
                time_str = datestr(alladj_data(sample_ct,2),'HH:MM');
                std_str  = sprintf(ODV_std_f, info.WMO, cast_num, 'C', ...
                    date_str, time_str, alladj_data(sample_ct,3), ...
                    alladj_data(sample_ct,4), alladj_data(sample_ct,5), alladj_data(sample_ct,6));
                std_str = regexprep(std_str,'NaN',MVI_str);
            end
            data_str = sprintf(ODV_adj_f, dummy_out(sample_ct,:));
            out_str = [std_str,data_str];
            % replace NaN' w/ missing value indicator
            out_str = regexprep(out_str,'NaN',MVI_str);  % replace NaN' w/ missing value indicator
            fprintf(fid_adj, '%s', out_str);
            line_ct = line_ct+1;
            clear out_str data_str date_str time_str
            if sample_ct == alladj_r
                fprintf(fid_adj,[std_str,ODV_space_f]); % add profile spacer line
                line_ct = line_ct+1;
            end
        end
        
        fclose(fid_adj);
        clear fid_raw cast_num sample_ct
        
        
        % MAKE CONFIG FILE
        %fid_adj  = fopen([dirs.FVlocal,'HRQC\', info.FloatViz_name,'_HRQC', '.CFG'],'w');
        fid_adj  = fopen([dirs.FVlocal,'HRQC\', info.FloatViz_name,...
            '_HRQC', '.CFG'],'W','n','UTF-8');
        
        fprintf(fid_adj,'//%0.0f\r\n',line_ct);
        fclose(fid_adj);
    end
    clear fid_raw line_ct
end
tf_odv = 1;

