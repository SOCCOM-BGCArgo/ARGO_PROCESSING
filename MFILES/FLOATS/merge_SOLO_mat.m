function d = merge_SOLO_mat(WMO_ID, dirs)
% This function merges all BGC sensors axes for a given cycle & then merges
% all cycles into one matrix in prep for printing to an ODV style *.TXT file

% CHANGES:
% 07/18/2022 JP - Updated  matix column assignments to deal with seansoe
%     being turned off. float param list determined from 1st cycle file,
%     then check field exisitence every time before assigning data

% *************************************************************************
% % TESTING
% % WMO_ID = '4903026'; % ss0001
% WMO_ID = '5906320'; % ua19191
% WMO_ID = '5906767'; % ss0004
% dirs   = [];

% fd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\5906320\';
% fn = '5906320.009.mat';

% *************************************************************************




% *************************************************************************
%  SET DATA DIRS AND PATHS
% *************************************************************************
RC.P    = [0 12000]; % from argo parameter list
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS. ONLY TWO ARE NEEDED FOR THIS FUNCTION TO OPERATE. ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end
list_fp  = [dirs.cal, 'MBARI_float_list.mat']; % MBARI master list file path

% *************************************************************************
% DO SOME PREP WORK
% *************************************************************************

%  BUILD LIST OF POSSIBLE BGC SOLO ARGO VARIABLES (UPDATE WHEN NEW VARIABLES ARISE)
%      *_QC, *_ADJUSTED, & *_ADJUSTED_QC BUILT FROM RAW PARAMETER NAME
raw_vars{1,1}  = 'PRES';  
raw_vars{2,1}  = 'TEMP';
raw_vars{3,1}  = 'PSAL'; 
raw_vars{4,1}  = 'DOXY';
raw_vars{5,1}  = 'DOXY2';
raw_vars{6,1}  = 'NITRATE'; 
raw_vars{7,1}  = 'CHLA';  
raw_vars{8,1}  = 'BBP700';
raw_vars{9,1}  = 'CDOM';     
raw_vars{10,1} = 'BBP532';
raw_vars{11,1} = 'PH_IN_SITU_TOTAL'; 
raw_vars{12,1} = 'DOWN_IRRADIANCE380'; 
raw_vars{13,1} = 'DOWN_IRRADIANCE412'; 
raw_vars{14,1} = 'DOWN_IRRADIANCE490'; 
raw_vars{15,1} = 'DOWNWELLING_PAR'; 

var_list       = raw_vars; % save a copy before qc params added

adj_vars = strcat(raw_vars,'_ADJUSTED');         % *_ADJUSTED names
raw_vars = [raw_vars'; strcat(raw_vars,'_QC')']; % add *_QC parameters
raw_vars = raw_vars(:);                % reshape so param & param qc pairs
adj_vars = [adj_vars'; strcat(adj_vars,'_QC')'];
adj_vars = adj_vars(:);

vars_ct = size(raw_vars,1); % raw & adjusted counts will be the same

% *************************************************************************
%   GET MASTER LIST, FIND FLOAT ON IT % GET *.MAT FILE LISTING FOR FLOAT
% *************************************************************************

% LOAD MBARI MASTER LIST
if exist(list_fp, 'file') % load list variables
    load(list_fp)
else
    disp('BUILDING FLOAT ID LIST ...')
    % leave dirs structure empty; the function will build it, ensuring all
    % fields are filled appropriately.
    d = MBARI_float_list([]); 
end
iWMO  = find(strcmp('WMO',d.hdr)             == 1); % col indexes
iMB   = find(strcmp('MBARI ID',d.hdr)        == 1);
iINST = find(strcmp('INST ID',d.hdr)         == 1);
ind   = find(strcmpi(WMO_ID, d.list(:,iWMO)) == 1); %row index for float

if isempty(ind)
    disp(['No WMO_ID match found in lookup table for ',WMO_ID])
    disp('NO MATCH FOUND - EXITING')
    d = []; % Return this if no data
    return
end
list = d.list(ind,:);
clear d

INST_ID_str  = list{iINST};
MBARI_ID_str = list{iMB};

disp(['MBARI ID for WMO ',WMO_ID, ' = ', MBARI_ID_str]);
disp(['Institution ID for ',WMO_ID, ' = ', INST_ID_str]);

% GET LIST OF *.MAT FILES IN THE WMO DIR
float_dir = fullfile(dirs.mat, WMO_ID,'\');
if isfolder(float_dir)
    tmp       = dir([float_dir,WMO_ID,'*mat']);
    file_list = {tmp.name}';
    tf        = cellfun(@isempty,regexp(file_list,'000\.mat$','once'));
    file_list = file_list(tf);
    r_list    = size(file_list,1);
    if r_list == 0 %directory exists but no .mat files exist yet.
        disp(['WMO directory found for ', MBARI_ID_str, ' (', WMO_ID, '), but no .mat files exist yet.']);
        d = []; % Return this if no data
        return
    end
else
    disp(['No WMO directory found for ', MBARI_ID_str, ' (', WMO_ID, ')']);
    d = []; % Return this if no data
    return
end


% ************************************************************************
%  PROCESS *.MAT FILES
% ************************************************************************
CTD_flag = 0; % Will go to 1 if CTD type & SN recovered
rdata    = [];
adata    = [];
hrrdata  = [];
disp(' ');
disp(['Merging float profiles for ' ,MBARI_ID_str,':'])

line_ct = 1;
for file_ct = 1:r_list
    fn      = file_list{file_ct};
    fp      = fullfile(float_dir, fn);

    % ALIGN ALL CYCLE DATA TO ONE PRESSURE AXIS
    sprof        = build_BSOLO_syn_prof(fp, var_list);
    INFO         = sprof.INFO;
    INFO.CTDtype = '';
    INFO.CTDsn   = '';
    fprintf('%0.0f ',INFO.cast)
    
    if CTD_flag == 0  % Look FOR CTD TYPE AND SN
        CTDtype = INFO.CTDtype;
        CTDsn   = INFO.CTDsn;
        if ~isempty(CTDtype) && ~isempty(CTDsn)
            CTD_flag = 1;
        end
    end
        
    % ONLY WANT 1 GPS POINT
    if size(INFO.gps,1) > 1
        gps = median(INFO.gps,1); % if multiple fixes take median
    else
        gps = INFO.gps(1,:);
    end  
    if gps(2) < 0, gps(2) = gps(2) + 360; end 
    
    % WE ARE NOT PULLING IN TRUE ADJUSTED PRES, TEMP OR PSAL YET SO 
    % ADJUSTED = RAW FOR NOW BUT DO A RANGE CHECK ON P FOR ANY FLIERS
    tbad = sprof.PRES.value < RC.P(1) | sprof.PRES.value > RC.P(2);
    sprof.PRES_QC.value(tbad)          = 4;
    sprof.PRES_ADJUSTED_QC.value(tbad) = 4;
    
    % ********************************************************************
    % BUILD HEADER ARRAY & DATA MATRICES 
    % this could be updated & cleaner but following flow of exisitng code
    % ********************************************************************
    if file_ct == 1
        rhdr = {}; % predim ARGO QC vars list
        ahdr = {};
        
        for i = 1: vars_ct % LOOP THROUGH MASTER LIST & FIND FLOAT VARIABLES
            if isfield(sprof, raw_vars{i,1})
                rhdr = [rhdr, raw_vars{i,1}];
                ahdr = [ahdr, adj_vars{i,1}];
            end
        end
        hdr_size = size(rhdr,2);
        rhdr2 =['Station' 'Matlab SDN' 'Lon [ºE]' 'Lat [ºN]' rhdr];
        ahdr2 =['Station' 'Matlab SDN' 'Lon [ºE]' 'Lat [ºN]' ahdr];
        
        % PREDIM OUTPUT DATA MATRICES
        rdata = ones(r_list*500, size(rhdr2,2)) * NaN;
        adata = rdata;
    end
    
    % BUILD CYCLE DATA MATRICES
    sample_rows = size(sprof.PRES.value,1);
    fill_0      = ones(sample_rows,1) *0; % for array building 
    
    rtmp = ones(sample_rows, hdr_size) * NaN;
    atmp = rtmp;

    
    % Float param list determined from 1st file, but sensor could be turned off
    % later in life! If so the param will no longere exist in the later
    % profile. Always check field existance before gettung data
    % Elsif block can be commented out, mostly for testing - jp 07/18/2022
    for i = 1:hdr_size
        if isfield(sprof, rhdr{i}) 
            rtmp(:,i) = sprof.(rhdr{i}).value;
            atmp(:,i) = sprof.(ahdr{i}).value;
        elseif ~isfield(sprof, rhdr{i}) && isempty(regexp(rhdr{i},'_QC','once'))
            fprintf(['%s sensor appears to be turned of for cycle ', ...
                '%0.0f\n'],rhdr{i},INFO.cast);
        end
    end
    rtmp = flipud(rtmp);
    atmp = flipud(atmp);

    % BUILD MATRIX OF CAST SDN LON and LAT
    tmp =[fill_0+INFO.cast, fill_0+INFO.sdn, fill_0+gps(2), ...
          fill_0+gps(3)]; %gps indexing changed after addition of sdn to gps vector (TM, 10/2020)
    
    rdata(line_ct:line_ct+sample_rows-1,:) = [tmp,rtmp];
    adata(line_ct:line_ct+sample_rows-1,:) = [tmp,atmp];
    line_ct = line_ct+sample_rows;
end
fprintf('\r\n ')
rdata = rdata(1:line_ct-1,:);
adata = adata(1:line_ct-1,:);

% UPDATE HEADERS TO INCLUDE META DATA COLS
rhdr  = rhdr2; % Update with cast sdn lat lon
ahdr  = ahdr2; % Update with cast sdn lat lon

% CHECK QC VALUES & SET QC = NaN to 99 & ASSOCIATED VALUE to 99999
tQC    = ones(line_ct-1, 1) * ~cellfun(@isempty, regexp(rhdr,'QC$','once'));
tD     = circshift(tQC,-1,2);
% Need to keep LAT as NaN if any for interp in argo2ODV_XXXX
tLAT   = ones(line_ct-1, 1) * strncmp(rhdr,'Lat',3); 
tr_nan = isnan(rdata);
ta_nan = isnan(adata);
rdata(tQC & tr_nan) = 99;
rdata(tD & tr_nan & ~tLAT)  = 99999;
adata(tQC & ta_nan) = 99;
adata(tD & ta_nan & ~tLAT)  = 99999;
clear tQQ tD tr_nan ta_nan tLAT

% FILL STRUCTURE AND CLEAN UP
d.INFO.INST_ID_num     = INST_ID_str;
d.INFO.WMO           = WMO_ID;
d.INFO.FloatViz_name = WMO_ID; %3/3/21 NOW USING WMO FOR FLOATVIZ FILE NAMING!!!  :-)
d.INFO.MBARI_ID_str  = MBARI_ID_str;
d.INFO.float_type    = INFO.float_type;
d.INFO.CTDtype       = CTDtype;
d.INFO.CTDsn         = CTDsn;

d.rhdr  = rhdr;
d.rdata = rdata;
d.ahdr  = ahdr;
d.adata = adata;

d.hrrdata = hrrdata;

clearvars -except d

   