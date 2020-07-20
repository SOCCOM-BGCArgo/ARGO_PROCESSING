function d = merge_ARGO_mat(MBARI_ID_str, dirs)
% ************************************************************************
% PURPOSE: 
%    This function merges all *.mat profile files for a given float.
%    This merged dataset is used to build the the ODV compatible text files
%    used by FloatViz and SOCOmViz
%
% USAGE:
%	d = merge_ARGO_mat(UW_ID_str, dirs)
%
% INPUTS:
%   MBARI_ID_str  = MBARI Float ID as a string
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
%   dirs.QCadj      = path to QC adjustment list for all floats
%   dirs.FVlocal   = path to Floatviz file made by matlab go here
%
% OUTPUTS:
%   A stucture of merged float data
%     d.INFO
%       .UD_ID_num
%       .WMO
%       .FloatViz_name
%       .CTDtype
%       .CTDsn
%
%    d.rhdr    = raw data header
%    d.rdata   = raw data
%    d.ahdr    = adjusted data header
%    d.adata   = adjusted data
%    d.hrrdata = high resoltion raw data
%
% ************************************************************************
% CHANGE HISTORY
% 4/21/2017 - added code to include S & T ADJUSTED & QC variables for LR
%    and HR data
% 02/11/2020 JP - Added PRES_QC, PRES_ADJUSTED & PRES_ADJUSTED_QC to the
%    parameter search list

% ************************************************************************
%TESTING
% MBARI_ID_str = '8501SOOCN';
%MBARI_ID_str = '0509SOOCN';
% MBARI_ID_str = '0412HAWAII';

% MBARI_ID_str = '9018SOOCN';
% MBARI_ID_str = '0568SOOCN';
% MBARI_ID_str = '0569SOOCN';
%MBARI_ID_str = '7593Hawaii';
% MBARI_ID_str = '0949STNP';
% MBARI_ID_str = '12888SOOCN';
% MBARI_ID_str = '18081SOOCN';
% MBARI_ID_str = '0566SOOCN';
% dirs =[];

% *************************************************************************
%  SET DATA DIRS AND PATHS
% *************************************************************************
RC.P    = [0 12000]; % from argo parameter list
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

%WMO_path  = [dirs.cal, 'float_WMO_list.mat']; % WMO file path
WMO_path  = [dirs.cal, 'MBARI_float_list.mat']; % WMO file path
% *************************************************************************
%  GET MBARI FLOAT NAME, FIND WMO DIR ASSOCIATED WITH FLOAT NAME,
%  AND GET FILE NAMES FOR CAST
% *************************************************************************
if exist(WMO_path, 'file') % load list variables
    load(WMO_path);
else
    disp('BUILDING FLOAT ID LIST ...')
    % build float name, UW_ID, WMO_ID lists
    %[float_names, UW_ID, WMO_ID] = get_MBARI_WMO_list(dirs);
    list = MBARI_float_list(dirs); %[float_name, UW_ID, WMO_ID, float type]
end
float_names = list(:,1);
UW_ID       = list(:,2);
WMO_ID      = list(:,3);
clear list

ind = find(strcmpi(MBARI_ID_str, float_names) == 1);
%pause
if isempty(ind)
    disp(['No MBARI_ID match found in lookup table for ',MBARI_ID_str])
    disp('NO MATCH FOUND - EXITING')
    d = []; % Return this if no data
    return
end

WMO       = WMO_ID{ind};
UW_ID_str = UW_ID{ind};

disp(['WMO# for ',MBARI_ID_str, ' = ', WMO]);
disp(['UW ID for ',MBARI_ID_str, ' = ', UW_ID_str]);


float_dir = [dirs.mat,WMO,'\'];
if exist(float_dir, 'dir') == 7 % list * .mat files
    file_list = ls([float_dir,WMO,'*.mat']);
    r_list    = size(file_list,1);
    if r_list == 0 %directory exists but no .mat files exist yet.
        disp(['WMO directory found for ', MBARI_ID_str, ' (', WMO, '), but no .mat files exist yet.']);
        d = []; % Return this if no data
        return
    end
else
    disp(['No WMO directory found for ', MBARI_ID_str, ' (', WMO, ')']);
    d = []; % Return this if no data
    return
end

% *************************************************************************
%                   LIST ALL POSSIBLE ARGO VARIABLES
%            THIS CAN BE UPDATED AS NEW VARIABLES COME ON LINE
% *************************************************************************
raw_vars(1,1)  = {'PRES'}; 
raw_vars(2,1)  = {'PRES_QC'}; 
raw_vars(3,1)  = {'TEMP'};
raw_vars(4,1)  = {'TEMP_QC'};
raw_vars(5,1)  = {'PSAL'}; 
raw_vars(6,1)  = {'PSAL_QC'};
raw_vars(7,1)  = {'DOXY'};
raw_vars(8,1)  = {'DOXY_QC'};
raw_vars(9,1)  = {'NITRATE'}; 
raw_vars(10,1)  = {'NITRATE_QC'}; 
raw_vars(11,1)  = {'CHLA'};  
raw_vars(12,1)  = {'CHLA_QC'};
raw_vars(13,1) = {'BBP700'};
raw_vars(14,1) = {'BBP700_QC'};
raw_vars(15,1) = {'CDOM'};   
raw_vars(16,1) = {'CDOM_QC'};  
raw_vars(17,1) = {'BBP532'};
raw_vars(18,1) = {'BBP532_QC'};
raw_vars(19,1) = {'PH_IN_SITU_TOTAL'}; 
raw_vars(20,1) = {'PH_IN_SITU_TOTAL_QC'}; 

adj_vars(1,1)  = {'PRES_ADJUSTED'}; 
adj_vars(2,1)  = {'PRES_ADJUSTED_QC'}; 
adj_vars(3,1)  = {'TEMP_ADJUSTED'};
adj_vars(4,1)  = {'TEMP_ADJUSTED_QC'};
adj_vars(5,1)  = {'PSAL_ADJUSTED'}; 
adj_vars(6,1)  = {'PSAL_ADJUSTED_QC'};
adj_vars(7,1)  = {'DOXY_ADJUSTED'}; % ADJUSTED VARIABLES
adj_vars(8,1)  = {'DOXY_ADJUSTED_QC'};
adj_vars(9,1)  = {'NITRATE_ADJUSTED'};
adj_vars(10,1)  = {'NITRATE_ADJUSTED_QC'};
adj_vars(11,1)  = {'CHLA_ADJUSTED'};
adj_vars(12,1)  = {'CHLA_ADJUSTED_QC'};
adj_vars(13,1) = {'BBP700_ADJUSTED'};
adj_vars(14,1) = {'BBP700_ADJUSTED_QC'};
adj_vars(15,1) = {'CDOM_ADJUSTED'};
adj_vars(16,1) = {'CDOM_ADJUSTED_QC'};
adj_vars(17,1) = {'BBP532_ADJUSTED'};
adj_vars(18,1) = {'BBP532_ADJUSTED_QC'};
adj_vars(19,1) = {'PH_IN_SITU_TOTAL_ADJUSTED'};
adj_vars(20,1) = {'PH_IN_SITU_TOTAL_ADJUSTED_QC'};  

raw_vars_ct = size(raw_vars,1);
adj_vars_ct = size(adj_vars,1);
if adj_vars_ct ~= raw_vars_ct
    dip('RAW and ADJ variable list should be the same size')
    return
end

% ************************************************************************
%  PROCESS *.MAT FILES
% ************************************************************************
CTD_flag = 0; % Will go to 1 if CTD type & SN recovered
rdata = [];
adata = [];
hrrdata = [];
disp(' ');
disp(['Merging float profiles for ' ,MBARI_ID_str,':'])
for file_ct = 1 : r_list
    float_file = strtrim(file_list(file_ct,:));
    load([float_dir,float_file])
    fprintf('%0.0f ',INFO.cast)
    
    if CTD_flag == 0 % Look FOR CTD TYPE AND SN
        CTDtype = INFO.CTDtype;
        CTDsn   = INFO.CTDsn;
        if ~isempty(CTDtype) && ~isempty(CTDsn)
            CTD_flag = 1;
        end
    end
        
    % ONLY WANT 1 GPS POINT
    if size(INFO.gps,1) > 1
        gps = median(INFO.gps); % if multiple fixes take median
    else
        gps = INFO.gps(1,:);
    end
        
    if gps(1) < 0, gps(1) = gps(1) + 360; end
    
    % RECENT ADDITION OF PRES, PRES_QC, PRES_ADJUSTED & PRES_ADJUSTED_QC
    % TO FLOAT PROCESSING (3/3/2020 JP). OLDER MAT FILES MAY NOT HAVE THEM
    % SO ADD TO LR & HR STRUCTURES IF NEED BE
    if exist('LR','var') == 1 && ~isfield(LR,'PRES_ADJUSTED')
        LR.PRES_QC = ones(size(LR.PRES));
        LR.PRES_QC(LR.PRES < RC.P(1) | LR.PRES > RC.P(2)) = 4;
        LR.PRES_QC(LR.PRES == 99999) = 99;
        LR.PRES_ADJUSTED = LR.PRES;
        LR.PRES_ADJUSTED_QC = LR.PRES_QC;
    end
    
    if exist('HR','var') == 1 && ~isfield(HR,'PRES_ADJUSTED')
        HR.PRES_QC = ones(size(HR.PRES));
        HR.PRES_QC(HR.PRES < RC.P(1) | HR.PRES > RC.P(2)) = 4;
        HR.PRES_QC(HR.PRES == 99999) = 99;
        HR.PRES_ADJUSTED = HR.PRES;
        HR.PRES_ADJUSTED_QC = HR.PRES_QC;
    end
    
    % ********************************************************************
    % BUILD HEADER ARRAY
    % ********************************************************************
    if file_ct == 1 
        rhdr = {}; % predim ARGO QC vars list
        ahdr = {};
        
        
        for i = 1: raw_vars_ct % LOOP THROUGH MASTER LIST & FIND FLOAT VARIABLES
            if isfield(LR, raw_vars{i,1})
                rhdr = [rhdr, raw_vars{i,1}];
                ahdr = [ahdr, adj_vars{i,1}];
            end
        end
        hdr_size = size(rhdr,2);
        rhdr2 =['Station' 'Matlab SDN' 'Lon [ºE]' 'Lat [ºN]' rhdr];
        ahdr2 =['Station' 'Matlab SDN' 'Lon [ºE]' 'Lat [ºN]' ahdr];
    end
   

    % ********************************************************************
    % IF NAVIS FLOAT MAKE HR NITRATE, MERGE DEEP LR WITH HR, & THEN SET 
    % LR fields=HR fields for ODV NITRATE IS REALLY THE ONLY LR FIELD
    % ********************************************************************
    % GET HR INDICES FOR LR MISSING VALUES
    if strcmp(INFO.float_type,'NAVIS') && isfield(LR,'NITRATE') && ...
            ~isempty(HR.PRES) 
        HR_fill0 = ones(size(HR.PRES))* 0;
        HR.NITRATE    = HR_fill0 + 99999;
        HR.NITRATE_QC = HR_fill0 + 99;
        HR.NITRATE_ADJUSTED    = HR_fill0 + 99999;
        HR.NITRATE_ADJUSTED_QC = HR_fill0 + 99;

        for i = 1: size(LR.PRES,1)
            if LR.DOXY(i) == 99999
                abs_diff = abs(HR.PRES - LR.PRES(i));
                min_diff = min(abs_diff);
                ind = find(abs_diff == min_diff); % 2 m max diff
                if min_diff < 2 % A MATCH!
                    HR.NITRATE(ind(1))    = LR.NITRATE(i);
                    HR.NITRATE_QC(ind(1)) = LR.NITRATE_QC(i);
                    HR.NITRATE_ADJUSTED(ind(1))    = LR.NITRATE_ADJUSTED(i);
                    HR.NITRATE_ADJUSTED_QC(ind(1)) = LR.NITRATE_ADJUSTED_QC(i);                    
                end
            end
        end
    end

    if strcmp(INFO.float_type,'NAVIS') && ~isempty(HR.PRES)
        
        % NOW MERGE DEEP LR WITH HR & THEN SET LR DATA = HR DATA FOR ODV
        % FILE PRINT OUT
        t1 = LR.DOXY ~= 99999; % discrete sample flags
        [~,LRHR_ind] = sort([LR.PRES(t1); HR.PRES],'ascend'); % sort ind's
        
        for i = 1: hdr_size
            if isfield(LR,rhdr{i}) && isfield(HR,rhdr{i})
                LR.(rhdr{i}) = [LR.(rhdr{i})(t1); HR.(rhdr{i})]; % deep LR + HR
                LR.(rhdr{i}) = LR.(rhdr{i})(LRHR_ind); % sorted
            end

           
            % DONT DO P TWICE  ADJUSTED & RAW HEADER NAME ARE EQUAL!
            % UPDATED 05/02/17 due to addind PSAl_ADJUSTED & TEMP_ADJUSTED
            % variables
            % Updated 2/26/20 PRES_ADJUSTED added as variable
            %if isempty(regexp(ahdr{i},'^PRES|^TEMP|^PSAL','once')) 
            %if isempty(regexp(ahdr{i},'^PRES','once')) % ADJUSTED DATA    
                if isfield(LR,ahdr{i}) && isfield(HR,ahdr{i})
                    LR.(ahdr{i}) = [LR.(ahdr{i})(t1); HR.(ahdr{i})]; % deep LR + HR
                    LR.(ahdr{i}) = LR.(ahdr{i})(LRHR_ind); % sorted
                end
            %end
        end
    end
 
    
    % ********************************************************************
    % MERGE LR APEX OR NAVIS FLOAT DATA
    
    sample_rows = size(LR.PRES,1);
    fill_0      = ones(sample_rows,1) *0; % for array building 
    
    rtmp = ones(sample_rows, hdr_size) * NaN;
    atmp = rtmp;
    
    ind = sample_rows:-1:1;
    for i = 1:hdr_size
        if isfield(LR,rhdr{i})
            rtmp(ind,i) = LR.(rhdr{i});
        end
        if isfield(LR,ahdr{i})
            atmp(ind,i) = LR.(ahdr{i});
        end
    end
    
    
    % BUILD MATRIX OF CAST SDN LON and LAT
    
    tmp =[fill_0+INFO.cast, fill_0+INFO.sdn, fill_0+gps(1), ...
          fill_0+gps(2)];
      
    rdata = [rdata; tmp,rtmp];
    adata = [adata; tmp,atmp];
    
   

    % ********************************************************************
    % NEXT MERGE HR APEX FLOAT DATA FOR HR+LR TXT FILES LATER (ALL DATA)
    % FOR APEX THIS IS ONLY PTS SO ADJ = RAW
    if strcmp(INFO.float_type,'APEX') && ~isempty(HR.PRES)
        hr_sample_rows = size(HR.PRES,1);
        hr_fill_0 = ones(hr_sample_rows,1) *0; % for array building
        
        hrrtmp = ones(hr_sample_rows, hdr_size) * NaN;
     
        ind = hr_sample_rows:-1:1;
        % There are no QC ajustments to S or T yet, but flags could be
        % altered in either ODV file in sirocco so get max QF value from
        % both
        for i = 1:hdr_size
            if isfield(HR, rhdr{i})
                if strcmp(rhdr{i},'PSAL_QC') % Man QC in raw or adj, get max
                    hrrtmp(ind,i) = max([HR.(rhdr{i}),HR.(ahdr{i})],[],2);
                elseif strcmp(rhdr{i},'TEMP_QC')
                    hrrtmp(ind,i) = max([HR.(rhdr{i}),HR.(ahdr{i})],[],2);
                else
                    hrrtmp(ind,i) = HR.(rhdr{i});
                end
            end            
        end           
        
        t_nan   = sum(~isnan(hrrtmp),1) == 0;
        hrrtmp(:,t_nan) = []; % remove cols after p,t,s
        
        tmp =[hr_fill_0+INFO.cast, hr_fill_0+INFO.sdn, ...
            hr_fill_0+gps(1), hr_fill_0+gps(2)];
        
        hrrdata = [hrrdata; tmp,hrrtmp];  
    end
    
end
fprintf('\r\n ')
rhdr = rhdr2; % Update with cast sdn lat lon
ahdr = ahdr2; % Update with cast sdn lat lon


% FILL STRUCTURE AND CLEAN UP
d.INFO.UW_ID_num     = UW_ID_str;
d.INFO.WMO           = WMO;
d.INFO.FloatViz_name = MBARI_ID_str;
d.INFO.float_type    = INFO.float_type;
d.INFO.CTDtype       = CTDtype;
d.INFO.CTDsn         = CTDsn;

d.rhdr  = rhdr;
d.rdata = rdata;
d.ahdr  = ahdr;
d.adata = adata;

d.hrrdata = hrrdata;

% clearvars -except d
    
    
