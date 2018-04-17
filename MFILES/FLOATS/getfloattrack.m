function d = getfloattrack(fname,varargin)
%
%
% INPUT:
%    fname    = float name with or without *.txt extension
%    varargin = QC_flag if present,
%               1         = use data flags to set flagged data = NaN
%               0 or none = no data adjustment
%
% OUTPUT:      returns a structure of float data
%    d.fname = float name
%    d.info  = various processing messages if triggered
%    d.raw   = matrix of float data [date, Sta, lat, lon, depth T S O2 NO3]   
%    d.track = a point per position [Matlab SDN, PROFILE #, Lon, Lat, Depth]
%    d.bnds  = time, lat & lon bounds for the float
%
% EXAPLES:
%   d = getfloattrack('5143StnP')
%   d = getfloattrack('5143StnP.txt')
%   d = getfloattrack('5143StnP',1)

% modified 5/19/2014 by jp ... deal with differing variable lengths
% modified 6/02/2014 by jp ... if QC file make float dir = QC dir
% modified 5/22/2015 by jp ... add pH and create a header cell array
% modified 12/2/2015 by jp ... add CHL and BSCAT, fix so *. or *.txt work
% modified 12/8/2015 by jp ... added a data QC_flag input

if isempty(varargin)
    QC_flag = 0;
    disp('QC FLAGS NOT USED')
elseif length(varargin) == 1 && varargin{1} == 0
    QC_flag = 0;
    disp('QC FLAGS NOT USED')    
elseif length(varargin) == 1 && varargin{1} == 1
    QC_flag = 1;
    disp('Flagged data will be set to NaN')
else
    disp('Too many input variables')
end

% TEST
%fname = '5143StnP';
%fname = '6963HOTpH';
%fname = '6975Bermuda';
%fname = '7593Hawaii'
%fname = '0276NoAtlantic';
%fname = '8514Hawaii'
%fname = '7672Hawaii'
%fname = '9092SoOcn'
%fname ='8514HAWAII.TXT'
% ***********************************************************************
% SET PATHS & FORMATS
% ***********************************************************************
%fname = regexprep(fname,'.txt','','ignorecase'); % Remove extension if present

data_source = 'internet';
data_source = 'network';

float_dir = '\\SIROCCO\wwwroot\lobo\Data\FloatVizData\';
% IF OPENING QC FILE MAKE FLOAT DIR QC FLOAT DIR
if regexp(fname,'QC$|QC\.')
    float_dir =[float_dir,'\QC\'];
end

float_url = 'http://www.mbari.org/lobo/Data/FloatVizData/';
local_dir   = 'C:\Users\jplant\Documents\MATLAB\TEMP\DATA';

% only want station, date, lat, lon, depth T S O2 NO3
% (*Cruise*) (Station) (*Type*) (mon/day/yr) (hh:mm) (Lon [°E])	(Lat [°N])
% (*Bot. Depth [m]*) .... will finish format string  once header is read

%var_exp      = 'Depth|Temp|Sal|Oxygen[|Nit|pHin'; % variables to keep

% variables to keep
var_exp      = 'Depth|Temp|Sal|Oxygen[|Nit|Chlor|BackSc|CDOM|pHin'; 

float_format = '%*s %f %*s %s %s %f %f %*f';
d_format     = 'mm/dd/yyyy HH:MM';

% ***********************************************************************
% EXTRACT TIME, LAT & LON FROM FLOAT DATA
% ***********************************************************************

% ***********************************************************************
% COPY FLOAT DATA TO LOCAL COMPUTER USING NETWORK OR INTERNET
% ***********************************************************************
if ~isempty(regexpi(fname,'\.txt','once')) % remove .txt if it exists
    fname = regexprep(fname,'\.txt','','ignorecase');
end

switch data_source
    case 'network' % GET TEXT FILE FROM SIROCCO AND STORE LOCALY
        s1 = [float_dir,fname,'.txt']; % build target string
        s2 = [local_dir,fname,'.txt'];   % build destination string
        copyfile(s1,s2);
        
    case 'internet'
        s1 = [float_url,fname,'.txt']; % build target string
        s2 = [local_dir,fname,'.txt'];   % build destination string
        
        
        [~,url_chk] = urlread(s1); % See if file exisit on the web
        if url_chk == 1
            f = urlwrite(s1,s2);
            disp(' ');
            disp(['Data for float ',fname,' retrieved from:  ',...
                float_url]);
            disp(['Saved as  ',f]);
        else
            disp('No file found!')
            return
        end
end

% ***********************************************************************
%    FINISH BUILDING FORMAT STRING BASED ON VARIABLE LIST FROM HEADER
%          STEP THROUGH COMMENT LINES TO GET TO HEADER LINE
%            DESIRED VARIABLES DETERMINED FROM "var_exp"
% ***********************************************************************
hdr_test = 1;
fid = fopen(s2);
while hdr_test == 1 % loop through until header line is found
   tline    = fgetl(fid);
   hdr_test = strncmp(tline,'//',2); % loop stops when no match found
end
float_vars  = textscan(tline,'%s','Delimiter','\t'); % get vars from header
float_vars  = float_vars{1,1}; % r x 1 cell array

saved_vars(1:4,1) = {'GMT', 'Station', 'Lon [°E]', 'Lat [°N]'};
ct = length(saved_vars);

ct1 = ct;
for i = 9:length(float_vars) % Finish building format string
    float_format = [float_format,' %f'];  % all remaining variables
    ct1 = ct1+1;
    saved_vars(ct1,1) = float_vars(i,1);
end
clear ct1

% ***********************************************************************
% READ TAB DELIMITED DATA INTO CELL ARRAY 
% ***********************************************************************
data = textscan(fid,float_format,'Delimiter','\t','CollectOutput',1,...
    'CommentStyle','//');
fclose(fid);
delete(s2)

%  Combine date & time strings but need to add space between the two,
%  use regexp & tokens
d_tmp = strcat(data{1,2}(:,1),regexprep(data{1,2}(:,2), '(\d+:\d+)',' $1'));
tmp   = [data{1,1}, data{1,3}]; % combine numbers

% ***********************************************************************
% CONDENSE DATA & QC
% ***********************************************************************

% TEXTSCAN SEEMS TO BE ADDING AN EXTRA ROW OF NaN's IN THE NUMBER ARRAYS SO REMOVE
t1 = isnan(tmp(:,1)); % Look in profile % col for NaN - should be at the end
tmp    = tmp(~t1,:); % Profile #, Lon, Lat, Depth
d_tmp  = d_tmp(~t1,:);

% REMOVE DUPLICATE SAMPLES IF THEY EXIST
 % Create a count & a unique identifier then look for duplicates
 % unique ID = cast # *10,000 + depth
 [r,~]  = size(tmp);
 tmp_srt  = [(1:r)',tmp(:,1)*1e4+ tmp(:,4),tmp(:,6)]; % ct ID S
 [B,IX] = sort(tmp_srt(:,2));
 tmp_srt = tmp_srt(IX,:); % sorted by unique ID
 df_ID  = [1;diff(tmp_srt(:,2))]; % add 1 to keep dim
 df_S  = [1;diff(tmp_srt(:,3))]; % add 1 to keep dim S
 tmp_srt    = [tmp_srt, df_ID, df_S];

 [B,IX] = sort(tmp_srt(:,1)); 
 tmp_srt = tmp_srt(IX,:); % sorted back by count
 t1 = tmp_srt(:,4) == 0 & tmp_srt(:,5) == 0; % diff in ID and Sal = 0
 out = tmp(t1,:);

 if ~isempty(out)
     tmp =tmp(~t1,:);
     d_tmp  = d_tmp(~t1,:);
     disp(['Duplicate data removed from cast #: ', ...
         sprintf('%3.0f ',out(:,1))]);
 end

% ONLY WANT ONE LINE PER PROFILE
[~,ind,~] = unique(tmp(:,1),'first'); % get index for start (bottom) of profile
tmp1       = tmp(ind,1:4); % Profile #, Lon, Lat, Depth
d_tmp1     = d_tmp(ind,:);

% CHECK THAT PROFILES HAVE DATA
t_size = ones(size(tmp1(:,1))); % predimension
for i = 1:length(tmp1(:,1))
    p_size = tmp(:,1) == tmp1(i,1); % # rows of raw data for profile
    if sum(p_size) < 4 % Not many samples
        t_size(i) = 0; % set flag to zero
        s3 = ['Removing profile ',num2str(tmp1(i,1)),' from ' fname,...
            '  -- only ',num2str(sum(p_size)),' sample(s)'];
        disp(s3);
        d.info(3) = {s3};
        tmp = tmp(~p_size,:);
        d_tmp = d_tmp(~p_size,:);
    end
end
t_size = logical(t_size);
tmp1   = tmp1(t_size,:);
d_tmp1 = d_tmp1(t_size) ;
        
% DO A QUICK CHECK FOR MISSING PROFILES
exp_profiles = min(tmp1(:,1)) : max(tmp1(:,1));
t1 = setdiff(exp_profiles,tmp1(:,1));
if ~isempty(t1)
    s4 = [fname,' is missing profile(s): ',num2str(t1)];
    disp(s4);
    d.info(2) = {s4};
end

% CHECK FOR MISSING POSITION - NO POSITION, NOT USABLE DATA !!!
t1    = isnan(tmp1(:,2)); % look for NaN's in  track longitude (no position)
t2    = isnan(tmp(:,2)) |isnan(tmp(:,5)); % look for NaN's in  all longitude & temp
if sum(t1) > 0
    s5 =[fname,' is missing position fix for profile(s): ',num2str(tmp1(t1,1)')];
    disp(s5)
    d.info(3) = {s5};
end

tmp1   = tmp1(~t1,:); % only want data with position fix for track data
d_tmp1 = d_tmp1(~t1,:);
tmp    = tmp(~t2,:); % only want data with position fix all data
d_tmp  = d_tmp(~t2,:);


% Matlab SDN, PROFILE #, Lon, Lat, Depth
dd     = [datenum(d_tmp1,d_format),tmp1]; % Convert & combine for track
d_all  = [datenum(d_tmp,d_format),tmp]; % Convert & combine for all data

% ***********************************************************************
% IF QC_FLAG = 1, SET FLAGGED DATA to NAN AND THEN REMOVE QC COLS
% OTHERWISE JUST REMOVE QC COLS
% ***********************************************************************
QC_cols = strcmp('QF',saved_vars);
QC_ind  = find(QC_cols == 1);
if QC_flag == 1
    tmp_data = d_all(:,QC_ind-1); % the data
    t1 = d_all(:,QC_ind) > 0; % find flagged data
    tmp_data(t1) = NaN; % set bad to NaN
    d_all(:,QC_ind-1) = tmp_data; % put adjusted data back in data set
end

% Now clear data not matching var exp
clear_ind = [];
for i = ct+1:length(saved_vars) 
    if isempty(regexpi(saved_vars{i,1},var_exp))
        clear_ind =[clear_ind,i];
    end
end
saved_vars(clear_ind) = [];
d_all(:,clear_ind)   = [];

% GET TIME, LON, LAT BOUNDS & DO A QUICK CHECK FOR MISSING PROFILES
time_bnds = [min(dd(:,1)) max(dd(:,1))];
lon_bnds  = [min(dd(:,3)) max(dd(:,3))];
lat_bnds  = [min(dd(:,4)) max(dd(:,4))];

% ***********************************************************************
% ASSIGN OUTPUT & CLEAN UP
% ***********************************************************************
d.fname      = fname;
d.raw        = d_all;
d.track      = dd;
d.bnds.t     = time_bnds;
d.bnds.lon   = lon_bnds;
d.bnds.lat   = lat_bnds;
d.hdr_vars   = float_vars;
d.saved_vars = saved_vars;
disp(['Float track extracted for ',fname]);
%clearvars -except d




 
 
