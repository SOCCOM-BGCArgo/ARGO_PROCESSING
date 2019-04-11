function list = MBARI_float_list(dirs)
% This function builds a cell array of float ID info:
%
% INPUT:
%   dirs is a structure of directories or []. If empty defaults are used
%
% OUTPUT:
%   An n x 4 cell array [float_name, UW_ID, WMO_ID, float type]
%       float names -   MBARI float ID's
%       UW_ID       -   Univ of Washington float ID #
%       WMO_ID      -   WMO ID #
%       float_type  -   APEX, NAVIS, or UNKNOWN at this point
%
%       A *.mat & *.txt file will also be created in dirs.cal
%
% EXAMPLES:
%       list = MBARI_float_list(dirs);
%       list = MBARI_float_list([]);
%
% CHANGE HISTORY
% 01/09/2017 - if UW site down function won't break now

% *************************************************************************
% SET OUTPUT FILE NAME AND DIRS STUCTURE
% *************************************************************************
fname = 'MBARI_float_list'; % save name
%   dirs.temp      = path to temporary working dir
%   dirs.msg       = path to float message file directories

if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.cal = [user_dir,'CAL\'];
    
    dirs.temp = 'C:\temp\';
    dirs.msg  = '\\atlas\ChemWebData\floats\';

elseif ~isstruct(dirs) 
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

% *************************************************************************
% DEFINE WORKING DIRS & URL ADDRESSES
% *************************************************************************

% MBARI MASTER FLOAT NAME LIST
FloatViz_names = ['\\sirocco\wwwroot\lobo\Data\FloatVizData\',...
                  'FloatVIZConfig.txt'];

% !!!! SIROCCO LINKS MAY NEED TO BE UPDATED SOON - jp 12/29/2016 !!!!             
              
% MOST WMO_ID #'s CAN BE PULLED FROM HERE
FloatViz_html = '\\Sirocco\wwwroot\chemsensor\FloatList.html'; 

% WMO ID #'s FROM SCRIPPS FOR SOCCOM
% http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html

% Univ of Wash list (UW float #'s , WMO_ID #'s here
UW_html = ['http://runt.ocean.washington.edu/argo/heterographs/', ...
           'rollcall.html'];

% 3#'s, at least one more, then non #'s - allow for #'s > 9999 coming soon       
float_exp    = '^\d{3}\d+\D+';        
exclude_expr = 'MTY|cor|\d{6}\d+'; % test deploys, data experiments, old HI

% *************************************************************************
% CHECK NETWORK CONNECTIONS IF NO GOOD EXIT
% *************************************************************************
if ~isdir('\\sirocco\wwwroot\lobo\Data\FloatVizData\')
    disp(['Could not access float list at: \\sirocco\wwwroot\lobo\', ...
          'Data\FloatVizData\'])
    list = {};
    return
end

if ~isdir(dirs.msg)
    disp(['Could not access float *.msg dirs to determine float type ',...
        'at: ', dirs.msg]);
    list = {};
    return
end

% *************************************************************************
% GET MASTER LIST OF MBARI FLOATS FROM FloatVIZConfig.txt
% *************************************************************************
if exist(FloatViz_names,'file') % check for exisiting file
    fid = fopen(FloatViz_names);
else
    disp(['Could not find: ', FloatViz_names])
    return
end

d = textscan(fid,'%s%s','CollectOutput',1,'Delimiter','.');
fclose(fid);

d = d{1,1}; % get rid of outer cell
d = d(:,1); % only want names, not suffix

% logical test - find names with#'s followed by txt
t1 = ~cellfun(@isempty,regexp(d(:,1),float_exp,'once')); %(floats)
t2 =  cellfun(@isempty,regexpi(d(:,1),exclude_expr,'once')); %exlude floats

float_names = d(t1&t2,:);
UW_ID       = regexp(float_names,'^\d+','once', 'match');
r           = size(float_names,1);

% SORT BY MBARI UW_ID
UW_ID_num = str2double(UW_ID);
[~,IX]    = sort(UW_ID_num);

float_names = float_names(IX);
UW_ID       = UW_ID(IX);
WMO_ID      = cell(r,1); % Predim for nextstep

clear d t1 t2 float_exp exclude_exp FloatViz_names IX UW_ID_num

% *************************************************************************
% GET ALL AVAILABLE WMO_ID #'s from FloatViz html page 1st
% *************************************************************************
if exist(FloatViz_html,'file')
    fid = fopen(FloatViz_html);
else
    disp(['Could not find: ', FloatViz_html])
    return
end

tline = fgetl(fid);
ct    = 0;
tmpID = cell(1000,3); % predim for [FLOAT NAME  UW_ID   WMO_ID]
while ischar(tline)
    % GET STR BLOCK WITH FLOAT NAME AND THEN WMO WITHIN ()
    str = regexp(tline,'\d{4}\w+\s\(\d+\)','once','match');
    if ~isempty(str)
        ct = ct+1;
        tmpID{ct,1} = regexp(str,'^\d+\w+','match','once'); % float name
        tmpID{ct,2} = regexp(str,'^\d+','match','once'); % UW
        tmpID{ct,3} = regexp(str,'(?<=\()\d+','match','once'); % WMO
    end
    tline = fgetl(fid);
end
tmpID = tmpID(1:ct,:);
fclose(fid);
clear ct tline str

% *************************************************************************
% TRY AND GET WMO #'S FOR EACH FLOAT NAME FROM MBARI DATA
% *************************************************************************
no_WMO = cell(100,3); % predim [float_name_index  float_name  UW_ID]
ct = 0;
for i = 1:r
    tf = strcmpi(float_names{i}, tmpID(:,1));
    if sum(tf) == 1
        WMO_ID(i) = tmpID(tf,3);
    else
        %disp(['WMO ID # not found for ', float_names{i}])
        ct = ct + 1;
        no_WMO{ct,1} = i; 
        no_WMO{ct,2} = float_names{i};  
        no_WMO{ct,3} = UW_ID{i};  
    end
end

if ct == 0 % ALL WMO ID's FOUND - END SCRIPT
    return
end

no_WMO = no_WMO(1:ct,:);
clear ct i r tf
r = size(no_WMO,1);

% *************************************************************************
% IF MISSING WMO #'s TRY UW LIST
% *************************************************************************
[~,status] = urlwrite(UW_html,[dirs.temp,'UW_rollcall.txt']);
WMO_found = ones(r,1)*0;
if ~status
    disp(['Could not access: ',UW_html])
else
    fid = fopen([dirs.temp,'UW_rollcall.txt']);
    
    for i = 1:r
        tline = '';
        while ischar(tline)
            if regexp(tline,['index/', no_WMO{i,3}],'once') % line w/ WMO ID
                str = textscan(tline,'%*s%*s%*s%*s%*s%s',1);
                str = char(str{1,1});
                if regexp(str,'\d+','once') % make sure a #
                    ind = no_WMO{i,1}; % index for missing WMO on orginal list
                    WMO_ID(ind) = {str};
                    WMO_found(i) = 1;
                    break
                end
            end
            tline = fgetl(fid);
        end
        frewind(fid); % back to top and look again - brute force
    end
    fclose(fid)
end
% *************************************************************************
% MAKE UP FAKE WMO # FOR FLOATS WITH NO WMO
% *************************************************************************
for i = 1: size(WMO_ID,1)
    if isempty(WMO_ID{i})
        WMO_ID{i} =['NO_WMO_',float_names{i}];
    end
end
   
for i = 1:r
    if WMO_found(i) == 0
        disp(['NO WMO_ID for ', no_WMO{i,2}])
    end
end
  
% *************************************************************************
% DETERMINE FLOAT TYPE AND ADD TO LIST
% *************************************************************************
float_type = cell(size(float_names));

for i = 1 : size(float_names,1)
    flt_name = float_names{i};
    % CHECK FOR FLOAT DIRECTORY
    if exist([dirs.msg,'f',UW_ID{i},'\'],'dir');
        msg_dir = [dirs.msg,'f',float_names{i},'\'];
        %disp('APEX float directory detected')
        FT = 'APEX';
        float_type{i} = FT;
    elseif exist([dirs.msg,'n',UW_ID{i},'\'],'dir');
        msg_dir = [dirs.msg,'n',UW_ID{i},'\'];
        FT = 'NAVIS';
        float_type{i} = FT;
    elseif strcmp(flt_name, '0412HAWAII')% SPECIAL CASE
            disp(['SPECIAL CASE: non standard float directory for : ', ...
                float_names{i}]);
            disp('Float type set to NAVIS in function')
            float_type{i} = 'NAVIS';
    elseif strcmp(flt_name, '6966HAWAII') % SPECIAL CASE
            disp(['SPECIAL CASE: non standard float directory for : ', ...
                float_names{i}]);
            disp('Float type set to APEX in function')
            float_type{i} = 'APEX';
    else
        disp(['could not find float directory for : ',float_names{i}]);
        float_type{i} = 'UNKNOWN';
        continue
    end
    clear flt_name msg_files file_path d type_test tf nrows j
end

% *************************************************************************
% BUILD FLOAT LIST CELL ARRAY
% *************************************************************************
list      = [float_names, UW_ID, WMO_ID, float_type];
clearvars -except list dirs fname

% *************************************************************************
% SAVE AS *.MAT & *.TXT
% *************************************************************************
save([dirs.cal,fname,'.mat'],'list');

r = size(list,1);
fid = fopen([dirs.cal,fname,'.txt'], 'w');
fprintf(fid,'%s\t%s\t%s\t%s\r\n','MBARI float name', 'UW ID#', ...
    'WMO ID #', 'float_type');
for i = 1:r
    fprintf(fid,'%s\t%s\t%s\t%s\r\n',list{i,1}, list{i,2}, list{i,3}, ...
        list{i,4});
end
fclose(fid);


