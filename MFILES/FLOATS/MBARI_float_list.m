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
% 01/09/2020 - JP, Complete over haul, removed obsolete tasks, condensed &
%              cleaned up
% 06/09/20 JP fixed  bugs in message string if no connection was made to UW

% TESTING
%dirs =[];

% *************************************************************************
% SET DIRS, FILENAMES, & FILTERS
% *************************************************************************
out_name = 'MBARI_float_list'; % save name

if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    
    dirs.cal     = [user_dir,'CAL\']; % save created lists here
    dirs.temp    = 'C:\temp\'; % path to temporary working dir
    dirs.msg     = '\\atlas\ChemWebData\floats\'; % to msg file dirs

elseif ~isstruct(dirs) 
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

dirs.FV_List = '\\sirocco\wwwroot\lobo\Data\FloatVizData\'; % MBRI master list
dirs.UW_html = ['http://runt.ocean.washington.edu/argo/', ...
    'heterographs/rollcall.html']; % UW list for WMO #'s

% FILTERS
float_exp    = '^\d{3}\d+\D+'; % 3#'s, at least one more, then non #'s             
exclude_expr = 'MTY|cor|\d{6}\d+'; % % Ignore these on Fv master list
no_wmo       = '6966|0948|0412'; % WMO's were never assigned

% *************************************************************************
% CHECK NETWORK CONNECTIONS IF NO GOOD EXIT
% *************************************************************************
if ~isdir(dirs.FV_List)
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
% GENERATE UW_ID FROM MBARI NAME
% *************************************************************************
if exist([dirs.FV_List,'FloatVIZConfig.txt'],'file') % check for exisiting file
    fid = fopen([dirs.FV_List,'FloatVIZConfig.txt']);
else
    disp(['Could not find: ', dirs.FV_List,'FloatVIZConfig.txt'])
    list = {};
    return
end

d = textscan(fid,'%s%s','CollectOutput',1,'Delimiter','.');
d = d{1,1}(:,1); % only want MBARI float names, not file extension
fclose(fid);

% LOGICAL TESTS TO SUBSET FLOAt LIST
t1 = ~cellfun(@isempty,regexp(d(:,1),float_exp,'once')); %(desired floats)
t2 =  cellfun(@isempty,regexpi(d(:,1),exclude_expr,'once')); %exluded floats

float_names = d(t1&t2,:);
UW_ID       = regexp(float_names,'^\d+','once', 'match');
r           = size(float_names,1);

% SORT BY MBARI UW_ID
UW_ID_num = str2double(UW_ID);
[~,IX]    = sort(UW_ID_num);

float_names = float_names(IX);
UW_ID       = UW_ID(IX);
WMO_ID      = cell(r,1); % Predim for nextstep
float_type  = cell(r,1); % i.e. APEX or NAVIS

clear ans fid d t1 t2 float_exp exclude_exp FloatViz_names IX UW_ID_num

% *************************************************************************
% GET WMO #'s & APEX (UW) ID's FROM UW ROLL CALL LIST
% Save roll call html page locally and then parse
% *************************************************************************
hdr_chk = 0;
ct      = 0;
[~,status] = urlwrite(dirs.UW_html,[dirs.temp,'UW_rollcall.txt']);

if ~status
    disp(['Could not access: ',dirs.UW_html,'. LIST WAS NOT UPDATED!!'])
    list = {};
    return
else
    fid = fopen([dirs.temp,'UW_rollcall.txt']);
    tline = '';
    
    % GET HEADER LINE & FLOAT LINES
    while ischar(tline) 
        % *** LOOK FOR HEADER ***
        if hdr_chk == 0 &&  ~isempty(regexp(tline,'^\s+Float', 'once'))
            hdr_chk = 1;
            rc_hdr = regexp(tline,'\s+','split');
            if isempty(rc_hdr{1,1})
                rc_hdr = rc_hdr(2:end);
            end
            % using space as splitting delimiter need to modify header a
            % bit
            I1       = find(strcmp(rc_hdr,'Date'));
            I2       = find(strcmp(rc_hdr,'Due-Date'));
            
            rc_hdr    = [rc_hdr(1:I1-1),'Month' 'Day' 'Year', ...
                         rc_hdr(I1+1:I2-1),'Due-Mon' 'Due-Day' 'Due-Yr', ...
                         rc_hdr(I2+1:end)];
            roll_call = cell(4000, size(rc_hdr, 2));
        % *** LOOK FOR FLOAT LINES ***   
        elseif regexp(tline,'index/\d+\.shtml','once','match') % float line
            ct = ct+1;
            tline = regexprep(tline,'\s+<.+index/',''); % remove some html encoding
            tline = regexprep(tline,'.shtml>\d+</a>',''); % remove some more                 
            tline = regexprep(tline,'*','NaN'); % replace
            dcell = regexp(tline,'\s+','split');
            roll_call(ct,1:size(dcell,2)) = dcell;
        end
        tline = fgetl(fid);
    end
end
roll_call = roll_call(1:ct,:);
fclose(fid);
if exist([dirs.temp,'UW_rollcall.txt'], 'file') == 2
    delete([dirs.temp,'UW_rollcall.txt']);
end
% For now just want UW_ID & WMO#. Later might want to utilize more info?
iUW       = find(strcmp(rc_hdr,'Float'));
iWMO      = find(strcmp(rc_hdr,'WmoId'));
roll_call = roll_call(:,[iUW,iWMO]);
rc_hdr    = rc_hdr([iUW,iWMO]); 

clear fid iUW iWMO I1 I2 str_fmt dcell tline hdr_chk status ct
    
% STEP THROUGH MBARI FLOAT NAMES & FILL IN WMO #'s
for i = 1:r
    fname = float_names{i};
    tUW = strcmp(roll_call(:,1), UW_ID{i});
    
    if strcmp(fname,'8514HAWAII') % SPECIAL CASE - deployed twice & dif WMO's
        WMO_ID{i} = '5904172';
    elseif strcmp(fname,'8501CALCURRENT') % SPECIAL CASE - deployed twice & dif WMO's
        WMO_ID{i} = '5904680';        
    elseif sum(tUW) == 1 % ALL GOOD -ONLY ONE MATCH
        WMO_ID{i} = roll_call{tUW,2};
    elseif sum(tUW) == 0 % NO MATCH, NOT NOTIFIED YET or WMO NEVER ASSIGNED
        WMO_ID{i} = ['NO_WMO_',float_names{i}];
        if regexp(UW_ID{i},no_wmo,'once')
            disp(['NO WMO_ID for ',float_names{i},'. WMO was never assigned'])
        else
            disp(['WMO_ID for ',float_names{i},' is not assigned yet'])
        end
    else % MULTIPLE MATCHES - SOMETHING NEEDS TO BE FIXED
        disp(['Multiple WMO IDs found for ',UW_ID{i},' - check records'])
    end
end
  
% *************************************************************************
% DETERMINE FLOAT TYPE AND ADD TO LIST BY CHECKING MSG DIR NAME
% *************************************************************************
for i = 1 : size(float_names,1)
    flt_name = float_names{i};

    if exist([dirs.msg,'f',UW_ID{i},'\'],'dir') == 7
        %msg_dir = [dirs.msg,'f',float_names{i},'\'];
%         FT = 'APEX';
%         float_type{i} = FT;
        float_type{i} = 'APEX';
    elseif exist([dirs.msg,'n',UW_ID{i},'\'],'dir') ==7
        %msg_dir = [dirs.msg,'n',UW_ID{i},'\'];
%         FT = 'NAVIS';
%         float_type{i} = FT;
        float_type{i} = 'NAVIS';
    elseif strcmp(flt_name, '0412HAWAII')% SPECIAL CASE
            disp([flt_name,' is a SPECIAL CASE: non standard float ', ...
                'directory. Float type set to NAVIS ']);
            float_type{i} = 'NAVIS';
    elseif strcmp(flt_name, '6966HAWAII') % SPECIAL CASE
            disp([flt_name,' is a SPECIAL CASE: non standard float ', ...
                'directory. Float type set to APEX ']);
            float_type{i} = 'APEX';        
%             disp(['SPECIAL CASE: non standard float directory for : ', ...
%                 float_names{i}]);
%             disp('Float type set to APEX in function')
%             float_type{i} = 'APEX';
    else
        disp(['could not find  ChemWebData float directory for : ',flt_name]);
        float_type{i} = 'UNKNOWN';
        continue
    end
    clear flt_name msg_files file_path d type_test tf nrows j
end

% *************************************************************************
% BUILD FLOAT LIST CELL ARRAY
% *************************************************************************
list      = [float_names, UW_ID, WMO_ID, float_type];
clearvars -except list dirs out_name

% *************************************************************************
% SAVE AS *.MAT & *.TXT
% *************************************************************************
save([dirs.cal,out_name,'.mat'],'list');

r = size(list,1);
fid = fopen([dirs.cal,out_name,'.txt'], 'w');
fprintf(fid,'%s\t%s\t%s\t%s\r\n','MBARI float name', 'UW ID#', ...
    'WMO ID #', 'float_type');
for i = 1:r
    fprintf(fid,'%s\t%s\t%s\t%s\r\n',list{i,1}, list{i,2}, list{i,3}, ...
        list{i,4});
end
fclose(fid);


