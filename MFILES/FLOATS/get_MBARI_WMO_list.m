function[float_names, UW_ID, WMO_ID] = get_MBARI_WMO_list(dirs)
% ************************************************************************
% PURPOSE: 
%   This function tries to build matching lists of UW_ID, WMO_IB, and
%   MBARI Float ID names by parsing MBARI and UW web page HTML code
%
% USAGE:
%   [float_names, UW_ID, WMO_ID] = get_MBARI_WMO_list(dirs)
%
% INPUTS:
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
%   Three cell arrays with ID names:
%      float_names = MBARI float names in FloatViz
%      UW_ID       = Univ. of Wash. float ID#
%      WMO_ID      = WMO ID# used by ARGO and others
%
%   Data are then saved as a "mat" and "txt" file

fclose all;
% *************************************************************************
%      SET DATA SOURCE PATHS AND REGULAR EXPRESSIONS
%      !!! MAY NEED TO CHANGE THESE FOR EACH USER !!! 
% ************************************************************************* 
fname     = 'float_WMO_list'; % Output file name

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    dirs.mat       = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\FLOATS\';
    dirs.cal       = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\Cal_files\';
    dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    dirs.NO3config = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\Cal_files\';
    dirs.temp      = 'C:\temp\';
    dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    dirs.QCadj      = '\\sirocco\wwwroot\lobo\Data\FloatVizData\QC\';
    dirs.FVlocal   = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\FloatViz\';
elseif ~isstruct(dirs) 
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

% MBARI MASTER FLOAT NAME LIST
FloatViz_names = ['\\sirocco\wwwroot\lobo\Data\FloatVizData\',...
                  'FloatVIZConfig.txt'];

% MOST WMO_ID #'s CAN BE PULLED FROM HERE
FloatViz_html = '\\Sirocco\wwwroot\chemsensor\FloatList.html'; 

% WMO ID #'s FROM SCRIPPS FO SOCCOM
% http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html

% Univ of Wash list (UW float #'s , WMO_ID #'s here
UW_html = ['http://runt.ocean.washington.edu/argo/heterographs/', ...
           'rollcall.html'];

float_exp    = '^\d{3}\d+\D+'; % 3#'s, at least one more, then non #'s       
exclude_expr = 'MTY|cor|\d{6}\d+'; % test deploys, data experiments, old HI

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

% % SORT BY UW ID
% UW_ID_num = str2double(UW_ID);
% if size(UW_ID_num,1) == size(float_names,1) % no empty cells lost in conversion
%     [~,IX] = sort(UW_ID_num);
% else
%     disp('UW #s not extracted correctly')
%     return
% end
% float_names = float_names(IX);
% UW_ID       = UW_ID(IX);
% WMO_ID      = cell(r,1); % Predim for nextstep

% SORT BY MBARI FLOAT NAME
%UW_ID_num = str2double(UW_ID);
[~,IX] = sort(float_names);

float_names = float_names(IX);
UW_ID       = UW_ID(IX);
WMO_ID      = cell(r,1); % Predim for nextstep

clear d t1 t2 float_exp exclude_exp FloatViz_names IX

% *************************************************************************
% GET ALL AVAILABLE WMO_ID #'s from FloatViz html page
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

% disp('NO WMO # found from MBARI list for:');
% str = '';
% for i = 1:r
%     str = [str, no_WMO{i,2},'  '];
% end
% disp(str)

% *************************************************************************
% IF MISSING WMO #'s TRY UW LIST
% *************************************************************************
[~,status] = urlwrite(UW_html,[dirs.temp,'UW_rollcall.txt']);
if ~status
    disp(['Could not access: ',UW_html])
    return
end

fid = fopen([dirs.temp,'UW_rollcall.txt']);
WMO_found = ones(r,1)*0;
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

% *************************************************************************
% MAKE UP FAKE WMO # FOR FLOATS WITH NO WMO
% *************************************************************************
for i = 1: size(WMO_ID,1)
    if isempty(WMO_ID{i})
        WMO_ID{i} =['NO_WMO_',UW_ID{i}];
    end
end

% PRINT TO SCREEN WMO ID #'s found at UW
% for i = 1:r
%     if WMO_found(i) == 1
%         disp(['WMO_ID for ', no_WMO{i,2}, ...
%             ' found on UW list: ', WMO_ID{no_WMO{i,1}}]);
%     end
% end
% PRINT TO SCREEN WMO ID #'s STLL NOT FOUND      
for i = 1:r
    if WMO_found(i) == 0
        disp(['NO WMO_ID for ', no_WMO{i,2}])
    end
end
    
fclose(fid);
delete([dirs.temp,'UW_rollcall.txt'])

%figure out what to do with 8514 soOcn and Hawaii from UW and MBARI list

% *************************************************************************
% SAVE THE INFO
% *************************************************************************
%float_list = [float_names, UW_ID, WMO_ID];
save([dirs.cal,fname,'.mat'],'float_names', 'UW_ID', 'WMO_ID');

r = size(float_names,1);
fid = fopen([dirs.cal,fname,'.txt'], 'w');
fprintf(fid,'%s\t%s\t%s\t\r\n','MBARI float name', 'UW ID#','WMO ID #');
for i = 1:r
    fprintf(fid,'%s\t%s\t%s\t\r\n',float_names{i}, UW_ID{i}, WMO_ID{i});
end
fclose(fid);

clear FloatViz_names FloatViz_html UW_html float_exp exclude_expr ans
clear status no_WMO fid i tline tmpID str r fname dirs UW_ID_num tf


