function list = get_msg_list(flt_name,dirs,varargin)
% PURPOSE: 
%   This function returns a structure with the regular and alternate
%   directory paths for the given float as well as msg file names within
%   each directory. The default message file  is '*.msg' but an optional
%   input string can change this. If files exist in the alternate
%   directory, duplicate file sizes are compared and the smaller files are
%   removed from a list. Designed for use with other functions / scripts
%   processing float data
%
% USAGE:
%	list = get_msg_list(flt_name,dirs,'msg')
%   list = get_msg_list(flt_name, [])
%
% INPUTS:
%   float name = a string with at least UW/MBARI float #
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
%   varargin = optional string defining msg type. Ex:
%                   'msg'  APEX message file (DEFAULT)
%                   'isus' ISUS message file
%
% OUTPUTS:
%   list = a structure of the data.
%       list.reg_dir  = main msg file diectory path
%       list.alt_dir  = alternate message file directory path if it exists
%       list.reg_list = main directory * msg file list
%       list.alt_list = alternate directory message file list
%       list.reg_sdn = main directory * msg file time stamp
%       list.alt_sdn = alternate directory message file time stamp
%
% EXAMPLES:
%   jp = get_msg_list('7601StnP',dirs);
%   jp = get_msg_list('7601',[]);
%   jp = get_msg_list('7601StnP','dirs,'msg');
%   jp = get_msg_list('7601StnP',[],'isus');

%
% Compile a list of float *.msg files from reg and alt directories
% use file size or date to choose which ones to keep

% %flt_name = '7601StnP';
% %flt_name ='9254SoOcn';
% flt_name ='6966Hawaii';
% %flt_name ='7593HAWAII';
% flt_name ='8514SoOcn';
% dirs =[];
% varargin =[];

% ************************************************************************
% CREATE REG & ALT DIR PATHS USING FLOAT NAME
% ************************************************************************
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\'];
    
    dirs.mat       = [user_dir,'DATA\FLOATS\'];
    dirs.cal       = [user_dir,'DATA\CAL\'];
    dirs.FVlocal   = [user_dir,'DATA\FLOATVIZ\'];
    dirs.FV        = [user_dir,'DATA\FLOATVIZ\'];
    %dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    
    %dirs.QCadj     = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\';
    dirs.QCadj     = [user_dir,'DATA\CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    %dirs.msg       = 'C:\temp\';
    %dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    
    dirs.log = [user_dir,'DATA\Processing_logs\'];
    dirs.bat = [user_dir,'batchfiles\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

duplicate_dir     = [dirs.msg,'duplicate\']; % A few older floats found here

flt_str  = regexp(flt_name,'^\d{3}\d+(?=\w*)','match'); %#'s chars may follow
if isempty(flt_str)
    flt_str = regexp(flt_name,'^\d{3}\d+','match'); % try just numbers next
end
flt_num  = str2num(flt_str{1,1});

% TEST FOR 'f' OR 'n' DIRS
if isdir([dirs.msg,'f',flt_str{1,1},'\'])      % 'f' dir is APEX UW/MBARI
    reg_dir = [dirs.msg,'f',flt_str{1,1},'\'];
elseif isdir([dirs.msg,'n',flt_str{1,1},'\'])
    reg_dir = [dirs.msg,'n',flt_str{1,1},'\']; % 'n' for NAVIS floats
elseif isdir([duplicate_dir,'f',flt_str{1,1},'\'])
    disp('NON STANDARD DATA DIR!!')
    disp(['DATA FOUND AT:  ',duplicate_dir,'f',flt_str{1,1},'\'])
    reg_dir = [duplicate_dir,'f',flt_str{1,1},'\']; % SOME OLD DATA HERE
else
    disp(['Could not find msg file directory for: ',flt_name])
    list = [];
    return
end
alt_dir = regexprep(reg_dir,'floats', 'floats\\alternate'); % need 2 \\

% ************************************************************************
% BUILD MESSSAGE FILE LISTS
% NEED TO LOOK AT REG AND ALT DIRS - CHOOSE LARGER FILE
% ************************************************************************

if isempty(varargin) % USE DEFAULTS
    end_str = 'msg';
elseif length(varargin) == 1 && ischar(varargin{1}) % USER DEFINED
    end_str = varargin{1};
else
    disp('Too many inputs')
    return   
end
msg_list1 = dir([reg_dir,'*.',end_str]);
msg_list2 = dir([alt_dir,'*.',end_str]);
    

if isempty(msg_list1) && isempty(msg_list2) % files to process?
    disp(['No  *.',end_str, ' files found for: ', flt_name])
    list = [];
    return
end

reg_list = {msg_list1.name}; % reg list names only
reg_sdn  = {msg_list1.datenum}; % reg list modification date only

if ~isempty(msg_list2) % COMPARE ALT FILE TO FILES IN MAIN DIR
%     disp(['Comparing *.',end_str, ' files in regular and ', ...
%         'alternate directories.......'])
     alt_list = {msg_list2.name}; % alt list names only
     alt_sdn  = {msg_list2.datenum}; % alt list modification date only
    
    % IF SAME FILE IN BOTH DIRS CHOOSE LARGER (COULD DO BY DATE INSTEAD?)
    reg_ct = []; alt_ct =[];
    for i = 1:length(alt_list)
        ind = find(strcmp(alt_list{i},reg_list)==1); %index in reg for alt file
        if isempty(ind) % alt file is not in reg dir so keep
            continue
        elseif msg_list1(ind).bytes < msg_list2(i).bytes % alt file bigger
            reg_ct = [reg_ct ,ind]; % index of files to remove from reg list
        else
            alt_ct = [alt_ct ,i];   % index of files to remove from alt list
        end
    end
    % REMOVE *.000.msg from lists if it exists
    reg_list(reg_ct) =[]; % Remove files specified by index
    reg_sdn(reg_ct) =[];
    ind = strcmp([flt_str{1,1},'.000.',end_str], reg_list);
    reg_list(ind)=[]; % Make sure no 000.msg in the list
    reg_sdn(ind) =[];
    ind = strcmp(['all.',end_str], reg_list);
    reg_list(ind)=[]; % Make sure no 000.msg in the list
    reg_sdn(ind) =[];
    
    alt_list(alt_ct) =[];
    alt_sdn(alt_ct) =[];
    ind = strcmp([flt_str{1,1},'.000.',end_str], alt_list);
    alt_list(ind)=[]; % Make sure no 000.msg in the list
    alt_sdn(ind)=[]; 
    ind = strcmp(['all.',end_str], alt_list);
    alt_list(ind)=[]; % Make sure no 000.msg in the list
    alt_sdn(ind) =[]; 
else
    disp(['No alternate directory found for: ',flt_name])
    alt_list ={};
    alt_sdn  =[];
end

% BUILD OUPUT STRUCTURE
list.reg_dir  = reg_dir;
list.alt_dir  = alt_dir;
list.reg_list = reg_list;
list.alt_list = alt_list;
list.reg_sdn = reg_sdn;
list.alt_sdn = alt_sdn;

%clearvars -except list
