function list = get_msg_list_sO2(flt_name,mydirs,varargin)
% PURPOSE: 
%   This function returns a structure with the regular and alternate
%   directory paths for the given float as well as msg file names within
%   each directory. The default message file  is '*.msg' but an optional
%   input string can change this. If files exist in the alternate
%   directory, duplicate file sizes are compared and the smaller files are
%   removed from a list. Designed for use with other functions / scripts
%   processing float data
%
%   If you are not working with an 'alternate' directory, and wish to just
%   retrieve a file list from a single working float directory, you may
%   set "mydirs" to a single string directory path.
%
% USAGE:
%	list = get_msg_list(flt_name,dirs,'msg')
%   list = get_msg_list(flt_name, [])
%
% INPUTS:
%   float name = a string with at least UW/MBARI float #
%   dirs       = Either a single directory path string, or a 2x1 cell
%                containing regular and alternate directory path strings.
%
%
%   varargin = optional string defining msg type. Ex:
%                   'msg'  APEX message file (DEFAULT)
%                   'isus' ISUS message file
% OUTPUTS:
%   list = a structure of the data.
%       list.reg_dir  = main msg file diectory path
%       list.alt_dir  = alternate message file directory path if it exists
%       list.reg_list = main directory * msg file list
%       list.alt_list = alternate directory message file list
% EXAMPLES:
%   jp = get_msg_list('7601StnP',dirs);
%   jp = get_msg_list('7601',[]);
%   jp = get_msg_list('7601StnP','dirs,'msg');
%   jp = get_msg_list('7601StnP',[],'isus');

%
% Compile a list of float *.msg files from reg and alt directories
% use file size or date to choose which ones to keep

%flt_name = '7601StnP';
%flt_name ='9254SoOcn';
%flt_name ='6966Hawaii';
%flt_name ='7593HAWAII';
%
%
% AUTHOR:  TANYA MAURER
% DATE:    01/04/17

% ************************************************************************
% CREATE REG & ALT DIR PATHS USING FLOAT NAME
% ************************************************************************
if isempty(mydirs);
    disp('Check "mydirs" input. Must be a single directory string, or a cell array of directory strings');
    return
end

if ~iscell(mydirs); %check if single entry or double
    dirs_msg{1} = mydirs;
    reg_dir = cell(1,1);
else
    dirs_msg{1} = mydirs{1};
    dirs_msg{2} = mydirs{2};    
    reg_dir = cell(2,1);
end
 
old_flt_dir     = '\\atlas\chemwebdata\floats\_oldftp\'; % A few old floats found here (***MBARI***)

flt_str  = regexp(flt_name,'^\d{3}\d+(?=\w*)','match'); %#'s but chars follow
% if isempty(flt_str)
%     flt_str = regexp(flt_name,'^\d{3}\d+','match'); % try just numbers next
% end
flt_num  = str2num(flt_str{1,1});

for i = 1:length(dirs_msg);
    % TEST FOR 'f' OR 'n' DIRS
    if isdir([dirs_msg{i},'f',flt_str{1,1},'\'])      % 'f' dir is APEX UW/MBARI
        reg_dir{i} = [dirs_msg{i},'f',flt_str{1,1},'\'];
    elseif isdir([dirs_msg{i},'n',flt_str{1,1},'\'])
        reg_dir{i} = [dirs_msg{i},'n',flt_str{1,1},'\']; % 'n' for NAVIS floats
    elseif isdir([old_flt_dir,flt_str{1,1},'\data\'])
        disp('NON STANDARD DATA DIR!!')
        disp(['DATA FOUND AT:  ',old_flt_dir,flt_str{1,1},'\data\'])
        reg_dir{i} = [old_flt_dir,flt_str{1,1},'\data\']; % SOME OLD DATA HERE
    else
        disp(['Could not find regular or alternate msg file directory for: ',flt_name])
%         list = [];
%         return
    end
end

if length(reg_dir) == 2; %alternate directory listed;
    alt_dir = reg_dir{2};
    reg_dir = reg_dir{1};
else
    alt_dir = [];
    reg_dir = char(reg_dir);
end

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
% REMOVE *.000.msg from lists if it exists
ind = strcmp([flt_str{1,1},'.000.',end_str], reg_list);
msg_list1(ind) = [];
reg_list(ind)=[]; % Make sure no 000.msg in the list
ind = strcmp(['all.',end_str], reg_list);
reg_list(ind)=[]; % Make sure no 000.msg in the list
msg_list1(ind) = [];

if ~isempty(msg_list2) % COMPARE ALT FILE TO FILES IN MAIN DIR
%     disp(['Comparing *.',end_str, ' files in regular and ', ...
%         'alternate directories.......'])
     alt_list = {msg_list2.name}; % alt list names only
    
    % IF SAME FILE IN BOTH DIRS CHOOSE LARGER (COULD DO BY DATE INSTEAD?)
    reg_ct = []; alt_ct =[];
    for i = 1:length(alt_list)
        log_i = strcmp(alt_list{i},reg_list);
        ind = find(log_i==1); %index in reg for alt file
        if isempty(ind) % alt file is not in reg dir so keep
            continue
        elseif msg_list1(ind).bytes < msg_list2(i).bytes % alt file bigger
            reg_ct = [reg_ct ,ind]; % index of files to remove from reg list
        else
            alt_ct = [alt_ct ,i];   % index of files to remove from alt list
        end
    end
    reg_list(reg_ct) =[]; % Remove files specified by index
    
    % REMOVE *.000.msg from alt_list if it exists
    alt_list(alt_ct) =[];
    ind = strcmp([flt_str{1,1},'.000.',end_str], alt_list);
    alt_list(ind)=[]; % Make sure no 000.msg in the list
    ind = strcmp(['all.',end_str], alt_list);
    alt_list(ind)=[]; % Make sure no 000.msg in the list  
else
    disp(['No alternate directory found for: ',flt_name])
    alt_list ={};
end

% BUILD OUPUT STRUCTURE
list.reg_dir  = reg_dir;
list.alt_dir  = alt_dir;
list.reg_list = reg_list;
list.alt_list = alt_list;

clearvars -except list