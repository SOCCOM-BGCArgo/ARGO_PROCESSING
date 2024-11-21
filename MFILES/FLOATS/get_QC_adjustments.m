function QC = get_QC_adjustments(WMO, dirs)
% ************************************************************************
% PURPOSE:
%   This function extracts quality control corrections from a master text
%   file that are applied to the raw float data to create adjusted ARGO
%   variables
%
% USAGE:
%   QC = get_QC_adjustments(WMO, dirs)
%
% INPUTS:
%   WMO        = WMO ID, as a string
%
%   dirs       = Either an empty variable or a structure with directory
%                strings where files are located. It must contain these
%                fields at a minimum or be empty. If dirs is empty default
%                paths will be used.
%
%   dirs.mat       = path to *.mat profile files for ARGO NETCDF
%   dirs.QCadj      = path to QC adjustment list for all floats
%
% OUTPUTS:
%    QC = a structure with correction coefficients for each QC'ed variable
%         coefficients = [SDN CYCLE GAIN OFFSET DRIFT]
%
% EXAMPLES:
%    QC = get_QC_adjustments('5905634',dirs);
%    QC = get_QC_adjustments('5905634',[]);
%    QC = get_QC_adjustments('NO_WMO_un0412',dirs);
%
% CHANGE LOG
%   02/01/2017 - added code ~line 84 to return QC = [] if no QC adjustment
%       for float
%   08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
%   02/5/2018 - Added code for including date of last QC in output structure.
%	This is used to assist in identifying whether a cycle is real-time or delayed mode, for BRtransfer purposes.
%
%   03/08/2021 TM Modifications to bring in line with the new MBARI master
%                 float list and switch to WMO for processed file names.
%
%  10/04/2024 JP changes made to accept changes to QCList files from
%                sageV2. Also some code clean up including eliminating the
%                need for the "get_QCstep_dates" function

% TESTING
% WMO = '4903591'
% WMO = '5906568'
% WMO = '5905102'; % ua12779 jp qclist has pH lines, VM qclist does not
% dirs = [];

% ************************************************************************
% SET PATHS AND FILE NAMES
% ************************************************************************
QC_list_file = sprintf('%s_FloatQCList.txt',WMO);
%QC_list_file = sprintf('%s_FloatQCList_JP.txt',WMO); % TESTING
QC           = []; % default output


% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
% need directory paths to point to QC list files & to * .mat files
if isempty(dirs)
    % i.e. userpath = 'C:\Users\jplant\Documents\MATLAB'
    user_dir    = fullfile(userpath,'ARGO_PROCESSING\DATA');    
    dirs.mat    = fullfile(user_dir,'FLOATS');
    dirs.QCadj = fullfile(user_dir,'CAL\QC_LISTS');
    %dirs.QCadj = fullfile('C:\temp\TEST\VM\QC_LISTS\'); %Testing
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

if ~isfield(dirs,'QCadj') && isfield(dirs, 'QClists')
    dirs.QCadj = dirs.QClists; % backwards compaitibility
end

% ************************************************************************
% GET QC ADJUSTMENTS FROM MASTER QC LIST
% ************************************************************************
qc_path = fullfile(dirs.QCadj, QC_list_file);
if ~isfile(qc_path)
    fprintf(['QC LIST FILE NOT FOUND : %s\n', ...
      'QC structure remains empty.'],  qc_path);
    return
end

% ************************************************************************
% READ QC LINES INTO A CELL ARRAY & SUBSET TO RECENT, GET QC DATE TOO
fid     = fopen(qc_path);
qc_cell = cell(300,1);
tline   = '';
ct      = 0;

while ischar(tline)
    if regexpi(tline,'PREVIOUS','once') % Done with parsing current QC
        break
    elseif ~isempty(regexp(tline,['^',WMO],'once')) % should be first line
        dstr = regexp(tline,'\d{2}/\d{2}/\d{2} \d{2}:\d{2}','match','once'); %old
        if ~isempty(dstr)
            QC.date = datenum(dstr, 'mm/dd/yy HH:MM');
        else
            QC.date = datenum(1900,01,01); % unrealistically old date
            fprintf('QC date stamp was not resolved for QC list file: %s\n',...
                QC_list_file)
        end
    else
        ct = ct+1;
        qc_cell{ct} = tline;
    end
    tline = fgetl(fid);
end
qc_cell = qc_cell(1:ct);
fclose(fid);

% ************************************************************************
% NOW BUILD CORRECTION STRUCTURES
params = {'Oxygen' 'Nitrate' 'pH';  % QC list file line ID
          'O'      'N'       'pH'}; % Fields for Process_APEX_float 

for ct = 1:size(params,2) % PREDIM BASED ON AVAILIBILITY
    if any(contains(qc_cell,params{1,ct}))
        QC.(params{2,ct}).steps = [];
        QC.(params{2,ct}).type = params{1,ct};
    end
end

% ADD QC CORRECTION MATRICES
for ct = 1: size(qc_cell,1)
    tmp = regexp(qc_cell{ct},',','split');

    if strcmp(tmp{1},'Oxygen')
        QC.O.steps  =[QC.O.steps; str2double(tmp(2:end))];

    elseif strcmp(tmp{1},'Nitrate') % cycle gain offset, drift
        QC.N.steps  =[QC.N.steps; str2double(tmp(2:end))];

    elseif regexp(qc_cell{ct},'pH.+offset','once') % string defining correction type
        QC.pH.pHpumpoffset = ...
            regexp(qc_cell{ct},'linear|poly|mixed','once','match');
        if isempty(QC.pH.pHpumpoffset)
            QC.pH.pHpumpoffset = 'none';
        end

    elseif strcmp(tmp{1},'pH')  % cycle, offset, drift  or  cycle gain offset, drift
        QC.pH.steps =[QC.pH.steps; str2double(tmp(2:end))];

    % Add QC_info if it exists. Down the line might be helpful for argo
    % sci meta info
    elseif strcmp(tmp{1},'QC parameter')
        info = tmp(3:end);
        QC.QC_info.(tmp{2}) = cell2struct( info(2:2:length(info)), ...
            info(1:2:length(info)-1) ,2);
    end
end

% Last step Backwards compatibility for older  pH QC lines (n x 3 with no gain col)
if isfield(QC,'pH') && size(QC.pH.steps,2) == 3 % old style. no gain col so add one
    QC.pH.steps = [QC.pH.steps(:,1), ones(size(QC.pH.steps(:,1))), ...
        QC.pH.steps(:,2:3)];
end

% ************************************************************************
% NOW GET PROFILE TIMES FOR EACH QC STEP AND ADD TO QC DATA. 
% TIME IS NEEDED FOR DRIFT CORRECTION TO DATA [sdn cycle gain offset drift]
% Use to use msg file, but if cycle exists in list file a mat file must
% exist so use it! Play some games with the path to the mat files depending
% on who is calling the function. If not running on the VM, default to chem
% otherwise use local source for mat files.
% DATES COULD BE SUPPLIED FROM SAGEv2 AT SOME POINT TO SIMPLIFY....

if contains(userpath,'bgcargovm')
    dirs.mat = dirs.mat; % No change
elseif isfolder('\\atlas\chem\ARGO_PROCESSING\DATA\FLOATS')
    dirs.mat = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATS';
end

for ct = 1:size(params,2)
    if ~isfield(QC, params{2,ct}) % no field so move on!
        continue
    end
    cgod = QC.(params{2,ct}).steps; %correction matrix
    cgod = [NaN(size(cgod(:,1))), cgod]; % add initial NaN col for time stamps
    fd   = fullfile(dirs.mat, WMO); % matfile dir

    for node_ct = 1:size(cgod,1) % step through node cycles
        fn = sprintf('%s.%03.0f.mat', WMO, cgod(node_ct,2));
        load(fullfile(fd,fn),'INFO') % load INFO structure from mat file
        if isfield(INFO, 'sdn')
            cgod(node_ct,1) = INFO.sdn; % profile termination time extracted
        else
            fprintf(['Cycle profile termination timed was not extracted', ...
                'from  %s\n'],fullfile(fd,fn));
        end
    end
    QC.(params{2,ct}).steps = cgod; % update data in QC structure
end

clearvars -except QC
