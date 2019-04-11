function QC = get_QC_adjustments(float_name, dirs)
% ************************************************************************
% PURPOSE: 
%   This function extracts quality control corrections from a master text
%   file that are applied to the raw float data to create adjusted ARGO
%   variables
%
% USAGE:
%   QC = get_QC_adjustments(float_name, dirs)
%
% INPUTS:
%   float_name = MBARI float ID name or path\MBARI float ID name.
%                 if a complete path is given it overrides the default
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
%    QC = a structure with correction coefficients for each QC'ed variable
%         coefficients = [SDN CYCLE GAIN OFFSET DRIFT]
%  
% EXAMPLES:
%    QC = get_QC_adjustments('9092SOOCN',dirs);
%    QC = get_QC_adjustments('9092SOOCN',[]);
%
% CHANGE LOG
%   02/01/2017 - added code ~line 84 to return QC = [] if no QC adjustment
%       for float
%   08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
% 02/5/2018 - Added code for including date of last QC in output structure.
%	This is used to assist in identifying whether a cycle is real-time or delayed mode, for BRtransfer purposes.

%float_name     = '9092SOOCN'; % TESTING
%float_name = '\\sirocco\wwwroot\lobo\Data\FloatVizData\QC\'FloatQCList.txt';
% ************************************************************************
% SET PATHS AND FILE NAMES
% ************************************************************************
flt_str        = regexpi(float_name,'^\d+','once','match');    
QC_adj_file    = [float_name,'_FloatQCList.txt'];

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
    
    %dirs.QCadj     = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\';
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
% GET QC ADJUSTMENTS FROM MASTER QC LIST
% ************************************************************************
if exist([dirs.QCadj,QC_adj_file],'file')
    fid = fopen([dirs.QCadj,QC_adj_file]);
    
    % FIND SECTION FOR SPECIFIC FLOAT IF IT EXISTS
    tline = ''; % prime engine
    while ischar(tline) && isempty(regexpi(tline,float_name,'once'))
        tline = fgetl(fid);
        tmp = textscan(tline,'%s%s%s');
        xx = [char(tmp{2}) ' ' char(tmp{3})];
        if strcmpi(char(tmp{1}),float_name)==1 && ~isempty(str2num(xx)) % no time of last QC entered
            thetimestamp = datenum(xx,'mm/dd/yy HH:MM');
            QC.date = thetimestamp;
        else
            QC.date = datenum(1900,01,01); % unrealistically old date
            break
        end
    end

    if ~ischar(tline) % -1 (end of file reached - no QC for float
        disp(['No QC adjustments for float: ',float_name])
        QC = [];
        return
        
    else % A match!
        file_pt = ftell(fid); % start of QC lines for given float
        tline   = fgetl(fid); % step to 1st line
        % PREDIMENSION BASED ON EXISTANCE
        while ischar(tline) && isempty(regexpi(tline,'^PREVIOUS','once'))
           if regexp(tline,'^Oxygen','once')
                QC.O.steps = [];
            elseif regexp(tline,'^Nitrate','once')
                QC.N.steps = [];
            elseif regexp(tline,'^pH','once')
                QC.pH.steps = [];
            elseif regexp(tline,'^CHL','once')
                QC.CHL.steps = [];
            elseif regexp(tline,'^BB','once')
                QC.BB.steps = [];
            elseif regexp(tline,'^CDOM','once')
                QC.CDOM.steps = [];
            end
            tline = fgetl(fid);
        end
        
        % NOW GO BACK AND ADD QC STEP DATA
        fseek(fid, file_pt, -1);
        tline   = fgetl(fid); % step to 1st line
        
        while ischar(tline) && isempty(regexpi(tline,'^PREVIOUS','once'))
      
			if regexp(tline,'^Oxygen','once') % only gain value & ONLY 1 LINE
                tmp = textscan(tline,'%s%f%f%f%f%s','Delimiter', ',');
                QC.O.steps  =[QC.O.steps;[tmp{1,2},tmp{1,3},tmp{1,4} tmp{1,5}]];
                QC.O.type   = 'Oxygen';
            elseif regexp(tline,'^Nitrate','once')% cycle gain offset, drift
                tmp = textscan(tline,'%s%f%f%f%f%s','Delimiter', ',');
                QC.N.steps  =[QC.N.steps; tmp{1,2},tmp{1,3},tmp{1,4}, ...
                    tmp{1,5}];
                QC.N.type   = 'Nitrate';
            elseif regexp(tline,'pH,\s+offset','once') % for pump ON V shift
                tmp = textscan(tline,'%s%s%f%s','Delimiter', ',');
                QC.pH.pHpumpoffset = tmp{3};
            elseif regexp(tline,'^pH','once') % cycle, offset, drift
                tmp = textscan(tline,'%s%f%f%f%s','Delimiter', ',');
                QC.pH.steps =[QC.pH.steps;[tmp{1,2}, 1, tmp{1,3}, tmp{1,4}]];
                QC.pH.type   = 'pH';
            elseif regexp(tline,'^CHL','once') % gain, offset & ONLY 1 LINE
                tmp = textscan(tline,'%s%f%f%s','Delimiter', ',');
                QC.CHL.steps =[QC.CHL.steps;[1, tmp{1,2}, tmp{1,3}, 0]];
                QC.CHL.type   = 'CHL';
            elseif regexp(tline,'^BB','once') % gain, offset & ONLY 1 LINE
                tmp = textscan(tline,'%s%f%f%s','Delimiter', ',');
                QC.BB.steps =[QC.BB.steps;[1, tmp{1,2}, tmp{1,3}, 0]];
                QC.BB.type   = 'BB';
            elseif regexp(tline,'^CDOM','once') % gain, offset & ONLY 1 LINE
                tmp = textscan(tline,'%s%f%f%s','Delimiter', ',');
                QC.CDOM.steps =[QC.CDOM.steps;[1, tmp{1,2}, tmp{1,3}, 0]];
                QC.CDOM.type   = 'CDOM';
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        
        % NOW GET PROFILE TIMES FOR EACH QC STEP FROM MESSAGE FILES AND ADD
        % TO QC DATA. TIME NEEDED FOR DRIFT CORRECTION
        % [sdn cycle gain offset drift]
        if isfield(QC,'O')
            sdn = get_QCstep_dates(flt_str,QC.O.steps,dirs);
            QC.O.steps = [sdn,QC.O.steps];
        end
        if isfield(QC,'pH')
            sdn = get_QCstep_dates(flt_str,QC.pH.steps,dirs);
            QC.pH.steps = [sdn,QC.pH.steps];
        end
        if isfield(QC,'N')
            sdn = get_QCstep_dates(flt_str,QC.N.steps,dirs);
            QC.N.steps = [sdn,QC.N.steps];
        end
        if isfield(QC,'CHL')
            sdn = get_QCstep_dates(flt_str,QC.CHL.steps,dirs);
            QC.CHL.steps = [sdn,QC.CHL.steps];
        end
        if isfield(QC,'BB')
            sdn = get_QCstep_dates(flt_str,QC.BB.steps,dirs);
            QC.BB.steps = [sdn,QC.BB.steps];
        end
        if isfield(QC,'CDOM')
            sdn = get_QCstep_dates(flt_str,QC.CDOM.steps,dirs);
            QC.CDOM.steps = [sdn,QC.CDOM.steps];
        end
    end
    
else
    disp(['No list of QC adjustments found: ', ...
        dirs.QCadj,QC_adj_file])
    disp(['No QC adjustments for float: ',float_name])
    disp('QC structure remains empty.')
    QC = [];
end
clearvars -except QC
