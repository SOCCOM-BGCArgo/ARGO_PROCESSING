function sdn = get_QCstep_dates(flt_name,QC_data,dirs)
% PURPOSE: 
%   This function parses an APEX biochemical float *.msg file(s) to get the
%   date(s) corresponding to the profile number(s) in the QC data. Only
%   used to build calibration structure. (Some overhead but used
%   infrequently)
%
% USAGE:
%	sdn = get_QCstep_dates(flt_str,cal_struct)
%
% INPUTS:
%	flt_name   = UW/MBARI Float ID # as a string
%   cal_struct = An n x 4 matrix of QC data [cast # Gain Offset Drift]
%   dirs       = Either an empty variable ora structure with directory
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
%       SDN =  Matlab sdn's corresponding to profile #'s
%               on data flag.
% EXAMPLES:
%   sdn = get_QCstep_dates(flt_str,cal_struct)

% TESTING
% flt_name = cal.info.UW_ID; % TESTING
% QC_data = QC.steps(:,2:end);       % TESTING
% QC_data = QC.O.steps;
% flt_name = '8501';
% ************************************************************************
% CREATE REG & ALT DIR PATHS USING FLOAT NAME
% ************************************************************************

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\'];
    
    dirs.mat       = [user_dir,'DATA\FLOATS\'];
    dirs.cal       = [user_dir,'DATA\CAL\'];
    %dirs.NO3config = [user_dir,'DATA\CAL\'];
    dirs.FVlocal   = [user_dir,'DATA\FLOATVIZ\'];
    dirs.FV        = [user_dir,'DATA\FLOATVIZ\'];
    %dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    
    dirs.QCadj     = [user_dir,'DATA\CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    %dirs.msg       = 'C:\temp\';
    dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    
    dirs.log = [user_dir,'DATA\Processing_logs\'];
    dirs.bat = [user_dir,'batchfiles\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

%#'s but chars follow
flt_num_str  = regexp(flt_name,'^\d{3}\d+(?=\w*)','match'); 
if isempty(flt_num_str) % try just numbers next
    flt_num_str = regexp(flt_name,'^\d{3}\d+','match'); 
end


% TEST FOR 'f' OR 'n' DIRS
if isdir([dirs.msg,'f',flt_num_str{1,1},'\'])      % 'f' dir is APEX UW/MBARI
    reg_dir = [dirs.msg,'f',flt_num_str{1,1},'\'];
elseif isdir([dirs.msg,'n',flt_num_str{1,1},'\'])
    reg_dir = [dirs.msg,'n',flt_num_str{1,1},'\']; % 'n' for NAVIS floats
else
    disp(['Could not find msg file directory for: ',flt_name])
    rows            = size(QC_data,1);
    sdn             = ones(rows,1)* NaN;
    return
end
alt_dir = regexprep(reg_dir,'floats', 'floats\\alternate');

% ************************************************************************
% OPEN MESSSAGE FILES CORRESPONDING TO CASTS IN QC_DATA & GET DATE
% ************************************************************************
msg_time_format = '%*s %*s %*s %*s %*s %s %s %s %s'; 
rows            = size(QC_data,1);
sdn             = ones(rows,1)* NaN;
for i = 1:rows
    tline = ' ';
    msg_file = [flt_num_str{1,1},'.',sprintf('%03.0f',QC_data(i,1)),'.msg'];
    % TRY REGULAR DIR FIRST
    if exist([reg_dir,msg_file],'file')
        fid = fopen([reg_dir,msg_file]);
        while ischar(tline)
            if regexp(tline,'^\$ Profile', 'once')
                tline = strtrim(tline);
                d_str = textscan(tline,msg_time_format,'CollectOutput',1);
                d_str = d_str{1}; % cells: mmm dd hh:mm:ss 2008
                % = Malab format - not sure why I did this but I think there was a
                % burp when I did it a simpler way (= Malab format)
                s1 = [d_str{2},'-',d_str{1},'-',d_str{4},' ',d_str{3}];
                sdn(i) = datenum(s1,0); % Profile end time as Matlab SDN
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        clear fid
    end
    
    % IF NO DATE RECOVERED TRY ALTERNATE DIR
    if isnan(sdn(i)) && exist([alt_dir,msg_file],'file') ==2
        fid = fopen([alt_dir,msg_file]);
        tline = ' ';
        while ischar(tline)
            if regexp(tline,'^\$ Profile', 'once')
                tline = strtrim(tline);
                d_str = textscan(tline,msg_time_format,'CollectOutput',1);
                d_str = d_str{1}; % cells: mmm dd hh:mm:ss 2008
                % = Malab format - not sure why I did this but I think there was a
                % burp when I did it a simpler way (= Malab format)
                s1 = [d_str{2},'-',d_str{1},'-',d_str{4},' ',d_str{3}];
                sdn(i) = datenum(s1,0); % Profile end time as Matlab SDN
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        clear fid
    end
    
    % STILL NO DATE RECOVERED - REPORT IT
    if isnan(sdn(i))
        disp(['No sdn recovered for profile ',num2str(QC_data(i,1))])
    end
    
end



