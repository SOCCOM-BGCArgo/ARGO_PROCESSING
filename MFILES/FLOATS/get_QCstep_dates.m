function sdn = get_QCstep_dates(WMO,QC_data,dirs)
% PURPOSE: 
%   This function parses an APEX biochemical float *.msg file(s) to get the
%   date(s) corresponding to the profile number(s) in the QC data. Only
%   used to build calibration structure. (Some overhead but used
%   infrequently)
%
% USAGE:
%	sdn = get_QCstep_dates(WMO,cal_struct)
%
% INPUTS:
%	WMO   = WMO ID # as a string
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
%   sdn = get_QCstep_dates(WMO,cal_struct)
% CHANGES
%   05/11/2020 JP Added code to correct gps weekday rollover bug if found
%      9634SOOCN is the only float we have with this problem
%
%   03/08/2021 TM Modifications to bring in line with the new MBARI master
%        float list and switch to WMO for processed file names.
%
%   06/15/2022 JP Updated SOLO time stamp extraction. Was using *.mat files
%        but if process mode is "all", all *.mat files are wiped before they
%        are loaded in this function producing an error
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

% LOAD FLOAT LIST FOR MSG FILE LOCATION
% LOAD MBARI FLOAT LIST
list_path  = [dirs.cal, 'MBARI_float_list.mat']; % file path
if exist(list_path, 'file') % load list variables
    load(list_path);
    FLOAT_LIST = d.list;
else
    disp('BUILDING FLOAT ID LIST ...')
    %[float_name, UW_ID, WMO_ID, float type]
    d = MBARI_float_list(dirs);
    FLOAT_LIST = d.list;
end
iMSG         = find(strcmp('msg dir', d.hdr)  == 1);
iWMO         = find(strcmp('WMO',d.hdr) == 1);
iMB          = find(strcmp('MBARI ID',d.hdr) == 1);
iINST        = find(strcmp('INST ID',d.hdr) == 1);
iFtype       = find(strcmp('float type',d.hdr) == 1);
IndexC       = strfind(d.list(:,iWMO),WMO);
Index        = find(not(cellfun('isempty',IndexC)));
MSGloc       = FLOAT_LIST{Index,iMSG};
MBARI_ID_str = FLOAT_LIST{Index,iMB};
INST_ID_str  = FLOAT_LIST{Index,iINST};
Ftype        = FLOAT_LIST{Index,iFtype};

if isfolder(MSGloc)      % from float list specification
    reg_dir = MSGloc;
else
    disp(['Could not find msg file directory for: ',WMO])
    rows            = size(QC_data,1);
    sdn             = ones(rows,1)* NaN;
    return
end

% CHECK FOR SECONDARY (ALTERNATE) DIRECTORY LISTING.
inst      = regexp(reg_dir,'(?<=floats\\)\w+','match','once'); % institute
alt_dir   = regexprep(reg_dir, inst,[inst,'\\alternate']);

% ************************************************************************
% OPEN MESSSAGE FILES CORRESPONDING TO CASTS IN QC_DATA & GET DATE
% ************************************************************************
msg_time_format = '%*s %*s %*s %*s %*s %s %s %s %s'; 
rows            = size(QC_data,1);
sdn             = ones(rows,1)* NaN;
for i = 1:rows
    if strcmp(Ftype,'SOLO')
        tmp = dir([reg_dir, sprintf('CTD\\*%03.0f.phy',QC_data(i,1))]);
        if ~isempty(tmp)
            fp    = fullfile(tmp(1).folder, tmp(1).name);
            tline = ' ';
            fid   = fopen(fp);
            while ischar(tline)
                if regexp(tline,'^MC600 TIME', 'once') % profile end time
                    dstr   = regexp(tline,'\d{4}\d+','once','match');
                    sdn(i) = datenum(dstr,'yyyymmddHHMMSS');
%                     disp(fp) % TESTING
%                     disp(tline)
%                     fprintf('cycle %0.0f %s\n',QC_data(i,1), datestr(sdn(i)))
                    break
                end
                tline = fgetl(fid);
            end
            fclose(fid);
            clear dstr fid tline fp tmp
        end

%     if strcmp(Ftype,'SOLO')
%         DD = load([dirs.mat,WMO,'\',WMO,'.',num2str(QC_data(i,1),'%03.f'),'.mat']);
%         sdn(i) = DD.INFO.sdn;
%         clear DD

    else
        tline = ' ';
        msg_file = [INST_ID_str,'.',sprintf('%03.0f',QC_data(i,1)),'.msg'];
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
    end

    % CHECK FOR BAD GPS TIME (week day roll over problem) 05/11/20 -jp
    % WE ARE GETTING TIME FROM PROFILE TERMINATION TIME NOT GPS FIX TIME!!!
    if isnan(sdn(i))
        disp(['No sdn recovered for profile ',num2str(QC_data(i,1))])
    else
        dvec = datevec(sdn(i));
        if sdn(i) > sdn(1) + 365*20 && dvec(1) == 2099 % 20 yrs from start?
            disp(['GPS time for this profile is > 20 years past start ', ...
                '- gps week day number bug?!!'])
            %dvec = datevec(INFO.sdn);
            dvec(1)  = 1999; % per aoml suggestion don't quite understand jump to 2099
            sdn(i) = datenum(dvec) + 1024*7;
        elseif sdn(i) < sdn(1) % bad gps time fix 10/30/19
            disp('GPS time for this profile is unreasonable - gps week day number bug!!')
            disp(['days since first profile = ',num2str((sdn(i) - ...
                sdn(1)),'%0.0f')]);
            sdn(i) = sdn(i) + 1024*7;
        end
    end
    
    
    
    % STILL NO DATE RECOVERED - REPORT IT
    %     if isnan(sdn(i))
    %         disp(['No sdn recovered for profile ',num2str(QC_data(i,1))])
    %     end
    
end



