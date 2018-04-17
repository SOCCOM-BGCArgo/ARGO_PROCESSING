function tf = get_argo_nc(WMO, cycle, site_flag, dirs)
% This function tries to copy an ARGO  'D' or 'R' NetCDF file  to your
% local computer from either IFREMER or USGODAE
%
% INPUTS:
%   WMO       - WMO ID as a number or char array
%   cycle     - float profile cycle as number or char array
%   site_flag - 1 = USGODAE, 0 = IFREMER
%   dirs      - a structure of directory paths, if '[]
%
% OUTPUTS:
%   tf.status - 1 = success, 0 = file was not copied
%   tf.path   - path to file
%   tf.name   - file name

% THESE WILL BE FUNCTION INPUTS EVENTUALLY
% 0506SOOCN	506	5904670	NAVIS 9313SOOCN	9313	5904474	APEX
% WMO       = 5904670; % TEST 0506SOOCN	506	5904670	NAVIS
% cycle     = 38;       % TEST
% site_flag = 1;       % 1 = usgodae, 0 = ifremer
% dirs      = [];


% ************************************************************************
% CHECK INPUTS AND BUILD FILE NAME
if isnumeric(WMO)
    WMO = num2str(WMO);
end
if isnumeric(cycle)
    cycle = sprintf('%03.0f',cycle);
end

argo_fname = ['D', WMO, '_', cycle, '.nc'];

tf.status = 0;
tf.path   = '';
tf.name   = '';

% SET DEFAULT DIR STRUCTURE
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.temp      = 'C:\temp\';
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end
    
% ************************************************************************
% OFTEN MATLAB FTP OBJECT NEEDS TO BE PUT IN PASIVE MODE
% Good hacks: 
%    http://blogs.mathworks.com/pick/2015/09/04/passive-mode-ftp/
%    http://undocumentedmatlab.com/blog/solving-an-mput-ftp-hang-problem
% ************************************************************************

% *************************************************************************
% CHOOSE DATA SITE
% *************************************************************************
if site_flag == 1 % United states
    ftp_target = 'usgodae.org';
    ftp_dir  = '/pub/outgoing/argo/dac/aoml/';
elseif site_flage == 0 % France
    ftp_target = 'ifremer.fr';
    ftp_dir  = '/ifremer/argo/dac/aoml/';
else
    disp('Unknown data site flag')
    return
end

% *************************************************************************
% TRY CONNECTING TO FTP SERVER
% *************************************************************************
try
    f = ftp(ftp_target); % Connect to FTP server
    binary(f)
catch
    disp(['Could not connect to ftp server at: ',ftp_target])
    disp('No files were obtained')
    return
end

% Enter passive mode by accessing the java object - this is the tricky part
cd(f);
sf = struct(f);
sf.jobject.enterLocalPassiveMode();

% GET LIST OF FILES IN WMO DIR
ftp_path = [ftp_dir, WMO, '/profiles/'];
cd(f, ftp_path);       
argo_dir   = dir(f);   % SOCCOM FLOAT DIR
argo_files = {argo_dir.name}'; % Cell array of ARGO netcdf files

% CHECK FOR D or R FILE
t1    = strcmp(argo_files, argo_fname);
if sum(t1) == 1 % D file exists
    fname = argo_files{t1}; %D file
elseif sum(t1) == 0 % No D file, look for for R file
    argo_fname = regexprep(argo_fname, 'D', 'R');
    t1         = strcmp(argo_files,argo_fname);
    if sum(t1) == 1
        fname = argo_files{t1}; % R file
    else
        disp(['No NetCDF file found for ', WMO,' cycle # ', cycle]);
        return % NO FILE!
    end
end

% TRY AND FTP COPY
try
    str = mget(f, fname, dirs.temp);
catch
    disp(['File found but ftp copy failed for ',WMO, ...
        ' cycle # ', cycle]);
    close(f);
    delete(f);
    return
end


tf.status = 1;
tf.path = str{1};
tf.name = fname;
close(f);
clearvars -except tf 







