function cal = get_float_cals(MBARI_ID_str, dirs)
% PURPOSE: 
%   This function tries to parse all the calibration data for the
%   biogoechemical sensors on a given float from a *FloatConfig.txt file
%   and from a *cal nitrate file if it exists. If a sensor is detected in 
%   the config or cal text files a stucture field with calibration data
%   will be created.
%   FUNCTION ASSUMES ONLY ONE SENSOR TYPE PER FLOAT (no good for 3xO2 floats)
%
% USAGE:
%	data = get_float_cals(MBARI_ID_str, dirs)
%
% INPUTS:
%   MBARI_ID_str  = MBARI Float ID as a string
%   dirs       = Either an empty variable or a structure with directory
%                strings where files are located. It must contain these
%                fields at a minimum or be empty. If dirs is empty default
%                paths will be used.
%
%   dirs.cal       = path to *.cal (nitrate) and cal####.mat (working float cal) files
%   dirs.temp      = path to temporary working dir
%
% OUTPUTS:
%       cal =   a structure of calibration info. Structure fields change
%               depending on the sensors detectd in the cal file and float
%               type (APEX or NAVIS).
%       
%       cal.O       Oxygen calibration coefficients
%                   (Aanderaa 3830&4330, SBE63)
%       cal.CHL     biooptics calibration coefficients (FLBBAP2 or MCOMS
%       cal.BB      biooptics calibration coefficients (FLBBAP2 or MCOMS
%       cal.CDOM    biooptics calibration coefficients (FLBBAP2 or MCOMS
%       cal.pH   	Durafet pH calibration coefficients
%       cal.N   	ISUS or SUNA calibration coefficients

%       cal.info    float information (WMO ID#, UW/MBARI ID#,
%                                      FILE CREATION DATE, ETC)
%
% REQUIRES: MBARI_float_list.m
%           
% EXAMPLES:
%   jp = get_float_cals('ua8501S', dirs)
%   jp = get_float_cals('ua8501S', [])
%
% Created 12/29/2015 by jp
% CHANGES
% 07/19/17 Added functionality to extract comments line at bottom of
%   *FloatConfig.txt files if they exist. These comments can then be
%   printed to the ODV txt files in argo2ODV*. m files
% 09/10/18 Added code to look for seconday pH pressure coefficients. At
%   this point this change onlly afects 0690. -jp
% 10/26/18 Added code to extract non zero Aanderaa 4330 SVUFoilCoef's. -jp
%
% 02/28/21 JP, Rebuilt for GO-BGC. Function is no longer float type
%              dependant. Now sensor type dependant.


% ************************************************************************
% FOR TESTING
%dirs =[];
% MBARI_ID_str = 'ua7614';
% MBARI_ID_str = 'un0412';
% MBARI_ID_str = 'ua12541';
% MBARI_ID_str = 'un0037';
% MBARI_ID_str = 'un0276';
% MBARI_ID_str = 'un0569';
%MBARI_ID_str = 'ua5143';
%MBARI_ID_str = 'un0565'; % has bbp532
% MBARI_ID_str = 'un0690';
% MBARI_ID_str = 'un1173';
%MBARI_ID_str = 'un0510';
%MBARI_ID_str = 'wn857';
%MBARI_ID_str = 'ua0068';
%MBARI_ID_str = 'wa760';
% ************************************************************************

% ************************************************************************
% SET DIRECTORY AND VARIABLES
% ************************************************************************
%INST_ID_str  = regexp(MBARI_ID_str,'^\d+','once','match');
config_file  = [MBARI_ID_str,'_FloatConfig.txt'];
NO3_calfile  = [MBARI_ID_str, '.cal'];
%master_list  = 'new_MBARI_float_list.mat';
master_list  = 'MBARI_float_list.mat';

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.cal        = [user_dir,'CAL\']; % save created lists here
    dirs.MasterDir  = [dirs.cal,'FLOAT_CONFIG\']; % MBARI master list source
    dirs.msg        = '\\seaecho.shore.mbari.org\floats\';
    dirs.bfile      = [getenv('USERPROFILE'),'\Documents\MATLAB\ARGO_MBARI2AOML\'];
    dirs.temp      = 'C:\temp\';   
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

info.WMO_ID     = '';
info.file_date  = datestr(now,'mm/dd/yyyy HH:MM');
info.float_type = '';

disp(['Building calibration file for ',MBARI_ID_str]);

% ************************************************************************
% CHECK IF MBARI ID IS ON THE MASTER LIST
% ************************************************************************
if exist([dirs.cal, master_list], 'file') == 2 % load master list
    load([dirs.cal, master_list]);
else
    disp('Master list not found, building a new list ....')
    d = MBARI_float_list(dirs); 
end

I.MBA  = strcmp(d.hdr,'MBARI ID'); % Define lookup indices
I.INST = strcmp(d.hdr,'INST ID');
I.WMO  = strcmp(d.hdr,'WMO');
I.TYP  = strcmp(d.hdr,'float type');
I.DIR  = strcmp(d.hdr,'msg dir');
I.PRG  = strcmp(d.hdr,'Program');
I.REG  = strcmp(d.hdr,'Region');
I.NC   = strcmp(d.hdr,'NC template');
I.BF   = strcmp(d.hdr,'tf Bfile');

t1     = strcmp(MBARI_ID_str, d.list(:,I.MBA));
if sum(t1) == 0
    disp([MBARI_ID_str,'Not on master list. Refreshing master list....'])
    d  = MBARI_float_list(dirs);
    t1 = strcmp(MBARI_ID_str, d.list(:,I.MBA));
    if sum(t1) == 0
        disp(['There is no master listing for ',MBARI_ID_str])
        tf = exist([dirs.cal,'FLOAT_CONFIG\', config_file],'file') == 2;
        if tf == 1
            disp(['But ', dirs.cal,'FLOAT_CONFIG\', config_file,' exists!']);
        else
            disp(['Please build ',config_file,' and try again'])
        end
        return
    end     
end

info.name       = MBARI_ID_str;
info.INST_ID    = d.list{t1, I.INST};
info.WMO_ID     = d.list{t1,I.WMO};
info.float_type = d.list{t1, I.TYP};
info.msg_dir    = d.list{t1, I.DIR};
info.tf_bfile   = d.list{t1, I.BF}; % WMO or not?
info.Program    = d.list{t1, I.PRG};
info.Region     = d.list{t1, I.REG};
info.NCtemplate = d.list{t1, I.NC};

if info.tf_bfile == 0 && isempty(regexp(info.WMO_ID,'\d{7}','once'))
    disp([MBARI_ID_str,' will never be assigned a WMO number'])
elseif isempty(regexp(info.WMO_ID,'\d{7}','once'))
    disp([MBARI_ID_str,' has not been assigned a WMO number yet'])
end

if info.float_type(1) ~= upper(info.name(2))
    disp(['WARNING: mbari ID & float type ID do not match (', ...
        MBARI_ID_str,' vs ',info.float_type,')'])
    disp('Please check config file name & info')
    return
end

% ************************************************************************
% CHECK IF FILES EXIST: FLOAT CONFIGURATION   & NO3 CAL FILES
% ************************************************************************
if exist([dirs.cal,'FLOAT_CONFIG\', config_file],'file') ~= 2
    disp(['Calibration *.txt file not found for ',MBARI_ID_str, ...
        ' - stopping script']);
    cal = [];
    return    
end

% CHECK FOR EXISTING ISUS CAL FILE 
info.isus_flag  = 1;
if exist([dirs.cal,'NO3_CAL\',NO3_calfile],'file') ~= 2
    disp([MBARI_ID_str,': ',NO3_calfile,' not found - nitrate spectra',...
        ' will not be re-processed']);
    info.isus_flag    = 0; % NO NO3 - add to ouput structure at end
end

% ************************************************************************
% PARSE CONFIG FILE TO GET OXYGEN CALIBRATION COEFFICIENTS
% Examples of OptodeSn output:
% 'OptodeSn = 737'                   (Aanderaa 3830 older)
% 'OptodeSn = 4330 1168'             (Aanderaa 4330 newer)
% 'OptodeSn = SBE63 rev J S/N 0455'  (Seabird 63)
% ************************************************************************
info.O2_flag    = 0;
tline = ' ';% initialize some variables
fid   = fopen([dirs.cal,'FLOAT_CONFIG\', config_file]);
if fid == -1
    disp(['Could not open *.txt config file for ',MBARI_ID_str])
    cal = [];
    return
end

while ischar(tline)   
    if regexp(tline,'OptodeSn =','once') % APEX
        optode_info = regexp(tline, '\d+', 'match'); % parse any numbers from SN line
        break
    elseif regexp(tline,'^O2 sensor','once') % NAVIS SBE63 or SBE83 
        optode_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
        break
    end
    tline = fgetl(fid);

end

if ~ischar(tline) % -1 (end of file reached - no O2)
    disp(['NO oxygen sensor detected for float ',MBARI_ID_str])
    
% NEWER STYLE AANDERAA OXYGEN CALIBRATION FORMAT (4330)    p02 = f(Phase)
elseif length(optode_info) == 2 && strcmp(optode_info{1},'4330')
    info.O2_flag = 1;
    O.type       = optode_info{1};
    O.SN         = optode_info{2};
    disp(['Aanderaa oxygen optode detected (', ...
        O.type,' SN: ',O.SN,')'])
    while tline ~= -1 % EXTRACT COEFFICIENTS
        % comments at end stop searching
        if ~isempty(regexp(tline,'^COMMENTS:', 'once'))
            break
        end
        
        if ~isempty(regexp(tline,'PhaseCoef','once'))
            ind = regexp(tline,'PhaseCoef\s+\d+\s+\d+\s+','once','end');
            O.PCoef = sscanf(tline(ind:end),'%f',4);
        end
        
        if ~isempty(regexp(tline,'FoilCoefA','once'))
            ind = regexp(tline,'FoilCoefA\s+\d+\s+\d+\s+','once','end');
            FCoefA = sscanf(tline(ind:end),'%f',14);
        end
        
        if ~isempty(regexp(tline,'FoilCoefB','once'))
            ind = regexp(tline,'FoilCoefB\s+\d+\s+\d+\s+','once','end');
            FCoefB = sscanf(tline(ind:end),'%f',14);
        end
        
        if ~isempty(regexp(tline,'FoilPolyDegT','once'))
            ind = regexp(tline,'FoilPolyDegT\s+\d+\s+\d+\s+','once','end');
            O.PolyDegT = sscanf(tline(ind:end),'%f',28);
        end
        
        if ~isempty(regexp(tline,'FoilPolyDegO','once'))
            ind = regexp(tline,'FoilPolyDegO\s+\d+\s+\d+\s+','once','end');
            O.PolyDegO = sscanf(tline(ind:end),'%f',28);
        end
        
        if ~isempty(regexp(tline,'SVUFoilCoef','once'))
            ind = regexp(tline,'SVUFoilCoef\s+\d+\s+\d+\s+','once','end');
            SVU = sscanf(tline(ind:end),'%f',7);
            if sum(SVU(:)) ~= 0
                O.SVUFoilCoef = SVU;
            end
        end
        
        if ~isempty(regexp(tline,'ConcCoef','once'))
            ind = regexp(tline,'ConcCoef\s+\d+\s+\d+\s+','once','end');
            O.ConcCoef = sscanf(tline(ind:end),'%f',2);
        end
        
        tline = fgetl(fid);
    end
    O.FCoef =[FCoefA;FCoefB];
    t1 = O.FCoef == 0; % Look for zero & remove to shorten calc, 0 * x = 0
    O.FCoef(t1)    = [];
    O.PolyDegT(t1) = [];
    O.PolyDegO(t1) = [];
    
    clear FCoefA FCoefB opt_info ind t1
    
% ************************************************************************
% NAVIS SBE63 / SBE83 OXYGEN CALIBRATION FORMAT
elseif length(optode_info) == 2 && strncmp(optode_info{1},'SBE63',3)
    info.O2_flag = 1;
    O.type       = optode_info{1};
    O.SN         = optode_info{2};
    disp(['SBE oxygen optode detected (', O.type,' SN: ',O.SN,')'])
    
    while ischar(tline)
        if regexp(tline,'^O2 sensor.+Temp coefficents', 'once')
            TC = regexp(tline,',','split');
            O.TempCoef = str2double(TC(2:end)); % TA0-3
        elseif regexp(tline,'^O2 sensor.+Phase coefficents', 'once')
            PC = regexp(tline,',','split');
            O.PhaseCoef = str2double(PC(2:end)); % TA0-3
        end
        tline =fgetl(fid);
    end
    clear TC PC
% ************************************************************************
% OLD STYLE AANDERAA OXYGEN CALIBRATION FORMAT (3830)    02 = f(Phase)
elseif length(optode_info) == 1 % NO senser type just S/N
info.O2_flag = 1;
 
    info.O2_flag = 1;
    O.type       = '3830';
    O.SN         = optode_info{1};
    disp(['Aanderaa oxygen optode detected (', O.type,' SN: ',O.SN,')'])
    coef = textscan(fid,'%s %*s %f'); %{1}= ID, {2}= value
    %BUILD COEFF ARRAYS for POLYVAL
    coef_srch = {'C0', 'C1', 'C2', 'C3', 'C4', 'P'};
    
    for i = 1:length(coef_srch) % build coeff arrays with eval
        % Aanderaa = ascending powers
        ind1 = find(strncmp(coef_srch{i},coef{1},length(coef_srch{i})));
        ind1 = sort(ind1,'descend')'; % Matlab needs descending powers
        O.(['p',coef_srch{i}]) = coef{1,2}(ind1);
    end
    
    clear coef_format coef_path coef_srch calfile ind1 i opt_info 
    
end

if info.O2_flag == 1
    cal.O = O;
else
    disp(['Can not determine oxygen optode type for ', MBARI_ID_str])
end

% ************************************************************************
% PARSE config FILE TO CHECK FOR FLBB/MCOMS CHL AND BACKSCATTER AND CDOM
% CAL COEFFICIENTS
% ************************************************************************
frewind(fid)
info.chl_flag   = 0;
tline = ' ';% initialize some variables
while ischar(tline)
    
    % FLBB BLOCK
    if regexp(tline,'FLBB','once') % FLBB SN
        CHL.type = regexp(tline,'FLBB\w+','match','once');
        CHL.SN   = regexp(tline,'\d+$','match','once');
        disp(['FLBB Bio-optics detected (',CHL.type,' ',CHL.SN,')'])
        info.chl_flag   = 1;
    elseif regexp(tline,'^\d+.+ChlDC','once') % FLBB CHL DC
        CHL.ChlDC = sscanf(tline,'%f',1);
    elseif regexp(tline,'^\d+.+ChlScale','once') % FLBB CHL DC
        CHL.ChlScale = sscanf(tline,'%f',1);
    elseif regexp(tline,'^\d+.+BetabDC','once') % FLBB BBP DC
        BB.type    = CHL.type;
        BB.SN      = CHL.SN;
        BB.BetabDC = sscanf(tline,'%f',1);
    elseif regexp(tline,'^\d+.+BetabScale','once') % FLBB CHL DC
        BB.BetabScale = sscanf(tline,'%f',1);
    
    % MCOMS BLOCK
    elseif regexp(tline,'^MCOM.+Chl','once') % MCOMS  CHL
        tmp          = regexp(tline,',','split');
        CHL.type     = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        CHL.SN       = regexp(tmp{1,1},'\d+(?=\))','once','match');
        CHL.ChlDC    = str2double(tmp{1,2});
        CHL.ChlScale = str2double(tmp{1,3});
        disp(['MCOMS Bio-optics detected (',CHL.type,' ',CHL.SN,')'])
        info.chl_flag   = 1;
    elseif regexp(tline,'^MCOM.+Backscatter700','once') % MCOMS BBP700
        tmp           = regexp(tline,',','split');
        BB.type       = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        BB.SN         = regexp(tmp{1,1},'\d+(?=\))','once','match');
        BB.BetabDC    = str2double(tmp{1,2});
        BB.BetabScale = str2double(tmp{1,3});        
    elseif regexp(tline,'^MCOM.+DOM|^MCOM.+Backscatter532','once')
        % MCOMS CDOM or BBP532
        tmp             = regexp(tline,',','split');
        CDOM.type       = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        CDOM.SN         = regexp(tmp{1,1},'\d+(?=\))','once','match');
        CDOM.CDOMDC    = str2double(tmp{1,2});
        CDOM.CDOMScale = str2double(tmp{1,3});         
    end
    
    tline = fgetl(fid);
end

if info.chl_flag == 0 
    disp(['No bio optics detected for float ',MBARI_ID_str])
else    % ADD CALIBRATION DATA TO STRUCTURE
    cal.CHL = CHL;
    cal.BB  = BB;
    if exist('CDOM', 'var') % NAVIS
        cal.CDOM = CDOM;
    end
end

clear tmp tline t1 BB CHL CDOM

% ************************************************************************
% PARSE Config FILE TO CHECK FOR pH COEFFICIENTS
% ************************************************************************
frewind(fid) % go back to top of file
tline = ' ';% initialize some variables
info.pH_flag = 0;

while ischar(tline)
    
    % MBARI pH CALIBRATION BLOCK
    if regexp(tline,'^Durafet SN','once') % MBARI Calibration
        pH.type = regexp(tline,'^\w+','match','once');
        pH.SN   = regexp(tline,'(?<=\=\s+)\w+','match','once'); % # after =
        disp(['MBARI pH sensor detected (',pH.type,' ',pH.SN,')'])
        info.pH_flag = 1;
    elseif regexp(tline,'^\d.+number\s+of\s+calibration','once')
        coef_ct = str2double(regexp(tline,'^\d','once','match'));
        tmp = ones(coef_ct,1)*NaN;
        for i = 1:coef_ct
            tline = fgetl(fid);
            tmp(i) = sscanf(tline,'%f',1);
        end
        pH.k0     = tmp(1);
        pH.k2     = tmp(2);
        pH.pcoefs = tmp(3:coef_ct);
    
    % SBE PRIMARY pH CALIBRATION BLOCK
    elseif regexp(tline,'^pH sensor coefficients','once') % SBE Calibration
        tmp       = regexp(tline,',','split');
        pH.type   = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        pH.SN     = regexp(tmp{1,1},'\w+(?=\))','once','match');
        pH.k0     = str2double(tmp(2));
        pH.k2     = str2double(tmp(3));
        pH.pcoefs = str2double(tmp(4:end))';
        disp(['SBE pH sensor detected (',pH.type,' ',pH.SN,')'])
        info.pH_flag = 1;
    % SBE SECONDARY pH CALIBRATION LINE (NAVIS 0690, others?)
    elseif regexp(tline,'^pH sensor secondary','once') % SBE Calibration 
        tmp = regexp(tline,',','split');
        pH.secondary_Zlimits = str2double(tmp(2:3));
        pH.secondary_pcoefs  = str2double(tmp(4:end))';
    end
    tline = fgetl(fid);
end

if info.pH_flag == 1 % ADD CALIBRATION DATA TO STRUCTURE
    cal.pH = pH;
end

% ************************************************************************
% LOOK FOR ANY SPECIAL NOTES / COMMENTS
% THESE MUST BE AT THE BOTTOM OF THE CONFIG FILE
% ************************************************************************
frewind(fid) % go back to top of file
ct        = 0;
note_flag = 0;
tline     = ' ';% initialize
while ischar(tline)
    if regexpi(tline,'^\s*COMMENTS:', 'once')
        note_flag = 1;
        tline = fgetl(fid); % get 1st comment line
    end
    if note_flag ==1 && size(tline,2) > 3 % some comments
        ct = ct+1;
        info.notes{ct,1} = tline;
    end
    tline = fgetl(fid);
end

fclose(fid); % ALL CHECKS DONE - CLOSE FILE
clear pH fid

% ************************************************************************
% PARSE *.cal TO CHECK FOR Nitrate calibration data
% ************************************************************************
if info.isus_flag == 1
    cal.N = parseNO3cal([dirs.cal,'NO3_CAL\',NO3_calfile]);
    disp(['Nitrate Sensor detected (',cal.N.type,' SN: ',cal.N.SN,')'])
end

% ************************************************************************
% Tidy up and save the data
% ************************************************************************
if exist('cal','var') ~= 1
    disp(['NO BGC SENSORS DETECTED FOR ',MBARI_ID_str, ' - CAL FILE ', ...
        'WAS NOT BUILT']);
    cal = [];
    return
end
    
cal.info = info;

s1 = [dirs.cal,'cal',MBARI_ID_str,'.mat'];
s2 = [dirs.temp,config_file];
s3 = [dirs.temp,NO3_calfile];

if exist(s1,'file')
    disp(['File already exist! (',s1,')'])
    str = input('Overwrite file? [y/n]','s');
    if strncmpi(str,'y',1)
        save(s1, 'cal');
    end
else
     save(s1, 'cal');
end

if exist(s2,'file')
    delete(s2);
end
if exist(s3,'file')
    delete(s3);
end
disp(' ')

clear info MBARI_ID_str ind ans fid tline
%clear s1 s2 s3
    
