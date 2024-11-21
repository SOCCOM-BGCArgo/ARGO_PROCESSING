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
%       cal.O       Oxygen calibration coefficients (Aanderaa 3830&4330, 
%                   SBE63,SBE83). cal.O2, cal.O3 if multiple O2 sensors on
%                   float
%       cal.CHL     biooptics calibration coefficients (FLBBAP2 or MCOMS
%       cal.BB      biooptics calibration coefficients (FLBBAP2 or MCOMS
%       cal.CDOM    biooptics calibration coefficients (FLBBAP2 or MCOMS
%       cal.pH   	pH calibration coefficients (MBARI durafet, SBE)
%       cal.N   	ISUS or SUNA calibration coefficients
%       cal.OCR

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
% 06/03/21 TM, Enhancement in support of OCR sensor (and SBE83 optode)
% 01/20/22 JP, Updated to pull in additional O2 calibrations for floats
%              with multiple O2 sensors (3XO2: un1173,un1342, 
%              2XO2: ua19298, ua19843)
% 01/20/22 JP, Updated OCR extraction to make more uniform across plotforms for
%              SOLO, NAVIS, APEX (ie SOLO 0001, ua19314, ua19191)
% 03/14/22 JP, Udpated to extract k2 f(P) if it exists in the config file
% 04/24/23 JP, Extended k2 f(P) to NAVIS float style pH cal - it was just APEX
% 05/23/23 JP, Extended bio-optics code to pull in CHL435 (FL2BB) if it
%              exists
% 07/25/23 JP, Minor update to DSD pH "coef_ct" number extraction to allow
%              coeficient count to go to double digits (i.e. > 9)

% ************************************************************************
% FOR TESTING
% MBARI_ID_str = 'un1341';
% MBARI_ID_str = 'ua7614';
% MBARI_ID_str = 'un0412';
% MBARI_ID_str = 'ua12541';
% MBARI_ID_str = 'un0037';
% MBARI_ID_str = 'un0276';
% MBARI_ID_str = 'un0569';
%MBARI_ID_str = 'ua5143';
%MBARI_ID_str = 'un0565'; % has bbp532
%MBARI_ID_str = 'un0690';
% 
%MBARI_ID_str = 'un0510';
%MBARI_ID_str = 'wn857';
%MBARI_ID_str = 'ua0068';


%MBARI_ID_str = 'un1173'; % 3X O2 Aanderaa, SBE63, SBE83
%MBARI_ID_str = 'ua5143'; 
%MBARI_ID_str = 'ua19314'; % OCR
%MBARI_ID_str = 'ua19191'; % OCR
%MBARI_ID_str = 'ua19298'; % %2XO2 Aanderaa
%MBARI_ID_str = 'ua19843'; % %2XO2 Aanderaa
%MBARI_ID_str = 'ss0001'; % %2XO2 Aanderaa
%MBARI_ID_str = 'wn1343'; % different O2 config format
%MBARI_ID_str = 'ua21844'; % MBARI GDF with k2_f(p) & 6th order f(P)
%MBARI_ID_str = 'ua20358'; % MBARI DSD with constant k2 & 6th order f(P)
%MBARI_ID_str = 'wn1475'; % SBE WHOI NAVIS with constant k2 & 9th order f(P)
% MBARI_ID_str = 'wn1487'; % SBE WHOI NAVIS with k2f(P) & 9th order f(P)
% MBARI_ID_str = 'un0949'; % SBE WHOI NAVIS with k2f(P) & 9th order f(P)
%MBARI_ID_str = 'ua21291'; % FLBBFL with CHL435
%MBARI_ID_str = 'ua21910'; % OCR 6 sensor
% MBARI_ID_str = 'un1512'; % 
% MBARI_ID_str = 'ua9630'; % 
%MBARI_ID_str = 'un0062'; % 
% MBARI_ID_str = 'un0061'; % 
% MBARI_ID_str = 'ua21286'; % 
% dirs =[];


% ************************************************************************

% ************************************************************************
% SET DIRECTORY AND VARIABLES
% ************************************************************************
config_file  = [MBARI_ID_str,'_FloatConfig.txt'];
NO3_calfile  = [MBARI_ID_str, '.cal'];
master_list  = 'MBARI_float_list.mat';

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.cal = [user_dir,'CAL\']; % save created lists here
    dirs.mat = [user_dir,'FLOATS\'];
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

info.O2_flag      = 0;
info.chl_flag     = 0;
info.chl435_flag  = 0;
info.pH_flag      = 0;
info.isus_flag    = 1; % default is 1
info.ocr_flag     = 0;

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
if exist([dirs.cal,'NO3_CAL\',NO3_calfile],'file') ~= 2
    disp([MBARI_ID_str,': ',NO3_calfile,' not found - nitrate spectra',...
        ' will not be re-processed']);
    info.isus_flag    = 0; % NO NO3 - add to ouput structure at end
end

% INITIALIZE CAL STRUCTURE
cal = struct;

% ************************************************************************
% PARSE CONFIG FILE TO GET OXYGEN CALIBRATION COEFFICIENTS
% Examples of OptodeSn output:
% 'OptodeSn = 737'                   (Aanderaa 3830 older)
% 'OptodeSn = 4330 1168'             (Aanderaa 4330 newer)
% 'OptodeSn = SBE63 rev J S/N 0455'  (Seabird 63)
% ************************************************************************
oxy_fnames   = {'O' 'O2' 'O3' 'O4'}; % field names for multiple O2 cals
tline        = ' ';% initialize some variables
fid          = fopen([dirs.cal,'FLOAT_CONFIG\', config_file]);
if fid == -1
    disp(['Could not open *.txt config file for ',MBARI_ID_str])
    cal = [];
    return
end

% CHECK FOR O2 CALS - MORE THAN 1??
optode_info = cell(10,3);  % pre over dimmension
cell_txt    = cell(100,1); % pre over dimmension
cal_ct      = 0;
line_ct     = 0;
while ischar(tline)
    if regexp(tline,'^\s*OptodeSn =','once') % APEX
        cal_ct = cal_ct +1;
        tmp = regexp(tline, '\d+', 'match'); % parse any numbers from SN line
        if size(tmp,2) == 2 % 4330
            optode_info(cal_ct,1:2) = tmp;
            optode_info{cal_ct,3}   = 'Aanderaa';
        elseif size(tmp,2) == 1 % 3830 no model number only SN
            optode_info{cal_ct,1} = '3830';
            optode_info{cal_ct,2} = tmp{1,1};
            optode_info{cal_ct,3} = 'Aanderaa';
        end
    elseif regexp(tline,'^\s*O2 sensor\s*(','once') % NAVIS SBE63 or SBE83
        cal_ct = cal_ct +1;
        %optode_info(cal_ct,1:2) = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
        % block preceded by begin paren or block followed by end paren
        % this assumes there will always be two blocks of chars
        filt_str = '(?<=\()[\w\.]*|[\w\.]*(?=\))'; 
        optode_info(cal_ct,1:2) = regexp(tline,filt_str,'match');
        optode_info{cal_ct,3}   = 'SBS';
    end
    
    line_ct           = line_ct + 1;
    cell_txt{line_ct} = tline; 
    tline             = fgetl(fid);
end

optode_info  = optode_info(1:cal_ct,:);
[~,ia,~]     = unique(strcat(optode_info(:,1),optode_info(:,2)), 'stable');
optode_info  = optode_info(ia,:);   % Unique calibrations
cell_txt     = cell_txt(1:line_ct); % cell array of config file lines
info.O2_flag = size(optode_info,1);

if info.O2_flag > 0
    for O2_ct = 1:info.O2_flag
        cal_field            = oxy_fnames{O2_ct};
        if isfield(cal, cal_field)
            cal = rmfield(cal, cal_field);
        end
        cal.(cal_field).type = optode_info{O2_ct,1};
        cal.(cal_field).SN   = optode_info{O2_ct,2};
        fprintf('%s oxygen optode detected  Model: %s SN: %s\n', ...
            optode_info{O2_ct,3},optode_info{O2_ct,1},optode_info{O2_ct,2});

        
        if regexp(optode_info{O2_ct,1},'^SBE|^4330') % SBE63, SBE83, 4330
            fmt  = [optode_info{O2_ct,1},'\s+',optode_info{O2_ct,2}];
            tCAL = ~cellfun(@isempty,regexp(cell_txt, fmt,'once'));
            tmp  = cell_txt(tCAL); % Config file O2 cal sub set as cell array

            for ct = 1:size(tmp,1)
                str = tmp{ct};
                
                % ****  Aanderaa 4330  ****
                if ~isempty(regexp(str,'PhaseCoef','once'))
                    ind = regexp(str,'PhaseCoef\s+\d+\s+\d+\s+','once','end');
                    cal.(cal_field).PCoef = sscanf(str(ind:end),'%f',4);
                    
                elseif ~isempty(regexp(str,'FoilCoefA','once'))
                    ind = regexp(str,'FoilCoefA\s+\d+\s+\d+\s+','once','end');
                    FCoefA = sscanf(str(ind:end),'%f',14);
                    
                elseif ~isempty(regexp(str,'FoilCoefB','once'))
                    ind = regexp(str,'FoilCoefB\s+\d+\s+\d+\s+','once','end');
                    FCoefB = sscanf(str(ind:end),'%f',14);
                    
                elseif ~isempty(regexp(str,'FoilPolyDegT','once'))
                    ind = regexp(str,'FoilPolyDegT\s+\d+\s+\d+\s+','once','end');
                    cal.(cal_field).PolyDegT = sscanf(str(ind:end),'%f',28);
                    
                elseif ~isempty(regexp(str,'FoilPolyDegO','once'))
                    ind = regexp(str,'FoilPolyDegO\s+\d+\s+\d+\s+','once','end');
                    cal.(cal_field).PolyDegO = sscanf(str(ind:end),'%f',28);
                    
                elseif ~isempty(regexp(str,'SVUFoilCoef','once'))
                    ind = regexp(str,'SVUFoilCoef\s+\d+\s+\d+\s+','once','end');
                    SVU = sscanf(str(ind:end),'%f',7);
                    if sum(SVU(:)) ~= 0
                        cal.(cal_field).SVUFoilCoef = SVU;
                    end
                    
                elseif ~isempty(regexp(str,'ConcCoef','once'))
                    ind = regexp(str,'ConcCoef\s+\d+\s+\d+\s+','once','end');
                    cal.(cal_field).ConcCoef = sscanf(str(ind:end),'%f',2);
                    
                    % ****  SBE 63 or SBE83 SEARCH  ****
                elseif regexp(str,'^O2 sensor.+Temp coefficents', 'once')
                    TC = regexp(str,',','split');
                    cal.(cal_field).TempCoef = str2double(TC(2:end)); % TA0-3
                    
                elseif regexp(str,'^O2 sensor.+Phase coefficents', 'once')
                    PC = regexp(str,',','split');
                    cal.(cal_field).PhaseCoef = str2double(PC(2:end)); % TA0-3
                end
            end
            
            if exist('FCoefA','var') && exist('FCoefB','var') % 4330 only
                cal.(cal_field).FCoef =[FCoefA;FCoefB];
                % Look for zero & remove to shorten calc, 0 * x = 0
                t1 = cal.(cal_field).FCoef == 0;
                cal.(cal_field).FCoef(t1)    = [];
                cal.(cal_field).PolyDegT(t1) = [];
                cal.(cal_field).PolyDegO(t1) = [];
            end
            clear FCoefA FCoefB  ind t1 PC TC SVU ct tmp tCAL fmt
            
        elseif regexp(optode_info{O2_ct,1},'^3830') % 3830
            fmt  = '^C[0-4][0-3]|^P[01]'; % 3830 coef line starts
            tCAL = ~cellfun(@isempty,regexp(cell_txt, fmt,'once'));
            tmp  = cell_txt(tCAL); % Config file 3830 O2 cal sub set as cell array
            
            coef_names = regexp(tmp,'^C\d+|^P\d+','match','once');
            coefs_vals = regexp(tmp,'(?<=\=\s+)[\d-+\.Ee]+','match','once');
            coefs_vals = str2double(coefs_vals);
            
            %BUILD COEFF ARRAYS for POLYVAL
            coef_srch = {'C0', 'C1', 'C2', 'C3', 'C4', 'P'};
            
            for i = 1:length(coef_srch) % build coeff arrays with eval
                t1 = strncmp(coef_srch{i}, coef_names, length(coef_srch{i}));
                % Aanderaa = ascending powers, Matlab needs descending powers
                cal.(cal_field).(['p',coef_srch{i}]) = flip(coefs_vals(t1),1);
            end
        end
        clear coef_names coef_vals tmp tCAL fmt ind1 i t1
    end
end
    
% ************************************************************************
% PARSE config FILE TO CHECK FOR FLBB/MCOMS CHL AND BACKSCATTER AND CDOM
% CAL COEFFICIENTS
% ************************************************************************
frewind(fid)
tline = ' ';% initialize some variables
while ischar(tline)
    % **************************************
    % FLBB / FLBBFL BLOCK
    if regexp(tline,'FLBBFL','once') % FLBBFL CHL435 SN
        CHL435.type = regexp(tline,'FLBB\w+','match','once');
        CHL435.SN   = regexp(tline,'\d+$','match','once');
        CHL.type    = CHL435.type;
        CHL.SN      = CHL435.SN;
        disp(['FLBBFL with CHL435 Bio-optics detected (', ...
            CHL435.type,' ',CHL435.SN,')'])
        info.chl435_flag   = 1;
        info.chl_flag   = 1;
    elseif regexp(tline,'FLBB','once') % FLBB SN
        CHL.type = regexp(tline,'FLBB\w+','match','once');
        CHL.SN   = regexp(tline,'\d+$','match','once');
        disp(['FLBB Bio-optics detected (',CHL.type,' ',CHL.SN,')'])
        info.chl_flag   = 1;

    elseif regexp(tline,'^\d+.+ChlDC','once') % FLBB CHL DC
        CHL.ChlDC = sscanf(tline,'%f',1);
    elseif regexp(tline,'^\d+.+ChlScale','once') % FLBB CHL DC
        CHL.ChlScale = sscanf(tline,'%f',1);

    elseif regexp(tline,'^\d+.+Chl435DC','once') % FLBB CHL DC
        CHL435.ChlDC = sscanf(tline,'%f',1);
    elseif regexp(tline,'^\d+.+Chl435Scale','once') % FLBB CHL DC
        CHL435.ChlScale = sscanf(tline,'%f',1);

    elseif regexp(tline,'^\d+.+BetabDC','once') % FLBB BBP DC
        BB.type    = CHL.type;
        BB.SN      = CHL.SN;
        BB.BetabDC = sscanf(tline,'%f',1);
    elseif regexp(tline,'^\d+.+BetabScale','once') % FLBB CHL DC
        BB.BetabScale = sscanf(tline,'%f',1);


        % **************************************
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
        
        % ECO TRIPLET BLOCK
    elseif regexp(tline,'^ECO.+Chl','once') % MCOMS  CHL
        tmp          = regexp(tline,',','split');
        CHL.type     = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        CHL.SN       = regexp(tmp{1,1},'\d+(?=\))','once','match');
        CHL.ChlDC    = str2double(tmp{1,2});
        CHL.ChlScale = str2double(tmp{1,3});
        disp(['O Bio-optics detected (',CHL.type,' ',CHL.SN,')'])
        info.chl_flag   = 1;
    elseif regexp(tline,'^ECO.+Backscatter700','once') % MCOMS BBP700
        tmp           = regexp(tline,',','split');
        BB.type       = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        BB.SN         = regexp(tmp{1,1},'\d+(?=\))','once','match');
        BB.BetabDC    = str2double(tmp{1,2});
        BB.BetabScale = str2double(tmp{1,3});
    elseif regexp(tline,'^ECO.+DOM','once')
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
    if exist('CHL435', 'var') % NAVIS
        cal.CHL435 = CHL435;
    end
end

clear tmp tline t1 BB CHL CDOM

% ************************************************************************
% PARSE Config FILE TO CHECK FOR pH COEFFICIENTS
% ************************************************************************
frewind(fid) % go back to top of file
tline = ' ';% initialize some variables
while ischar(tline)
    % MBARI pH CALIBRATION BLOCK
    if regexp(tline,'^Durafet SN','once') % MBARI Calibration
        pH.type = regexp(tline,'^\w+','match','once');
        pH.SN   = regexp(tline,'(?<=\=\s+)\w+','match','once'); % # after =
        disp(['MBARI pH sensor detected (',pH.type,' ',pH.SN,')'])
        info.pH_flag = 1;
    elseif regexp(tline,'^\d.+number\s+of\s+calibration','once')% # MBARI pH cal lines
        coef_ct = str2double(regexp(tline,'^\d+','once','match'));
        %tmp = ones(coef_ct-2,1)*NaN;
        pH.pcoefs = ones(coef_ct-2,1)*NaN; % -2 for k0 & k2
        for i = 1:coef_ct
            tline = fgetl(fid);
            if ~isempty(regexp(tline,'\=k0','once')) % k0 line
                pH.k0 = sscanf(tline,'%f',1);
            elseif ~isempty(regexp(tline,'\=k2','once')) % k2 line
                k2_tmp = regexp(tline,',','split');
                pH.k2  = str2double(k2_tmp(1:end-1))';
                if size(pH.k2,1)>1
                    fprintf('k2 as a function of pressure detected\n')
                end
            else
                pH.pcoefs(i-2) = sscanf(tline,'%f',1);
            end
        end

        
    % SBE PRIMARY pH CALIBRATION BLOCK
    elseif regexp(tline,'^pH sensor coefficients','once') % SBE Calibration
        tmp         = regexp(tline,',','split');
        pH.type     = regexp(tmp{1,1},'(?<=\()\w+','once','match');
        pH.SN       = regexp(tmp{1,1},'\w+(?=\))','once','match');

        coef_id_str = regexp(tmp{1,1},'(?<=\[).+(?=\])','match','once');
        coef_id     = regexp(coef_id_str,'[\w\^]+','match');
        tK0         = strcmp(coef_id,'k0');
        tK2         = strncmp(coef_id,'k2',2);
        tFP         = ~(tK0|tK2);
        coef_vals   = str2double(tmp(2:end));

        % SANITY CHECK: ID & VAL SIZES MUST BE EQUAL
        if size(coef_id,2) ~=  size(coef_vals,2)
            fprintf(['WARNING: pH coefficient ID count (%0.0f) does not match ', ...
                'coefficent value count (%0.0f) - please check %s config file ', ...
                'for typos\n'],size(coef_id,2), size(coef_vals,2), MBARI_ID_str)
        end

        pH.k0       = coef_vals(tK0);
        pH.k2       = coef_vals(tK2)';
        pH.pcoefs   = coef_vals(tFP)';

        disp(['SBE pH sensor detected (',pH.type,' ',pH.SN,')'])
        info.pH_flag = 1;


%         pH.k0     = str2double(tmp(2));
%         pH.k2     = str2double(tmp(3));
%         pH.pcoefs = str2double(tmp(4:end))';
%         disp(['SBE pH sensor detected (',pH.type,' ',pH.SN,')'])
%         info.pH_flag = 1;

        % SBE SECONDARY pH CALIBRATION LINE (NAVIS 0690, others?)
        % this code may require updating if new formats get secondary
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
% CHECK FOR OCR SENSOR CAL COEFFICIENTS - USE CELL ARRAY OF TEXT LINES
%                             JP 01/20/2021
% ************************************************************************
fmt  = '^\s*OCR';
tCAL = ~cellfun(@isempty,regexp(cell_txt, fmt,'once'));
if sum(tCAL,1) > 0
    if isfield(cal, 'OCR')
        cal = rmfield(cal, 'OCR');
    end
    tmp  = cell_txt(tCAL); % Config file OCR cal sub set as cell array
    wls = {};
    for i = 1:size(tmp,1) % Step through OCR Channels
        str = tmp{i};
        model_sn = regexp(str,'(?<=\()[\w\s]+(?=\))','match', 'once');
        if i == 1
            cal.OCR.type = regexp(model_sn,'OCR\d+','match','once');
            cal.OCR.SN   = regexp(model_sn,'\d+$','match','once');
        end 
        chan_str = ['CH',regexp(str,'(?<=CHANNEL\s*)\d+','match', 'once')];
        WL = regexp(str,'(?<=CHANNEL\s+\d{2}\s+)\w+','match','once');
        wls = [wls;WL];
        
        cal.OCR.(chan_str).WL = WL;
        coef_str = regexp(str,'(?<=\],).+','match','once');
        ocr_coefs = str2double(regexp(coef_str,',','split'));
        cal.OCR.(chan_str).a0 = ocr_coefs(1);
        cal.OCR.(chan_str).a1 = ocr_coefs(2);
        cal.OCR.(chan_str).im = ocr_coefs(3);
    end
    wl_str = sprintf('%s ',wls{:});
    fprintf('%s %s radiometer detected with %0.0f channels: %s\n', ...
        cal.OCR.type, cal.OCR.SN, size(tmp,1),wl_str);
    info.ocr_flag = 1;
end
clear tmp tline chan_str coef_str i ocr_coefs WL wls wl_str

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
if exist('cal','var') ~= 1 || isempty(fieldnames(cal))
    disp(['NO BGC SENSORS DETECTED FOR ',MBARI_ID_str, ' - CAL FILE ', ...
        'WAS NOT BUILT']);
    cal = [];
    return
end

cal.info = info;

s1 = [dirs.cal,'cal',MBARI_ID_str,'.mat'];
% s1 = [dirs.cal,'TEST\cal',MBARI_ID_str,'.mat']; %  **** TESTING ****
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

