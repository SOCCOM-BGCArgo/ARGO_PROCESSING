function cal = get_float_cals(MBARI_ID_str, dirs)
% PURPOSE: 
%   This function tries to recover all the calibration data for the
%   biogoechemical sensors on a given float. If a sensor is detected in
%   the config or cal text files a sturcture field with calibration data
%   will be created
%
% USAGE:
%	data = get_float_cals(UW_ID_str, dirs)
%
% INPUTS:
%   MBARI_ID_str  = MBARI Float ID as a string
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
% REQUIRES: get_QCstep_dates
%           get_MBARI_WMO_list
%           
% EXAMPLES:
%   jp = get_float_cals('8501SOOCN', dirs)
%   jp = get_float_cals('8514SOOCN', [])
%
% Created 12/29/2015 by jp
% CHANGES
% 07/19/17 Added functionality to extract comments line at bottom of
%   *FloatConfig.txt files if they exist. These comments can then be
%   printed to the ODV txt files in argo2ODV*. m files
% 09/10/18 Added code to look for seconday pH pressure coefficients. At
%   this point this change onlly afects 0690. -jp
% 10/26/18 Added code to extract non zero Aanderaa 4330 SVUFoilCoef's. -jp


% ************************************************************************
% FOR TESTING
% dirs =[];
% MBARI_ID_str = '9125SoOcn';
%MBARI_ID_str = '8501SoOcn';
% MBARI_ID_str = '12747SoOcn';
% MBARI_ID_str = '0690SOOCN';
% %MBARI_ID_str = '0569SoOcn';
% MBARI_ID_str = '7614SoOcn';

% ************************************************************************

% ************************************************************************
% SET DIRECTORY AND VARIABLES (Once rubust make some function inputs)
% ************************************************************************
UW_ID_str = regexp(MBARI_ID_str,'^\d+','once','match');
config_file  = [MBARI_ID_str,'_FloatConfig.txt'];
% APEX_config_file  = [MBARI_ID_str,'_FloatConfig.txt'];
% NAVIS_config_file = [MBARI_ID_str,'_FloatConfig.txt'];
NO3_calfile       = [MBARI_ID_str, '.cal'];

% APEX_config_file  = 'FloatConfig.txt';
% NAVIS_config_file = ['NAVISConfig_',UW_ID_str,'.txt'];
% NO3_calfile       = [UW_ID_str, '.cal'];

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    %dirs.NO3config = [user_dir,'CAL\'];
    dirs.FVlocal   = [user_dir,'FLOATVIZ\'];
    dirs.FV        = [user_dir,'FLOATVIZ\'];
    %dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    
    dirs.QCadj     = [user_dir,'CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    %dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    
    dirs.log = [user_dir,'Processing_logs\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

info.isus_flag  = 1; % 0 = no isus cal file
info.chl_flag   = 1;
info.pH_flag    = 1;
info.O2_flag    = 1;

info.WMO_ID     = '';
info.file_date  = datestr(now,'mm/dd/yyyy HH:MM');
info.float_type = 'APEX';

clear cal % clear any pre-exisiting structure

clear tmp r choose_dir
disp(['Building calibration file for ',MBARI_ID_str]);
% ************************************************************************
% CHECK IF FILES EXIST: FLOAT CONFIGURATION   & NO3 CAL FILES
% ************************************************************************
if exist([dirs.cal,'FLOAT_CONFIG\', config_file],'file') ~= 2
    disp(['Calibration *.txt file not found for ',MBARI_ID_str, ...
        ' - stopping script']);
    cal = [];
    return
end

% GET ISUS CAL FILE LOCALLY FOR NOW - THIS IS A PLCEHOLDER
if exist([dirs.cal,'NO3_CAL\',NO3_calfile],'file') ~= 2 % copy file
    disp([MBARI_ID_str,': ',NO3_calfile,' not found - nitrate spectra',...
        ' will not be re-processed']);
    info.isus_flag    = 0; % add to ouput structure at end
end

% ************************************************************************
% GET WMO # FROM KEN's WEB PAGE HTML FILE(IS THERE A BETTER WAY?)
% ADD UW ID TOO and MBARI FLOATVIZ NAME
% COULD ALSO GET FOR UW WEB PAGE
% ************************************************************************

if exist([dirs.cal,'MBARI_float_list.mat'], 'file') %load list variables
    load([dirs.cal,'MBARI_float_list.mat']);
else
    disp('BUILDING FLOAT ID LIST ...')
    list = MBARI_float_list(dirs); %float_names, UW_ID, WMO_ID, float_type;
end
float_names = list(:,1);
UW_ID       = list(:,2);
WMO_ID      = list(:,3);

ind = find(strcmpi(MBARI_ID_str,float_names) == 1);
if ~isempty(ind)
    info.WMO_ID = WMO_ID{ind};
else
    disp(['NO WMO # FOUND FOR',MBARI_ID_str,' UPDATINIG FLOAT ID LIST....'])
    %[float_names, UW_ID, WMO_ID] = get_MBARI_WMO_list(dirs);
    list = MBARI_float_list(dirs); %float_names, UW_ID, WMO_ID, float_type;
    float_names = list(:,1);
    UW_ID       = list(:,2);
    WMO_ID      = list(:,3);

    ind = find(strcmp(MBARI_ID_str,float_names) == 1);
    if  ~isempty(ind)
        info.WMO_ID = WMO_ID{ind};
    else
        disp('NO WMO # FOUND!')
    end
end

info.name = MBARI_ID_str;
info.UW_ID = UW_ID_str;

clear fid tline list

% ************************************************************************
% PARSE CONFIG FILE TO GET OXYGEN CALIBRATION COEFFICIENTS
% Examples of OptodeSn output:
% 'OptodeSn = 737'                   (Aanderaa 3830 older)
% 'OptodeSn = 4330 1168'             (Aanderaa 4330 newer)
% 'OptodeSn = SBE63 rev J S/N 0455'  (Seabird 63)
% ************************************************************************
tline = ' ';% initialize some variables
fid   = fopen([dirs.cal,'FLOAT_CONFIG\', config_file]);
if fid == -1
    disp(['Could not open *.txt config file for ',MBARI_ID_str])
    cal = [];
    return
end

while ischar(tline)
    if regexp(tline,'^NAVIS')
        info.float_type = 'NAVIS';
    end
        
    if regexp(tline,'OptodeSn =','once')
        optope_info = regexp(tline, '\d+', 'match'); % parse any numbers from SN line
        break
    elseif regexp(tline,'^O2 sensor','once') % NAVIS SBE63- get MODEL and SN
        optope_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
        break
    end
    tline = fgetl(fid);
end

if ~ischar(tline) % -1 (end of file reached - no O2)
    disp(['NO oxygen detected for float ',UW_ID_str])
    info.O2_flag = 0; % add to ouput structure at end
end

% ************************************************************************
% NEWER STYLE AANDERAA OXYGEN CALIBRATION FORMAT (4330)    p02 = f(Phase)

if length(optope_info) == 2 && ...
        ~isempty(regexp(optope_info{1},'4330','once'))
    O.type = optope_info{1};
    O.SN   = optope_info{2};
    disp(['Oxygen optode detected on APEX float (',O.type,' SN: ',O.SN,')'])
    while tline ~= -1 % EXTRACT COEFFICIENTS
        % comments at end stop searching
        if ~isempty(regexp(tline,'^COMMENTS:', 'once')) 
            break
        end 
        if ~isempty(regexp(tline,'PhaseCoef','once'))
            ind = regexp(tline,'PhaseCoef','once')+ 19;
            O.PCoef = cell2mat((textscan(tline(ind:end),'%f',4)));
        end
        if ~isempty(regexp(tline,'FoilCoefA','once'))
            ind = regexp(tline,'FoilCoefA','once')+ 19;
            FCoefA = cell2mat((textscan(tline(ind:end),'%f',14)));
        end
        if ~isempty(regexp(tline,'FoilCoefB','once'))
            ind = regexp(tline,'FoilCoefB','once')+ 19;
            FCoefB = cell2mat((textscan(tline(ind:end),'%f',14)));
        end
        if ~isempty(regexp(tline,'FoilPolyDegT','once'))
            ind = regexp(tline,'FoilPolyDegT','once')+ 22;
            O.PolyDegT = cell2mat((textscan(tline(ind:end),'%f',28)));
        end
        if ~isempty(regexp(tline,'FoilPolyDegO','once'))
            ind = regexp(tline,'FoilPolyDegO','once')+ 22;
            O.PolyDegO = cell2mat((textscan(tline(ind:end),'%f',28)));
        end
        
        if ~isempty(regexp(tline,'SVUFoilCoef','once'))
            ind = regexp(tline,'SVUFoilCoef \d+ \d+ ','once','end')+1;
            SVU = sscanf(tline(ind:end),'%f');
            if sum(SVU(:)) ~= 0
                O.SVUFoilCoef = SVU;
            end
        end
        
        if ~isempty(regexp(tline,'ConcCoef','once'))
            ind = regexp(tline,'ConcCoef','once')+ 19;
            O.ConcCoef = cell2mat((textscan(tline(ind:end),'%f',2)));
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
% APEX SBE63 OXYGEN CALIBRATION FORMAT
    
elseif strcmp(info.float_type, 'APEX') && length(optope_info) == 2 &&  ...
        ~isempty(regexp(tline,'SBE63','once'))
    O.type = ['SBE',optope_info{1}];
    O.SN   = optope_info{2};
    disp(['Oxygen optode detected on APEX float (',O.type,' SN: ',O.SN,')'])
    
    A =[]; B =[]; C =[]; TA =[]; SOLB =[];
    
    while tline ~= -1 % EXTRACT COEFFICIENTS
        % comments at end stop searching
        if ~isempty(regexp(tline,'^COMMENTS:', 'once'))
            break
        end
        ind = regexp(tline,'Sbe63LogCoef\(\)','once','end');
        if ~isempty(ind) % coef lines
            coef = textscan(tline(ind+1:end),'%s %*s %f');
            if strncmp('A',coef{1},1)
                A      = [A,coef{2}];
            elseif strncmp('B',coef{1},1)
                B      = [B,coef{2}];
            elseif strncmp('C',coef{1},1)
                C      = [C,coef{2}];
            elseif strncmp('E',coef{1},1)
                E      = coef{2};
            elseif strncmp('TA',coef{1},2)
                TA     = [TA,coef{2}];
            elseif strncmp('SOLB',coef{1},4)
                SOLB   = [SOLB,coef{2}];
            elseif strncmp('SOLC',coef{1},4)
                SOLC   = coef{2};
            elseif strncmp('RefSal',coef{1},6)
                RefSal = coef{2};
            elseif strncmp('RefP',coef{1},4)
                RefP   = coef{2};
            end
        end
        tline = fgetl(fid);
    end
    O.TempCoef = TA;
    O.PhaseCoef = [A,B,C]; 
    % Can add more to 'O' later if desired
    clear A B C E TA SOLB SOLC RefSal Refp
    clear ind coef  
    
% ************************************************************************
% OLD STYLE AANDERAA OXYGEN CALIBRATION FORMAT (3830)    02 = f(Phase)  
%elseif length(optope_info) == 1 && ~isempty(regexp(tline,'3830','once'))
elseif length(optope_info) == 1 
    O.type = '3830';
    O.SN   = optope_info{1};
    disp(['Oxygen optode detected on APEX float(',O.type,' SN: ',O.SN,')'])
    coef = textscan(fid,'%s %*s %f'); %{1}= ID, {2}= value
    %BUILD COEFF ARRAYS for POLYVAL
    coef_srch = {'C0', 'C1', 'C2', 'C3', 'C4', 'P'};
    
    for i = 1:length(coef_srch) % build coeff arrays with eval
        % Aanderaa = ascending powers
        ind1 = find(strncmp(coef_srch{i},coef{1},length(coef_srch{i})));
        ind1 = sort(ind1,'descend')'; % Matlab needs descending powers
        eval(['O.p',coef_srch{i},'=coef{2}(ind1);']);
    end;
    clear coef_format coef_path coef_srch calfile ind1 i opt_info  
    
% ************************************************************************
% NAVIS SBE63 OXYGEN CALIBRATION FORMAT
elseif strcmp(info.float_type, 'NAVIS') && strncmp(optope_info{1},'SBE63',5)
    O.type = optope_info{1};
    O.SN   = optope_info{2};
    disp(['Oxygen optode detected on NAVIS float(',O.type,' SN: ',O.SN,')']) 
    
    while regexp(tline,'^O2 sensor', 'once')
        if regexp(tline,'Temp coefficents', 'once')
            TC = textscan(tline,'%*s%f%f%f%f', 'Delimiter', ',', ...
                'collectoutput',1);
            O.TempCoef = TC{1}; % TA0-3
        elseif regexp(tline,'Phase coefficents', 'once')
            PC = textscan(tline,'%*s%f%f%f%f%f%f%f%f', 'Delimiter', ',',...
                'collectoutput',1);
            O.PhaseCoef = PC{1}; % A0-2 B0-1 C0-2
        end
        tline =fgetl(fid);
    end
    clear TC PC
else
    disp(['Can not determine oxygen optode type for ', UW_ID_str])
    info.O2_flag = 0;
end

if info.O2_flag == 1
    cal.O = O;
end
%clear O optode_info
% ************************************************************************
% PARSE config FILE TO CHECK FOR APEX CHL AND BACKSCATTER AND
% CDOM CAL COEFFICIENTS
% ************************************************************************
frewind(fid)
tline = ' ';% initialize some variables

if strcmp(info.float_type,'APEX')
    while ischar(tline) && isempty(regexp(tline,'CHLFLUOR','once'))
        tline = fgetl(fid);
    end
    
    if ~ischar(tline) % -1 (end of file reached - no CHL)
        disp(['No chlorophyll fluorescence or backscatter detected for', ...
            ' float ',UW_ID_str])
        info.chl_flag = 0; % add to ouput structure at end
    else
        ind = regexp(tline,'FLBB|MCOMS', 'once');
        chl_info = regexp(tline(ind:end), '\w+', 'match'); % model & SN
        OPT_type = chl_info{1};
        OPT_SN   = chl_info{2};
        disp(['Bio-optics detected (',tline(ind:end),')'])
    end
    
    if info.chl_flag == 1
        while ischar(tline)
            tline = strtrim(tline);
            % comments at end stop searching
            if ~isempty(regexp(tline,'^COMMENTS:', 'once'))
                break
            end
            if ~isempty(regexp(tline,'ChlDC','once'))
                CHL.ChlDC = str2double(strtok(tline));
                CHL.type  = OPT_type;
                CHL.SN    = OPT_SN;
            end
            if ~isempty(regexp(tline,'ChlScale','once'))
                CHL.ChlScale = str2double(strtok(tline));
            end
            if ~isempty(regexp(tline,'BetabDC','once'))
                BB.BetabDC = str2double(strtok(tline));
                BB.type    = OPT_type;
                BB.SN      = OPT_SN;
            end
            if ~isempty(regexp(tline,'BetabScale','once'))
                BB.BetabScale = str2double(strtok(tline));
            end
            if ~isempty(regexp(tline,'CDOMDC','once'))
                CDOM.CDOMDC = str2double(strtok(tline));
                CDOM.type   = OPT_type;
                CDOM.SN     = OPT_SN;
            end
            if ~isempty(regexp(tline,'CDOMScale','once'))
                CDOM.CDOMScale = str2double(strtok(tline));
            end
            
            tline = fgetl(fid);
        end
    end
    
elseif strcmp(info.float_type,'NAVIS')
    while ischar(tline) && isempty(regexp(tline,'^MCOMS','once'))
        % comments at end stop searching
        if ~isempty(regexp(tline,'^COMMENTS:', 'once'))
            tline = -1; % set to numeric so it fails next test
            break
        end
        tline = fgetl(fid);
    end

    if ~ischar(tline) % -1 (end of file reached - no CHL)
        disp(['No chlorophyll fluorescence or backscatter detected for', ...
            ' float ',UW_ID_str])
        info.chl_flag = 0; % add to ouput structure at end
    elseif info.chl_flag == 1;
        chl_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
        CHL.type = chl_info{1};
        CHL.SN   = chl_info{2};
        
        while regexp(tline,'^MCOMS', 'once')
            if regexp(tline,'Chl', 'once')
                disp(['CHL Bio-optics detected (',CHL.type,' ',CHL.SN,')'])
                chl = textscan(tline,'%*s%f%f', 'Delimiter', ',', ...
                    'collectoutput',1);
                chl = chl{1};
                CHL.ChlDC    = chl(1);
                CHL.ChlScale = chl(2);
            end
            if regexp(tline,'Backscatter700', 'once')
                if ~exist('BB','var') == 1
                    bb_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
                    BB.type = bb_info{1};
                    BB.SN   = bb_info{2};
                    disp(['Backscatter700 Bio-optics detected (', ...
                        BB.type,' ',BB.SN,')'])
                end
                bscat = textscan(tline,'%*s%f%f', 'Delimiter', ',', ...
                        'collectoutput',1);                  
                BB.BetabDC = bscat{1}(1);
                BB.BetabScale = bscat{1}(2);
            end
            % 3rd MCOMS channel can be CDOM or BBP532
            if regexp(tline,'DOM', 'once')
                if ~exist('CDOM','var') == 1
                    cdom_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
                    CDOM.type = cdom_info{1};
                    CDOM.SN   = cdom_info{2};
                    disp(['CDOM Bio-optics detected (', ...
                        CDOM.type,' ',CDOM.SN,')'])
                end
                cdom = textscan(tline,'%*s%f%f', 'Delimiter', ',', ...
                        'collectoutput',1);                  
                CDOM.CDOMDC    = cdom{1}(1);
                CDOM.CDOMScale = cdom{1}(2);
            end 
            if regexp(tline,'Backscatter532', 'once') % !!!! TEMPORARY FIX
                if ~exist('CDOM','var') == 1
                    cdom_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
                    CDOM.type = cdom_info{1};
                    CDOM.SN   = cdom_info{2};
                    disp(['Backscatter532 Bio-optics detected (', ...
                        CDOM.type,' ',CDOM.SN,')'])
                end
                cdom = textscan(tline,'%*s%f%f', 'Delimiter', ',', ...
                        'collectoutput',1);                  
                CDOM.CDOMDC    = cdom{1}(1);
                CDOM.CDOMScale = cdom{1}(2);
            end 

             tline =fgetl(fid);
             if ~ischar(tline)
                 break
             end

        end
    end
        clear chl_info chl bb_info bscat cdom_info cdom
end
 
% ADD CAL DATA TO THE STRUCTURE (IF IT EXISTS)
if info.chl_flag == 1;   
    % ADD CALIBRATION DATA TO STRUCTURE
    cal.CHL = CHL;
    cal.BB = BB;
    if exist('CDOM', 'var') % NAVIS
        cal.CDOM = CDOM;
    end
end

clear  chl_info OPT_type OPT_SN ind BCHL BB CDOM tline bb_info bscat

% ************************************************************************
% PARSE Config FILE TO CHECK FOR pH COEFFICIENTS
% ************************************************************************
frewind(fid) % go back to top of file
tline = ' ';% initialize some variables

if strcmp(info.float_type,'APEX')
    while ischar(tline) && isempty(regexp(tline,'Durafet','once'))
        % comments at end stop searching
        if ~isempty(regexp(tline,'^COMMENTS:', 'once'))
            tline = -1; % set to numeric so it fails next test
            break
        end
        tline = fgetl(fid);
    end
    
    if ~ischar(tline) % -1 (end of file reached - no pH)
        disp(['No pH sensor detected in cal file for float ',UW_ID_str])
        info.pH_flag = 0; % add to ouput structure at end
    else
        pH.type = 'Durafet';
        pH.SN   = regexp(tline,'(?<=Durafet\s*SN\s*=\s*)\d+', ...
                         'match', 'once');
        %ph_info  = textscan(tline,'%s %*s %*s %s','CollectOutput',1);
        %pH.type  = ph_info{1}{1};
        %pH.SN    = ph_info{1}{2};
        
        disp(['pH detected (',pH.type,' ',pH.SN,')'])
    end
    
    if info.pH_flag == 1
        tline = fgetl(fid); %this line should give coeff count
        count_str = regexp(tline,'\d+(?=,\s+number)','once','match');
        if ~isempty(count_str)
            coef_ct = str2double(count_str);
            pH.pcoefs =[];
            for i = 1:coef_ct
                tline = strtrim(fgetl(fid));
                str_ct = regexp(tline,',','once'); % look for first comma
                if ~isempty(str_ct)
                    num_str = tline(1:str_ct);
                else % no coefficent label
                    num_str = tline;
                end
                if i == 1
                    pH.k0 = str2double(num_str);
                elseif i == 2
                    pH.k2 = str2double(num_str);
                else
                    pH.pcoefs = [pH.pcoefs; str2double(num_str)];
                end
            end
        else
            disp(['Could not determine coefficient count for float ',UW_ID_str])
            info.pH_flag = 0; % add to ouput structure at end
        end
        clear count_str coef_ct str_ct num_str i
    end
    
elseif strcmp(info.float_type,'NAVIS')
    while ischar(tline) && isempty(regexp(tline,'^pH sensor','once'))
        % comments at end stop searching
        if ~isempty(regexp(tline,'^COMMENTS:', 'once'))
            tline = -1; % set to numeric so it fails next test
            break
        end
        tline = fgetl(fid);
    end
    
    if ~ischar(tline) % -1 (end of file reached - no pH)
        disp(['No pH sensor detected in cal file for float ',UW_ID_str])
        info.pH_flag = 0; % add to ouput structure at end
    elseif info.pH_flag == 1;
        ph_info = regexp(tline,'(?<=\()\w+|\w+(?=\))','match');
        pH.type = ph_info{1};
        pH.SN   = ph_info{2};
        disp(['pH detected (',pH.type,' ',pH.SN,')'])
        
        % Find 1st comma (data starts after that)
        comma_ct = regexp(tline,',','once');
        ph_tmp   = textscan(tline(comma_ct+1:end),'%f','Delimiter',',', ...
            'CollectOutput',1);
        ph_tmp    = ph_tmp{1}; % cell to number array
        pH.k0     = ph_tmp(1);
        pH.k2     = ph_tmp(2);
        pH.pcoefs = ph_tmp(3:end);
        
        % OK NOW CHECK FOR SCONDARY PCOEF LINE - will be next line
        tline = fgetl(fid);
        if ischar(tline) && ~isempty(regexp(tline,'^pH sensor secondary','once'))
            comma_ct = regexp(tline,',','once');
            ph_tmp   = textscan(tline(comma_ct+1:end),'%f','Delimiter',',', ...
                'CollectOutput',1);
            ph_tmp     = ph_tmp{1}; % cell to number array
            pH.secondary_Zlimits = [ph_tmp(1) ph_tmp(2)];
            pH.secondary_pcoefs = ph_tmp(3:end);
        end
    end
    clear ph_info comma_ct ph_tmp
end

if info.pH_flag == 1; % ADD CALIBRATION DATA TO STRUCTURE
    cal.pH = pH;
end

% ************************************************************************
% LOOK FOR ANY SPECIAL NOTES / COMMENTS AT BOTTOM OF CAL FILE
% ************************************************************************
frewind(fid) % go back to top of file
ct = 0;
note_flag = 0;
tline = ' ';% initialize some variables
while ischar(tline)
    tline = strtrim(tline);
    if regexpi(tline,'^COMMENTS:', 'once')
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
cal.info = info;


% s1 = [dirs.cal,'cal',cal.info.UW_ID,'.mat'];
% s2 = [dirs.temp,APEX_config_file];
% s3 = [dirs.temp,NO3_calfile];

s1 = [dirs.cal,'cal',MBARI_ID_str,'.mat'];
s2 = [dirs.temp,config_file];
s3 = [dirs.temp,NO3_calfile];

if exist(s1,'file')
    disp(['File already exist! (',s1,')'])
    str = input('Overwrite file? [y/n]','s');
    if strncmpi(str,'y',1)
        save(s1, 'cal')
    end
else
    save(s1, 'cal')
end

if exist(s2,'file')
    delete(s2);
end
if exist(s3,'file')
    delete(s3);
end

clear info MBARI_ID_str UW_ID_str ind ans fid tline
%clear s1 s2 s3
    
