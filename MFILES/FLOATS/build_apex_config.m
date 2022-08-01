function tf = build_apex_config(mbari_fn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT TO BUILD APEX FLOAT CONFIGS:
%
% A: OBTAIN ISUS CALIBRATION FILE
%     (1) Get MSC # for float from Dana's Inventory list using APF ID # as index
%         'http://runt.ocean.washington.edu/swift/Argo.....'
%     (2) Use MSC# to find nitrate cal file on CHEM
%         ISUSFloatCalibrationListing.xlsx
%     (3) Rename nitrate cal file & add needed header lines
%
% B: O2 CAL FILE
%     (1) cal lines from text file in float msg directory
%
% C: EXCEL PH CAL EXTRACTION
%     (1) Get pH cal data from excel spread sheet using MSC# as index
%     (2) pHlog.xlsx
%
% D: PDF FLBB CAL EXTRACTION 
%     (1) Check for FFLBB pdf cal files in float msg directory
%     (2) Extract FLBB cal info
%         parse_cal_pdf.m
%
% E: Generate float config file
%     (1) Write to local staging folder
%
% EMILY CLARK, TANYA MAURER, JOSH PLANT
% MBARI
%
% CHANGE HISTORY:
% 08/27/2021 EC changed section D to parse_cal_pdf function (JP)
% 08/30/2021 TM updated ISUS inventory loop to cycle through program types
% 09/01/2021 EC updated header builder/drop down menus
% 09/02/2021 JP updated section A to account for floats without NO3 cals
% 10/07/2021 EC updated ISUS directory to Mari's new version
% 11/03/2021 EC incorporated new JP parse_runt_inventory function
% 11/05/2021 JP changed directory to send config to 'staging' folder
% 12/06/2021 EC Added 'EQPac' as a basin option in drop down list
% 03/25/2021 JP udates for 8th order f(P) and 6th order k2 f(p) & bio
%               optics fixes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% add print to staging folder
% change xls read to readcell
% add optional input for MSC which would circumvent the UW search for the
% MSC

% % TESTING
%mbari_fn = 'ua12783'; % 2018 Tpos
%mbari_fn = 'ua19875'; % 2020 Tpos
%mbari_fn = 'ua12786'; % 2019
% mbari_fn = 'ua19970'; % 2 pdfs
%mbari_fn = 'ua19443'; % 1 pdf 
%mbari_fn = 'ua18861'; % 0 PDF
%mbari_fn = 'ua18601' % no nitrate cal
%mbari_fn = 'ua19751' % Resing 09/21
%mbari_fn = 'ua20043' % Resing 09/21
%mbari_fn = 'ua20704' % Resing 09/21
%mbari_fn = 'ua8497'

%% ************************************************************************
% CHOOSE FLOAT
 %mbari_fn = 'ua20002'
 %mbari_fn = 'ua20675'
 %mbari_fn = 'ua20148'

%% ************************************************************************
% DEFINE DIRS STRUCTURE
user_dir     = getenv('USERPROFILE'); %returns user path,i.e.'C:\Users\eclark'
user_dir     = [user_dir, '\Documents\MATLAB\']; % sets local file path
dirs.cal     = [user_dir,'ARGO_PROCESSING\DATA\CAL\']; % sets local cal file path
dirs.temp    = 'C:\temp\'; % sets temp file
dirs.msg     = '\\seaecho.shore.mbari.org\floats\UW\'; % sets location of float msg file
%dirs.ph      = '\\atlas\Chem\DuraFET\APEX+pH Calibrations\pHlogFiles\'; % sets location of pH cal
dirs.ph      = '\\atlas\chem\DuraFET\CALIBRATION\';
%dirs.ph      = 'C:\Users\jplant\Documents\MATLAB\Apps\pH_Calibration\test data\';
dirs.isus    = '\\atlas\Chem\ISUS\mari\'; % sets location of nitrate cal NEW
dirs.ph_fn   = 'pHlog_jp.xlsx';

% DO SOME PREP WORK
APEX_ID_str = regexp(mbari_fn,'\d+','once','match'); % isolates float # w/o letters
config_fn   = [mbari_fn,'_FloatConfig.txt']; % creates a text file name with float #
ncal_fn     = [mbari_fn,'.cal']; % creates the cal file name with the float #

dvec_now    = datevec(now); % sets exact time and date of run
yr          = dvec_now(1); % use this to limit file search to yr & yr-1
yr_rewind   = -3; % Number of previous yr files or directories to check for matches

O2_cal = {}; % creating a tbd space(cell array) for the O2 cal

%% ************************************************************************
% GET MSC # FOR A GIVEN APF ID # FROM UW
% CHECK ISUS INVENTORY FOR PREVIOUS & PRESENT YEARS

INVtype = {'IsusInventory','GobgcInventory','TposInventory'};
tfg     = 0; % flag to break out of 1st loop (inventory type)

for inv_type = 1:size(INVtype, 2)
    INV = INVtype{inv_type};
    
    for yr_ct  = yr_rewind:0 % for previous and present years
        yr_str = num2str(yr+yr_ct,'%0.0f');
        
        url = ['http://runt.ocean.washington.edu/swift/Argo',yr_str, ...
            'Logistics/',INV,'.mbari'];
        disp(['Checking ',url, ' for APF # = ',APEX_ID_str])
        
        d = parse_runt_inventory(dirs.temp, url);

        if isempty(d)
        %if size(d) == [0 0]
            disp('This combination of year and inventory probably does not exist')
            continue
        else
            tf = strcmp(d.data(:,2), APEX_ID_str);
            if sum(tf)== 1
            %if sum(tf)>0
                tfg = 1;
                % extract MSC
                MSC_str = d.data{tf, 5};
                MSC = str2double(MSC_str);
                if isnan(MSC)
                    disp(['APF found but NO MSC # present on ',APEX_ID_str])
                    return
                end
            elseif sum(tf) > 1 % multiple MSC's something wrong!
                disp(['Multiple MSC matches for ',APEX_ID_str, ...
                    ' something is fishy!'])
                return
            else
                disp('APF not found...but I will keep looking :)')
                continue
            end
        end
        
        if tfg == 1
            disp([APEX_ID_str, ' found in ',url])
            break
        end
    end
    if tfg == 1
        break
    end
end

if tfg == 0
    disp('WARNING: Float was not found on any of the runt urls searched')
    disp('WARNING: No MSC extracted')
end

% clean up variables that are no longer needed:
clear i d url dvec_now ans chk fid fixed_w tline str line_ct tmp tmp_str yr_ct

%TESTING
%  MSC = 916;
% MSC_str = '916';


%% ************************************************************************
% NOW TRY AND FIND ISUS CALIBRATION FILE
% LOOK IN MASTER LIST FIRST, IF NOT FOUND THEN HUNT WITH MSC# TO FIND
% DIRECTORY & THEN CALIBRATION FILE

tf_no3  = 1; % switch to 0 if NO3 cal file not found

%mbari_ncals_fn = 'ISUSFloatCalibrationListing.xlsx'; % pulls up excel directory of nitrate cal filenames
mbari_ncals_fn = 'ISUS_Summary2.xlsx'; % pulls up excel directory of nitrate cal filenames NEW

[status,msg] = copyfile([dirs.isus, mbari_ncals_fn], dirs.temp); % copy to local

if status==0
    disp('WARNING: Could not copy ISUSFloatCalibrationListing.xlsx to local')
else
    % read the file names in the xls into a new cell array
    [~,~,ncals] = xlsread([dirs.temp, mbari_ncals_fn], ...
        'ISUSCalFilenames','C14:E500'); % AH500 # may need to be extended
    
    
    % create an array to hold pertinent info for finding float of choice
    ncals_hdr = {'MSC' 'dir path' 'file name'};
    ncals = ncals(~cellfun(@isnan,ncals(:,1)),:); % remove NaN MSC rows
    tMSC = cell2mat(ncals(:,1)) == MSC;  % create regular array to suit needs
    
    % ALL GOOD BUT CHECK THAT FILE EXISTS TOO
    if sum(tMSC) == 1  && exist([ncals{tMSC,2},'\',ncals{tMSC,3}],'file')%if the file name is present in the preceeding columns...
        disp(['ISUS calibration file found in master list for  MSC # ',MSC_str, ...
            ' and APEX ID # ',APEX_ID_str,':'])
        cal_fn = ncals{tMSC,3}; % separate filename
        fp = [ncals{tMSC,2},'\',cal_fn]; % path to ncal file
    else
        disp(['Nitrate cal for MSC# ',MSC_str,' not found in ',mbari_ncals_fn])
        disp('Searching directory structure for a MSC match ...')
        
        ncal_dir ={};
        for yr_ct = yr_rewind:0 % look in current & previous year for file
            yr_str      = num2str(yr+yr_ct,'%0.0f');
            apex_dir = ['APEX',yr_str ];
            list = dir([dirs.isus,apex_dir]);
            dir_names = {list.name}';
            tf = ~cellfun(@isempty,regexp(dir_names,[MSC_str,'(?=\_Spec)'], ...
                'once','match'));
            if sum(tf,1) > 0
                ncal_dir = dir_names(tf);
                break
            end
        end
        if isempty(ncal_dir)
            disp(['ISUS calibration directory could not be found for MSC # ',...
                MSC_str])
            tf_no3  = 0;
            %return
            
        else  % look for isus cal files
            n_tmp = cell(size(ncal_dir,1),3); % predim
            tg = ones(size(ncal_dir,1),1);    % predim
            for n_ct = 1:size(ncal_dir,1)
                n_tmp{n_ct,1} = ncal_dir{n_ct,1};
                
                nd = dir([dirs.isus, apex_dir,'\', ncal_dir{n_ct,1},'\*cal']);
                if isempty(nd)
                    tg(n_ct) = 0;
                    continue
                end
                
                n_tmp{n_ct,2} = nd.name;
                n_tmp{n_ct,3} = nd.datenum;
            end
            n_tmp = n_tmp(logical(tg),:); % remove empty file name rows
            [~, ind] = sort(cell2mat(n_tmp(:,3)), 'descend');
            n_tmp = n_tmp(ind,:); % sorted by cal file time stamp, newest first
            
            if size(n_tmp,1) > 1
                disp(['Multiple ISUS calibration files found for MSC # ',...
                    MSC_str,'. Choosing most recent file:'])
            end
            cal_path = [dirs.isus, apex_dir,'\',n_tmp{1,1},'\'];
            cal_fn   = n_tmp{1,2};
            fp       = [cal_path, cal_fn];
            disp(['ISUS calibration file found for  MSC # ',MSC_str, ...
                ' and APEX ID # ',APEX_ID_str,':'])
            disp(fp)
            clear n_tmp tg n_ct nd
        end
    end
    
   % **********************************************************************
    % NO3 cal file found copy to local, add header info & rename
    
    tf      ='Y';
    if tf_no3 == 1 && exist([dirs.cal,'NO3_CAL\',ncal_fn],'file')
        str = [dirs.cal,'NO3_CAL\',ncal_fn,' already exists!'];
        disp(str)
        tf = input(['Do you want to replace the existing nitrate ', ...
            'calibration file (Y/N)?'],'s'); % input argument displays the text
        % in the prompt (before 1st comma)
        % DOES NOT EVALUATE THE RESPONSE
    end
    
    if tf_no3 == 1 && ~isempty(regexpi(tf,'^Y','once')) %adding headers
        add_hdr{1,1} = ['H,Original source file: ',fp,',,,,'];
        add_hdr{2,1} = 'H,Pixel base,1,,,';
        add_hdr{3,1} = 'H,Sensor Depth offset,0,,,';
        add_hdr{4,1} = 'H,Br wavelength offset,210,,,';
        add_hdr{5,1} = 'H,Min fit wavelength,,,,';
        add_hdr{6,1} = 'H,Max fit wavelength,,,,';
        add_hdr{7,1} = 'H,Use seawater dark current,No,,,';
        add_hdr{8,1} = 'H,Pressure coef,0.0265,,,';
        
        [status,msg] = copyfile(fp,dirs.temp);
        if status==0
            disp(['WARNING: Could not copy original source file: ',fp, 'to local'])
        else
            
            fid     = fopen([dirs.temp,cal_fn]);
            fid_new = fopen([dirs.cal,'NO3_CAL\',ncal_fn],'w');
            tline = ' ';
            while ischar(tline)
                if regexp(tline,'^H,CalTemp','once')
                    fprintf(fid_new,'%s\r\n', tline);
                    for i = 1:size(add_hdr,1)
                        fprintf(fid_new,'%s\r\n', add_hdr{i,1});
                    end
                elseif regexp(tline,'^E|^H','once')
                    fprintf(fid_new,'%s\r\n', tline);
                end
                tline = fgetl(fid);
            end
            fclose(fid);
            fclose(fid_new);
            clear fid fid_new tline i
        end
    elseif tf_no3 == 1
        disp([dirs.cal,'NO3_CAL\',ncal_fn,' already exists and will not be updated'])
    else
        disp(['NO ISUS CALIBRATION FILE WAS FOUND FOR MSC # ',MSC_str])
        disp([mbari_fn,' May not have a nitrate sensor - check to confirm!!'])
    end
end

% ************************************************************************
% NOW BUILD CONFIG.TXT FILE
tf ='Y';
if exist([dirs.cal,'\FLOAT_CONFIG\',config_fn],'file') % see if it already exists
    str = [dirs.cal,'\FLOAT_CONFIG\',config_fn,' already exists!'];
    disp(str)
    tf = input(['Do you want to replace the existing config ', ...
        'calibration file (Y/N)?'],'s');
    if ~strncmpi(tf,'Y',1)
        disp([dirs.cal,'\FLOAT_CONFIG\',config_fn,' already exists and will not be updated'])
        return
    end
end

% CHECK FOR OXYGEN CALIBRATION FILE
msg_path =[dirs.msg,'f',APEX_ID_str,'\']; % set path to float in seaecho
O2cal_fn = ls([msg_path,'ox*calibration*']); % separate O2 filename

% if, by chance, there is more than 1 O2 cal, choose the right one
if size(O2cal_fn,1) > 1
    O2cal_fn = uigetfile([msg_path,'ox*calibration*'],'Please choose O2 cal file');
end

if size(O2cal_fn,1) == 1 % only one file found
    fid = fopen([msg_path, strtrim(O2cal_fn)]);
    tline = ' ';
    O2_chk = 0;
    O2_ct  = 0;
    srch_str = ['PhaseCoef|FoilID|FoilCoefA|FoilCoefB|FoilPolyDegT|', ...
        'FoilPolyDegO|SVUFoilCoef|ConcCoef'];
    tf_coef_repeat = 0; % deal with optode cals like 20148 with repeated cal block
    while ischar(tline)
        if regexp(tline,'PhaseCoef','once') % build SN line
            tmp = regexp(tline(regexp(tline,'PhaseCoef\s+','end','once')...
                :end),'\d+','match');
            O2_ct = O2_ct+1;
            str = ['OptodeSn = ',tmp{1},' ',tmp{2}];
            tf_coef_repeat = tf_coef_repeat +1;
            if tf_coef_repeat < 2
                O2_cal{O2_ct,1} = str;
            else
                break
            end
            clear tmp str
        end
        if regexp(tline, srch_str,'once')
            O2_ct = O2_ct+1;
            O2_cal{O2_ct,1} = tline;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    disp('O2 cal file found!')
elseif size(O2cal_fn,1) > 1
    disp('More than 1 O2 cal file found - please choose file')
    [file,path] = uigetfile([msg_path,'oxygen*'],'Please choose O2 cal file');
else
    disp('No O2 cal file found - no cal data extracted')
end

%% ************************************************************************
% TRY AND GET PH CAL FROM excel SHEET
disp('Trying to extract pH calibration data ...')
[status,msg] = copyfile([dirs.ph, dirs.ph_fn], dirs.temp); % copy to local temp file
d = get_ph_version();

% GET INDICES
iDF     = find(strcmp(d.hdr,'DF ID') == 1, 1);
iAPX    = find(strcmp(d.hdr,'UW ID') == 1);
iMSC    = find(strcmp(d.hdr,'MSC') == 1);
iK2     = find(strncmp(d.hdr,'k2',2) == 1);
iFP     = find(strncmp(d.hdr,'f(P)',4) == 1);
iK0_ON  = find(strcmp(d.hdr,'k0 pON') == 1);
iK0_OFF = find(strcmp(d.hdr,'k0 pOFF') == 1);
iK0_HCL = find(strcmp(d.hdr,'k0 HCl') == 1);
iVER    = find(strcmp(d.hdr,'Build notes') == 1);

ph_ct    = 0;
tf_pH_continue = 1;
ph_cal   = {}; %predim empty incase no pH sensor
APEX_num = str2double(APEX_ID_str);

% TRY & LOCATE MSC # & APEX ID # IN pHlog.xlsx FILE
tMSC = cell2mat(d.data(:,iMSC))  == MSC;
tAPX = cell2mat(d.data(:,iAPX))  == APEX_num;
tpH  = tMSC;

if sum(tMSC) == 0 && sum(tAPX) == 0
    fprintf('WARNING: NEITHER MSC OR APEX ID FOUND in %s FOR %s\n', dirs.ph_fn, mbari_fn);
    fprintf('NO pH SENSOR OR CAL DATA HAS NOT BEEN ENTERED IN %s',dirs.ph_fn);
    tf_pH_continue = 0;
elseif sum(tMSC,1) == 0 && sum(tAPX,1) == 1
    fprintf('WARNING: APEX ID EXISTS BUT NO MSC YET IN %s FOR %s\n', dirs.ph_fn, mbari_fn);
    tpH = tAPX;
    tf_pH_continue = 0;
elseif sum(tMSC) > 1
    fprintf('WARNING: MORE THAN 1 MSC # IDENTIFIED IN %s FOR %s\n', dirs.ph_fn, mbari_fn);
    return
end

if tf_pH_continue == 0 % CHECK TO SEE IF YOU WANT TO CONTINUE
    tf ='N';
    tf = input('No pH calibration data found or MSC# missing. Continue anyway? (Y/N)?','s');
    if ~strncmpi(tf,'Y',1)
        disp('EXITING! Config file will not be built');
        return
    end
    tf_pH_continue = 1;
end

% NOW TRY AND GET COEFFICIENTS
if tf_pH_continue == 1
    DF      = d.data{tpH,iDF};
    if isnumeric(DF)
        DF = sprintf('%0.0f',DF); % make a string because some are
    end
    APEX_ID = d.data{tpH,iAPX};
    K0_ON   = d.data{tpH,iK0_ON};
    K0_OFF  = d.data{tpH,iK0_OFF};
    K0_HCL  = d.data{tpH,iK0_HCL};
    K2      = cell2mat(d.data(tpH, iK2));
    FP      = cell2mat(d.data(tpH, iFP));
    
    % DO SOME CHECKS
    if isnan(APEX_ID)
        fprintf('WARNING: No APEX ID FOR pH SN %s in %s. Building cal anyway\n',...
            DF, dirs.ph_fn); 
    elseif APEX_ID ~= APEX_num
        fprintf('APEX ID(%0.0f) does not match ISUS inventory(%0.0f)\n', ...
            APEX_ID, APEX_num);
        tf_pH_continue = 0;
    end
    if isnan(K0_ON)
        fprintf('No valid k0 pump ON value (%0.0f) for %s\n', ...
            K0_ON, mbari_fn);
        tf_pH_continue = 0;
    end
    if all(isnan(K2))
        fprintf('No valid k2 value(s) for %s\n',mbari_fn);
        tf_pH_continue = 0;
    end
    if all(isnan(FP))
        fprintf('No valid f(P) value(s) for %s\n',mbari_fn);
        tf_pH_continue = 0;
    end
end

% BUILD CELL ARRAY FOR PRINTING
if tf_pH_continue == 0
    fprintf('No pH calibration will be built for %s\n',mbari_fn);
else
    K2      = K2(~isnan(K2)); % Remove any NaN's
    FP      = FP(~isnan(FP)); % Remove any NaN's
    ph_cal{1,1} = sprintf('Durafet SN = %s, APEX# %s',DF, APEX_ID_str);
    
    % CHECK FOR LAB KO (PUMP ON IS DEFAULT VALUE TO USE)
    if ~isnan(K0_ON)
        k0_str = sprintf(['%0.5f, =k0 Pon; %0.5f =k0 Poff; ',...
            '%0.5f =k0 HCl'],K0_ON, K0_OFF, K0_HCL);
    elseif isnan(K0_ON) && ~isnan(K0_HCL)
        k0_str = sprintf(['%0.5f, =k0 HCl; WARNING: PROVISONAL ',...
            'No lab k0 was found'],K0_HCL);
    end
    ph_cal{3,1} = k0_str;
    if size(K2,2) == 1 % K2 constant & not a function of P
        ph_cal{4,1} = sprintf('%0.5e, =k2*T',K2);
    else
        ph_cal{4,1} = sprintf([repmat('%0.5e,',1,size(K2,2)),' =k2(fP)*T'],flip(K2));
    end
    
    % CHECK & BUILD Pcoef strings
    FP(end) = []; % Remove constant coef so p(0) = 0
    P = (flip(FP))';
    ph_cal{2,1} = sprintf(['%0.0f, number of calibration ', ...
        'coefficients'],size(P,1)+2);
    for i = 1:size(P,1)
        ph_cal{4+i,1} = sprintf(['%0.4E, =k',num2str(i+2),'*P^',...
            num2str(i)],P(i));
    end
end

clear DF_num  DF_txt DF_APX K0_ON K0_OFF K0_HCL
clear iDF iAPX iMSC iK6 iK2T iK0_HCL iK0_ON iK0_OFF ph_ct


%% ************************************************************************
% OK NOW CHECK FOR FLBB PDF CAL SHEETS & TRY & EXTRACT SCALE & DARK COUNTS
% This uses a function from the files exchange by Derek Wood:
% https://www.mathworks.com/matlabcentral/fileexchange/
% 63615-read-text-from-a-pdf-document
% Before pdf reader can work add a couple lines to startup.m to
% switch to the dir where the mfile is found and execute the following
% line: javaaddpath('iText-4.2.0-com.itextpdf.jar')
% it only needs to be executed once upon the start of a matlab session

clear SBE SBEbbp bbp_dc bbp_scale SBEChl chl_dc chl_scale SNchar SN

disp(['Looking for bio-optical calibration data in: ',msg_path]);
optics_cal  = {};

% ANY PDF FILES IN THE MSG DIR
tmp = dir([msg_path,'*CharSheet?.pdf']);

if isempty(tmp) % no pdfs with "Char sheets"
    disp(['No PDF files found in: ',msg_path]);
    disp('No bio-optical calibration data will be extracted');
    optics_cal = {};
else
    % loop through 1+ pdf Char Sheets
    t_optics = 1; chl_str = ''; bbp_str = '';
    for i = 1:size(tmp,1)
        fn = tmp(i);
        pdf_fn = fullfile(fn.folder, fn.name);
        SBE = parse_cal_pdf(pdf_fn); % call JP pdf parse cal function
        if isfield(SBE,'Chl')
            optics_cal{1,1} = sprintf('CHLFLUOR SN %s-%s', SBE.Chl.Model, SBE.Chl.SN);
            optics_cal{2,1} = sprintf('%0.0f ChlDC', SBE.Chl.DC);
            optics_cal{3,1} = sprintf('%0.4f ChlScale', SBE.Chl.Scale);
            tg_chl = SBE.Chl.DC > 1 & SBE.Chl.DC < 100 & ...
                SBE.Chl.Scale > 1e-3 & SBE.Chl.Scale < 0.01;
            if ~tg_chl
                t_optics = 0;
            end
        end

        if isfield(SBE,'bbp700')
            optics_cal{4,1} = sprintf('%0.0f BetabDC', SBE.bbp700.DC);
            optics_cal{5,1} = sprintf('%0.4e BetabScale', SBE.bbp700.Scale);
            tg_bbp = SBE.bbp700.DC > 1 & SBE.bbp700.DC < 100 & ...
                SBE.bbp700.Scale > 1e-7 & SBE.bbp700.Scale < 3e-6;
            if ~tg_bbp
                t_optics = 0;
            end
        end
        % SANITY CHECK BIO OPTICS CAL COEFFICIENTS
        if t_optics == 0
            disp('PDF bioptical extraction did not seem reasonable- add manually')
            fprintf(['Chl DC = %0.0f   Chl Scale = %0.4f\n',...
                'bbp DC = %0.0f   bbp Scale = %0.4f\n'], SBE.chl.DC, ...
                SBE.chl.Scale, SBE.bbp700.DC, SBE.bbp700.Scale);
        end
    end
end

            
            
        
%         
%         
%             if chl_dc > 1 && chl_dc < 99 && chl_scale < 0.02 && ...
%             chl_scale > 0.003 && bbp_dc < 99 && bbp_dc > 1 && ...
%             bbp_scale < 10e-6 && bbp_scale > 0.5e-6
%         
%         
%         
%         
%         SBE = parse_cal_pdf(pdf_fn); % call JP pdf parse cal function
%         % look for bbp and chl in one or both pdfs
%         if isfield(SBE,'bbp700')
%             SBEbbp = struct2cell(SBE.bbp700);
%             
%             bbp_dc    = SBEbbp{3,1};
%             bbp_scale = SBEbbp{2,1};
%         end
%         if isfield(SBE,'Chl')
%             SBEChl = struct2cell(SBE.Chl);
%             SNchar = regexp(fn.name, '^FLBB\w+', 'once', 'match'); % find acronym for SN
%             chl_dc    = SBEChl{3,1};
%             chl_scale = SBEChl{2,1};
%         end
%     end
%     % in cases where the FLBB model is not in the file name, default to AP2
%     if isempty(SNchar)
%         SN = append('FLBBAP2-', SBEChl{1,1});
%         disp( 'Warning: FLBBAP2 label added by default')
%     else
%         SN = append(SNchar,'-', SBEChl{1,1});
%     end
%     
%     % MANUALLY ENTER IF PARSER NOT WORKING:
%     SN = 'FLBBAP2-6329'
%     chl_dc = 47
%     chl_scale = 0.0073 
%     bbp_dc = 46
%     bbp_scale = .000001762
%     
%     % sanity check and build optics_cal
%     if chl_dc > 1 && chl_dc < 99 && chl_scale < 0.02 && ...
%             chl_scale > 0.003 && bbp_dc < 99 && bbp_dc > 1 && ...
%             bbp_scale < 10e-6 && bbp_scale > 0.5e-6
%         
%         optics_cal{1,1} = sprintf('CHLFLUOR SN %s', SN);
%         optics_cal{2,1} = sprintf('%s ChlDC', string(chl_dc));
%         optics_cal{3,1} = sprintf('%s ChlScale', chl_scale);
%         optics_cal{4,1} = sprintf('%s BetabDC', string(bbp_dc));
%         optics_cal{5,1} = sprintf('%s BetabScale', bbp_scale);
%         disp('Bio optical calibration coefficients have been extracted')
%     else
%         disp ...
%             ('PDF bioptical extraction did not seem reasonable- add manually')
%         pause
%     end
%end

%% ************************************************************************

% CHOOSE PROGRAM AFFILIATION
ProgramList = {'GO-BGC' 'SOCCOM' 'UW/MBARI BGC-Argo' 'TPOS' 'EXPORTS', ...
     'SBE83 O2 TEST' 'MBARI GDF-TEST' 'NOT DEFINED'};
 [ind,tf] = listdlg('ListString', ProgramList, 'SelectionMode','single',...
     'Name','Choose Program Affiliation','ListSize',[300,200], ...
     'PromptString',['Choose for: ', mbari_fn]);

 if tf == 0
    disp(['No Program afiliation was choosen! Setting affiliation to: ',...
        'UW/MBARI BGC-Argo']);
    Program = 'UW/MBARI BGC-Argo';
else
    Program = ProgramList{ind};
end

% CHOOSE OCEAN REGION
RegionList = {'SO' 'NPac' 'SPac' 'EQPac' 'NAtl' 'SAtl' 'IO' 'ARC' 'NOT DEFINED'};
[ind,tf] = listdlg('ListString', RegionList, 'SelectionMode','single',...
    'Name','Choose Region','ListSize',[300,200], ...
    'PromptString',['Choose for: ', mbari_fn]);
if tf == 0
    disp(['No ocean region was choosen! Setting affiliation to: ',...
        'NOT DEFINED']);
    Region = 'NOT DEFINED';
else
    Region = RegionList{ind};
end

%% ************************************************************************
% IF YOU GET TO HERE O2 CAL & PH CAL HAVE BEEN SUCCESFULLY EXTRACTED
% PRINT TO FILE

fid = fopen([dirs.cal,'\FLOAT_CONFIG\staging\',config_fn], 'w'); 
%fid = fopen([dirs.cal,'\FLOAT_CONFIG\',config_fn], 'w'); 
fprintf(fid,'%s\r\n',['Institution ID: ', APEX_ID_str]);
fprintf(fid,'%s\r\n',['MBARI ID: ', mbari_fn]);
fprintf(fid,'%s\r\n',['Program: ', Program]);
fprintf(fid,'%s\r\n',['Region: ', Region]);
fprintf(fid,'%s\r\n', msg_path);

out = [O2_cal; ph_cal; optics_cal]; % print O2, pH, and optics to file

for i = 1:size(out,1)
    fprintf(fid,'%s\r\n',out{i,1});
end
fclose(fid);

% SPOOL TO SCREEN AS CHECK
fid = fopen([dirs.cal,'\FLOAT_CONFIG\staging\',config_fn]);
%fid = fopen([dirs.cal,'\FLOAT_CONFIG\',config_fn]);

tline = ' ';
while ischar(tline)
    disp(tline)
    tline = fgetl(fid);
end
fclose(fid);
