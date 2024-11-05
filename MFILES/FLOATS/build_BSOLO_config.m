%function tf = build_BSOLO_config(mbari_id)

%  01/29/24: TM, Modifications for more cal flavors!  ss4003 pending
%  deployment...
%
% ***********************************************************************
% TESTING
%mbari_id = 'ss0001';
%mbari_id = 'ss0002';
%mbari_id = 'ss0003';

mbari_id = 'ss4004' %not processed with code  - missing O2 pdf cab get 
%data from *.meta but not incorperated yet


% ***********************************************************************
% DO SOME PREP WORK &  FILE TESTS
% ***********************************************************************
disp(['Building config file for: ', mbari_id])

% SET UP DIRS
user_dir    = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir    = [user_dir, '\Documents\MATLAB\'];
config_dir  = [user_dir, 'ARGO_PROCESSING\DATA\CAL\FLOAT_CONFIG\staging\'];
no3_dir     = config_dir;
%no3_dir     = [user_dir, 'ARGO_PROCESSING\DATA\CAL\NO3_CAL\'];
msg_dir     = '\\seaecho.shore.mbari.org\floats\SIO\';
temp_dir    = 'C:\temp\';


% DO SOME PREP WORK
SOLO_ID_str = regexp(mbari_id,'\d+','once','match'); % isolates float # w/o letters
config_fn   = [mbari_id,'_FloatConfig.txt']; % creates a text file name with float #
ncal_fn     = [mbari_id,'.cal']; % creates the cal file name with the float #
Program     = '';
Region      = '';

%TEST FOR SIO SOLO FLOAT 
if strncmp(mbari_id,'ss',2) % SIO SOLO float
    float_type = 'BGCSOLO';
else
    disp(['"',mfilename,'" is expecting a BGCSOLO float!'])
    disp(['mbari id = "',mbari_id, '" suggests a different float type'])
    return
end

%TEST FOR FLOAT MSG & CAL DIRS
msg_fp = fullfile(msg_dir, mbari_id,'\');
if ~isfolder(msg_fp)
    disp(['msg dir for ', mbari_id,' not found at: ',msg_fp]);
    return
end

cal_fp = fullfile(msg_fp,'cals\');
if ~isfolder(cal_fp)
    disp(['cal dir for ', mbari_id,' not found at: ',cal_fp]);
    return
else % get cal dir listing
    tmp = dir(cal_fp);
    tf_dir = cell2mat({tmp.isdir}');
    cal_fnames = {tmp.name}';
    cal_fnames = cal_fnames(~tf_dir);
    clear tmp tf_dir
end

% ***********************************************************************
%                     CHOOSE PROGRAM AFFILIATION
% ***********************************************************************
ProgramList = {'GO-BGC' 'SOCCOM' 'UW/MBARI BGC-Argo' 'TPOS' 'EXPORTS', ...
     'SBE83 O2 TEST' 'MBARI GDF-TEST' 'NOT DEFINED'};
 [ind,tf] = listdlg('ListString', ProgramList, 'SelectionMode','single',...
     'Name','Choose Program Affiliation','ListSize',[300,200], ...
     'PromptString',['Choose for: ', mbari_id]);

 if tf == 0
    disp(['No Program afiliation was choosen! Setting affiliation to: ',...
        'UW/MBARI BGC-Argo']);
    Program = 'UW/MBARI BGC-Argo';
else
    Program = ProgramList{ind};
 end

% ***********************************************************************
%                       CHOOSE OCEAN REGION
% ***********************************************************************
RegionList = {'SO' 'NPac' 'SPac' 'EQPac' 'NAtl' 'SAtl' 'IO' 'ARC' 'NOT DEFINED'};
[ind,tf] = listdlg('ListString', RegionList, 'SelectionMode','single',...
    'Name','Choose Region','ListSize',[300,200], ...
    'PromptString',['Choose for: ', mbari_id]);
if tf == 0
    disp(['No ocean region was choosen! Setting affiliation to: ',...
        'NOT DEFINED']);
    Region = 'NOT DEFINED';
else
    Region = RegionList{ind};
end

% ***********************************************************************
% ***********************************************************************
%                    TRY AND GET OXYGEN CALIBRATION DATA
% ***********************************************************************
% ***********************************************************************
tf_O2C = ~cellfun(@isempty, regexp(cal_fnames,'Oxygen.pdf'));
tf_O2T = ~cellfun(@isempty, regexp(cal_fnames,'Temperature','once'));
% tf_O2  = tf_O2C & ~tf_O2T; %TM I'm not understanding this logic... for
% 4002 this need to be modified...argh.
tf_O2  = tf_O2C;

if sum(tf_O2) == 1
    SBE = parse_cal_pdf(fullfile(cal_fp, cal_fnames{tf_O2}));
elseif sum(tf_O2) > 1
    O2cal_fn = uigetfile([cal_fp,'*Oxygen*.pdf'],'Please choose O2 cal file');
    SBE = parse_cal_pdf(fullfile(cal_fp, O2cal_fn));
end

if sum(tf_O2) == 0 % NO PDF file found
    fprintf('No pdf O2 calibration file found for %s at : %s\n', ...
        mbari_id, cal_fp);
    fprintf('Looking for meta file to extract O2 cal info ....\n');
    tf_m = ~cellfun(@isempty, regexp(cal_fnames,'meta$'));
    if sum(tf_m) ~= 1
        disp(['No O2 phase calibration found at: ',cal_fp,' for float ',mbari_id])
        return
    end
    SBE = parse_BSOLO_meta(fullfile(cal_fp, cal_fnames{tf_m}));
end

% REMOVE SPACE IN O2 MODEL NAME if it exists (i.e. "SBE 63" => "SBE63")
SBE.O2.Model = regexprep(SBE.O2.Model,' ','');
if regexp(mbari_id,'ss000[123]','once')
    SBE.O2.Model = regexprep(SBE.O2.Model,'63','83'); % SBE defines incorrectly
end

O2 = SBE.O2;
if ~isfield(O2,'TA') % Do temperature coeffiecients exist yet?
    tf_O2 = ~cellfun(@isempty, regexp(cal_fnames,'Temperature.pdf'));
    if sum(tf_O2) == 1
        SBE = parse_cal_pdf(fullfile(cal_fp, cal_fnames{tf_O2}));
    elseif sum(tf_O2) > 1
        O2cal_fn = uigetfile([cal_fp,'*Temperature*pdf'], ...
            'Please choose O2 cal file');
        SBE = parse_cal_pdf(fullfile(cal_fp, O2cal_fn));
    else
        disp(['No O2 T calibration found at: ',cal_fp,' for float ',mbari_id])
        return
    end
    O2.TA = SBE.O2.TA;
end

% O2 SANITY CHECKS
if size(O2.A,2) ~= 3
    disp(['Incorrect # of O2 "A" coeffcients! Please Check calibration ', ...
        'file & pdf parser ']);
    return
end
if size(O2.B,2) ~= 2
    disp(['Incorrect # of O2 "B" coeffcients! Please Check calibration ', ...
        'file & pdf parser ']);
    return
end
if size(O2.C,2) ~= 3
    disp(['Incorrect # of O2 "C" coeffcients! Please Check calibration ', ...
        'file & pdf parser ']);
    return
end
if size(O2.TA,2) ~= 4
    disp(['Incorrect # of O2 "TA" coeffcients! Please Check calibration ', ...
        'file & pdf parser ']);
    return
end

% USE TEMPLATES FOR CONFIG FILE O2 LINES
O2_temp_str = ['O2 sensor (%s %s) Temp coefficents [TA0 TA1 TA2 TA3],', ...
    '%0.6e,%0.6e,%0.6e,%0.6e'];
O2_phase_str = ['O2 sensor (%s %s) Phase coefficents [A0 A1 A2 ',...
    'B0 B1 C0 C1 C2],%0.6e,%0.6e,%0.6e,%0.6e,%0.6e,%0.6e,%0.6e,%0.6e'];

% *************************************************************************
% *************************************************************************
% IF YOU GET HERE AT LEAST O2 CAL DATA EXIST
% *************************************************************************
% *************************************************************************
out = cell(30,1);

out{1} = ['Institution ID: ', SOLO_ID_str];
out{2} = ['MBARI ID: ', mbari_id];
out{3} = ['Program: ', Program];
out{4} = ['Region: ', Region];
out{5} = '';
out{6} = msg_fp;
out{7} = sprintf(O2_temp_str, O2.Model, O2.SN, O2.TA);
out{8} = sprintf(O2_phase_str, O2.Model, O2.SN, O2.A, O2.B, O2.C);

% *************************************************************************
%        NOW TRY AND GET ECO TRIPLET CALIBRATIONS
% *************************************************************************
disp(['Looking for bio-optical calibration files(s) in: ',cal_fp]);
tf_eco   = ~cellfun(@isempty, regexp(cal_fnames,'CharSheet.\.pdf'));
if sum(tf_eco)==0
    tf_eco   = ~cellfun(@isempty, regexpi(cal_fnames,'Char Sheet.\.pdf')); % What the heck?
end
eco_cals = cal_fnames(tf_eco);
chl_str  = '%s Chl fluorescence (%s %s) [ChlDC ChlScale],%0.0f,%0.4e';
bbp_str  = '%s Backscatter700 (%s %s) [Betab700DC Betab700Scale],%0.0f,%0.4e';
fdom_str = '%s Fluorescent DOM (%s %s) [FDOMDC FDOMScale],%0.0f,%0.4e';

if isempty(eco_cals) % no pdfs with "Char sheets"
    disp(['No PDF files found in: ',cal_fp]);
    disp('No bio-optical calibration data will be extracted');
    optics_cal = {};
else
    for eco_ct = 1: size(eco_cals,2) % loop through 1+ pdf Char Sheets
        fn = eco_cals{eco_ct};
        pdf_fn = fullfile(cal_fp, fn);
        SBE = parse_cal_pdf(pdf_fn); % call JP pdf parse cal function
        
        if isfield(SBE,'Chl')
            cal = SBE.Chl;
            out{9} = sprintf(chl_str, cal.Model, cal.Model, cal.SN, ...
                cal.DC, cal.Scale);
        end
        if isfield(SBE,'bbp700')
            cal = SBE.bbp700;
            out{10} = sprintf(bbp_str, cal.Model, cal.Model, cal.SN, ...
                cal.DC, cal.Scale);
        end            
        if isfield(SBE,'FDOM')
            cal = SBE.FDOM;
            out{11} = sprintf(fdom_str, cal.Model, cal.Model, cal.SN, ...
                cal.DC, cal.Scale);
        end  
    end
end


% *************************************************************************
%                     SBE pH CALIBRATIONS
% *************************************************************************
disp(['Looking for pH calibration data file(s): ',cal_fp]);

tf_pH   = ~cellfun(@isempty, regexpi(cal_fnames,'CalSheet\.pdf'));
if sum(tf_pH) == 1
    SBE = parse_cal_pdf(fullfile(cal_fp, cal_fnames{tf_pH}));
elseif sum(tf_pH) > 1
    pHcal_fn = uigetfile([cal_fp,'*CalSheet.pdf'],'Please choose pH cal file');
    SBE = parse_cal_pdf(fullfile(cal_fp, pHcal_fn));
else
    disp(['No pH calibration found at: ',cal_fp,' for float ',mbari_id])
    skip_pH = 1;
end

if isfield(SBE,'pH')
    k0   = SBE.pH.k0;
    k2   = SBE.pH.k2;
    fP   = SBE.pH.fP;
    % SBE reports F0 for f(P) which is not used so remove 1st coefficient
    fP(1)=[];

    str1 = sprintf('pH sensor coefficients (SBE %s) [k0 ',SBE.pH.SN);
    if size(k2,2) == 1
        meta_k2 = 'k2 ';
    else
        meta_k2 = sprintf(repmat('k2%0.0f ',1,size(k2,2)),(1:size(k2,2))-1);
    end
    meta_fP  = sprintf(repmat('F%0.0f ',1,size(fP,2)),(1:size(fP,2)));
    pH_coefs = [k2, fP];
    pH_str   = sprintf(['%s%s%s],%0.5f', repmat(',%0.6E',1,size(pH_coefs,2))], ...
               str1, meta_k2, strtrim(meta_fP), k0, pH_coefs);
%     config_ct = config_ct+1;
%     config_cell{config_ct} = pH_str;
    cal = SBE.pH;
    out{12} = sprintf(pH_str, cal.SN, cal.k0, cal.k2, flip(cal.fP));
%     keyboard
    clear  O2Pstr Pcoef Pmeta t1 str2
end

% TM Can delete this section....now that we have more complex options for
% cal coef use the pH extraction section that was more recently written from build_NAVIS_float
%--------------------------------------------------------------------------
% % % skip_pH = 0;
% % % pH_meta_cell = {'k0' 'k2' 'k3xP' 'k4xP^2' 'k5xP^3' 'k6xP^4', ...
% % %                 'k7xP^5' 'k8xP^6'};
% % % pH_fmt_cell  = {'%0.4f' '%0.4e' '%0.4e' '%0.4e' '%0.4e' '%0.4e', ...
% % %                 '%0.4e' '%0.4e'};
            
% % % pH_str = ['pH sensor coefficients (SBE %s) ', ...
% % %     '[k0 k2 k3xP k4xP^2 k5xP^3 k6xP^4 k7xP^5 k8xP^6],', ...
% % %     '%0.4f,%0.4e,%0.4e,%0.4e,%0.4e,%0.4e,%0.4e,%0.4e'];
% % % 
% % % jp = sprintf('%s ',pH_meta_cell{:})


% % % if sum(tf_pH) == 1
% % %     SBE = parse_cal_pdf(fullfile(cal_fp, cal_fnames{tf_pH}));
% % % elseif sum(tf_pH) > 1
% % %     pHcal_fn = uigetfile([cal_fp,'*CalSheet.pdf'],'Please choose pH cal file');
% % %     SBE = parse_cal_pdf(fullfile(cal_fp, pHcal_fn));
% % % else
% % %     disp(['No pH calibration found at: ',cal_fp,' for float ',mbari_id])
% % %     skip_pH = 1;
% % % end
% keyboard
% % % if skip_pH == 0 && isfield(SBE,'pH')
% % %     cols_fP  = size(SBE.pH.fP,2); % of poly coefs (excuding constant = 0)
% % %     pH_meta_cell = pH_meta_cell(1:cols_fP+2);
% % %     pH_meta_str  = strtrim(sprintf('%s ',pH_meta_cell{:}));
% % %     
% % %     pH_fmt_cell = pH_fmt_cell(1:cols_fP+2);
% % %     pH_fmt_str  = strtrim(sprintf('%s,',pH_fmt_cell{:}));
% % %     pH_fmt_str  = pH_fmt_str(1:end-1); % loose extra comma at end
% % %     
% % %     pH_str = ['pH sensor coefficients (SBE %s) [', pH_meta_str,'],', ...
% % %         pH_fmt_str];
% % %     
% % %     
% % %     
% % %     cal = SBE.pH;
% % %     out{12} = sprintf(pH_str, cal.SN, cal.k0, cal.k2, flip(cal.fP));
% % % end
%--------------------------------------------------------------------------




% *************************************************************************
%                     OCR CALIBRATION - a txt file
%          4 channels & a0, a1, im vales for each channel
% *************************************************************************
disp(['Looking for OCR504 calibration data file(s): ',cal_fp]);
tf_OCR   = ~cellfun(@isempty, regexpi(cal_fnames,'OCR504'));
skip_OCR = 0;
default_wls = {'380' '412' '490' 'PAR'}; % This may need to be a user input at some point
%ocr_str  = 'OCR504 %s CHANNEL %02.0f %s [a0 a1 im],%s,%s,%s';
ocr_str  = 'OCR CHANNEL %02.0f %s (OCR504 %s) [a0 a1 im],%s,%s,%s';
a0 = cell(1,4);
a1 = a0;
im = a0;
ct = 1;

if sum(tf_OCR) == 1
    OCRcal_fn = cal_fnames{tf_OCR};
else
    OCRcal_fn = uigetfile([cal_fp],'Please choose OCR cal file');
% else
%     disp(['No OCR504 calibration found at: ',cal_fp,' for float ',mbari_id])
%     skip_OCR = 1;
end

if skip_OCR == 0
    copyfile(fullfile(cal_fp, OCRcal_fn), temp_dir);
    fid = fopen(fullfile(temp_dir, OCRcal_fn));
    tline = ' ';
    if contains(OCRcal_fn,'ocr') %ocr log file with specific coeff structure
        while ischar(tline)
            if regexp(tline,'^Serial Number','once')
                SN = regexp(tline,'\d+','match','once');
            end

            if regexp(tline,'^\ta0','once')
                a0{ct} = regexp(tline,'(?<=a0:)[\d\.e-]+','match','once');
            end

            if regexp(tline,'^\ta1','once')
                a1{ct} = regexp(tline,'(?<=a1:)[\d\.e-]+','match','once');
            end

            if regexp(tline,'^\tim','once')
                im{ct} = regexp(tline,'(?<=im:)[\d\.e-]+','match','once');
                ct = ct+1;
            end

            if ct == 5
                break
            end
            tline = fgetl(fid);
        end
    else %alternate format
%         keyboard
        while ischar(tline)
            if regexp(tline,'^SN','once')
                SN = regexp(tline,'\d+','match','once');
            end

            if regexp(tline,'OPTIC2')
                tline = fgetl(fid);
                thecoeffs = str2num(tline);
                a0{ct} = thecoeffs(1);
                a1{ct} = thecoeffs(2);
                im{ct} = thecoeffs(3);
                ct = ct+1;
            end

            if ct == 5
                break
            end
            tline = fgetl(fid);
        end 
    end
    fclose(fid);
    clear ct
    for i = 1:size(a0,2)
        cell_ct = 12+i;
        out{cell_ct} = sprintf(ocr_str, i, default_wls{i},  SN, a0{i}, ...
            a1{i}, im{i});
    end
end

% ************************************************************************
% PRINT OUT CELL ARRAY TO FILE

% First spool to screen for visual check
out = out(1:cell_ct);
D = out';
fprintf('\n')
fprintf('%s\n', D{:});
fprintf('\n')

config_fp = fullfile(config_dir, config_fn);

% Next check if config file already exists
% out_fn   = ['test',mbari_id,'_FloatConfig.txt'];
test_str = 'Y';
fchk     = ls(config_fp);

if ~isempty(fchk)
    disp(['Config file already exists: ',config_fp]);
    test_str = input(['Are you sure you want to overwrite this config ',...
        'file(Y/N)'],'s');
end

% Print to file
if regexpi(test_str,'^Y','once')
    fid = fopen(config_fp, 'w');
    fprintf(fid, '%s\r\n', D{:}); %
    fclose(fid);
end

% ************************************************************************
% LAST STEP GRAB SUNA CAL FILE, FORMAT AND BRING TO MBARI SIDE 
% ************************************************************************
disp(['Looking for SUNA nitrate calibration file(s): ',cal_fp]);
tf_NO3   = ~cellfun(@isempty, regexpi(cal_fnames,'SUNA.+CAL'));
if sum(tf_NO3)==0
    tf_NO3   = ~cellfun(@isempty, regexpi(cal_fnames,'SNA.+CAL')); %arg...
end

skip_NO3 = 0;

if sum(tf_NO3) == 1
    NO3cal_fn = cal_fnames{tf_NO3};
elseif sum(tf_NO3) > 1
    NO3cal_fn = uigetfile([cal_fp,'SUNA*CAL'],'Please choose SUNA cal file');
else
    disp(['No SUNA nitrate calibration found at: ',cal_fp,' for float ',mbari_id])
    skip_NO3 = 1;
end

if skip_NO3 == 0 % FILE EXISTS IN SIO FLOAT CAL DIR
    copyfile(fullfile(cal_fp, NO3cal_fn), temp_dir);
    temp_cal_fp = fullfile(temp_dir, NO3cal_fn);
    dest_cal_fp = fullfile(no3_dir, ncal_fn);
    
    % CHECK IF SUNA CAL FILE ALREADY EXISTS
    test_str = 'Y';
    fchk     = ls(dest_cal_fp);
    if ~isempty(fchk)
        disp(['NO3 cal file already exists: ', dest_cal_fp]);
        test_str = input(['Are you sure you want to overwrite this NO3 ',...
            'cal file(Y/N)'],'s');
    end

    if regexpi(test_str,'^Y','once')
        fid1 = fopen(temp_cal_fp); % FILE FROM SIO CAL DIR
        fid  = fopen([temp_dir, ncal_fn],'w'); % FILE FOR MBARI PROCESSING
        tline =' ';
        while ischar(tline)
            tline = fgetl(fid1);
            if ischar(tline) && contains(tline,'CORRECTABLE')
                fprintf(fid, '%s\r\n', tline) ;
                fprintf(fid, 'H,Original source file: %s\n', [cal_fp, NO3cal_fn]);
                fprintf(fid, 'H,Pixel base,1,,\n') ;
                fprintf(fid, 'H,Sensor Depth offset,0,,,\n') ;
                fprintf(fid, 'H,Br wavelength offset,210,,,\n') ;
                fprintf(fid, 'H,Min fit wavelength,,,,\n') ;
                fprintf(fid, 'H,Max fit wavelength,,,,\n') ;
                fprintf(fid, 'H,Use seawater dark current,No,,,\n') ;
                fprintf(fid, 'H,Pressure coef,0.0265,,,\n') ;
            else
                fprintf(fid, '%s\r\n', tline) ;
            end
        end
        fclose(fid);
        fclose(fid1);
        
        [successfinal,~,~] = copyfile([temp_dir, ncal_fn],no3_dir);
        if successfinal == 1
            disp(['SNA file successfully generated for float ', ...
                mbari_id,' and copied to ',no3_dir,'.'])
        end
    end
end





