%function tf = build_UW_NAVIS_config(mbari_fn)

% *************************************************************************
% *************************************************************************
% % To exract information from pdf files This function uses a function from
% the files exchange by Derek Wood:
% https://www.mathworks.com/matlabcentral/fileexchange/
% 63615-read-text-from-a-pdf-document
% Before pdf reader can work add a couple lines to startup.m to
% switch to the dir where the mfile is found and execute the following
% line: javaaddpath('iText-4.2.0-com.itextpdf.jar')
% it only needs to be executed once upon the start of a matlab session
%
% CHANGE HISTORY:

% 2024-07-10 - Nicola added a check for k0 and fP(0) comparisons to flag if
%   anything could be out of range and incorrect in the pH cal files.

% ************************************************************************
% CHOOSE FLOAT

% mbari_id = 'un1521';
mbari_id = 'un1538';


% ************************************************************************
% DO SOME PREP WORK
INST_ID     = regexp(mbari_id,'\d+','once','match'); % isolates float # w/o letters
config_fn   = [mbari_id,'_FloatConfig.txt']; % creates a text file name with float #
new_ncal_fn = [mbari_id,'.cal']; % creates the cal file name with the float #

% DEFINE DIRS STRUCTURE
user_dir     = getenv('USERPROFILE'); %returns user path,i.e.'C:\Users\eclark'
user_dir     = [user_dir, '\Documents\MATLAB\']; % sets local file path
dirs.cal     = [user_dir,'ARGO_PROCESSING\DATA\CAL\']; % sets local cal file path
dirs.temp    = 'C:\temp\'; % sets temp file
dirs.config  = [dirs.cal, 'FLOAT_CONFIG\staging\'];
dirs.msg     = fullfile('\\seaecho.shore.mbari.org\floats\UW\',... % sets location of float msg file
               regexprep(mbari_id,'u',''),'\');


dvec_now    = datevec(now); % sets exact time and date of run
% yr          = dvec_now(1); % use this to limit file search to yr & yr-1
% yr_rewind   = -3; % Number of previous yr files or directories to check for matches

O2_cal = {}; % creating a tbd space(cell array) for the O2 cal




% ************************************************************************
%            TRY AND FIND SUNA NO3 CALIBRATION FILE IN MSG FILE DIR
% ************************************************************************
tf_no3      = 1; % switch to 0 if NO3 cal file not found
no3_dir_str = '^NTR|^SUNA';
mdir     = dir(dirs.msg); % float msg file listing (files & dirs)
fnames   = {mdir.name}';
tf_dir   = cell2mat({mdir.isdir})';
tg_dirs  = ~cellfun(@isempty, regexpi(fnames, no3_dir_str,'once')) & tf_dir;
no3_dirs = fnames(tg_dirs); % root & subdirs for any potential suna cal dirs

% CHECK TO SEE IF PROCESSING NO3 CAL FILE ALREADY EXISTS
if tf_no3 == 1 && exist([dirs.config, new_ncal_fn],'file')
    fprintf('%s\n%s\n','SUNA cal file already exists!',[dirs.config, new_ncal_fn]);
    tf = input(['Do you want to replace the existing nitrate ', ...
        'calibration file (Y/N)?'],'s'); 
    if isempty(regexpi(tf,'^Y','once'))
        tf_no3 = 0;
        disp([dirs.config, new_ncal_fn,' will not be updated'])
    end
end

% STEP THROUGH DIRS & LOOK FOR SUNA CAL FILES
if tf_no3 == 1 % TRY AND COPY TO LOCAL
    no3_cals = cell(4,1);
    cal_ct   = 1;
    for ct = 1:size(no3_dirs,1)+1 % step through potential dirs
        if ct == 1
            fp = fullfile(dirs.msg);
        else
            fp = fullfile(dirs.msg, no3_dirs{ct-1},'\');
        end
        tmp = dir([fp,'SNA*cal']);
        if ~isempty(tmp)
            nct = size(tmp,1);
            no3_cals(cal_ct:cal_ct+nct-1) = strcat({tmp.folder}',{'\'},{tmp.name}');
            cal_ct = cal_ct+nct;
        end
    end
    no3_cals = no3_cals(1:cal_ct-1); % SUNA cal file listings
    clear cal_ct tf_dir tg_dirs no3_dirs ct fp

    if size(no3_cals,1) == 1
        ncal_fp = no3_cals{1};
    elseif size(no3_cals,1) >1 % multilple options
        [indx,~] = listdlg('PromptString',{['WARNING: ',...
            'mulitple SUNA NO3 cal files detected'], ...
            'Please choose correct one.',''},...
            'SelectionMode','single','ListString',no3_cals','ListSize',[450,80], ...
            'Name','CHOOSE SUNA CAL FILE');
        if isempty(indx) % operation was cancelled
            fprintf(['Multiple SUNA cal files detected but none were selected\n',...
                'No NO3 cal file will be produced for this float\n'])
            tf_no3 = 0;
        else
            ncal_fp = no3_cals{indx};
        end
    else % no SUNA cal files were found
        fprintf(['No SUNA cal files were detected. No NO3 cal file will be ',...
            'produced for this float\n'])
        tf_no3 = 0;
    end
end

if tf_no3 == 1 % TRY AND COPY TO LOCAL
    ncal_fn = regexpi(ncal_fp,'SNA.+cal','match','once');
    dest_fp = [dirs.temp, ncal_fn];
    [status,~] = copyfile(ncal_fp, dest_fp); % copy to local
    if status== 0
        fprintf(['WARNING: Could not copy SUNA cal file to local. No NO3 cal ', ...
            'file will be produced for this float\n'])
        tf_NO3 = 0;
    end
end

% FILE FOUND & LOCAL COPY SUCCESSFULL
% LAST BUT NOT LEAST - GENERATE NO3 CAL FILE FOR MBARI PROCESSING
if tf_no3 == 1
    add_hdr{1,1} = ['H,Original source file: ',ncal_fp,',,,,'];
    add_hdr{2,1} = 'H,Pixel base,1,,,';
    add_hdr{3,1} = 'H,Sensor Depth offset,0,,,';
    add_hdr{4,1} = 'H,Br wavelength offset,210,,,';
    add_hdr{5,1} = 'H,Min fit wavelength,,,,';
    add_hdr{6,1} = 'H,Max fit wavelength,,,,';
    add_hdr{7,1} = 'H,Use seawater dark current,No,,,';
    add_hdr{8,1} = 'H,Pressure coef,0.0265,,,';

    fid     = fopen(dest_fp);
    fid_new = fopen([dirs.config, new_ncal_fn],'w');
    tline = ' ';
    TCal_chk = 0;
    while ischar(tline)
        if regexp(tline,'^H,T_CAL','once') & TCal_chk == 0
            fprintf(fid_new,'%s\r\n', tline);
            for i = 1:size(add_hdr,1)
                fprintf(fid_new,'%s\r\n', add_hdr{i,1});
            end
            TCal_chk = 1;
        elseif regexp(tline,'^E|^H','once')
            fprintf(fid_new,'%s\r\n', tline);
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    fclose(fid_new);
    clear fid fid_new tline i add_hdr TCal_chk
end

% ************************************************************************
% ************************************************************************
%                     NOW BUILD CONFIG.TXT FILE
% ************************************************************************
% ************************************************************************

if exist([dirs.config,config_fn],'file') % see if it already exists
    fprintf('\n%s\n%s\n','Config file already exists!',[dirs.config,config_fn]);
    tf = input(['Do you want to replace the existing config ', ...
        'calibration file (Y/N)?'],'s');
    if ~strncmpi(tf,'Y',1)
        disp('Config file will not be updated')
        return
    end
end

config_cell = cell(20,1);
config_ct   = 0;
config_cell{config_ct+1} = sprintf('Institution ID: %s',INST_ID);
config_cell{config_ct+2} = sprintf('MBARI ID: %s',mbari_id);
config_ct = config_ct+2;


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
config_ct = config_ct+1;
config_cell{config_ct} = sprintf('Program: %s', Program);


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
config_ct = config_ct+1;
config_cell{config_ct} = sprintf('Region: %s', Region);


% ***********************************************************************
%                        ADD MSG FILE DIR
% ***********************************************************************
config_ct = config_ct+1;
config_cell{config_ct} = sprintf('%s',dirs.msg);


% ***********************************************************************
%    CHECK FOR OXYGEN CALIBRATION FILE (Always in main NAVIS msg dir?)
% ***********************************************************************
O2cal_fp  = [dirs.msg,'ox*calibration*'];
O2cal_fn = ls(O2cal_fp); % separate O2 filename

% if, by chance, there is more than 1 O2 cal, choose the right one
if size(O2cal_fn,1) > 1
    O2cal_fn = uigetfile(O2cal_fp,'Please choose correct O2 cal file');
elseif size(O2cal_fn,1) == 0 % NO O@ cal file found
    disp('No O2 cal file found - no cal data extracted')
end

if size(O2cal_fn,1) == 1 % only one file found now, size = zero means no O2 found
    fid   = fopen([dirs.msg, strtrim(O2cal_fn)]);
    d     = textscan(fid,'%s','Delimiter','\n');
    O2tmp = d{1,1}; % text file lines now in a cell array
    fclose(fid);

    % TEST CELL ARRAY FOR NEEDED INFO
    % MODEL
    t1       = ~cellfun(@isempty, regexp(O2tmp,'Model\s+\=','once'));
    O2_model = regexp(O2tmp{t1},'(?<=Model\s+\=\s+)\w+.+','match','once');

    % SN
    t1       = ~cellfun(@isempty, regexp(O2tmp,'Config\()\s+Serial\#\s+\=','once'));
    O2_SN    = regexp(O2tmp{t1},'(?<=Serial\#\s+\=\s+)\d+','match','once');

    str1 = sprintf('O2 sensor (%s %s) Temp coefficents ',regexprep(O2_model,' ',''),O2_SN);
    clear t1 O2_model O2_SN

    % TEMPERATURE COEFFICIENTS
    t1       = ~cellfun(@isempty, regexp(O2tmp,'TA\d{1}\s+\=','once'));
    Tmeta    = regexp(O2tmp(t1),'TA\d{1}','match','once');
    Tcoef    = regexp(O2tmp(t1),'(?<=TA\d{1}\s+\=\s+\+*)[\d\.\-e\+]+','match','once'); 
    str2     = strtrim(sprintf('%s ',Tmeta{:}));
    O2Tstr   = sprintf(['%s[%s]',repmat(',%0.7E',1,size(Tcoef,1))],str1, ...
               str2,str2double(Tcoef));
    config_ct = config_ct+1;
    config_cell{config_ct} = O2Tstr;
    clear  O2Tstr Tcoef Tmeta t1 str2

     % PHASE COEFFICIENTS
    t1       = ~cellfun(@isempty, regexp(O2tmp,'\s+[ABC]\d{1}\s+\=','once'));
    Pmeta    = regexp(O2tmp(t1),'[ABC]\d{1}','match','once');
    Pcoef    = regexp(O2tmp(t1),'(?<=\s+[ABC]\d{1}\s+\=\s+\+*)[\d\.\-e\+]+','match','once');
    str2     = strtrim(sprintf('%s ',Pmeta{:}));
    O2Pstr   = sprintf(['%s[%s]',repmat(',%0.7E',1,size(Pcoef,1))],...
        regexprep(str1,'Temp','Phase'), str2,str2double(Pcoef));
    config_ct = config_ct+1;
    config_cell{config_ct} = O2Pstr;

    clear  O2Pstr Pcoef Pmeta t1 str2
end


% ************************************************************************
%  TRY AND FIND pH CAL PDF & EXTRACT DATA (Always in main  msg dir?)
%                            NO! i.e un1116
% ************************************************************************
% STEP THROUGH DIRS & LOOK FOR SBE pH CAL FILES
tf_pH       = 1; % switch to 0 if NO3 cal file not found
pH_dir_str  = '^PH'; % sub dir regexp filter
pH_fn_str   = '^SBE\s+[A-Z]\d+.+pdf|^\d+-\d+-\d+.+pdf|^\d+-\d+.+pdf';
mdir        = dir(dirs.msg); % float msg file listing (files & dirs)
fnames      = {mdir.name}';
tf_dir      = cell2mat({mdir.isdir})';
tg_dirs     = ~cellfun(@isempty, regexpi(fnames, pH_dir_str,'once')) & tf_dir;
pH_dirs     = fnames(tg_dirs); % root & subdirs for any potential pH cal dirs

pH_cals = cell(10,1);
cal_ct   = 1;
for ct = 1:size(pH_dirs,1)+1 % step through potential dirs
    if ct == 1
        fp = fullfile(dirs.msg);
    else
        fp = fullfile(dirs.msg, pH_dirs{ct-1},'\');
    end

    tmp         = dir([fp,'*pdf']); % pdf file listing in msg dir
    t1          = ~cellfun(@isempty, regexp({tmp.name}',pH_fn_str)); %pH pdf logical
    pdf_listing = strcat({tmp.folder}',{'\'},{tmp.name}'); % pdf file paths
    pdf_listing = pdf_listing(t1); % pH pdf file paths


    if isempty(pdf_listing)
        fprintf('No pH pdf files match ph pdf file pattern for %s\n',mbari_id);
        fprintf('%s\n',fp);
    else
        nct = size(pdf_listing,1);
        pH_cals(cal_ct:cal_ct+nct-1) = pdf_listing;
        cal_ct = cal_ct+nct;
    end
end
t1 = ~cellfun(@isempty,pH_cals);
pH_cals = pH_cals(t1); % SUNA cal file listings
clear cal_ct tf_dir tg_dirs pH_dirs ct pdf_listing

if sum(t1) == 0
    fprintf('No pH pdf cal files found for %s\n',mbari_id);
    tf_pH = 0;
end

% IF YOU GET HERE A PDF CAL FILE EXISTS FOR pH
if tf_pH == 1
    if size(pH_cals,1) >1
        [indx,~] = listdlg('PromptString',{['WARNING: ',...
            'mulitple pH pdf cal files detected'], ...
            'Please choose correct one.',''},...
            'SelectionMode','single','ListString',pH_cals,'ListSize',[450,80], ...
            'Name','CHOOSE pH CAL FILE');
        if isempty(indx) % operation was cancelled
            fprintf(['Multiple pH cal files detected but none were selected\n',...
                'No pH cal data will be produced for this float\n'])
            tf_pH = 0;
        else
            ph_cal_fp = pH_cals{indx};
        end
    else
        ph_cal_fp = pH_cals{1};
    end
end

if tf_pH == 1
    SBE  = parse_cal_pdf(ph_cal_fp);
    k0   = SBE.pH.k0;
    k2   = SBE.pH.k2;
    fP   = SBE.pH.fP;

    % fP(0) / k0 range check
    if fP(1)-k0 >= 0.015 
        fP_diff = fP(1)-k0;
    elseif k0-fP(1) >= 0.015
        k0_diff = k0-fP(1);
    end

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
    config_ct = config_ct+1;
    config_cell{config_ct} = pH_str;
    clear  O2Pstr Pcoef Pmeta t1 str2


end

% ************************************************************************
%  TRY AND FIND MCOMS CAL PDF & EXTRACT DATA (Always in main  msg dir?)
%                            NO! i.e un1116
% ************************************************************************
% STEP THROUGH DIRS & LOOK FOR SBE MCOMS CAL FILES
tf_MCOMS       = 1; % switch to 0 if NO3 cal file not found
MCOMS_dir_str  = '^MCOMS'; % sub dir regexp filter
MCOMS_fn_str   = '^MCOMS';
mdir        = dir(dirs.msg); % float msg file listing (files & dirs)
fnames      = {mdir.name}';
tf_dir      = cell2mat({mdir.isdir})';
tg_dirs     = ~cellfun(@isempty, regexpi(fnames, MCOMS_dir_str,'once')) & tf_dir;
MCOMS_dirs  = fnames(tg_dirs); % root & subdirs for any potential pH cal dirs

MCOMS_cals = cell(10,1);
cal_ct   = 1;
for ct = 1:size(MCOMS_dirs,1)+1 % step through potential dirs
    if ct == 1
        fp = fullfile(dirs.msg);
    else
        fp = fullfile(dirs.msg, MCOMS_dirs{ct-1},'\');
    end

    tmp         = dir([fp,'*pdf']); % pdf file listing in msg dir
    t1          = ~cellfun(@isempty, regexp({tmp.name}',MCOMS_fn_str)); %pH pdf logical
    pdf_listing = strcat({tmp.folder}',{'\'},{tmp.name}'); % pdf file paths
    pdf_listing = pdf_listing(t1); % pH pdf file paths

    if isempty(pdf_listing)
        fprintf('No MCOMS pdf files match MCOMS pdf file pattern for %s\n',mbari_id);
        fprintf('%s\n',fp);
    else
        nct = size(pdf_listing,1);
        MCOMS_cals(cal_ct:cal_ct+nct-1) = pdf_listing;
        cal_ct = cal_ct+nct;
    end
end
t1 = ~cellfun(@isempty,MCOMS_cals);
MCOMS_cals = MCOMS_cals(t1); % MCOMS cal file listings
clear cal_ct tf_dir tg_dirs MCOMS_dirs ct pdf_listing

if sum(t1) == 0
    fprintf('No MCOMS calibration files were detected for %s\n',mbari_id);
    tf_MCOMS = 0;
end

if tf_MCOMS == 1
    param_str = 'CHL|\s+700|\(700|FDOM|CDOM';
    t1     = ~cellfun(@isempty, regexp(MCOMS_cals, param_str,'once'));
    if sum(t1) == 0 % try generic flavor name (1 file all cals)
        %generic_pat ='MCOMSC*-\d+\_Char\s*Sheet';
        generic_pat ='MCOMSC*-\d+\_*\s*Char\s*Sheet';
        t1     = ~cellfun(@isempty, regexp(MCOMS_cals, generic_pat,'once'));
    end

    if sum(t1) == 0
        fprintf(['No MCOMS calibration files match any calibration ', ...
            'file name patterens %s\n'], mbari_id);
        tf_MCOMS = 0;
    else
        MCOMS_cals = MCOMS_cals(t1); % Triplet (chl,bbp,fdom) cal file listings
    end
end

% ALL GOOD SOME APPROPRIATE CAL FILES EXIST - GET COEFFICIENTS
if tf_MCOMS == 1
    param_pat{1}  = 'MCOMS Chl fluorescence (%s %s) [ChlDC ChlScale],%0.0f,%0.4E';
    param_pat{2}  = 'MCOMS Backscatter700 (%s %s) [Betab700DC Betab700Scale],%0.0f,%0.4E';
    param_pat{3}  = 'MCOMS Fluorescent DOM (%s %s) [FDOMDC FDOMScale],%0.0f,%0.4E';

    % LOOP THROUGH TRIPLET PARAMS
    param_cell  = {'CHL' '\s+700|\(700' 'FDOM|CDOM'}; % to fimd file name
    param_field = {'Chl' 'bbp700' 'FDOM'}; % to find pdf parse cal field name
    for ct = 1:size(param_cell,2)
        t1 = ~cellfun(@isempty, regexp(MCOMS_cals, param_cell{ct}, 'once'));
        if sum(t1) > 1
            Pcell = {sprintf('WARNING: mulitple %s pdf cal files detected.',...
                param_field{ct}), 'Please choose the correct file.'};
            tmp_cals = MCOMS_cals(t1);
            [indx,~] = listdlg('PromptString', Pcell, ...
                'SelectionMode','single','ListString',tmp_cals,'ListSize', ...
                [450,80],'Name','CHOOSE CAL FILE');
            fp = tmp_cals{indx};
            if isempty(indx) % operation was cancelled
                fprintf(['Multiple %s cal files detected but none were selected\n',...
                    'No %s cal data will be produced for this float\n'], ...
                    param_field{ct}, param_field{ct})
                continue
            end
        elseif sum(t1) == 1
            fp = MCOMS_cals{t1};
        elseif sum(t1) == 0 && size(MCOMS_cals,1) == 1 % generic multipage pdf
            fp = MCOMS_cals{1};
        else
            fprintf(['No %s cal files detected!\n',...
                'No %s cal data will be produced for this float\n'], ...
                param_field{ct}, param_field{ct})
            continue
        end

        SBE  = parse_cal_pdf(fp);
        if ~isempty(SBE.(param_field{ct}).DC) & ~isempty(SBE.(param_field{ct}).Scale)
            MCOMS_str = sprintf(param_pat{ct}, SBE.(param_field{ct}).Model, ...
                SBE.(param_field{ct}).SN, SBE.(param_field{ct}).DC, ...
                SBE.(param_field{ct}).Scale);
            config_ct = config_ct+1;
            config_cell{config_ct} = MCOMS_str;

        else
            fprintf(['WARNING %s coeficients were not properly extracted ',...
                'from pdf. Check config file for errors\n'],param_field{ct})
        end
    end
end
config_cell = config_cell(1:config_ct);

% ************************************************************************
% ************************************************************************
%                      PRINT CONFIG FILE
% ************************************************************************
% ************************************************************************

fp  = [dirs.config,config_fn];
fid = fopen(fp,'w');
fprintf(fid,'%s\r\n',config_cell{:});
fclose(fid);

if exist('fP_diff')
    fprintf(['!!! CHECK PH CALS !!! fP(0) is larger by %.2f mV !!! \n'], (fP_diff*1000))
elseif exist('k0_diff')
    fprintf(['!!! CHECK PH CALS !!! k0 is larger by %.2f mV !!! \n'], (k0_diff*1000))
end




















