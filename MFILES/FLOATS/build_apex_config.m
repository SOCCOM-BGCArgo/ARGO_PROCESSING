%function tf = build_apex_config(mbari_fn)
% ad optional input for MSC which would circumvent the UW search for the
% MSC

% 1) Get MSC # for float from Dana's Inventory list using APF ID # as index
% 2) Use MSC# to find nitrate cal file on CHEM
% 3) Rename nitrate cal file & add needed header lines
% 4) Get O2 cal lines from text file in float msg directory
% 5) Get pH cal data from excel spread sheet using MSC# as index
% 6) Check for FFLBB bdf cal files in float msg directory
% 7) Extract FLBB cal info
% 8) Generate float config file

% TESTING
%mbari_fn = '18771SOOCN';
%mbari_fn = '18320SOOCN';
%mbari_fn = '18721SOOCN';
%mbari_fn = '18545SOOCN';
%mbari_fn = '18098SOOCN';

%mbari_fn = '18299SOOCN';
%mbari_fn = '18994SOOCN';
%mbari_fn = '18821SOOCN';
%mbari_fn = '18739SOOCN';
%mbari_fn = '17898SOOCN';
%mbari_fn = '18082SOOCN';
%mbari_fn = '18013SOOCN';
%mbari_fn = '18864SOOCN';

%mbari_fn = '18081SOOCN';
%mbari_fn = '18796SOOCN';
%mbari_fn = '18110SOOCN';
%mbari_fn = '18852SOOCN';
%mbari_fn = '18862SOOCN';
%mbari_fn = '18417SOOCN';

mbari_fn = '18643SOOCN';
mbari_fn = '18861SOOCN';
mbari_fn = '18829SOOCN';
mbari_fn = '18169SOOCN';
mbari_fn = '18097SOOCN';



% ************************************************************************
% DEFINE DIRS STRUCTURE
user_dir     = getenv('USERPROFILE'); %returns user path,i.e.'C:\Users\jplant'
user_dir     = [user_dir, '\Documents\MATLAB\'];
%dirs.cal     = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
dirs.cal     = 'C:\temp\'; % TESTING
dirs.temp    = 'C:\temp\';
dirs.msg     = '\\atlas\ChemWebData\floats\';
dirs.ph      = '\\atlas\Chem\DuraFET\APEX+pH Calibrations\pHlogFiles\';
dirs.isus    = '\\atlas\Chem\ISUS\Carole\';


% ************************************************************************
% DO SOME PREP WORK
APEX_ID_str = regexp(mbari_fn,'^\d+','once','match');
config_fn   = [mbari_fn,'_FloatConfig.txt'];
ncal_fn     = [mbari_fn,'.cal'];

dvec_now    = datevec(now);
yr          = dvec_now(1); % use this to limit file search to yr & yr-1
isus_inv_hdr = {'WrcId' 'ApfId' 'SbeId' 'IsusId' 'OptodeId' 'FlbbId' , ...
    'DeploymentOpportunity' 'action' 'Specification' 'Deploy location'};
isus_inv = cell(300,size(isus_inv_hdr,2));
iAPF = find(strcmp('ApfId', isus_inv_hdr) == 1);
iMSC = find(strcmp('IsusId', isus_inv_hdr) == 1);

O2_cal ={};

% ************************************************************************
% ************************************************************************
% GET MSC # FOR A GIVEN APF ID # FROM UW
% CHECK ISUS INVENTORY FOR PREVIOUS & PRESENT YEARS
for yr_ct =-1:0
    yr_str      = num2str(yr+yr_ct,'%0.0f');
    
    url = ['http://runt.ocean.washington.edu/swift/Argo',yr_str, ...
           'Logistics/IsusInventory.mbari'];
    url_fn = ['UW_IsusInventory_',yr_str];
    websave(['C:\temp\',url_fn,'.txt'],url);
    disp(['Checking ',url, ' for APF # = ',APEX_ID_str])
    
    % PARSE URL FILE
    fid     = fopen([dirs.temp,url_fn,'.txt']);
    tline   = ' ';
    chk     = 0;
    line_ct = 1;
    fixed_w = [];
    
    % DATA IS MOSTLY FIXED WIDTH - FIGURE OUT WIDTHS BY SCANNING FOR
    % LARGEST FIXED WIDTH COLUMN COUNT THEN GO THROUGH AGAIN & PARSE
    % chk = 1 Start of APEX float lines, chk = 2 end float lines rewind
    % chk = 3, start of float lines again
    while ischar(tline)
        if chk == 0 && ~isempty(regexp(tline,'^\s+#---','once'))
            chk = 1; % entering list of APEX floats
        elseif chk == 1 && ~isempty(regexp(tline,'\S+','once'))
            % start scanning for largest fixed with column count
            tmp = regexp(tline,'\S+', 'end');
            if size(tmp,2) > size(fixed_w,2)
                fixed_w = tmp;
                tmp_str = tline;
            end
        elseif chk == 1 && isempty(regexp(tline,'\S+','once'))
            chk = 2;
            fixed_w =[0, fixed_w]; % add for indexing
            frewind(fid);
        elseif chk == 2 && ~isempty(regexp(tline,'^\s+#---','once'))
            chk = 3;
        elseif chk == 3 && ~isempty(regexp(tline,'\S+','once'))
            for i = 1:size(isus_inv_hdr,2)
                if i == size(isus_inv_hdr,2)
                    str = strtrim(tline(fixed_w(i)+1:end));
                else
                    str = strtrim(tline(fixed_w(i)+1:fixed_w(i+1)));
                end
                if isempty(str) && i ==1 % very incomplete  so skip
                    line_ct = line_ct -1;
                    break
                elseif isempty(str)
                    str = 'NaN';
                end
                isus_inv{line_ct,i} = str;
            end
            line_ct = line_ct +1;
        elseif chk == 3 && isempty(regexp(tline,'\S+','once'))
            break
        end
        tline = fgetl(fid);
    end
    fclose(fid);

    % NOW CHECK FOR APEX ID #
    tf = strcmp(isus_inv(:,iAPF), APEX_ID_str);
    if sum(tf) > 0
        disp([APEX_ID_str, ' found in ',url_fn])
        break
    end
end
isus_inv = isus_inv(1:line_ct-1,:);
clear i dvec_now ans chk fid fixed_w tline str line_ct tmp tmp_str yr_ct

MSC_str = isus_inv{tf,iMSC};
MSC = str2double(MSC_str);
if isnan(MSC)
    disp(['NO MSC # found for ',APEX_ID_str, ' on runt'])
    return
end

%TESTING
% MSC = 807;
% MSC_str = '807';
% ************************************************************************
% ************************************************************************
% NOW TRY AND FIND ISUS CALIBRATION FILE
% LOOK IN MASTER LIST FIRST, IF NOT FOUNDTHEN HUNT WITH MSC# TO FIND
% DIECTORY & THEN CALIBRATION FILE
mbari_ncals_fn = 'ISUSFloatCalibrationListing.xlsx';
[status,msg] = copyfile([dirs.isus, mbari_ncals_fn], dirs.temp); % copy to local

[~,~,ncals] = xlsread([dirs.temp, mbari_ncals_fn], ...
     'ISUSCalFilenames','A2:C200'); % AH500 # may need to be extended
ncals_hdr = {'MSC' 'dir path' 'file name'};
ncals = ncals(~cellfun(@isnan,ncals(:,1)),:); % remove NaN MSC rows
tMSC = cell2mat(ncals(:,1)) == MSC;

% ALL GOOD BUT CHECK THAT FILE EXISTS TOO
if sum(tMSC) == 1  && exist([ncals{tMSC,2},'\',ncals{tMSC,3}],'file')
    disp(['ISUS calibration file found in master list for  MSC # ',MSC_str, ...
        ' and APEX ID # ',APEX_ID_str,':'])
    cal_fn = ncals{tMSC,3};
    fp = [ncals{tMSC,2},'\',cal_fn]; % path to ncal file
else
    disp(['Nitrate cal for MSC# ',MSC_str,'not found in ',mbari_ncals_fn])
    disp('Searching directory structure for a MSC match ...')

    ncal_dir ={};
    for yr_ct = -1:0 % look in current & previous year for file
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
        return
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
    
% ************************************************************************
% NO3 cal file found copy to local, add header info & rename

tf ='Y';
if exist([dirs.cal,'NO3_CAL\',ncal_fn],'file')
    str = [dirs.cal,'NO3_CAL\',ncal_fn,' already exists!'];
    disp(str)
    tf = input(['Do you want to replace the existing nitrate ', ...
                'calibration file (Y/N)?'],'s');
end

if regexpi(tf,'^Y','once')
    add_hdr{1,1} = ['H,Original source file: ',fp,',,,,'];
    add_hdr{2,1} = 'H,Pixel base,1,,,';
    add_hdr{3,1} = 'H,Sensor Depth offset,0,,,';
    add_hdr{4,1} = 'H,Br wavelength offset,210,,,';
    add_hdr{5,1} = 'H,Min fit wavelength,,,,';
    add_hdr{6,1} = 'H,Max fit wavelength,,,,';
    add_hdr{7,1} = 'H,Use seawater dark current,No,,,';
    add_hdr{8,1} = 'H,Pressure coef,0.0265,,,';
    
    [status,msg] = copyfile(fp,dirs.temp);

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
else
    disp([dirs.cal,'NO3_CAL\',ncal_fn,' already exists and will not be updated'])
end

% ************************************************************************
% ************************************************************************
% NOW BUILD CONFIG.TXT FILE
tf ='Y';
if exist([dirs.cal,config_fn],'file') % see if it already exists
    str = [dirs.cal,config_fn,' already exists!'];
    disp(str)
    tf = input(['Do you want to replace the existing config ', ...
                'calibration file (Y/N)?'],'s');
    if ~strncmpi(tf,'Y',1)
        disp([dirs.cal,config_fn,' already exists and will not be updated'])
        return
    end      
end

% ************************************************************************
% CHECK FOR OXYGEN CALIBRATION FILE
msg_path =[dirs.msg,'f',APEX_ID_str,'\'];
O2cal_fn = ls([msg_path,'ox*calibration*']);

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
    while ischar(tline)
        if regexp(tline,'PhaseCoef','once') % build SN line
            tmp = regexp(tline(regexp(tline,'PhaseCoef\s+','end','once')...
                  :end),'\d+','match');
            O2_ct = O2_ct+1;
            str = ['OptodeSn = ',tmp{1},' ',tmp{2}];
            O2_cal{O2_ct,1} = str;
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

% ************************************************************************
% ************************************************************************
% TRY AND GET PH CAL FROM excell SHEET
disp('Trying to extract pH calibration data ...')
[status,msg] = copyfile([dirs.ph, 'pHlog.xlsx'], dirs.temp); % copy to local
PH_hdr ={'DF#' 'APEX#' 'MSC#' 'k2 (T)' 'k0 (HCl)' 'k6' 'k5' 'k4' 'k3', ...
         'k2' 'k1', 'k0' 'R^2' 'k0 (pOFF)' 'k0 (pON)'};
iDF     = find(strcmp('DF#',       PH_hdr) == 1);
iAPX    = find(strcmp('APEX#',     PH_hdr) == 1);
iMSC    = find(strcmp('MSC#',      PH_hdr) == 1);
iK6     = find(strcmp('k6',         PH_hdr) == 1);
iK2T    = find(strcmp('k2 (T)',    PH_hdr) == 1);
iK0_HCL = find(strcmp('k0 (HCl)',  PH_hdr) == 1);
iK0_ON  = find(strcmp('k0 (pON)',  PH_hdr) == 1);
iK0_OFF = find(strcmp('k0 (pOFF)', PH_hdr) == 1);

ph_ct = 0;

% GET data from CALIBRATION SUMMARY sheet in pHlog.xlsx
[PH_num,PH_txt,PH_raw] = xlsread([dirs.temp, 'pHlog.xlsx'], ...
    'CALIBRATION SUMMARY','A5:AH500'); % AH500 # may need to be extended
PH_num = PH_num(:,[2,4,6,16,17,19:26,32,34]);
PH_txt = PH_txt(:,2:3); % Just to Grab any non numeric Durafet ID's (SBE0
tf = isnan(PH_num(:,iMSC));
PH_num(tf,:) =[]; % remove all rows with no MSC values
PH_txt(tf,:) =[]; % remove all rows with no MSC values

% NOW try and find MSC value
tMSC = PH_num(:,iMSC) == MSC;
if sum(tMSC,1) == 1 % ALL good only 1 line recovered
    MSC_ph_num = PH_num(tMSC,:);
    MSC_ph_txt = PH_txt(tMSC,:);
    DF_num = MSC_ph_num(iDF);
    DF_txt = MSC_ph_txt{iDF};
    DF_APX = MSC_ph_txt{iAPX};
    
    % TRY AND GET DF ID
    if isnan(DF_num) && isempty(DF_txt) 
        disp(['NO DF ID found for MSC# ',MSC_str])
        return
    elseif isnan(DF_num)
        DF_ID = DF_txt;
    else
        DF_ID = num2str(DF_num,'%0.0f');
    end
    
    
    % CONFIRM IDENTICAL APEX ID NUMBER'S FOR MSC #
    if isnan(MSC_ph_num(iAPX))
        disp(['WARNING: NO APEX ID # found in pHlog.xlsx for MSC# ',...
            MSC_str]);
    elseif MSC_ph_num(iAPX) ~= str2num(APEX_ID_str)
        disp(['APEX ID in pHlog.xlsx (',num2str(MSC_ph_num(iAPX)),') ',...
            'does not match ISUS Inventory (',APEX_ID_str,') for ',...
            'MSC# ',MSC_str])
        return
    end
    
    ph_cal{1,1} = sprintf('Durafet SN = %s, APEX# %s',DF_ID, APEX_ID_str);
    
    % CHECK FOR LAB KO (PUMP ON IS DEFAULT VALUE TO USE)
    K0_ON  = MSC_ph_num(iK0_ON);
    K0_OFF = MSC_ph_num(iK0_OFF);
    K0_HCL = MSC_ph_num(iK0_HCL);
    if ~isnan(K0_ON)
        k0_str = sprintf(['%0.5f, =k0 Pon; %0.5f =k0 Poff; ',...
            '%0.5f =k0 HCl'],K0_ON, K0_OFF, K0_HCL);
    elseif isnan(K0_ON) && ~isnan(K0_HCL)
        tf = input(['No lab k0 was found! Do you want to build pH ',...
            'calibration using HCL k0 Y/N)?'],'s');
        if ~strncmpi(tf,'y',1)
            disp('No config file built')
            return
        end
        k0_str = sprintf(['%0.5f, =k0 HCl; WARNING: PROVISONAL ',...
            'No lab k0 was found'],K0_HCL);
    else
        disp(['No k0 could be found for pH calibration - add data to ',...
            'pHlog.xlsx']);
    end
    ph_cal{3,1} = k0_str;
    ph_cal{4,1} = sprintf('%0.5f, =k2*T',MSC_ph_num(iK2T));
    
    % CHECK & BUILD Pcoef strings
    P = (flip(MSC_ph_num(iK6:iK6+5)))';
    if sum(~isnan(P),1) == 0
        disp(['No P Coefs found in pHlog.xlsx for MSC# ',MSC_str])
        return
    end
    ph_cal{2,1} = sprintf(['%0.0f, number of calibration ', ...
        'coefficients'],size(P,1)+2);
    for i = 1:size(P,1)
        ph_cal{4+i,1} = sprintf(['%0.4E, =k',num2str(i+2),'*P^',...
            num2str(i)],P(i));
    end
elseif sum(tMSC,1) > 1
    disp(['More than one pH cal row found for MSC# ',MSC_str])
    return
else
    disp(['NO pH cal row found for MSC# ',MSC_str])
    return
end
clear DF_num  DF_txt DF_APX K0_ON K0_OFF K0_HCL
clear iDF iAPX iMSC iK6 iK2T iK0_HCL iK0_ON iK0_OFF ph_ct

% ************************************************************************
% OK NOW CHECK FOR FLBB PDF CAL SHEETS & TRY & EXTRACT SCALE & DARK COUNTS
% This uses a function from the files exchane by Derek Wood:
% https://www.mathworks.com/matlabcentral/fileexchange/
% 63615-read-text-from-a-pdf-document
% Before pdf reader can work andd a cuople lines to startup.m to
% switch to the dir where the mfile is found and execute the following
% line: javaaddpath('iText-4.2.0-com.itextpdf.jar')
% it only needs to be executed once upon the start of a matlab session
disp(['Looking for bio-optical calibration data in: ',msg_path]);
clear optics_list
optics_cal  = {};

% ANY PDF FILES IN THE MSG DIR
tmp = dir([msg_path,'*pdf']);
if isempty(tmp)
    disp(['No PDF files found in: ',msg_path]);
    disp('No bio-optical calibration data will be extracted');
elseif size(tmp,1) == 2 % 2 pdf files, probably FLBB cal files
    optics_list = {tmp.name}';
    tf_chl = ~cellfun(@isempty,regexp(optics_list,'CHL','once'));
    tf_bbp = ~cellfun(@isempty,regexp(optics_list,'700','once'));
else
    disp(['More than 2 pdf files ',msg_path]);
    disp(optics_list)
    disp('Please manually add optics cal data to config file');
end

% TRY & GET CHL CAL
if exist('optics_list','var') && sum(tf_chl,1) > 0
    chl_fn  = optics_list{tf_chl};
    chl_txt = pdfRead([msg_path, chl_fn]); % get text from pdf file
    chl_txt = chl_txt{1,1};
    % THE NEXT 2 LINES ARE IMPORTANT! The pdf reader pulls in some odd
    % hidden characters. The reg exp line only grabs common visible
    % characters & leaves the hidden ones behind via indexing
    inds = regexp(chl_txt,'[\w.-() /=]'); % this removes hidden chars in str
    chl_txt = chl_txt(inds);
    chl_dc_str  = regexp(chl_txt,'\d+(?=\s+counts)','match','once');
    chl_scale_str = regexp(chl_txt,'\d\.\d+(?=\s+.g/l/count)','match','once');
    chl_sn = regexp(chl_txt,'FLBB\w+\-\d+','match','once');
    % DO SOME SANITY CHECKS ON EXTRACTED VALUS
    chl_dc    = str2double(chl_dc_str);
    chl_scale = str2double(chl_scale_str);
    tbad_chl = isnan(chl_dc) || isnan(chl_scale) || chl_dc < 1 || ...
        chl_dc > 99 || chl_scale < 0.003 || chl_scale > 0.02;
end

% TRY & GET BBP CAL
if exist('optics_list','var') && sum(tf_bbp,1) > 0
    bbp_fn = optics_list{tf_bbp};
    bbp_txt = pdfRead([msg_path, bbp_fn]);
    bbp_txt = bbp_txt{1,1};
    inds = regexp(bbp_txt,'[\w.-() /=]'); % this removes hidden chars in str
    bbp_txt = bbp_txt(inds);
    bbp_dc_str  = regexp(bbp_txt,'\d+(?=\s+counts)','match','once');
    bbp_scale_str = regexp(bbp_txt,'\d\.\d+E-\d+(?=\(m-1sr-1)','match','once');
    bbp_sn = regexp(chl_txt,'FLBB\w+\-\d+','match','once');  
    % DO SOME SANITY CHECKS ON EXTRACTED VALUS
    bbp_dc    = str2double(bbp_dc_str);
    bbp_scale = str2double(bbp_scale_str);
    tbad_bbp = isnan(bbp_dc) || isnan(bbp_scale) || bbp_dc < 1 || ...
        bbp_dc > 99 || bbp_scale < 0.5e-6 || bbp_scale > 10e-6;
end

% BUILD BIO OPTICAL CAL STRINGS
if exist('optics_list','var') && ~tbad_bbp && ~tbad_chl % valid numbers extraxted
    optics_cal{1,1} = sprintf('CHLFLUOR SN %s',chl_sn);
    optics_cal{2,1} = sprintf('%s ChlDC',chl_dc_str);
    optics_cal{3,1} = sprintf('%s ChlScale',chl_scale_str);
    optics_cal{4,1} = sprintf('%s BetabDC',bbp_dc_str);
    optics_cal{5,1} = sprintf('%s BetabScale',bbp_scale_str);
    disp('Bio optical calibration coefficients have been extracted')
elseif  exist('tbad_bbp','var') || exist('tbad_chl','var')
    disp('PDF bioptical extraction did not seem reasonable- add mannually')
end
    



% ************************************************************************
% IF YOU GET TO HERE O2 CAL & PH CAL HAVE BEEN SUCCESFULLY EXTRACTED
% PRINT TO FILE
fid = fopen([dirs.cal,config_fn], 'w');
fprintf(fid,'%s\r\n',APEX_ID_str);
fprintf(fid,'%s\r\n',mbari_fn);
fprintf(fid,'%s\r\n','\\Atlas\chemwebdata\floats\f');
fprintf(fid,'%s\r\n','\\Sirocco\wwwroot\lobo\Data\FloatVizData\');
fprintf(fid,'0\r\n');
out = [O2_cal; ph_cal; optics_cal];
for i = 1:size(out,1)
    fprintf(fid,'%s\r\n',out{i,1});
end
fclose(fid);

% SPOOL TO SCREEN AS CHECK
fid = fopen([dirs.cal,config_fn]);
tline = ' ';
while ischar(tline)
    disp(tline)
    tline = fgetl(fid);
end
fclose(fid);
    


        

            



