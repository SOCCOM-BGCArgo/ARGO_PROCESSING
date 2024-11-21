%function tf = build_WHOI_config(mbari_id)


% MODIFICATIONS:
% 03/12/2023 - JP - updated code to deal with variable length f(P) from API
%    call. Added warning block to detect if k2 is a function of P but no
%    code block yet to deal with it as not certain how it will be served by
%    API. Cleaned up NO3 cal file part to look in 3 different dir options
%    for souce file
% 04/24/23 - JP - updated code to deal with k2 f(P)
% 07/10/24 - Nicola added a check for k0 and fP(0) comparisons to flag if
%   anything could be out of range and incorrect in the pH cal files.

% ***********************************************************************
% TESTING
mbari_id = 'wn1567'; %1565 1567


% ***********************************************************************
disp(['Building config file for: ', mbari_id])

% SET UP DIRS
user_dir    = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir    = [user_dir, '\Documents\MATLAB\'];
config_dir  = [user_dir, 'ARGO_PROCESSING\DATA\CAL\FLOAT_CONFIG\staging\'];
no3_dir     = [user_dir, 'ARGO_PROCESSING\DATA\CAL\NO3_CAL\'];
msg_dir     = '\\seaecho.shore.mbari.org\floats\WHOI\';
%temp_dir    = 'C:\Users\eclark\Documents\temp\';
temp_dir    = 'C:\temp\';

% SET SOME DEFAULT INFO
WHOI_domain  = 'https://db.whoifloatgroup.org/';
api_call_str = '%sapi/%s?FLOAT_SERIAL_NO=%s&PLATFORM_TYPE=%s';

%GET INSTITUTE ID FROM MBARI ID
inst_id = regexp(mbari_id,'\d+','once','match');
if mbari_id(2) == 'n'
    float_type = 'NAVIS_EBR';
elseif mbari_id(2) == 'a'
    float_type = 'APEX'; % CHECK IF MOR TO STRING?
end

% ************************************************************************
% TRY AND GET WMO USING INSTITUTE ID & FLOAT TYPE USING API CALL
api_call = sprintf(api_call_str, WHOI_domain,'wmo', inst_id, float_type);
try
    info = webread(api_call);
    wmo = info.WMO;
catch
    disp(['Could not find WMO for ', mbari_id,'. Setting temporary WMO #'])
    wmo = ['NO_WMO_', mbari_id];
end


% NEXT TRY & GET CALIBRATION INFO
api_call = sprintf(api_call_str, WHOI_domain,'cal', inst_id, float_type);
try
    info = webread(api_call);
catch
    disp(['Could not find calibration data for ', mbari_id,'. EXITING!'])
    return
end

% *************************************************************************
% *************************************************************************
% IF YOU GET HERE CAL DATA EXIST
% *************************************************************************
% *************************************************************************
out = cell(30,1);
% CHOOSE PROGRAM AFFILIATION
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


% CHOOSE OCEAN REGION
% RegionList = {'SO' 'NPac' 'SPac' 'NAtl' 'IO' 'ARC' 'NOT DEFINED'};
RegionList = {'SO' 'NPac' 'EQPac' 'SPac' 'NAtl' 'SAtl' 'IO' 'ARC' 'NOT DEFINED'}; % Nicola added SAtl 12/27/2024
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

out{1} = ['Institution ID: ', inst_id];
out{2} = ['MBARI ID: ', mbari_id];
out{3} = ['Program: ', Program];
out{4} = ['Region: ', Region];
out{5} = '';

% GET MSG FILE DIRECTORY
d         = dir(msg_dir);
whoi_dirs = regexpi({d.name}','^w[ans]\d+','once','match');
whoi_dirs = whoi_dirs(~cellfun(@isempty, whoi_dirs)); % remove empty cells
t1 = strcmp(whoi_dirs, mbari_id);
if sum(t1) == 1
    out{6} = [msg_dir, mbari_id, '\'];
else
    out{6} = 'Msg file directory not resolved!';
end
out_ct = 7;

% *************************************************************************
% STEP THROUGH SENSORS
if isfield(info,'OPTODE_DOXY')
    MODEL = info.OPTODE_DOXY.SENSOR_MODEL;
    SN    = info.OPTODE_DOXY.SENSOR_SERIAL_NO;
    
    if regexp(MODEL,'^SBE63|^SBE83','once')
        O2 = info.OPTODE_DOXY.PREDEPLOYMENT_CALIB_COEFFICIENT;
        O2_temp_str = ['O2 sensor (%s %s) Temp coefficents ', ...
            '[TA0 TA1 TA2 TA3],%0.6e,%0.6e,%0.6e,%0.6e'];
        out{out_ct} = sprintf(O2_temp_str, MODEL, SN, O2.TA0, O2.TA1, ...
            O2.TA2, O2.TA3);
        out_ct = out_ct +1;
        
        O2_phase_str = ['O2 sensor (%s %s) Phase coefficents [A0 A1 A2 ',...
            'B0 B1 C0 C1 C2],%0.6e,%0.6e,%0.6e,%0.6e,%0.6e,%0.6e,%0.6e,%0.6e'];
        out{out_ct} = sprintf(O2_phase_str, MODEL, SN, O2.A0, O2.A1, ...
            O2.A2, O2.B0, O2.B1, O2.C0, O2.C1, O2.C2);
        out_ct = out_ct +1;
    end
end

if isfield(info,'FLUOROMETER_CHLA')
    MODEL = info.FLUOROMETER_CHLA.SENSOR_MODEL;
    SN    = info.FLUOROMETER_CHLA.SENSOR_SERIAL_NO;
    
    if regexp(MODEL,'^MCOMS','once')
        MODEL = 'MCOMS';
        cal   = info.FLUOROMETER_CHLA.PREDEPLOYMENT_CALIB_COEFFICIENT;
        tmp_str = ['MCOMS Chl fluorescence (%s %s) [ChlDC ChlScale],' ...
            '%0.0f,%0.4e'];
        out{out_ct} = sprintf(tmp_str, MODEL, SN, cal.DARK_CHLA, ...
            cal.SCALE_CHLA);
        out_ct = out_ct +1;
    end
end

if isfield(info,'BACKSCATTERINGMETER_BBP700')
    MODEL = info.BACKSCATTERINGMETER_BBP700.SENSOR_MODEL;
    SN    = info.BACKSCATTERINGMETER_BBP700.SENSOR_SERIAL_NO;
    
    if regexp(MODEL,'^MCOMS','once')
        MODEL = 'MCOMS';
        cal   = info.BACKSCATTERINGMETER_BBP700.PREDEPLOYMENT_CALIB_COEFFICIENT;
        tmp_str = ['MCOMS Backscatter700 (%s %s) [Betab700DC Betab700Scale],' ...
            '%0.0f,%0.4e'];
        out{out_ct} = sprintf(tmp_str, MODEL, SN, cal.DARK_BACKSCATTERING700, ...
            cal.SCALE_BACKSCATTERING700);
        out_ct = out_ct +1;
    end
end

if isfield(info,'FLUOROMETER_CDOM')
    MODEL = info.FLUOROMETER_CDOM.SENSOR_MODEL;
    SN    = info.FLUOROMETER_CDOM.SENSOR_SERIAL_NO;
    
    if regexp(MODEL,'^MCOMS','once')
        MODEL = 'MCOMS';
        cal   = info.FLUOROMETER_CDOM.PREDEPLOYMENT_CALIB_COEFFICIENT;
        tmp_str = ['MCOMS Fluorescent DOM (%s %s) [FDOMDC FDOMScale],' ...
            '%0.0f,%0.4e'];
        out{out_ct} = sprintf(tmp_str, MODEL, SN, cal.DARK_CDOM, ...
            cal.SCALE_CDOM);
        out_ct = out_ct +1;
    end
end



if isfield(info,'TRANSISTOR_PH')
    MODEL = info.TRANSISTOR_PH.SENSOR_MODEL;
    SN    = info.TRANSISTOR_PH.SENSOR_SERIAL_NO;
    
    if regexp(MODEL,'^SEAFET','once')
        MODEL = 'SBE';
        cal = info.TRANSISTOR_PH.PREDEPLOYMENT_CALIB_COEFFICIENT;
        cal_fields   = fieldnames(cal);
        cal_vals     = cell2mat(struct2cell(cal)); % struct > cell > values
        
        % Store f0 for calibration check later
        f0_tf        = strncmp(cal_fields,'F0',2);
        f0 = cal_vals(f0_tf);

        % excuded F0 coef
        tF0          = ~strcmp(cal_fields,'F0');
        cal_fields   = cal_fields(tF0); 
        cal_vals     = cal_vals(tF0);
        
        coef_ct      = size(cal_fields,1); %
        pcoef_tf     = strncmp(cal_fields,'F',1); % f(P) logical
        pcoef_fields = cal_fields(pcoef_tf);
        k2_tf        = strncmp(cal_fields,'K2',2); % k2 f(P) logical
        k0_tf        = strncmp(cal_fields,'K0',2); % k0 logical

%         if sum(k2_tf) ~=1 || max(size(cal.K2)) ~=1
%             fprintf(['WARNING: MULTIPLE K2 coefficients found - ', ...
%                 'add k2 f(p) code block!\n']);
%             return
%         end        

        % BUILD pH cal sting
        SNk0_str = sprintf('pH sensor coefficients (%s %s) [k0 ',MODEL, SN);
        % depending how WHOI's k2 f(P) is served could use k2 field names instead
        k2_str   = sprintf(repmat('k2%0.0f ',1 , sum(k2_tf)),(1:sum(k2_tf))-1); % start coef index at 0
        fP_str   = sprintf(repmat('%s ',1 , sum(pcoef_tf)), pcoef_fields{:}); 
        fP_str   = regexprep(fP_str,' $',']');
        fmt_str  = [',%0.5f',repmat(',%0.4e',1,coef_ct-1)]; % value format string

        out{out_ct} = sprintf([SNk0_str, k2_str, fP_str, fmt_str], ...
            cal_vals(k0_tf), cal_vals(k2_tf), cal_vals(pcoef_tf));
        out_ct = out_ct +1;
        

        % fP(0) / k0 range check
        if f0-cal_vals(k0_tf) >= 0.015
            fP_diff = f0-cal_vals(k0_tf);
        elseif cal_vals(k0_tf)-f0 >= 0.015
            k0_diff = cal_vals(k0_tf)-f0;
        end

%         str1 = 'pH sensor coefficients (%s %s) [k0 k2 k3xP ';
%         str2 = sprintf('k%0.0fxP^%0.0f ',[4:coef_ct;2:sum(pcoef_tf)]);
%         str2 = regexprep(str2,' $',']'); % replace last trailing space w/ bracket
%         str3 = [',%0.4f',repmat(',%0.4e',1,coef_ct-1)];
%         tmp_str = [str1, str2,str3];

%         out{out_ct} = sprintf(tmp_str, MODEL, SN, cal.K0, cal.K2, ...
%             cal_vals(pcoef_tf));
%         out_ct = out_ct +1;


       % pause % TESTING

        clear cal cal_fields cal_vals tF0 coef_ct pcoef_tf pcoef_fields
        clear k2_tf  str1 str2 str3 tmp_str
        clear  SNk0_str k2_str fP_str fmt_str k0_tf 
    end
end


% ************************************************************************
% PRINT OUT CELL ARRAY TO FILE

% First spool to screen for visual check
out = out(1: out_ct-1);
D = out';
fprintf('\n')
fprintf('%s\n', D{:});
fprintf('\n')

% Next check if config file already exists
out_fn   = [mbari_id,'_FloatConfig.txt'];
test_str = 'Y';
fchk     = ls([config_dir, out_fn]);

if ~isempty(fchk)
    disp(['Config file already exists: ',config_dir, out_fn]);
    test_str = input(['Are you sure you want to overwrite this config ',...
        'file(Y/N)'],'s');
end

% Print to file
if regexpi(test_str,'^Y','once')
    fid = fopen([config_dir,out_fn], 'w');
    fprintf(fid, '%s\r\n', D{:}); %
    fclose(fid);
end

% ************************************************************************
% NOW GRAB SUNA CAL FILE, FORMAT AND BRING TO MBARI SIDE %%%
% ************************************************************************
%rawdir      = out{6}; %defined above
% rawdir = [out{6}, 'calsheets\']
SNA_dest_fn = [mbari_id,'.cal'];
SNA_dest_fp = [config_dir, SNA_dest_fn];

% Double check that the msg dir exists and you have access
if exist(out{6},'dir')~=7
    disp(['No raw float directory identified for float ',mbari_id, ...
        '.  Cannot grab  nitrate SNA cal file.  Exiting.'])
    tf = 0;
    return
end

% CHECK IF SUNA CAL FILE ALREADY EXISTS
test_str = 'Y';
fchk     = ls([no3_dir, SNA_dest_fn]);
if ~isempty(fchk)
    disp(['NO3 cal file already exists: ', no3_dir, SNA_dest_fn]);
    test_str = input(['Are you sure you want to overwrite this NO3 ',...
        'cal file(Y/N)'],'s');
    if isempty(regexpi(test_str,'^Y','once'))
        fprintf('Nitrate cal file was not updated!\n')
        return
    end
end

% TRY AND FIND NO3 CAL FILE(S) IN WHOI FLOAT MSG DIR
NO3_dirs = {[out{6},'calsheets\SUNA\'], 'primary'; ...
            [out{6},'calsheets\SUNA\Calibration\'], 'secondary'; ...
            [out{6},'calsheets\'], 'final'};

for nct = 1: size(NO3_dirs,1)
    dir_info = dir([NO3_dirs{nct,1},'SNA*.cal']); 
    if isempty(dir_info)
        if nct == size(NO3_dirs,1)
            fprintf('NO3 cal not found in %s dir: %s - exiting\n', ...
                NO3_dirs{nct,2}, NO3_dirs{nct,1});
            return
        else
            fprintf('NO3 cal not found in %s dir: %s\n', ...
                NO3_dirs{nct,2}, NO3_dirs{nct,1});
            continue
        end
    else
        fprintf('NO3 cal found in %s dir: %s\n', ...
            NO3_dirs{nct,2}, NO3_dirs{nct,1});
        break
    end
end

% ARE THERE MULTIPLE NO3 CAL FILES TO CHOOSE FROM?
no3_fn = {dir_info.name}';
if size(no3_fn,1) > 1
    disp(['Multiple SNA files exist in ',dir_info(1).folder, ...
        ' for float ',mbari_id,'!'])
    % How often will this happen (perhaps if an updated cal)?  Maybe never. Could
    % add a block with dropdown for user to pick file, or leave as is.
    return
end

SUNA_sourcefile = fullfile(dir_info(1).folder, no3_fn{1});
[sunasuccess,~,~] = copyfile(SUNA_sourcefile, temp_dir);
fid1 = fopen([temp_dir, no3_fn{1}]);
fid = fopen(SNA_dest_fp,'w');
k=0;
while ~feof(fid1)
    k=k+1;
    % line = fgetl(fid1); % switched "line" str name to "tline". line is a
    % function name
    tline = fgetl(fid1);
    if contains(tline,'CORRECTABLE')
        fprintf(fid, '%s\r\n', tline) ;
        fprintf(fid, 'H,Original source file: %s\n', SUNA_sourcefile);
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

fprintf('SUNA cal file successfully generated for float %s\n',mbari_id)

%[successfinal,~,~] = copyfile(SNA_dest_fp,config_dir);
% if successfinal == 1
%     disp(['SNA file successfully generated for float ',mbari_id,' and copied to ',no3_dir,'.'])
% end

% Print pH cal warning if necessary
if exist('fP_diff')
    fprintf(['!!! CHECK PH CALS !!! fp(0) is larger by: %.2f mV !!! \n'], (fP_diff*1000))
elseif exist('k0_diff')
    fprintf(['!!! CHECK PH CALS !!! k0 is larger by: %.2f mV !!! \n'], (k0_diff*1000))
end

