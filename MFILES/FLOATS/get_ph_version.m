function d = get_ph_version()

% SCRIPT TO BUILD LIST OF PH SENSOR SN & VERSION NUMBER
d =[];

% UPDATED BY JP ON 08/07/2021 - 
%     1) Switched excel reading from xlsread() to readcell() since xlsread
%        will be phased out by matlab.
%     2) Used  to require a valid UW ID for float to be included in output
%        table, but this canbe empty now, all that is required is a valid
%        DF SN & MSC number. pH cal GUI uses this function and UW ID may
%        not always be available at calibration time
%     3) Better clean up and accounting of duplicate entries, invalid
%        entries, wrong variable time (lots of mixed types!!!
%     4) Added stem versions 5 & 6
% 09/24/2021 - JP - Pointed to new location for pHlog.xlsx
% 03/22/2022 - JP - Updated parse to new pHlog.xlsx format. Up to 8th order
%        f(P) function & k2 now up to 6th order function of pressure too

% ************************************************************************
% SET DIR STRUCTURE
% ************************************************************************
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

dirs.temp = 'C:\temp\';
%dirs.ph  = '\\atlas\Chem\DuraFET\APEX+pH Calibrations\pHlogFiles\';
dirs.ph   = '\\atlas\chem\DuraFET\CALIBRATION\';
%dirs.ph   = 'C:\Users\jplant\Documents\MATLAB\Apps\pH_Calibration\test data\'; % TESTING
dirs.fn   = 'pHlog.xlsx'; 
%dirs.save = user_dir;
dirs.save = [user_dir,'\PH_DIAGNOSTICS\'];

% ************************************************************************
% LOAD AND PARSE pH_log file
% This is a messy excel sheet, condensing down to useful info for APEX
% floats and then adding version defination and descriptions
% ************************************************************************

disp('Trying to extract info from pHlog.xlsx calibration file ...')
[status,msg] = copyfile([dirs.ph, dirs.fn], dirs.temp); % copy to local

%ph_ver_hdr = {'DF ID' 'UW ID' 'MSC' 'K2' 'K0_pON' 'Build notes'};
ph_ver_hdr = {'DF ID' 'UW ID' 'MSC', ...
    'k2 C6' 'k2 C5' 'k2 C4' 'k2 C3' 'k2 C2' 'k2 C1' 'k2 C0', ...
    'f(P) C8' 'f(P) C7' 'f(P) C6' 'f(P) C5' 'f(P) C4' 'f(P) C3' 'f(P) C2' 'f(P) C1' 'f(P) C0', ...
    'k0 pON' 'k0 pOFF' 'k0 HCl' 'Build notes'}; % output header

C    = readcell([dirs.temp, dirs.fn],'Sheet','CALIBRATION SUMMARY'); 
% REMOVE ROWS & COLS WITH ALL <missing>

% ***********************************************************************
% FIND HEADER ROW  GET INDICES & DO 1st ROUND OF CLEANING /CONDENSING
t_hdr      = find(strcmp(C(:,2),'DF#') == 1); % index not a logical
excel_hdr  = C(t_hdr,:);

iDF     = find(strcmp(excel_hdr,'DF#') == 1, 1);
iAPX    = find(strcmp(excel_hdr,'APEX#') == 1);
iMSC    = find(strcmp(excel_hdr,'MSC#') == 1);
iK2     = find(strncmp(excel_hdr,'k2 f(P)',7) == 1);
iFP     = find(strncmp(excel_hdr,'f(P)',4) == 1);
iK0_ON  = find(strncmp(excel_hdr,'k0 (pON)',8) == 1);
iK0_OFF = find(strncmp(excel_hdr,'k0 (pOFF)',9) == 1);
iK0_HCL = find(strncmp(excel_hdr,'k0 (HCl)',8) == 1);
iVER    = find(strcmp(excel_hdr,'Sensor version') == 1);


keep_inds = [iDF iAPX iMSC iK2 iFP iK0_ON iK0_OFF iK0_HCL iVER];
C = C(t_hdr+1:end,keep_inds);

tf_char  = cellfun(@ischar,C); % find characters
tf_num   = cellfun(@isnumeric,C); % find numbers
tg       = tf_char | tf_num; % either a number or a charcter
tg_rows   = sum(tg,2) > 0;


C = C(tg_rows,:);
clear tf_char tf_num tg tg_cols tg_rows t_hdr keep_inds


% ***********************************************************************
% REASSIGN INDICES & 2nd ROUND OF CLEANING
iDF     = find(strcmp(ph_ver_hdr,'DF ID') == 1, 1);
iAPX    = find(strcmp(ph_ver_hdr,'UW ID') == 1);
iMSC    = find(strcmp(ph_ver_hdr,'MSC') == 1);
iK2     = find(strncmp(ph_ver_hdr,'k2',2) == 1);
iFP     = find(strncmp(ph_ver_hdr,'f(P)',4) == 1);
iK0_ON  = find(strcmp(ph_ver_hdr,'k0 pON') == 1);
iK0_OFF = find(strcmp(ph_ver_hdr,'k0 pOFF') == 1);
iK0_HCL = find(strcmp(ph_ver_hdr,'k0 HCl') == 1);
iVER    = find(strcmp(ph_ver_hdr,'Build notes') == 1);

% VALID ENTRIES MUST HAVE MSC & APEX ID's
tf_DF  = cellfun(@ischar,C(:,iDF)) | cellfun(@isnumeric,C(:,iDF));
tf_MSC  = cellfun(@isnumeric,C(:,iMSC)); % find numbers only
tg = tf_DF & tf_MSC;
C = C(tg,:);
clear tf_DF tf_MSC tg

rC   = size(C,1);
tSN  = ones(rC,1) *0; % predim
tMSC = tSN;

% ***********************************************************************
% !!!! 3rd ROUND OD CLEANING MESSY DATA !!!!
% STEP THROUGH ROWS & MAKE LOGICAL ARRAYS. DEAL WITH <missing> CELLS &
% OTHER FIXES DUE TO CRAZINESS OF pHlog.xlsx!!!
for i = 1:rC
    try
        % DF SN FILTERS:
        % Is DF SN a number if so  OK
        if isnumeric(C{i,iDF})
            tSN(i) = 1;
        % Does DF SN start with a "C" or "F" followed bt a #, early SBE stems: OK
        elseif ischar(C{i,iDF}) && ~isempty(regexp(C{i,1},'^[CF]*\d+','once'))
            tSN(i) = 1;
        end
        
        % MSC FILTER: Is MSC a number? If so OK
        if isnumeric(C{i,iMSC})
            tMSC(i) = 1;
        end
        
        % UW ID FILTERS:
        % Is UW ID a <missing> value? If so set = NaN
        if ismissing(C{i,iAPX})
            C{i,iAPX} = NaN;
        % FIND UW ID WITH "CF" suffix (for "carbon fiber") & convert to #
        elseif ischar(C{i,iAPX}) && ~isempty(regexp(C{i,iAPX},'CF','once'))
            C{i,iAPX} =  str2double(regexp(C{i,iAPX},'\d+','once','match'));
        % remove line if char string
        elseif ischar(C{i,iAPX})
            tSN(i) = 0;
        % multiple UW ID 8514 lines no associated K0 value
        elseif C{i,iAPX} == 8514 && ismissing(C{i,iK0_ON})
            tSN(i) = 0;
        end
        
        % K2 FILTERS
        tmp = C(i,iK2); % get k2 data
        tg  = cellfun(@isnumeric,tmp);
        tmp(~tg) = {NaN};
        C(i,iK2) = tmp;
        
        % F(P) FILTERS
        tmp = C(i,iFP); % get k2 data
        tg  = cellfun(@isnumeric,tmp);
        tmp(~tg) = {NaN};
        C(i,iFP) = tmp;

        % K0 FILTERS:
        tmp = C(i,[iK0_ON, iK0_OFF, iK0_HCL]); % get k2 data
        tg  = cellfun(@isnumeric,tmp);
        tmp(~tg) = {NaN};
        C(i,[iK0_ON, iK0_OFF, iK0_HCL]) = tmp;
        
        if isnumeric(C{i,iK0_ON}) && C{i,iK0_ON} > 10
            tSN(i) = 0;
        end
        
        % VERSION COMMENT FILTER: Not a string? make it one & set to 'no info'
        if ~ischar(C{i,iVER})
            C{i,iVER} = 'no info';
        end
    catch
        continue
    end
end

ph_ver_data = C(tSN & tMSC, :);

% NOW DO A SECOND LOOP AND ASSIGN Version type
build_type = cell(size(ph_ver_data,1),2);
iDF = find(strcmp(ph_ver_hdr,'DF ID'));

for i = 1:size(ph_ver_data,1)
    DF      = ph_ver_data{i,iDF};
    ver_str = ph_ver_data{i,iVER};
    
    if isnan(DF)
        continue
    elseif ~isempty(regexp(ver_str,'AE GDF','once'))
        build_type{i,1} = 8.0;
        build_type{i,2} = 'AE GDF';  
    elseif ~isempty(regexp(ver_str,'GDF','once'))
        build_type{i,1} = 7.0;
        build_type{i,2} = 'MBARI GDF';
    elseif ischar(DF) || DF > 10000
        build_type{i,1} = 0.0;
        build_type{i,2} = 'SBE stem';
    elseif DF < 128
        build_type{i,1} = 1.0;
        build_type{i,2} = 'Ag wire, thinner ISFET cover';
    elseif DF < 164
        build_type{i,1} = 1.1;
        build_type{i,2} = 'Ag wire, thicker ISFET cover';
    elseif DF < 175
        build_type{i,1} = 1.2;
        build_type{i,2} = 'Pt wire, thicker ISFET cover';
    elseif DF < 198
        build_type{i,1} = 2.0;
        build_type{i,2} = 'Double ring, Pt wire, thicker ISFET cover';
    elseif DF < 251 % not certain about this divission
        build_type{i,1} = 3.0;
        build_type{i,2} = 'Smaller PT feed holes + V2.0';     
     elseif DF < 279 % not certain about this divission
        build_type{i,1} = 4.0;
        build_type{i,2} = 'Roll pin, Smaller PT feed holes + V2.0';
    elseif DF < 289 % not certain about this divission
        build_type{i,1} = 5.0;
        build_type{i,2} = 'V4 + Pellet Oring changed to 6.5 x 1.5 Buna 70A';
    elseif DF >= 289 % not certain about this divission
        build_type{i,1} = 6.0;
        build_type{i,2} = 'V5 + bus wire in CE crimp'; 
    end
end

%FINALIZE
ph_ver_hdr = [ph_ver_hdr(1:iVER-1), 'sensor version', ...
    'verion notes', ph_ver_hdr(iVER)];
ph_ver_data = [ph_ver_data(:,1:iVER-1), build_type, ph_ver_data(:,iVER)];

d.hdr = ph_ver_hdr;
d.data = ph_ver_data;

save([dirs.save,'MBARI_pHsensor_versions.mat'],'d')
clearvars -except d

