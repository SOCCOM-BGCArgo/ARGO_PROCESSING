function d = get_shipboard_data(file_path)
% PURPOSE: 
%   This function parses shipboard data to get calibration cast data
%   for SOCCOM float deployments
%
% USAGE:
%   data = get_shipboard_data(path\data_file)
%
% INPUTS:
%   path\data_file = Path to file

%
% OUTPUTS: a structure
%             .hdr    = cell array of column headers
%             .data   = a matrix of the calibration data data
%
% EXAMPLES:
%    get_shipboard_data(['C:\Users\jplant\Documents\MATLAB\ARGO\DATA\', ...
%                        'Shipboard\320620140320_hy1.csv'])
% CHANGE HISTORY
%   12/6/2016 Added code to estimate alkalinity if not measured (LIARv2)
%       and to not break if no pH data exists. -jp
%   03/02/2017 Added code to fix bad HHMM data format in file & to not use
%       DIC or ALK values < 0 when estimating pH. -jp
%   08/09/2017 Added code to look for PH on seawater scale (ie P18). If
%       found (and ph Total not found) change pHscalein input in CO2SYSSOCCOM
%       TO 2 (SEAWATER SCALE). Output used changed col 18 to col 37.
%   08/15,2017 Added code to calculate depth if NaN's in depth col (HOT
%       ALOHA data)
%   01/09/2018 Added code to capture chlorophyll data fron files. This
%       incudes two new parameters & their QC fields: 'CHLORA' &
%       'TOT_CHL_A'. The latter is HPLC derived while the former is
%       presumabley fluorometrically derived. -jp
%   08/0/2018 Added code to estimate missing Si & PO4 if NO3 exists using
%       Redfield approximation. Si & PO4 are used for inputs into CO2SYS to
%       estimate pH in situ. The nutrient estimates are not returned to the
%       data set for output. This is added around line 319. SR1B cruise
%       (floats 9652 9655 9657 9662) has no Silicate data.  - jp
%   08/28/2018 Fixed QF screening.  Now only allows "2 = No problems noted".
%   03/05/2019 Added CTDOXY to wanted vars list

% TESTING             
% file_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%              'SHIPBOARD\096U20160426.exc.csv'];

% file_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%              'SHIPBOARD\33RO20161119.exc.csv'];

%file_path = ' C:\temp\MR16090320170303220727_forSOCCOM.csv';

%file_path = ' C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\SHIPBOARD\096U20160426.exc.csv';
%file_path = ' C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\SHIPBOARD\096U20160314.exc.csv';

%file_path = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\SHIPBOARD\74JC20151217_hy1.csv';

%file_path = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\SHIPBOARD\06AQ20141202_hy1.csv';

% ************************************************************************
% LIST OF DESIRED VARIBLES AND FORMAT STRING
% ADD / REMOVE WHAT YOU WANT
wanted_vars ={'SECT_ID'          '%s'; ...
              'STNNBR'           '%f'; ...
              'CASTNO'           '%f'; ...
              'DATE'             '%s'; ...
              'TIME'             '%s'; ...
              'LATITUDE'         '%f';...
              'LONGITUDE'        '%f'; ...
              'DEPTH'            '%f'; ...
              'CTDPRS'           '%f'; ...
              'CTDTMP'           '%f'; ...
              'CTDSAL'           '%f'; ...    %'SALNITY' could use bottle or ctd
              'CTDSAL_FLAG_W'    '%f'; ...
              'CTDOXY'           '%f'; ...    %'SALNITY' could use bottle or ctd
              'CTDOXY_FLAG_W'    '%f'; ... 
              'OXYGEN'           '%f'; ...
              'OXYGEN_FLAG_W'    '%f'; ...
              'SILCAT'           '%f'; ...
              'SILCAT_FLAG_W'    '%f'; ...
              'NITRAT'           '%f'; ...
              'NITRAT_FLAG_W'    '%f'; ...
              'PHSPHT'           '%f'; ...   
              'PHSPHT_FLAG_W'    '%f'; ... 
              'TCARBN'           '%f'; ...
              'TCARBN_FLAG_W'    '%f'; ...
              'ALKALI'           '%f'; ...
              'ALKALI_FLAG_W'    '%f'; ...
              'PH_TOT'           '%f'; ...
              'PH_TOT_FLAG_W'    '%f'; ...
              'PH_SWS'           '%f'; ... % NOAA PMEL OFTEN REPORT DATA IN SW SCALE 
              'PH_SWS_FLAG_W'    '%f'; ...           
              'PH_TMP'           '%f';...
              'CHLORA'           '%f';...
              'CHLORA_FLAG_W'    '%f';...
              'TOT_CHL_A'        '%f';
              'TOT_CHL_A_FLAG_W' '%f'};
              
% *************************************************************************
% OPEN FILE, STEP DOWN TO 1st HEADER LINE, BUILD FORMAT STRING
% *************************************************************************

fid = fopen(file_path);
if fid == -1
    d= [];
    return
end
tline = ' ';
while ischar(tline)
    if regexp(tline,'^EXPOCODE,\w+', 'once') % stop at 1st header line
        break
    end
    tline = fgetl(fid);
end

if ~ischar(tline)
    disp('No header line found')
    return
end

% CHECK HEADER VARIANTS - RENAME IF NEEDED
tline = regexprep(tline,'SECT,','SECT_ID,');

tmp_hdr    = regexp(tline,',','split'); % CELL ARRAY OF HEADER VARIABLES 
tline = fgetl(fid); % STEP TO SECOND HEADER LINE (UNITS)
tmp_units  = regexp(tline,',','split'); % CELL ARRAY OF UNITS 
hdr_cols   = size(tmp_hdr,2);
format_str = '';
hdr        = {};


for i = 1:hdr_cols
    tf = strcmp(tmp_hdr{i}, wanted_vars(:,1));
    if sum(tf) > 0 % A data match!
        format_str = [format_str,wanted_vars{tf,2}];
        hdr   = [hdr,tmp_hdr{i}];        
    else
        format_str = [format_str,'%*s'];
    end
end

% *************************************************************************
% PARSE THE DATA
% *************************************************************************
tline = fgetl(fid); % STEP PAST SECOND HEADER LINE

%d     = textscan(fid,format_str,'Delimiter',',','CollectOutput',1);
d     = textscan(fid,format_str,'Delimiter',',','CommentStyle', 'END');
fclose(fid);

clear fid format_str hdr_rows tmp_hdr tf i

% *************************************************************************
% CONDENSE THE DATA
% *************************************************************************
cruise_ID = d{1,1}{1,1};

iDAY  = strcmp(hdr,'DATE');
iTIME = strcmp(hdr,'TIME');
iSECT = strcmp(hdr,'SECT_ID');

if regexp(d{1,iDAY}{1,1},'/')
    DAY = datenum(d{1,iDAY},'mm/dd/yyyy');
else
    DAY = datenum(d{1,iDAY},'yyyymmdd');
end

% THE HOUR FORMAT CAN VARY IN SOME OF THESE FILES. THE LEADING ZEROS CAN BE
% MISSING  example: '660' instead of '0660'
% USE REGULAR EXPRESSIONS TO FILL IN MISSING 0's
hr_tmp = d{1,iTIME};
% Look only 3,2, or 1 #'s, capture as token and replace with appropriate 
% number of 0's added to token
hr_tmp = regexprep(hr_tmp,'\:',''  ,'once'); % remove ":" in time if it exists
hr_tmp = regexprep(hr_tmp,'(^\d{3}$)','0$1'  ,'once'); %only 3 #'s
hr_tmp = regexprep(hr_tmp,'(^\d{2}$)','00$1' ,'once'); %only 2 #'s
hr_tmp = regexprep(hr_tmp,'(^\d{1}$)','000$1','once'); %only 1 #'s

HR  = datenum(hr_tmp,'HHMM');
sdn = DAY + HR - fix(HR);
hdr(iTIME) = []; % remove "Time" from header
hdr(iSECT) = []; % remove "SECT_ID" from header
clear hr_tmp



data = [d{1,2}, d{1,3}, sdn]; % BUILD MATRIX
for i = 6:size(d,2)
    data = [data,d{1,i}];
end
clear d

% SET MISSING VALUES = NaN
data(data == -999) = NaN;

% ************************************************************************
% CHECK QUALITY FLAGS & SET BAD DATA TO NaN 10/30/2017 JP

ind  = regexp(hdr,'_FLAG_W');
tfQC = ~cellfun(@isempty,ind); % quality flag cols
for i = 1:length(tfQC)
    if tfQC(i) == 1 %QF flag
        QCtmp = find(data(:,i)~=2);
        if ~isempty(QCtmp)
            data(QCtmp,i-1) = nan;
        end
    end
end

% tfD  = logical(tfQC*0); % predim for data cols
% tfD(1:end-1) = tfQC(2:end); % shifted by -1 for data index assoc with QC flag
% % QCtmp = data(:,tfQC) == 4; % Flag bad data for removal, 4 = "did not trip correctly"
% QCtmp = data(:,tfQC) ~= 2; % Flag bad data for removal.  2 = "no problems noted"
% Dtmp = data(:,tfD );
% Dtmp(QCtmp) = NaN; % bad values to nan
% data(:,tfD) = Dtmp; % put back into data matrix
data(:,tfQC) = []; % remove QC cols -don't need any more
hdr(tfQC)    = []; % remove QC cols -don't need any more
clear ind tfQC tfD QCtmp Dtmp

% *************************************************************************
% NOW ESTIMATE PH AT IN SITU TEMPERATURE
% *************************************************************************
% GET SHIPBOARD DATA INDICES INDICES
iP    = find(strcmp('CTDPRS', hdr) == 1); % dbar
iT    = find(strcmp('CTDTMP', hdr) == 1);
iS    = find(strcmp('CTDSAL',hdr)  == 1);
iDIC  = find(strcmp('TCARBN',hdr)  == 1);
iALK  = find(strcmp('ALKALI',hdr)  == 1);
iPH   = find(strcmp('PH_TOT',hdr)  == 1);
iPHSW = find(strcmp('PH_SWS',hdr)  == 1);
iPHT  = find(strcmp('PH_TMP',hdr)  == 1);
iSI   = find(strcmp('SILCAT',hdr)  == 1);
iPO4  = find(strcmp('PHSPHT',hdr)  == 1);
iNO3  = find(strcmp('NITRAT',hdr)  == 1);


iLAT  = find(strcmp('LATITUDE',hdr)  == 1); % THESE ARE NEEDED FOR LIAR
iLON  = find(strcmp('LONGITUDE',hdr) == 1); % TO ESTIMATE ALKALINITY IF
iO    = find(strcmp('OXYGEN',hdr)    == 1); % IT IS NOT PRESENT IN THE 
iSDN  = find(strcmp('DATE',hdr)      == 1); % DATASET
iZ    = find(strcmp('DEPTH',hdr)     == 1);

% CALCULATE DEPTH IF NAN's IN DEPTH COL
nan_Z = isnan(data(:,iZ));
data(nan_Z, iZ) = sw_dpth(data(nan_Z,iP),data(nan_Z,iLAT)); 


if isempty(iALK) % NO ALKALINITY MEASUREMENT ESTIMATE WITH LIAR
    LIAR_pos = [data(:,iLAT), data(:,iLON), data(:,iZ)];
    ptemp    = theta(data(:,iP), data(:,iT), data(:,iS), 0);
    Equations    = 7; % S, Theta, AOU
    
    
    %       µmol/kg unless MolalityTF is set to 0
    % MeasIDVec Parameter Key:                  Equation Key:
    % 1. Salinity                               5.  S, Theta, N, AOU
    % 2. Potential temperature                  6.  S, Theta, N
    % 3. Nitrate                                7.  S, Theta, AOU
    % 4. AOU                                    8.  S, Theta
    % 5. Silicate                               13. S, N, AOU
    % 6. O2                                     14. S, N
    % 7. Temperature                            15. S, AOU
    %                                           16. S

    
%     MeasIDVec    = [1 2 6]; % PSAL, Pot_TEMP, DOXY_ADJ,
%     Measurements = [data(:,iS), ptemp, data(:,iO)];
    
    MeasIDVec    = [1 6 7]; % PSAL, DOXY_ADJ, TEMP, % updated 8/11/17 jp
    Measurements = [data(:,iS), data(:,iO), data(:,iT) ];
    
    % using default MeasUncerts
    disp('Creating Alkalinity data column using LIAR_v2........')
%     [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos, Measurements, MeasIDVec, ...
%         Equations, [], 1); % update 04/25/17 
    
    [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos, ...
            Measurements, MeasIDVec,'Equations', Equations); % update 08/11/17
    data = [data, Alk_Est];
    hdr  = [hdr, 'ALKALI'];
    iALK = find(strcmp('ALKALI',hdr)  == 1); % % find index again
    clear LIAR_pos ptemp MeasIDVec Equations Alk_Est Uncert_Est
end

% Alkalinity column exist but contains NaN - added 05/5/17 
if ~isempty(iALK) && any(isnan(data(:,iALK)))
    disp('Adding missing Alkalinty values using LIAR_v2........')
    tnan_alk = isnan(data(:,iALK));
    
    LIAR_pos = [data(tnan_alk,iLAT), data(tnan_alk,iLON), ...
                data(tnan_alk,iZ)];
    ptemp    = theta(data(tnan_alk,iP), data(tnan_alk,iT), ...
               data(tnan_alk,iS), 0);
    
%     MeasIDVec    = [1 2 6]; % PSAL, Pot_TEMP, DOXY_ADJ,
%     Measurements = [data(tnan_alk,iS), ptemp, data(tnan_alk,iO)];
%     Equations    = 7; % S, Theta, AOU
    
    MeasIDVec    = [1 6 7]; % PSAL, DOXY_ADJ, TEMP, % updated 8/11/17 jp
    Measurements = [data(tnan_alk,iS), data(tnan_alk,iO), data(tnan_alk,iT)];
    Equations    = 7; % S, Theta, AOU

%     [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos, Measurements, ...
%              MeasIDVec, Equations, [], 1); 
    [Alk_Est, Uncert_Est, MinUncert_Equ] = LIAR(LIAR_pos, ...
            Measurements, MeasIDVec,'Equations', Equations); % update 08/11/17
    data(tnan_alk,iALK) = Alk_Est;
    clear LIAR_pos ptemp MeasIDVec Equations Alk_Est Uncert_Est tnan_alk
end

% SET UP CO2 SYSTEM PARAMETERS TO GET IN SITU pH ON THE TOTAL SCALE
ph_flag1 = 1;
ph_flag2 = 1;
if isempty(iPH) && ~isempty(iPHSW) % CHECK FOR PH ON SW SCALE
    iPH = iPHSW;
    pHSCALEIN     = 2;  % SEAWATER SCALE
else
    pHSCALEIN     = 1;  % TOTAL SCALE
end

K1K2CONSTANTS = 10; % Lueker et al, 2000
%KSO4CONSTANTS = 1;  % KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED) 
KSO4CONSTANTS = 3;  % KSO4 of Dickson & TB of Lee 2010 (USE THIS ONE !!!)

SAL     = data(:,iS);
TEMPOUT = data(:,iT);

PRESOUT = data(:,iP);

if ~isempty(iSI) && ~isempty(iPO4) %prelim datafiles often dont contain these
    SI      = data(:,iSI);
    PO4     = data(:,iPO4);

    %testing
    % t_nan = isnan(data(:,iSI));
    % data(:,iSI) = 20;
    % SI = data(:,iSI);

    % ***********************
    % SOME TIMES NEED TO ESTIMATE SILICATE OR PHOSPHATE BOTTLE DATA
    % USE REDFILED RATIO TO APROXIMATE IF GOOD NITRATE EXISTS, OTHERWISE SET = 0
    % JP FIX 08/08/2018 (ie SR1B cruise & floats 9652 9655 9657 9662) 
    nan_SI  = isnan(SI)  & ~isnan(data(:,iNO3));
    nan_PO4 = isnan(PO4) & ~isnan(data(:,iNO3));
    if sum(nan_SI) > 0 || sum(nan_PO4) > 0
        disp('NO3 data exists but some complimentary Si or PO4 data is missing')
        disp('Estimating missing data with Redfield ratio aproximation')

        SI(nan_SI)   = data(nan_SI, iNO3) * 2.5;  % APROX WITH REDFIELD
        PO4(nan_PO4) = data(nan_PO4, iPO4) / 16;
    end
    clear nan_SI nan_PO4
end
% ***********************

if ~isempty(iALK) && ~isempty(iPH) % check first: Alkalinity & pH exist
    PAR1 = data(:,iALK);
    PAR1(PAR1 < 0) = NaN;
    PAR1TYPE   = 1;
    PAR2 = data(:,iPH);
    PAR2(PAR2 < 0) = NaN;
    PAR2TYPE   = 3;
    TEMPIN = data(:,iPHT);
    PRESIN  = 0;
    
    [OUT1] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
            TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
            pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);
        
    %pH_tot_iAiPH = OUT1(:,37);
    % Col 37 in OUT is pH output total scale regardless of input scale
else 
    ph_flag1 = 0; % O if no DIC, ALk or pH
end
     
if ~isempty(iALK) && ~isempty(iDIC) %check 2nd: Alk & DIC
    PAR1 = data(:,iALK);
    PAR1(PAR1 < 0) = NaN;
    PAR1TYPE   = 1;
    PAR2 = data(:,iDIC);
    PAR2(PAR2 < 0) = NaN;
    PAR2TYPE   = 2;
    TEMPIN = data(:,iT);
    PRESIN = data(:,iP);
    
    [OUT2] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
        TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
        pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);
    
    %pH_tot_iAiDIC = OUT2(:,37);
    % Col 37 in OUT is pH output total scale regardless of input scale
else
    ph_flag2 = 0; % O if no DIC, ALk or pH
end

if ph_flag1 ==1
    % ADD TO DATA AND HEADER
    %data = [data, OUT(:,18)]; % Col 18 in OUT is pH output
    % Col 37 in OUT is pH output total scale regardless of input scale
    data = [data, OUT1(:,37)]; 
    hdr  = [hdr, 'PH_TOT_INSITU'];
end
if ph_flag2 ==1
    % ADD TO DATA AND HEADER
    %data = [data, OUT(:,18)]; % Col 18 in OUT is pH output
    % Col 37 in OUT is pH output total scale regardless of input scale
    data = [data, OUT2(:,37)]; 
    hdr  = [hdr, 'PH_TOT_INSITU_ALKDIC'];
end


% % MAKE SOME PLOTS FOR TESTING
% iSDN = find(strcmp('DATE',  hdr) == 1);
% iZ   = find(strcmp('DEPTH', hdr) == 1);
% iP   = find(strcmp('CTDPRS',hdr) == 1);
% iT   = find(strcmp('CTDTMP',hdr) == 1);
% iS   = find(strcmp('CTDSAL',hdr) == 1);
% iO   = find(strcmp('OXYGEN',hdr) == 1);
% iN   = find(strcmp('NITRAT',hdr) == 1);
% iALK  = find(strcmp('TCARBN',hdr) == 1); % pH & ALK
% iDIC  = find(strcmp('ALKALI',hdr) == 1); % pH & ALK
% iPH0  = find(strcmp('PH_TOT',hdr) == 1); % pH & ALK
% iPH1  = find(strcmp('PH_TOT_INSITU',hdr) == 1); % pH & ALK
% iPH2  = find(strcmp('PH_TOT_INSITU_ALKDIC',hdr) == 1); % ALK & DIC
% 
% tgood = ~isnan(data(:,iPH1)) & ~isnan(data(:,iPH2));
% tbad  = data(:,1) == 2 & data(:,2) == 1;
% 
% F1 = figure(1);
% F1.Position = [273 141.8000 1.0056e+03 716];
% y_lim = [0 2000];
% 
% subplot(2,3,1)
% plot(data(tgood,iALK),data(tgood,iP), 'bo', 'MarkerSize', 3)
% set(gca,'Ydir', 'Reverse')
% ylim(y_lim);
% xlabel(hdr{iALK});
% ylabel(hdr{iP});
% hold on
% plot(data(tbad,iALK),data(tbad,iP), 'r*', 'MarkerSize', 3)
% hold off
% set(gca,'FontSize', 14);
% 
% 
% subplot(2,3,2)
% plot(data(tgood,iDIC),data(tgood,iP), 'bo', 'MarkerSize', 3)
% set(gca,'Ydir', 'Reverse')
% ylim(y_lim);
% xlabel(hdr{iDIC});
% %ylabel(hdr{iP});
% hold on
% plot(data(tbad,iDIC),data(tbad,iP), 'r*', 'MarkerSize', 3)
% hold off
% set(gca,'FontSize', 14);
% 
% subplot(2,3,3)
% plot(data(tgood,iPH0),data(tgood,iP), 'bo', 'MarkerSize', 3)
% set(gca,'Ydir', 'Reverse')
% ylim(y_lim);
% xlabel([hdr{iPH0}, ' 25C'],'Interpreter', 'none');
% %ylabel(hdr{iP});
% hold on
% plot(data(tbad,iPH0),data(tbad,iP), 'r*', 'MarkerSize', 3)
% hold off
% set(gca,'FontSize', 14);
% 
% subplot(2,3,4)
% plot(data(tgood,iPH1),data(tgood,iP), 'bo', 'MarkerSize', 3)
% set(gca,'Ydir', 'Reverse')
% ylim(y_lim);
% xlabel(hdr{iPH1},'Interpreter', 'none');
% ylabel(hdr{iP});
% hold on
% plot(data(tbad,iPH1),data(tbad,iP), 'r*', 'MarkerSize', 3)
% hold off
% set(gca,'FontSize', 14);
% 
% subplot(2,3,5)
% plot(data(tgood,iPH2),data(tgood,iP), 'bo', 'MarkerSize', 3)
% set(gca,'Ydir', 'Reverse')
% ylim(y_lim);
% xlabel(hdr{iPH2},'Interpreter', 'none');
% %ylabel(hdr{iP});
% hold on
% plot(data(tbad,iPH2),data(tbad,iP), 'r*', 'MarkerSize', 3)
% hold off
% set(gca,'FontSize', 14);
% 
% subplot(2,3,6)
% plot(data(tgood,iPH1)-data(tgood,iPH2), data(tgood,iP), 'bo', 'MarkerSize', 3)
% set(gca,'Ydir', 'Reverse')
% ylim(y_lim);
% xlabel({[hdr{iPH1}, '  -  '],hdr{iPH2}},'Interpreter', 'none');
% %ylabel(hdr{iP});
% hold on
% plot(data(tbad,iPH1)-data(tbad,iPH2),data(tbad,iP), 'r*', 'MarkerSize', 3)
% hold off
% set(gca,'FontSize', 14);

% *************************************************************************
% ASSIGN TO STRUCTURE AND CLEAN UP
% *************************************************************************
d.hdr    = hdr;
d.data   = data;
d.cruise = cruise_ID;
%d.units  = units;

%clearvars -except d
