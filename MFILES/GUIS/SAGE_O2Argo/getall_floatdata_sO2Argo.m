function DATA = getall_floatdata_sO2Argo(thedatadir,floatID)
% ************************************************************************
% getall_floatdata_sO2Argo.m
% ************************************************************************
%
% Function to loop through all data from External (non-MBARI) Argo float,
% utilizing the BRtrj and ODV(from Mprof.nc) files.
% 
%
% INPUTS:  thedatadir  = input directory path structure
%          floatID = WMO ID
%
% OUTPUTS: 
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 10/25/17
% UPDATES:
% NOTES: Modified from MBARI's internal SageO2 function
% getall_floatdata_sO2 which works off of float msg files.
% ************************************************************************
% ************************************************************************

DATA.floatNAME = floatID;
% datadir = [dirs.Argo,floatID,'\'];

%__________________________________________________________________________
%__________________________________________________________________________
% CHECK THAT NECESSARY FILES EXIST
BRtraj = [thedatadir,floatID,'_BRtraj.nc'];
metaF = [thedatadir,floatID,'_meta.nc'];
Mprof = [thedatadir,floatID,'_Mprof.nc'];
Sprof = [thedatadir,floatID,'_Sprof.nc'];
odvfile = [thedatadir,'ODV',floatID,'.TXT'];
% if exist(BRtraj,'file')==0 || exist(metaF,'file')==0 || exist(Mprof,'file')==0 % if missing any of the 3
% if exist(metaF,'file')==0 || exist(Mprof,'file')==0 % if missing any of the 2;  if missing BRtraj is ok, just cal O2 to WOA
%     errorstring = [floatID,'_meta.nc, or ',floatID,'_Mprof.nc is missing from ',thedatadir,'.  Download these files from the GDAC before proceeding.'];
if exist(odvfile,'file')==0  % if missing ODV file;  if missing BRtraj is ok, just cal O2 to WOA; meta may not be needed anymore
    errorstring = ['ODV',floatID,'.TXT is missing from ',thedatadir,'.  Must be present before proceeding.  Use Sprof converter code to generate ODV file.'];
    h = errordlg(errorstring,'Missing File');
end
 
%__________________________________________________________________________
%__________________________________________________________________________
% PARSE CORIOLIS ODV FILE TO GET PROFILE DATA
D = get_FloatViz_data(odvfile);
% PARSE HEADER INFORMATION TO GET COLUMN INDICES
iST   = find(strcmp('Station',      D.hdr) == 1); % station
iSDN   = find(strcmp('SDN',      D.hdr) == 1); % datetime
iLAT   = find(strcmp('Lat [°N]',      D.hdr) == 1); % lat
iLON   = find(strcmp('Lon [°E]',      D.hdr) == 1); % lon
iS   = find(strcmp('Salinity[pss]',      D.hdr) == 1); % CTD S
iT  = find(strcmp('Temperature[°C]',   D.hdr) == 1); % CTD T
iP = find(strcmp('Pressure[dbar]', D.hdr) == 1); % CTD P
iO = find(strcmp('Oxygen[µmol/kg]', D.hdr) == 1); % Oxygen (umol/kg)
iOsat = find(strcmp('OxygenSat[%]', D.hdr) == 1); % OxygenSat (%)
% all quality flags columns are (iX)+1 

%__________________________________________________________________________
%__________________________________________________________________________
% GET FLOAT TRACK INFO
MPROF = D.data;
myC = unique(MPROF(:,iST));
DATA.Ncycles = length(myC);
DATA.cycles = myC;
%get track info
for jk = 1:length(myC)
    myLAT(jk,1) = nanmean(MPROF(MPROF(:,iST)==myC(jk),4));
    myLON(jk,1) = nanmean(MPROF(MPROF(:,iST)==myC(jk),3));
    mySDN(jk,1) = nanmean(MPROF(MPROF(:,iST)==myC(jk),1));
end
DATA.track = [mySDN myC myLON myLAT];

%__________________________________________________________________________
%__________________________________________________________________________
% ORGANIZE DATA AND PERFORM SOME QC SCREENING
% remove all bad quality flags, for oxygen keep ODV flag 4 (Argo flag 3,
% probably good but uncorrected)
% ODV quality flags:
%   0 = good
%   1 = missing or uninspected
%   4 = questionable (uncorrected)
%   8 = bad

%station and time (no existing qc flags)
SDN = MPROF(:,iSDN);
ST = MPROF(:,iST);
LAT = MPROF(:,iLAT);
LON = MPROF(:,iLON);

%pres, keep only 1s and 0s
P = MPROF(:,iP);
P(MPROF(:,iP+1)>1,:)=nan;

%sal, keep only 1s and 0s
S = MPROF(:,iS);
S(MPROF(:,iS+1)>1,:)=nan;

%temp, keep only 1s and 0s
T = MPROF(:,iT);
T(MPROF(:,iT+1)>1,:)=nan;

%oxygen, keep only qf<=4 (ODV 4 = questionable, Argo QF 3 equivalent)
O = MPROF(:,iO);
O(MPROF(:,iO+1)>4,:)=nan;
OSAT = MPROF(:,iOsat);
OSAT(MPROF(:,iOsat+1)>4,:)=nan;

%__________________________________________________________________________
%__________________________________________________________________________
% GENERATE VARIABLE TO HOLD PTS PROFILE DATA
%[sdn, cast, S, LAT, LON, P, T, S]; 
TMPpts = [SDN ST LAT LON P T S];
%double check for missing values
[tmp1, tmp2] = find(TMPpts < - 100000000);
TMPpts(tmp1,tmp2) = nan;
DATA.PTSdata = TMPpts;
DATA.PTSdata_hdr = {'sdn' 'cast' 'lat' 'lon' 'p' 't' 's'};

%__________________________________________________________________________
%__________________________________________________________________________
% GENERATE VARIABLE TO HOLD OXYGEN PROFILE DATA
% Contains many same variables as TMPpts, but different order, to maintain
% consistency with needed GUI variables downstream.
% Oxygen Phase data is empty, not needed unless you want to recalculate oxygen. Column structure is maintained.
%
%[sdn, cast, S, P, T, Phase,  O2]; % umol /L
TMP = [SDN ST S P T nan(size(SDN,1),1) O];
%double check for missing values
[tmp1, tmp2] = find(TMP < - 100000000);
TMP(tmp1,tmp2) = nan;
DATA.O2data{1} = TMP; 
% ----------------------------------------
% CALCULATE OXYGEN PERCENT SATURATION %
% % % O2sol = oxy_sol(TMP(:,5),TMP(:,3),0); % Benson & Krause coeffs used, see oxy_sol.m in MFILES/MISC
% % % profph2o = exp(52.57 -(6690.9./TK) - 4.681 .* log(TK));
% % % DATA.O2data{1}(:,11) = (TMP(:,end)./O2sol).*100; % percent saturation
%------------------------------------------
% Just grab O2 %sat from MPROF file.  Has already been calculated.
% Calculations are also written above (also see oxy_sol.m) for additional reference in case
% future code mods required due to change in MPROF structure.  
%double check for missing values
[xsat1, xsat2] = find(OSAT < - 100000000);  %remove fill values
OSAT(xsat1,xsat2) = nan;
whos OSAT
DATA.O2data{1}(:,11) = OSAT;
%DATA.O2data columns are now: 
%[sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat]; % umol /L & % sat
% O2sol, pO2, pH20 columns are empty for profile data, since we grabbed O2Sat from the Mprof
% file directly, however, these columns are kept as placeholders.

%__________________________________________________________________________
%__________________________________________________________________________
% EXTRACT FLOAT TYPE from Josh's ODV file (this way don't need meta file.)
fid = fopen(odvfile);
tline = ' ';
while ischar(tline)
    if regexp(tline,'type', 'once') % stop at header line
        break
    end
    tline = fgetl(fid);
end
exprs = '\:';
splitStr = regexp(tline,exprs,'split');
dft = strtrim(splitStr(2));
dft = cell2mat(dft);
xi = regexp(dft,'_');
if isempty(xi)
    DATA.floatTYPE = dft;
else
    DATA.floatTYPE = [dft(1:xi-1),' ',dft(xi+1:end)];
end

%__________________________________________________________________________
%__________________________________________________________________________
% PARSE TRAJ FILE TO GET AIR-O2 DATA
% USER SHOULD REVIEW THIS CODE!!
% SOME "IN-AIR" DATA FROM CERTAIN FLOATS/DACS IS UNRELIABLE.  YOU MAY WANT
% TO REVIEW THE PRESSURE THRESHOLD (PRESthres) DEFINED BELOW, OR IMPOSE
% ALTERNATIVE METHOD OF EXTRACTING IN-AIR DATA BASED ON KNOWLEDGE OF YOUR
% OWN FLOAT'S ACTIVITY/MISSION.
if exist(BRtraj,'file')==0 %NO BRtraj data!!  Cal to WOA & maintain empty O2air structures.
    DATA.O2air = cell(1,2);
    DATA.O2air{1} = []; %empty, no air cals
    DATA.O2air{2} = []; 
else
    mytarget = [thedatadir,floatID,'_BRtraj.nc'];   
    ds  = ncdataset(mytarget); 
    ga = ds.attributes;       % Global Attributes
    VARs = ds.variables;
    % check for needed variables.
    Index = strfind(VARs, 'PPOX_DOXY');
    Index_ppox = find(not(cellfun('isempty', Index)));
    Index = strfind(VARs, 'PPOX_DOXY_QC');
    Index_ppoxqc = find(not(cellfun('isempty', Index)));
    Index = strfind(VARs, 'MEASUREMENT_CODE');
    Index_mc = find(not(cellfun('isempty', Index)));
    Index = strfind(VARs, 'PRES');
    Index_pres = find(not(cellfun('isempty', Index)));
    Index = strfind(VARs, 'TEMP_DOXY');
    Index_tdoxy = find(not(cellfun('isempty', Index)));
    Index = strfind(VARs, 'TEMP_DOXY_QC');
    Index_tdoxyqc = find(not(cellfun('isempty', Index)));
    Index = strfind(VARs, 'CYCLE_NUMBER');
    Index_cycle = find(not(cellfun('isempty', Index)));
% NOT ALL DACS ARE PROPERLY POPULATING JULD, LATITUDE AND LONGITUDE WITHIN
% BRTRAJ FILES.  WE CAN USE LAT,LON,JULD WITHIN THE MPROF FILES FOR NOW,
% BUT KEEP THIS HERE JUST IN CASE FOR FUTURE.
%     Index = strfind(VARs, 'JULD');
%     Index_jday = find(not(cellfun('isempty', Index)));
%     Index = strfind(VARs, 'LATITUDE');
%     Index_lat = find(not(cellfun('isempty', Index)));
%     Index = strfind(VARs, 'LONGITUDE');
%     Index_lon = find(not(cellfun('isempty', Index)));
    if isempty(Index_ppox) || isempty(Index_mc) || isempty(Index_pres) || isempty(Index_tdoxy) || isempty(Index_cycle) ...
           || isempty(Index_ppoxqc) || isempty(Index_tdoxyqc) % NO Air-cal data
        f = msgbox({'BRtraj file does not include necessary variables for aircal.','Calibrate to WOA2013.'});
        DATA.O2air = cell(1,2);
        DATA.O2air{1} = []; %empty, no air cals
        DATA.O2air{2} = []; 
    else % Air-cal data; populate structures:
        % Keep only QC flags <4 for ppox_doxy 
        surfO = double(ds.data('PPOX_DOXY')); %air-cal samples;  Units of millibar (hPa equivalent)
        surfO_qc = ds.data('PPOX_DOXY_QC'); %you really only want good data, however, I don't think DACs aren't performing QC on air-samples as of now, so QCflags=0 (Argo) 
        tmp_qcO = ones(size(surfO_qc))*NaN; % predim QC matrix
        for j = 1:size(surfO_qc,1)
            qcvals = sscanf(surfO_qc(j,:),'%1f')';
            tmp_qcO(j,1:size(qcvals,2)) = qcvals;
        end
        surfOQC = tmp_qcO;
        clear qcvals tmp_qcO
        oi = find(surfOQC>3);
        surfO(oi)=nan;
        
        Tdoxy = double(ds.data('TEMP_DOXY')); %This variable is needed for aircal!  But not all DACs are storing it...
        Tdoxy_qc = ds.data('TEMP_DOXY_QC'); 
        tmp_Td = ones(size(Tdoxy_qc))*NaN; % predim QC matrix
        for j = 1:size(Tdoxy_qc,1)
            qcvals = sscanf(Tdoxy_qc(j,:),'%1f')';
            tmp_Td(j,1:size(qcvals,2)) = qcvals;
        end
        TdoxyQC = tmp_Td;
        clear qcvals tmp_Td
        ti = find(TdoxyQC>3);
        Tdoxy(ti)=nan;
       
%%%%% NOTE:: USER TO CHANGE MEASUREMENT-CODE SPECIFICATIONS BASED ON KNOWLEDGE OF THEIR OWN FLOAT PLATFORM'S ACQUISITION OF IN-AIR SAMPLES %%%%%     
% SEE THE MEASUREMENT CODE SPECIFICATIONS FOR IN-AIR DAT IN SECTION 2.5.2 OF THE ARGO QUALITY CONTROL MANUAL FOR DISSOLVED OXYGEN   
% Note that prior to ADMT18, the measurement code specifications for in-air oxygen samples were MC = 1090 (1100-10) or 1099 (1100-1) (sequence or single obs with telemetry, respectively)
        mscode = ds.data('MEASUREMENT_CODE'); %
        pres = ds.data('PRES'); 
        cycnum = ds.data('CYCLE_NUMBER');
%         jday = ds.data('JULD');
%         Lat = ds.data('LATITUDE');
%         Lon = ds.data('LONGITUDE');
%         cycnumi = ds.data('CYCLE_NUMBER_INDEX'); %keep, in case of future need
%         doxy = ds.data('DOXY'); %keep, in case of future need
        % XX = [double(cycnum) jday Lat Lon double(pres) double(Tdoxy) double(rphase) double(tphase) double(surfO)];
        XX = [double(cycnum) double(pres) double(Tdoxy) double(surfO)];
        X = XX(mscode==711,:); % MC = 711; X+11 for in-air samples, part of surface sequence(X = 700 (TST))
        T = X(:,3); %temp doxy
        S = zeros(length(T),1);
%         P = zeros(length(T),1);
        P = X(:,2);
        % P2 = X(:,2);
        % P2(isnan(P2))=0;
        O2 = X(:,end);

%__________________________________________________________________________
%__________________________________________________________________________
        % GET OTHER O2 PARAMETERS, THIS CODE SNIPPET CAME FROM "calc_O2_4ARGO.m"
        % The BRtraj file already gives us pO2 in millibars, so all we need
        % is pH20.
        %
        % -------------------------------------------------------------------------
        % % O2 params--------------------------------------------------------------
        % % Benson & Krause (cm^3/dm^3) *** USE THESE COEFFICIENTS!!!! ***
        % pA = [3.88767 -0.256847 4.94457 4.05010 3.22014 2.00907]; % Temperature Coeff
        % pB = [-8.17083e-3 -1.03410e-2 -7.37614e-3 -6.24523e-3]; % Salinity coff
        % Co = -4.88682e-7;
        % O2_volume =  22.3916;
        % % ************************************************************************
        % % %Aanderaa default coeff (cm^3/dm^3) (Combined fit Garcia & Gordon 1992)
        % % pA = [1.71069 0.978188 4.80299 3.99063 3.22400 2.00856]; % Temperature Coeff
        % % pB = [-4.29155e-3 -6.90358e-3 -6.93498e-3 -6.24097e-3]; % Salinity coff
        % % Co = -3.11680e-7;
        % % O2_volume =  22.414; % Volume ideal gas
        % % ************************************************************************
        TK     = T + 273.15;
        % if TEMP_DOXY is empty or all bad, use ctd temperature at surface
        DATA.PTSdata_hdr = {'sdn' 'cast' 'lat' 'lon' 'p' 't' 's'};
        
        
        % Ts     = log((298.15-T)./ TK);
        % part1  = polyval(pB,Ts);
        % S_corr = S.*part1 + S.^2 * Co;
        % L      = polyval(pA,Ts) + S_corr;
        % O2sol  = (1000/O2_volume) * exp(L); % Oxygen solubility real gas, mmol/m^3 =uM/L
        % 
        % % Calculate vapor pressure h20
        ph2o = exp(52.57 -(6690.9./TK) - 4.681 .* log(TK));
        % pO2 = (O2 ./ O2sol) .* ((1013.25- pH2O) * 0.20946);
        % -------------------------------------------------------------------------
        % -------------------------------------------------------------------------
        
%__________________________________________________________________________
%__________________________________________________________________________
        % STORE PO2 AND PH20 IN ORDER TO CALCULATE GAIN
        % Ok, now have pO2 and pH20 for all air cal samples.  For cycles with multiple
        % telemetry air samples, do some averaging.
        mycycs = X(:,1);
        cycn = unique(mycycs); %includes "0" cycle?
        [~,IA,IB] = intersect(myC,cycn); %get lat lon sdn associated with each
        % aircal sample.  Use Mprof for this (not all DACs are storing Lat,
        % Lon, Juld in the BRtraj files!!!)
        LAT = myLAT(IA);
        LON = myLON(IA);
        SDN = mySDN(IA);
        cycn = cycn(IB); % eliminate "0" cycle
        for i = 1:length(cycn)
            myO2 = O2(mycycs==cycn(i));
            myph2o = ph2o(mycycs==cycn(i));
            myp = P(mycycs==cycn(i));
            %find only air-cal sample!!!  Can use pressure to home in on
            %"in-air" samples, although perhaps better method?  This is
            %conservative, some floats have optodes mounted on a stick, so
            %checking "CONFIG_OptodeVerticalPressureOffset_dbar" in meta
            %file could give a more precise measure of optode location -->
            %more datapoints.
            PRESthresh = 0.10;
            bp = find(myp<PRESthresh); 
            myO2 = myO2(bp);
            myph2o = myph2o(bp);
            PO2(i) = nanmean(myO2);
            PH2O(i) = nanmean(myph2o);
            PO2st(i) = nanstd(myO2);
            PH2Ost(i) = nanstd(myph2o);
%             Is dedicated air-cal sequence.  For
%             certain PROVOR floats, in-air samples have a sequence of 15
%             samples.  First 5 are "in-water", last 10 are "in-air".  But
%             this does not apply to all floats so do not keep in main
%            primary software version.
%             if length(myO2) == 15 
%                 PO2(i) = nanmean(myO2(6:end));
%                 PH2O(i) = nanmean(myph2o(6:end));
%                 PO2st(i) = nanstd(myO2(6:end));
%                 PH2Ost(i) = nanstd(myph2o(6:end));
%             else
%                 PO2(i) = nanmean(myO2);
%                 PH2O(i) = nanmean(myph2o);
%                 PO2st(i) = nanstd(myO2);
%                 PH2Ost(i) = nanstd(myph2o);
%             end
        end

%__________________________________________________________________________
%__________________________________________________________________________
        % FOR DATA.O2air FIELD, MAINTAIN SAME STRUCTURE AS DONE IN SAGE_O2 FOR MBARI, IN
        % ORDER TO UTILIZE SIMILAR SUPPORTING FUNCTIONS
        %[sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat]; % umol /L & % sat
        DATA.O2air = cell(1,2);
        DATA.O2air{1}{1} = [SDN cycn nan(length(cycn),6) PO2' PH2O' nan(length(cycn),1)]; %average for each cycle
        DATA.O2air{1}{3} = []; %maintain compatability with pre-existing code-base
        DATA.O2air{2}{1} = [SDN cycn nan(length(cycn),6) PO2st' PH2Ost' nan(length(cycn),1)]; %std for each cycle
        DATA.O2air{2}{3} = []; %maintain compatability with pre-existing code-base
    end
end

% END
