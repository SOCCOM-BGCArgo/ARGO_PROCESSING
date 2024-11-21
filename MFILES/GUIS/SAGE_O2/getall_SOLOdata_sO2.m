function DATA = getall_SOLOdata_sO2(dirs,MBARI_ID)
% ************************************************************************
% getall_SOLOdata_sO2.m
% ************************************************************************
%
% Function to loop through all BSOLO data files for a particular float in
% specified directory.
%
%
% INPUTS:  pn      = input directory path structure
%          floatID = MBARI float ID (ie ss0001)
%
% OUTPUTS:
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 03/24/2022
% UPDATES:

%LG and JP 9/13/2024: SPECIFIC CONDITION FOR FLOAT 3902533 (ss4018):
%   Cycle 12 has fill values of -999 in the aircal data leading to imaginary numbers when calculating air ref data.
%   Added a test to check for missing values in the air cal data (-999), and then a loop to remove this missing data
%   from the dataset while calculating pO2 references. Fixes the issue for ss4018, hopefully doesn't pop up for more.
%
% NOTES:
%
% *************************************************************************
% TESTING
%
% dirs = SetupDirs_sO2;
% MBARI_ID = 'ss4018';
% *************************************************************************
%
% ************************************************************************
% ************************************************************************
DATA.floatTYPE = 'SOLO';
% GET FLOAT CALS___________________________________________________________
float_cal_path = [dirs.cal,'cal',MBARI_ID,'.mat'];
if exist(float_cal_path,'file')
    disp(['Loading existing calibration file: ',float_cal_path])
    load(float_cal_path);
else
    cal = get_float_cals(MBARI_ID, dirs);
end
cal_info = cal.info;

if isempty(cal)
    disp(['NO CALIBRATION DATA FOR ',MBARI_ID])
    return
end

%DEFINE DATA ROOT LOCATIONS & COUNT N FILES_________________________________
Dpath = [dirs.mat,cal_info.WMO_ID,'\'];
F = ls([Dpath,cal.info.WMO_ID,'*.mat']);
FVf = [dirs.FV,cal.info.WMO_ID,'.TXT'];
%remove the 000.mat file, if exists (should always be first listed, so...):
if contains(F(1,:),'000') %~isempty(SF)
    F(1,:)=[]; %remove
end
DATA.matfiles = F;
DATA.Nprof = size(F,1);

%EXTRACT & ORGANIZE DATA FOR A GIVEN FLOAT_________________________________
FVd = get_FloatViz_data(FVf);
iSDN = find(strcmp('SDN',FVd.hdr)              == 1);
iSTN = find(strcmp('Station',FVd.hdr)              == 1);
iLAT = find(strcmp('Lat [°N]',FVd.hdr)              == 1);
iLON = find(strcmp('Lon [°E]',FVd.hdr)              == 1);
iP = find(strcmp('Pressure[dbar]',FVd.hdr)              == 1);
iT = find(strcmp('Temperature[°C]',FVd.hdr)              == 1);
iS = find(strcmp('Salinity[pss]',FVd.hdr)              == 1);
iO = find(strcmp('Oxygen[µmol/kg]',FVd.hdr)              == 1);
iOsat = find(strcmp('OxygenSat[%]',FVd.hdr)              == 1);
% Lots of fill values in the O2 data from the axes merge.  Replace with NaN
% here.
FVd.data(FVd.data(:,iO)<-1000,iO)=NaN;
FVd.data(FVd.data(:,iOsat)<-1000,iOsat)=NaN;
if ~isempty(F)
    TRACK = zeros(size(F,1),4);
    DATA.PTSdata_hdr = {'sdn'  'cast'  'lat'  'lon'  'p'  't'  's'};
    DATA.PTSdata = [FVd.data(:,iSDN) FVd.data(:,iSTN) FVd.data(:,iLAT) FVd.data(:,iLON) FVd.data(:,iP) FVd.data(:,iT) FVd.data(:,iS)];
    % organize O2 data (profile phase and computed o2).
    %--------------------
    % (1) PROFILE: Fill profile data first.  Don't need to carry
    % through phase data for profile O2...
    DATA.O2phase{1} = NaN; %don't need to carry phase data for profile category.  Not used explicitly in sageO2.
    DATA.O2data{1} = nan(size(FVd.data,1),11); %to be in line with getall_floatdata_sO2, but only cols 1-5,7,11 are needed (?) [%profile data: sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat]
    DATA.O2data{1}(:,1) = FVd.data(:,iSDN);
    DATA.O2data{1}(:,2) = FVd.data(:,iSTN);
    DATA.O2data{1}(:,3) = FVd.data(:,iS);
    DATA.O2data{1}(:,4) = FVd.data(:,iP);
    DATA.O2data{1}(:,5) = FVd.data(:,iT);
    DATA.O2data{1}(:,7) = FVd.data(:,iO);
    DATA.O2data{1}(:,11) = FVd.data(:,iOsat);
    %--------------------
    % (2) telemetry air data -- Not available for SOLOs!
    DATA.O2phase{2} = NaN;
    DATA.O2data{2} = NaN;
    % (3) initialize in-air sqnc (nr sfc)
    DATA.O2phase{3} = nan(50000,11);
    DATA.O2data{3} = nan(50000,11);
    % (3) initialize in-air sqnc (true in air)
    DATA.O2phase{4} = nan(50000,11);
    DATA.O2data{4} = nan(50000,11);
    
    DATA.O2air{1}{1} = []; %No avg inair with telemetry, only sequence data.
    DATA.O2air{2}{1} = []; %No std inair with telemetry, only sequence data.
    DATA.O2air{1}{2} = nan(DATA.Nprof,11); %avg of near-sfc in-air sequence
    DATA.O2air{1}{3} = nan(DATA.Nprof,11); %avg of true air in-air sequence
    DATA.O2air{2}{2} = nan(DATA.Nprof,11); %std of near-sfc in-air sequence
    DATA.O2air{2}{3} = nan(DATA.Nprof,11); %std of true air in-air sequence
    k1 = 1;
    k2 = 1;
    k3 = 1; %index o2air only for cycles that have in-air data (not every cycle taken, previously was indexing using iF)
    for iF = 1:DATA.Nprof
        tmpfile = [Dpath,F(iF,:)];
        D = load(tmpfile);
        
        %LG and JP 9/13/2024: SPECIFIC CONDITION FOR FLOAT 3902533 (ss4018 :[ ):
        %Cycle 12 has fill values of -999 in the aircal data leading to imaginary numbers when calculating air ref data

        %air_flds retrieves the names of objects in D.TRAJ.InAir as a cell array
        air_flds = fieldnames(D.TRAJ.InAir);

        %"missingAirCal" tests whether any of the data in the air cal section of the SOLO .dox file (last few lines of data)
        % are missing, which pertains to -999, in columns 6 (phase delay doxy) or column 8 (temp doxy). If a -999 is present
        % in these columns, it screws up the NCEP reference calculation as imaginary numbers.
        missingAirCal = D.TRAJ.InAir.RAW(:,6) == -999 | D.TRAJ.InAir.RAW(:,8) == -999;

        %Checks if there are any missing air cal values (-999)
        if any(missingAirCal)
            for ct = 1:size(air_flds,1) %Loops through all sub structures in D.TRAJ.InAir
                %Subset the data within the sub structure so that only data that IS NOT missing is included. Missing data is dropped.
                D.TRAJ.InAir.(air_flds{ct})=D.TRAJ.InAir.(air_flds{ct})(~missingAirCal,:);
            end
        end
        %Continue with the code as normal.

        %take median of gps if multiple.
        gps = median(D.INFO.gps,1,'omitnan');
        trck = [gps(1) D.INFO.cast gps(2:3)];
        TRACK(iF,:) = trck;
        % organize O2 data (in air sequence data)
        %--------------------
        % (3) in-air sequence near-sfc data -- population of this field is
        % subject to change...depending on how we identify "usable" in air
        % samples...per John!  (SOLO is a bit different than standard APEX
        % sample scheme)
        nrsfc = find(D.TRAJ.InAir.IN_AIR_LOGICAL==0); %using 'find' to get around the "-1" first sample..cludgy
        truesfc = find(D.TRAJ.InAir.IN_AIR_LOGICAL==1);
        if isempty(truesfc) %no true in-air data
            continue
        end
        % [sdn, cast, s, p, t, phase, (Tvolt)]
        DATA.O2phase{3}(k1:k1+length(nrsfc)-1,1) = repmat(D.INFO.sdn,length(nrsfc),1);
        DATA.O2phase{3}(k1:k1+length(nrsfc)-1,2) = repmat(D.INFO.cast,length(nrsfc),1);
        DATA.O2phase{3}(k1:k1+length(nrsfc)-1,3) = zeros(length(nrsfc),1); %not carried through in phy files
        DATA.O2phase{3}(k1:k1+length(nrsfc)-1,4) = zeros(length(nrsfc),1); %not carried through in phy files
        DATA.O2phase{3}(k1:k1+length(nrsfc)-1,5) = D.TRAJ.InAir.TEMP_DOXY(nrsfc);
        DATA.O2phase{3}(k1:k1+length(nrsfc)-1,6) = D.TRAJ.InAir.PHASE_DELAY_DOXY(nrsfc);
        %[sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat]; % umol /L & % sat
%         DATA.O2data{3} = nan(length(nrsfc),11);
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,1) = repmat(D.INFO.sdn,length(nrsfc),1);
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,2) = repmat(D.INFO.cast,length(nrsfc),1);
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,3) = zeros(length(nrsfc),1); %not carried through in phy files
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,4) = zeros(length(nrsfc),1); %not carried through in phy files
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,5) = D.TRAJ.InAir.TEMP_DOXY(nrsfc);
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,6) = D.TRAJ.InAir.PHASE_DELAY_DOXY(nrsfc);
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,9) = D.TRAJ.InAir.PPOX_DOXY(nrsfc);
        DATA.O2data{3}(k1:k1+length(nrsfc)-1,10) = D.TRAJ.InAir.PH2O(nrsfc);
        
        %--------------------
        % (4) in-air sequence true in-air data -- population of this field is
        % subject to change...depending on how we identify "usable" in air
        % samples...per John!  (SOLO is a bit different than standard APEX
        % sample scheme)
        % [sdn, cast, s, p, t, phase, (Tvolt)]
        DATA.O2phase{4}(k2:k2+length(truesfc)-1,1) = repmat(D.INFO.sdn,length(truesfc),1);
        DATA.O2phase{4}(k2:k2+length(truesfc)-1,2) = repmat(D.INFO.cast,length(truesfc),1);
        DATA.O2phase{4}(k2:k2+length(truesfc)-1,3) = zeros(length(truesfc),1); %not carried through in phy files
        DATA.O2phase{4}(k2:k2+length(truesfc)-1,4) = zeros(length(truesfc),1); %not carried through in phy files
        DATA.O2phase{4}(k2:k2+length(truesfc)-1,5) = D.TRAJ.InAir.TEMP_DOXY(truesfc);
        DATA.O2phase{4}(k2:k2+length(truesfc)-1,6) = D.TRAJ.InAir.PHASE_DELAY_DOXY(truesfc);
        %[sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat]; % umol /L & % sat
        DATA.O2data{4}(k2:k2+length(truesfc)-1,1) = repmat(D.INFO.sdn,length(truesfc),1);
        DATA.O2data{4}(k2:k2+length(truesfc)-1,2) = repmat(D.INFO.cast,length(truesfc),1);
        DATA.O2data{4}(k2:k2+length(truesfc)-1,3) = zeros(length(truesfc),1); %not carried through in phy files
        DATA.O2data{4}(k2:k2+length(truesfc)-1,4) = zeros(length(truesfc),1); %not carried through in phy files
        DATA.O2data{4}(k2:k2+length(truesfc)-1,5) = D.TRAJ.InAir.TEMP_DOXY(truesfc);
        DATA.O2data{4}(k2:k2+length(truesfc)-1,6) = D.TRAJ.InAir.PHASE_DELAY_DOXY(truesfc);
        DATA.O2data{4}(k2:k2+length(truesfc)-1,9) = D.TRAJ.InAir.PPOX_DOXY(truesfc);
        DATA.O2data{4}(k2:k2+length(truesfc)-1,10) = D.TRAJ.InAir.PH2O(truesfc);
        
        %Now compute avg and std in-air for near sfc and true in air.
        DATA.O2air{1}{2}(k3,:) = mean(DATA.O2data{3}(k1:k1+length(nrsfc)-1,:),1,'omitnan'); %avg of near-sfc in-air sequence
        DATA.O2air{1}{3}(k3,:) = mean(DATA.O2data{4}(k2:k2+length(truesfc)-1,:),1,'omitnan'); %avg of true air in-air sequence
        DATA.O2air{2}{2}(k3,:) = std(DATA.O2data{3}(k1:k1+length(nrsfc)-1,:),1,'omitnan'); %std of near-sfc in-air sequence
        DATA.O2air{2}{3}(k3,:) = std(DATA.O2data{4}(k2:k2+length(truesfc)-1,:),1,'omitnan'); %std of true air in-air sequence
        k1 = k1+length(nrsfc);
        k2 = k2+length(truesfc);
        k3 = k3+1;
    end
    DATA.track = TRACK;
    %trim arrays
    DATA.O2air{1}{2}(k3:end,:)=[]; %trim
    DATA.O2air{1}{3}(k3:end,:)=[]; %trim
    DATA.O2air{2}{2}(k3:end,:)=[]; %trim
    DATA.O2air{2}{3}(k3:end,:)=[]; %trim
    DATA.O2phase{3}(k1:end,:)=[];
    DATA.O2data{3}(k1:end,:)=[];
    DATA.O2phase{4}(k2:end,:)=[];
    DATA.O2data{4}(k2:end,:)=[];
else  % if ~isempty(F); else
    warning('There are no .mat files to extract');
    msgbox('WARNING: There are no .mat files to extract.','WARNING')
end %end if ~isempty(F)
% close(WB)
%end % end getall_SOLOdata_sO2


