function adj_data = apply_QC_corr(data, cast_sdn, QC)
% PURPOSE:
%   This function applies the quality control adjustments to nitrate, pH,
%   and other(?) float data. Adjustments are made based on gain, offset
%   and drift corrections found in FloatQCList.txt file.
%
% USAGE:
%	adj_data = apply_QC_corr(data, cast_SDN, QC)
%              PTS are needed for pH QC
%
% INPUTS:
%	data      = matrix [P, T, S, data array to be adjusted]
%               P & T needed for pH adjustment, but carried S along - maybe
%               can be used later(?)
%   cast_sdn  = profile termination time as a Matlab SDN
%   QC        = QC structure for given sensor structure includes QC data
%              [SDN CAST# GAIN OFFSET DRIFT]
%
% OUTPUTS:
%   adj_data =  array of adjusted data
%
% SPECIAL NOTES:A GAIN value appears in every QC_data row but to date only
%               one gain value is applied to the whole float data set.
%
%               pH QC data may have a pump offset term to correct for a
%               slight VRS change when the pump turns on continuously for
%               CP mode. This is incorperated in ALL corrects, just set to
%               zero when not used
%
% Created 1/11/2016 by jp
%
% CHANGES:
% 9/10/2018 TM, change correction scheme to better match MATLAB ischange.  
%    Now corrections are derived and applied on a per-cycle basis
%    (segments between change points are treated as discontinuous with independent
%    drifts and offsets), instead of cumulative across cycles (previously,
%    complex code, offsets were removed then added back in, so offset showing
%    in the correction matrix was cumulative. This resulted in essentially
%    same result as current scheme, but math/code was much more confusing....).
%
% 07/13/2020 TM, previously application of pH adjustments assumed reference
%    temperature of 2degC, modified to use temperature at 1500m (if exist)
%
% 10/10/2024 JP, Clean up for SAGEv2 update QClist file output and for
%    cycle by cycle pH pump offset correction which is now done in
%    Process_APEX_float and its allies The pH ref temperature correction
%    still needs to be updated to match the ref P range which is now
%    available from the QC list file

% **** TEST ****
% cast_sdn  = INFO.sdn;
% %data      = [LR.PRES, LR.TEMP, LR.PSAL, phtot];
% QC_data   = cal.pH.QC;
% data      = [LR.PRES, LR.TEMP, LR.PSAL, LR.NITRATE];
% %QC_data   = cal.N.QC;
% cal       = cal.N;

% data = QCD;
% cast_sdn = d.sdn;
% cal = cal.pH;
% **************

% ************************************************************************
fv.bio    = 99999; % bio argo data fill value

adj_data = data(:,1) * 0 + fv.bio; % defualt if NO QC
if isempty('QC')
    disp('NO QC STUCTURE FOUND')
    return
end

% ************************************************************************
% DETERMINE IF pH SENSOR
if strcmpi(QC.type,'pH') % pH data!!
    % pump_offset = QC.pHpumpoffset; % for pH only
    P           = data(:,1);
    T           = data(:,2);
    %S           = data(:,3); % NOT USED FOR NOW
    data        = data(:,4);
    % Find temperature at 1500m, otherwise use T = 2C (if shallow)
    pres_tol = 500; % if shallow, take T up to 1000m depth
    p1500   = abs(P- 1500);
    min1500 = min(p1500);
    if min1500 < pres_tol % look for 1500m sample first
        ind = find(p1500 == min1500,1);
        if ~isempty(ind)
            disp('Calculating reference temp at 1500m for application of pH correction to k0...')
            TREF = T(ind);
            disp(['T = ',num2str(TREF),' degC'])
            disp(['refdepth = ',num2str(P(ind))])
            TCOR        = (TREF + 273.15)./(T + 273.15); % drift + offset is f(T) too, brings correction to k0 space
        else
            disp('Could not calculate reference temp at 1500m for application of pH correction to k0.  Using default reference temp of 2degC.')
            TCOR        = (2 + 273.15)./(T + 273.15); % drift + offset is f(T) too????
        end
    else
        disp('Could not calculate reference temp at 1500m for application of pH correction to k0.  Using default reference temp of 2degC.')
        TCOR        = (2 + 273.15)./(T + 273.15); % drift + offset is f(T) too????
    end
    clear P T S tP
    
else % ALL OTHER DATA
    data = data(:,4);
end

% ************************************************************************
t_fv = data == fv.bio; % flag fill values in data

% IF QC LEG DOES NOT START ON PROFILE #1 SET QC DATA = NAN
if QC.steps(1,1) > cast_sdn % 1st QC step does not start at profile 1
    fprintf(['WARNING: QC node does not start on 1st float cycle. Cycles ', ...
        'Adj data for cycles < 1st node wil be set = NaN\n'])
    adj_data = data * 0 + fv.bio; %jp 9/23/19 9634 cycle 138- clock bad cuasing problems
    return
end

% REMOVE STEPS NOT NEEDED FOR DATA CORRECTION CALCULATION BECAUSE FLOAT
% HAS NOT GOT THERE YET
t1 = QC.steps(:,1) > cast_sdn;
QC.steps(t1,:) = []; % present cast is before step so remove step (not needed)
% [sdn cycle gain offset drift]

% ADJUSTMENTS ARE DISCONTINUOUS ACROSS SEGMENTS AS OF 9/2018 %%%
% since we removed all steps with date later than current cycle, the
% drift and gain of interest for current cycle are in the last step
G_at_last_step = QC.steps(end,3); %this is the gain at the last step
D_at_last_step = QC.steps(end,5); %this is the drift, over time since the last step
O_at_last_step = QC.steps(end,4); %this is the offset at last step

%now find time difference between last step and the current cycle.
time_since_last_step = cast_sdn - QC.steps(end,1);

% *************************************************************************
% OXYGEN: FIRST CALCULATE WHAT GAIN SHOULD BE FOR THE INPUT CYCLE USING THE
% QC MATRIX PROVIDED.  THEN APPLY GAIN AS SIMPLE MULTIPLIER.
% [SDN CAST# GAIN OFFSET DRIFT]
% *************************************************************************
if strcmpi(QC.type,'Oxygen')
    %now use time difference between last step and the current cycle.  If
    %they are the same, then the gain is simply "G_at_last_step", otherwise
    %it is "G_at_last_step" + "D_at_last_step"*time_since_last_step
    time_since_last_step = cast_sdn - QC.steps(end,1);
    if time_since_last_step == 0
        usethisG = G_at_last_step;
    else
        usethisG = D_at_last_step.*time_since_last_step./365+G_at_last_step;
    end
    adj_data = data.*usethisG;
end

% *************************************************************************
% FOR NO3, PH, APPLY OFFSET AND DRIFT TO DATA AND THEN APPLY GAIN
% GAIN DOES NOT VARY WITH TIME
% [SDN CAST# GAIN OFFSET DRIFT]
%
% AS OF 9/2018 WE ARE NOW TREATING EACH SEGMENT INDEPENDENTLY (ADJUSTMENTS
% ARE NOT CUMULATIVE ACROSS TIME SERIES).  THIS IS MUCH MORE
% STRAIGHTFORWARD FOR THE USER, AND REMOVES THE COMPLEXITY OF THE
% CALCULATIONS.  RESULTS ARE THE SAME. CODE IS JUST MORE STRAIGHTFORWARD.
% *************************************************************************

% *************************************************************************
% NITRATE QC
% *************************************************************************
if strcmpi(QC.type,'Nitrate')
    % CALCULATE OFFSET AND DRIFT
    dt  = (cast_sdn - QC.steps(end,1)) ./ 365; % elapsed years in QC leg
    cor = O_at_last_step + D_at_last_step .* dt; % TM (as of 9/2018)
    adj_data       = (data - cor) ./ G_at_last_step;
    adj_data(t_fv) = fv.bio; % Put fill values back in
end

% *************************************************************************
% pH QC
% pH ASLO HAS A PUMP OFFSET TERM
% *************************************************************************
if strcmpi(QC.type,'pH'); %DO pH QC
    dt  = (cast_sdn - QC.steps(end,1)) ./ 365; % elapsed years in QC leg
    cor = O_at_last_step + D_at_last_step .* dt; % TM (as of 9/2018)
    %adj_data       = data + (pump_offset - cor) .* TCOR;
    adj_data       = data  - cor.* TCOR;
    adj_data(t_fv) = fv.bio; % Put fill values back in
end

% % *************************************************************************
% % CHL QC
% % *************************************************************************
% % ONLY ONE LINE OF CORRECTION DATA (GAIN AND OFFSET)
% if strcmpi(QC.type,'CHL')
%     adj_data       = data * QC.steps(1,3) - QC.steps(1,4);
%     adj_data(t_fv) = fv.bio; % Put fill values back in
% end
% 
% % *************************************************************************
% % BACKSSCATTER QC
% % *************************************************************************
% % ONLY ONE LINE OF CORRECTION DATA (GAIN AND OFFSET)
% % SO FAR GAIN IS 1 so basically a DC correction after the fact but should
% % the subtraction be done at the cts level not the VSF
% if strcmpi(QC.type,'BB')
%     adj_data       = data * QC.steps(1,3) - QC.steps(1,4);
%     adj_data(t_fv) = fv.bio; % Put fill values back in
% end
% 
% % *************************************************************************
% % CDOM QC
% % *************************************************************************
% % ONLY ONE LINE OF CORRECTION DATA (GAIN AND OFFSET)
% if strcmpi(QC.type,'CDOM')
%     adj_data       = data * QC.steps(1,3) - QC.steps(1,4);
%     adj_data(t_fv) = fv.bio; % Put fill values back in
% end


