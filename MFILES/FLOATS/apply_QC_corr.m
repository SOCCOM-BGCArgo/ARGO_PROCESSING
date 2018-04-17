function adj_data = apply_QC_corr(data, cast_sdn, QC)
% PURPOSE: 
%   This function applies the quality control adjustments to nitrate, pH,
%   and other(?) float data. Adjustments are made based on gain, offset
%   and drift corrections found in FloatQCList.txt file.
%
% USAGE:
%	adj_data = apply_QC_corr(data, cast_SDN, cal)
%              PTS are included for pH QC
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
fv.QC     = 99;
adj_data = data(:,1) * 0 + fv.bio; % defualt if NO QC
if isempty('QC')
    disp('NO QC STUCTURE FOUND')
    return
end


%DETERMINE IF pH SENSOR
if strcmpi(QC.type,'pH') % pH data!!
    pump_offset = QC.pHpumpoffset; % for pH only
    P           = data(:,1);
    tP          = P < 980; % logical flag for pump on calc
    pump_offset = pump_offset * tP; % array of 0's and offset
    T           = data(:,2);
    %S           = data(:,3); % NOT USED FOR NOW
    data        = data(:,4);
    TCOR        = (2 + 273.15)./(T + 273.15); % drift + offset is f(T) too????
    clear P T S tP
    
else % ALL OTHER DATA
    data = data(:,4); 
end
t_fv = data == fv.bio; % flag fill values in data

% IF QC LEG DOES NOT START ON PROFILE #1 SET QC DATA = NAN
if QC.steps(1,1) > cast_sdn % 1st QC step does not start at profile 1
    adj_data = data * 0 + fv.QC;
    return
end  
      
% REMOVE STEPS NOT NEEDED FOR DATA CORRECTION CALCULATION BECAUSE FLOAT
% HAS NOT GOT THERE YET    
t1 = QC.steps(:,1) > cast_sdn;
QC.steps(t1,:) = []; % present cast before step so remove step (not needed)
% [sdn cycle gain offset drift]

rr = size(QC.steps,1); % # of QC steps
cor = 0; % start with zero corection

% *************************************************************************
% APPLY OFFSET AND DRIFT TO DATA AND THEN APPLY GAIN
% GAIN DOES NOT VARY WITH TIME
% [SDN CAST# GAIN OFFSET DRIFT]
% *************************************************************************

% *************************************************************************
% NITRATE QC
% *************************************************************************
if strcmpi(QC.type,'Nitrate');
% CALCULATE OFFSET AND DRIFT  ***  JP STYLE  ***  
    for i = 1:rr
        if i == rr % final QC row, dt varies with time
            dt  = (cast_sdn - QC.steps(i,1)) ./ 365; % elapsed years in QC leg
        else       % dt for previous QC rows is constant
            dt  = (QC.steps(i+1,1) - QC.steps(i,1)) ./ 365;
        end
        cor = cor + QC.steps(i,4) + QC.steps(i,5) .* dt; % JP & KJ
    end
    adj_data       = (data - cor) ./ QC.steps(i,3); 
    adj_data(t_fv) = fv.bio; % Put fill values back in
end

% *************************************************************************
% pH QC
% pH ASLO HAS A PUMP OFFSET TERM
% *************************************************************************
if strcmpi(QC.type,'pH'); %DO pH QC
    for i = 1:rr
        if i == rr % final QC row, dt varies with time
            dt  = (cast_sdn - QC.steps(i,1)) ./ 365; % elapsed years in QC leg 
        else       % dt for previous QC rows is constant
            dt  = (QC.steps(i+1,1) - QC.steps(i,1)) ./ 365;  
        end
        cor = cor + QC.steps(i,4) + QC.steps(i,5) .* dt; % JP & KJ
    end
    adj_data       = data + (pump_offset - cor) .* TCOR;
    adj_data(t_fv) = fv.bio; % Put fill values back in
end
    
% *************************************************************************
% CHL QC 
% *************************************************************************
% ONLY ONE LINE OF CORRECTION DATA (GAIN AND OFFSET)
if strcmpi(QC.type,'CHL')
    adj_data       = data * QC.steps(1,3) - QC.steps(1,4); 
    adj_data(t_fv) = fv.bio; % Put fill values back in
end  

% *************************************************************************
% BACKSSCATTER QC 
% *************************************************************************
% ONLY ONE LINE OF CORRECTION DATA (GAIN AND OFFSET)
% SO FAR GAIN IS 1 so basically a DC correction after the fact but should
% the subtraction be done at the cts level not the VSF
if strcmpi(QC.type,'BB')
    adj_data       = data * QC.steps(1,3) - QC.steps(1,4); 
    adj_data(t_fv) = fv.bio; % Put fill values back in
end

% *************************************************************************
% CDOM QC 
% *************************************************************************
% ONLY ONE LINE OF CORRECTION DATA (GAIN AND OFFSET)
if strcmpi(QC.type,'CDOM')
    adj_data       = data * QC.steps(1,3) - QC.steps(1,4); 
    adj_data(t_fv) = fv.bio; % Put fill values back in
end
    
        