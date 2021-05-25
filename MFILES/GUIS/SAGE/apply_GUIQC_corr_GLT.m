function new_qc_data = apply_GUIQC_corr_GLT(handles,DATA)
% PURPOSE:
%   This function applies the quality control adjustments to the
%   "handles.qc_data" matrix built by the float_qc.m GUI. Only used for N03
%   O2 and pH
%
%   These changes only occur locally in the handles structure
%
%   Adjustments are made based on gain, offset and drift corrections
%
% INPUTS:
%	handles     The structure that is passed to all the functions in the
%
%
% SPECIAL NOTES:A GAIN value appears in every QC_data row but to date only
%               one gain value is applied to the whole float data set.
%
%               pH QC data may have a pump offset term to correct for a
%               slight VRS change when the pump turns on continuously for
%               CP mode. This is incorperated in ALL corrects, just set to
%               zero when not used
%
% APPLY OFFSET AND DRIFT TO DATA AND THEN APPLY GAIN
% GAIN DOES NOT VARY WITH TIME
% [CYCLE GAIN OFFSET DRIFT]
% *************************************************************************
%
% Created 9/28/2016 by jp
% Updated 07/13/2020 by TM, previously application of pH adjustments assumed reference temperature of 2degC, modified to use temperature at 1500m (if exist)

% GET SOME DATA AND INFORMATION
raw_data  = handles.raw_data;
qc_data   = handles.qc_data;
QCA       = handles.QCA;
qca       = QCA.(DATA.paramtag);% Sensor specific
rr        = size(qca,1);
%[~,ia,~]  = unique(raw_data(:,2)); % unique casts
%casts     = raw_data(ia,1:2);% SDN, CAST
cor       = ones(size(raw_data.data(:,1)));

% GET INDICES
iP    = find(strcmp('Pressure[dbar]',raw_data.hdr)   == 1);
iT    = find(strcmp('Temperature[°C]',raw_data.hdr)  == 1);
iS    = find(strcmp('Salinity[pss]',raw_data.hdr)    == 1);
iZ    = find(strcmp('Depth[m]',raw_data.hdr)         == 1);
iO    = find(strcmp('Oxygen[µmol/kg]',raw_data.hdr)  == 1);
iN    = find(strcmp('Nitrate[µmol/kg]',raw_data.hdr) == 1);
iPh   = find(strcmp('pHinsitu[Total]',raw_data.hdr)  == 1);
iCYC   = find(strcmp('Station',raw_data.hdr)  == 1);


% IF QC LEG DOES NOT START ON PROFILE #1 SET DATA = NAN

if qca(1,1) > min(raw_data.data(:,2)) % 1st QC step does not start at profile 1
    t1 = raw_data.data(:,2) < qca(1,1);
    raw_data.data(t1,:) = NaN;
end

% *************************************************************************
% pH QC
% pH ALSO HAS A PUMP OFFSET TERM
% NO GAIN CORRECTION FOR PH
% *************************************************************************
if strcmp(DATA.paramtag, 'PH')
    IND = iPh;
    if isempty(QCA.PH_OFFSET)
        pump_offset = 0;
    else
        pump_offset = QCA.PH_OFFSET;
    end
    tP          = raw_data.data(:,iP) < 980; % logical flag for pump on calc
    pump_offset = pump_offset * tP; % array of 0's and offset
    
    % Find temperature at 1500m for each profile, otherwise use T = 2C (if shallow)
    thecycs = unique(raw_data.data(:,iCYC));
    TCOR = []; %initialize
    for icyc = 1:length(thecycs)
        proftmp = raw_data.data(raw_data.data(:,iCYC)==thecycs(icyc),:);
        pres_tol = 500; %if shallow, take T up to 1000m depth
        p1500   = abs(proftmp(:,iP)- 1500);
        min1500 = min(p1500);
        if min1500 < pres_tol; % look for 1500m sample first
            ind = find(p1500 == min1500,1);
            if ~isempty(ind)
                disp('Calculating reference temp at 1500m for application of pH correction to k0...')
                TREF = proftmp(ind,iT);
                disp(['T = ',num2str(TREF),' degC'])
                disp(['refdepth = ',num2str(proftmp(ind,iP))])
                TCORtmp        = (TREF + 273.15)./(proftmp(:,iT) + 273.15); % drift + offset is f(T) too, brings correction to k0 space
            else
                disp('Could not calculate reference temp at 1500m for application of pH correction to k0.  Using default reference temp of 2degC.')
                TCORtmp        = (2 + 273.15)./(proftmp(:,iT) + 273.15); % drift + offset is f(T) too, brings correction to k0 space
            end
        else
            disp('Could not calculate reference temp at 1500m for application of pH correction to k0.  Using default reference temp of 2degC.')
            TCORtmp        = (2 + 273.15)./(proftmp(:,iT) + 273.15); % drift + offset is f(T) too, brings correction to k0 space
        end
        TCOR = [TCOR;TCORtmp];
    end
    if length(TCOR)~=size(raw_data.data,1)
        disp('WARNING, TCOR VECTOR LENGTH MISMATCH!')
    end
    
    
    last_cor = 0;
    for i = 1:rr
        if i < rr
            t1 = raw_data.data(:,2) >= qca(i,1) &  ...
                raw_data.data(:,2) < qca(i+1,1); %get block
        else
            t1 = raw_data.data(:,2) >= qca(i,1); %get block
        end
        dt = (raw_data.data(t1,1) - min(raw_data.data(t1,1))) ./ 365; % elapsed years in block
        tmp_cor  = last_cor + qca(i,3) + qca(i,4) .* dt; %offset + drift
        %         last_cor = tmp_cor(end);
        cor(t1)  = tmp_cor;
    end
    %     qc_data.data(:,IND) = raw_data.data(:,IND) + cor + ...
    %                            (pump_offset - cor) .* TCOR;
    qc_data.data(:,IND) = raw_data.data(:,IND) + ...
        (pump_offset - cor) .* TCOR;
    
    % SET BAD QC DATA TO NaN
    t1 = qc_data.data(:,IND+1) == 8;
    qc_data.data(t1,IND) = NaN;
end

% *************************************************************************
% NITRATE QC
% CYCLE, GAIN, OFFSET, DRIFT
% *************************************************************************
if strcmp(DATA.paramtag, 'NO3')
    IND = iN;
    last_cor = 0;
    for i = 1:rr
        if i < rr
            t1 = raw_data.data(:,2) >= qca(i,1) & ...
                raw_data.data(:,2) < qca(i+1,1); %get data block
        else
            t1 = raw_data.data(:,2) >= qca(i,1); %get data block
        end
        dt = (raw_data.data(t1,1) - min(raw_data.data(t1,1))) ./ 365; % elapsed years in block
        tmp_cor  = last_cor + qca(i,3) + qca(i,4) .* dt; %offset + drift
        %         last_cor = tmp_cor(end);
        cor(t1)  = tmp_cor;
    end
    % %     figure
    % %     plot(cor,'r.-')
    qc_data.data(:,IND)   = (raw_data.data(:,IND) - cor) ./ qca(i,2);
    
    % SET BAD QC DATA TO NaN
    t1 = qc_data.data(:,IND+1) == 8;
    qc_data.data(t1,IND) = NaN;
end

% % % % *************************************************************************
% % % % OXYGEN QC
% % % % GAIN ONLY AND ONLY ONE GAIN
% % % % THIS WILL PROBALY CHANGE IN THE NEAR FUTURE
% % % % *************************************************************************
% % % % if strcmp(handles.info.data_type, 'O2')
% % % %     IND = iO;
% % % %     qc_data.data(:,IND)  = raw_data.data(:,IND) .* qca(i,2);
% % % % end;
% % %
% % % if strcmp(DATA.paramtag, 'O2')
% % %     IND = iO;
% % %     last_cor = 0;
% % %     for i = 1:rr
% % %         if i < rr
% % %             t1 = raw_data.data(:,2) >= qca(i,1) & ...
% % %                 raw_data.data(:,2) <= qca(i+1,1); %get block
% % %         else
% % %             t1 = raw_data.data(:,2) >= qca(i,1); %get block
% % %         end
% % %         dt = (raw_data.data(t1,1) - min(raw_data.data(t1,1))) ./ 365; % elapsed years in block
% % %         tmp_cor  = last_cor + qca(i,3) + qca(i,4) .* dt; %offset + drift
% % %         last_cor = tmp_cor(end);
% % %         cor(t1)  = tmp_cor;
% % %     end
% % %     qc_data.data(:,IND)   = (raw_data.data(:,IND) - cor) .* qca(i,2);
% % %
% % %     % SET BAD QC DATA TO NaN
% % %     t1 = qc_data.data(:,IND+1) == 8;
% % %     qc_data.data(t1,IND) = NaN;
% % % end


if isfield(handles,'new_qc_data') == 1 %if other variables have been updated, preserve them!
    handles.new_qc_data.data(:,IND) = qc_data.data(:,IND);
    handles.new_qc_data.data(:,IND+1) = qc_data.data(:,IND+1);
    new_qc_data = handles.new_qc_data;
else
    new_qc_data = qc_data;
end


