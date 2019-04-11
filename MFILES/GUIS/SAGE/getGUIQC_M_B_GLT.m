function CGOD = getGUIQC_M_B_GLT(DATA,handles)
%
% getGUIQC_M_B is a SAGE helper function.
% This function calculates quality control adjustments: drift(M) and
% offset(B) between cycle breaks in the QC table([Cycle Gain Offset Drift])
% for N & pH. If the cycle difference between steps = 1 the slope is
% set = 0 and only an offset is determined.

% CHANGE LOG
% 03/23/2017 - Updated to include Ken's global MLR for pH -jp
% 08/14/2017 -  Added code to include nitrate gain correction -jp 
%

CGOD = []; % default output if no claculations performed

% ************************************************************************
% GET SOME REQUIRED INFO ABOUT THE DATA TO USE
% END FUNCTION IF NOT NO3 or pH FOR NOW
data_type    = DATA.paramtag; % NO3 PH or O2
compdata_str = DATA.reftag; % comp data type

if regexp(data_type, 'O2|S|T') % NO OXYGEN DRIFT YET, NO S & T corr applied 
    return
end

% ************************************************************************
% GET SELECTED FLOAT DATA (RAW & QC)
% SET MISSING VALUES TO NaN, GET COLUMN INDICES
% ************************************************************************
% dRAW = DATA.datasub;
dRAW = DATA.rawsub;
% dRAW(dRAW== -1e10) = NaN; % missing values

% IF NO QC DATA THIS WILL EQUAL RAW
dQC  = DATA.qcsub; % O in QC data needed for MLR calcs
dQC(dQC == -1e10) = NaN; % missing values

% GET FLOAT DATA INDICES
iP    = find(strcmp('Pressure[dbar]',handles.raw_data.hdr)   == 1);
iT    = find(strcmp('Temperature[°C]',handles.raw_data.hdr)  == 1);
iS    = find(strcmp('Salinity[pss]',handles.raw_data.hdr)    == 1);
iO    = find(strcmp('Oxygen[µmol/kg]',handles.raw_data.hdr)  == 1);
iN    = find(strcmp('Nitrate[µmol/kg]',handles.raw_data.hdr) == 1);
iPH   = find(strcmp('pHinsitu[Total]',handles.raw_data.hdr)  == 1);

if strcmp(data_type, 'NO3') % GET PROPER FLOAT DATA INDEX
    IX = iN;
elseif strcmp(data_type, 'PH')
    IX = iPH;
else
    IX = [];
end

% ONLY USE NON NAN DATA
%tnan = isnan(dQC.data(:,IX));
% tnan = isnan(dQC(:,IX))| isnan(dRAW(:,IX));
% dRAW(tnan,:) =[];
% dQC(tnan,:) =[];

% ************************************************************************
% GET SELECTED COMPARISON DATA
% ************************************************************************
compdata_str = DATA.reftag;
comp_data = DATA.refsub;
comp_data(dRAW== -1e10) = NaN;
dRAW(dRAW== -1e10) = NaN; % missing values
% diff_data = DATA.diffsub;
raw_data = dRAW;
diff_data = raw_data(:,IX)-comp_data;

% % COMP DATA IS WOA 2013 NITRATE
% if strcmp(compdata_str,'WOA') && strcmp(data_type, 'NO3')
%     comp_data = handles.WOA;
%     comp_data(tnan) = [];
%     
% % COMP DATA IS CANYON NEURAL NETWORK ESTIMATE    
% elseif strcmp(compdata_str,'CANYON') 
%     C = handles.canyon;
%     if isempty(C) % NO CANYON DATA SO END PROGRAM
%         return
%     end;
%     if strcmp(data_type, 'NO3')
%         IND = find(strcmp('canyon_no3',C.hdr) == 1);
%         comp_data  = C.data(~tnan,IND);
%     elseif strcmp(data_type, 'PH')
%         IND = find(strcmp('canyon_ph',C.hdr)  == 1);
%         comp_data  = C.data(~tnan,IND);
%     else
%         comp_data  = [];
%     end
%     clear C IND
%     
% % CHECK FOR WILLIAMS MLRs - USE QC DATA FOR THIS!!!!    
% elseif strncmp(compdata_str,'MLR W50to80',8) 
%     MLR  = handles.WRK_MLR;
%     potT = theta(dQC.data(:,iP), dQC.data(:,iT), dQC.data(:,iS),0);
%     sig_theta  = density(dQC.data(:,iS), potT)-1000; %p =0, t= pot temp
%     
%     comp_data = MLR.cC + dQC.data(:,iO)*MLR.cO + dQC.data(:,iS)*MLR.cS + ...
%         dQC.data(:,iT)*MLR.cT + sig_theta*MLR.cST + dQC.data(:,iP)*MLR.cP;
% 
% elseif strcmp(compdata_str,'LIR') % LINR, LIPHR
%     L = handles.LIR;
%     if isempty(L) % NO LIR DATA SO END PROGRAM
%         return
%     end;
%     if strcmp(data_type, 'NO3')
%         IND = find(strcmp('LIR_no3',L.hdr) == 1);
%         comp_data  = L.data(~tnan,IND);
%     elseif strcmp(data_type, 'PH')
%         IND = find(strcmp('LIR_ph',L.hdr)  == 1);
%         comp_data  = L.data(~tnan,IND);
%     else
%         comp_data  = [];
%     end
%     clear L IND
% 
% else
%     disp('Could not determine comparison data source')
%     return
% end
% 
% clear potT sig_theta MLR  

% ************************************************************************
% NOW SUBSET DATA OVER DEPTH RANGES & GET DATA-COMP DIFF &
% CLEAN UP A BIT
% ************************************************************************
% t1 = dRAW.data(:,iP) >= handles.depth_min.UserData & ...
%      dRAW.data(:,iP) <= handles.depth_max.UserData; 
% 
% % t2 = dRAW.data(:,2) >= handles.profile_min.UserData & ...
% %      dRAW.data(:,2) <= handles.profile_max.UserData;
%  
% raw_data  = dRAW.data(t1,:); % Raw data subsetted over depth range
% comp_data = comp_data(t1);   % Comparison data subsetted over dapth range
% diff_data = raw_data(:,IX) - comp_data;
% clear dRAW dQC t1 t2 iN iO iP iPH iS iT  
%[raw_data(:,IX) comp_data diff_data]
% ************************************************************************
% BRING IN TABLE DATA & GET FLOAT TRACK & DO SOME CHECKS
% ************************************************************************
track = DATA.track;           % [SDN CYCLE LON LAT]
QC    = DATA.tableDATA;   % [CYLE GAIN OFFSET DRIFT]

% CHECK FOR GAIN IN NITRATE - if EXISTS REDO DIFF CALC 08/14/17 jp
% Divide comp data by gain before doing regression
% ideal: Ncomp = sum(Mi*Nflt + B) / G   so... Ncomp * G = sum(Mi*Nflt + B)
if strcmp(data_type, 'NO3') && QC(1,2) ~= 1
    diff_data = raw_data(:,IX) - comp_data.* QC(1,2);
end

% qca's should not end with a slope(drift) = 0, if not 0 add a row
if QC(end,4) ~= 0
    QC = [QC; QC(end,1)+1, QC(1,2), 0, 0];
end

% MAX TABLE CYCLE CAN NOT BE > MAX GOOD DATA CYCLE - RESET IF NEEDED
t_nan = isnan(raw_data(:,IX));
if QC(end,1) > max(raw_data(~t_nan,2))
    QC(end,1) = max(raw_data(~t_nan,2));
end

rows  = size(QC,1);
% sdn   = ones(rows,1) * NaN;
% dt    = diff(sdn);

% for i = 1:size(QC,1) % Get time stamp for each step in qc table
%     t1 = track(:,2) == QC(i,1);
%     sdn(i) = track(t1,1);
% end
    
dPQC = diff(QC(:,1)); % profile # differences

% ************************************************************************
% DATA PREP DONE NOW STEP THROUH QC ADJUSTMENT BREAKS AND GET M & B
% [Cycle Gain Offset Drift]
% ************************************************************************
K = 0; % Coefficient count for correction model = used for AIC & BIC estimate
QCnew = QC; % predimension new QC table data
QCnew(:,3:4) = NaN;
reg = ones(size(QC,1),5)*NaN;
t_good    = ~isnan(diff_data); % Non NaN data
for i = 1:rows
    if i == rows || dPQC(i) == 1 % step = 1 or last line of QC
        K  = K+1; % intercept (offset) only
        t1 = raw_data(:,2) == QCnew(i,1);
        % all data for cast
        if sum(t1) == 0
            reg(i,1:2) = [0 0]; % NO profile!
        else
            reg(i,1:2) = [0 nanmean(diff_data(t1 & t_good))];
        end
    else
        K  = K+2; % intercept and slope
        t1 = raw_data(:,2) >= QCnew(i,1) & raw_data(:,2) < QCnew(i+1,1);
        reg_x = raw_data(t1&t_good,1); % Float regression subset
        if isempty(reg_x)
            continue % could have wrong depth range so no data returned
        end  
        reg_x = reg_x - reg_x(1);
        reg_y = diff_data(t1&t_good);

        [m,b,r,sm,sb] = lsqfity(reg_x, reg_y); %
                    
%       SEE IF SLOPE IS SIGNIFICANT BEFORE USING!
%       SSDX = sqrt(sum((reg_x - mean(reg_x)).^2)); %
%       RSS  = sum((reg_y - (m .* reg_x + b)).^2); % sum squared residuals
%       SEE  = sqrt(RSS ./ (size(reg_x,1)-2)); % std error est RSS/(n-2)
%       T_stat = (m-0) ./ (SEE ./ SSDX)
%       FOR NOW JUST LOOK AT R*R, IF < 0.1 JUST USE OFFSET
        if r*r < 0.1
            reg(i,1:2) = [0 nanmean(diff_data(t1 & t_good))];
        else
            reg(i,:)  = [m,b,r.*r,sm,sb];
        end
      
    end
   
end

% dt    = diff(sdn); % fractions of year
% node_end = ones(size(reg(:,1)))*0;
% node_end(2:end) = dt.*reg(1:end-1,1) + reg(1:end-1,2);
% cum_offset = reg(:,2) - node_end; % To be compatible with Ken's scheme
QCnew(:,3) = reg(:,2); % Offset
QCnew(:,4) = reg(:,1) .* 365; % Drift /yr
% if ~any(isnan(QCnew))
%     CGOD = QCnew;
% else
%     disp(['NaN''s found in QCnew - check function or depth range: ', ...
%           mfilename])
% end
CGOD = QCnew;
if any(isnan(QCnew))
    disp(['NaN''s found in QCnew - check function or depth range: ', ...
          mfilename])
end

    










