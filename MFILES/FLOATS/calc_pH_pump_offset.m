function out = calc_pH_pump_offset(data, CpAct, BadQF)
%function out = calc_pH_pump_offset(data, CpAct, QFtype, method)
% This function attempts to calculate the pH pump offset for a a given pH
% profile using a 1st order polynomial approach as well as a slightly more
% complicated 3rd order polynomial iterative approach.
% function can probably use some refinement
%
%
% INPUTS:
%   data - an Nx2 or Nx3 profile data matrix [PRES PH PHQC]
%          PRES & PH are mandatory, PHQC can be an Argo QF array, 
%          ODV QF array, an array of NaN's or omited
%
%   CpAct - the pressure (dbar) at which the CTD pump switches from spot
%           sampling to continuously on (continuous profiling mode)
%
%   BadQF - bad data quality flag value (ie 4 for Argo & 8 for ODv
%
% OUTPUT:
%   out - a structure with the following fields
%     hdr   -  {'spline offset' 'spline SSR' 'resid std' 'Zish', ...
%               'lin offset' 'lin diff' 'lin_resid std' 'lin Zish'}
%     data  - row array, columns defined by hdr
%     info - nested structure with the following fields
%       CpAct = CpAct P
%       P     = Pressure subset for iterative solution
%       PH    = pH subset for iterative solution

%
%
% ************************************************************************
% ************************************************************************

% % ******* START BUILD TEST DATA *************
% fn = '5906033.TXT'; % 12878
% cycle_num = 42;
% % fn = '5906493.TXT'; % 19378
% % cycle_num = 1;
% % fn = '5906540.TXT'; % un1452
% % cycle_num = 5;
% 
% 
% fd     = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\';
% fp     = fullfile(fd,fn);
% d      = get_FloatViz_data(fp);
% f.STA  = find(strcmp('Station',d.hdr)         == 1);
% f.P    = find(strcmp('Pressure[dbar]',d.hdr)  == 1);
% f.PH   = find(strcmp('pHinsitu[Total]',d.hdr)   == 1); % Raw pH!!!!
% tCYCLE = d.data(:,f.STA) == cycle_num;
% 
% data   = d.data(tCYCLE,[f.P,f.PH,f.PH+1]);
% data(data == -1e10) = NaN;
% CpAct  = 985;
% %Pwin   = 350;
% BadQF  = 8;
% method = 2;
% % ******* END BUILD TEST DATA *************

% ************************************************************************
% ************************************************************************
% SET SOME DEFAULTS
out  = [];
Pwin = 250;
CpAct_buffer = 15; % dP offset (safety gap for hi res sampling
% ************************************************************************

% ************************************************************************
% DO SOME SANITY CHECKS  - ROUND 1

% CHECK SIZE OF DATA MATRIX - IF QF's EXIST USE TO SUBSET DATA
if size(data,2) == 3
    tbad = data(:,3) == BadQF;
    if sum(tbad) > 0
        fprintf('WARNING: %d bad data points detected and removed\n', ...
            sum(tbad));
        data = data(~tbad,:);
    end
    if isempty(data)
        fprintf('WARNING: No data left after quality flag check!\n')
        return
    end
end

% DOES THE FLOAT GO DEEP ENOUGH TO WORRY ABOUT THE PUMP OFFSET?
if max(data(:,1),[],1,'omitnan') <= CpAct
    fprintf(['Max depth(%0.1fm) is <= Cp Activation depth - ',...
        'Can not calculate pump offset\n'], max(data(:,1),[],1,'omitnan'));
    return
end

% ************************************************************************
%     ARRANGE DATA & CREATE SUB SAMPLE INDICES TO DEAL WITH DIFFERENT
%                      FLOAT DP'S & SORT ASCENDING
% ************************************************************************

% FIND DATA ON BOTH SIDES OF CP ACTVATION - may not need really
tCP   = data(:,1) < CpAct; % Cp part of profile
tDS   = data(:,1) > CpAct; % discrete sampling part of profile

% NOW SET UP SAMPLE PRESSURE RANGE TARGET ARRAYS: 
% USE 100 DBAR SPACING & SUBSAMPLE & work away from CpAct P
CP_targets = (CpAct - CpAct_buffer : -100 : CpAct - CpAct_buffer - Pwin)';
CP_targets = flip(CP_targets); % shallow to deep
DS_targets = (CpAct + CpAct_buffer: 100 : CpAct + CpAct_buffer + Pwin)';
CP_inds    = ones(size(CP_targets,1),2)*NaN; % [index, dbar] predim
DS_inds    = ones(size(CP_targets,1),2)*NaN;

% *** CP lookup index & pressure ***
for ct = 1:size(CP_targets,1)
    tmp = abs(data(:,1) - CP_targets(ct));
    [~, ind] = min(tmp,[],1,'omitnan');
    CP_inds(ct,:) = [ind(1),data(ind(1))]; % specify index in case duplicates
end
[~,ia]  = sort(CP_inds(:,2),1,'ascend'); % want CP shllow to deep
CP_inds = CP_inds(ia,:); %indices for offset calcs
tg      = CP_inds(:,2) < CpAct; % double check all P < CpAct
CP_inds = CP_inds(tg,:);

% *** DS Lookup index & pressure ***
for ct = 1:size(DS_targets,1)
    tmp = abs(data(:,1) - DS_targets(ct));
    [~, ind] = min(tmp,[],1,'omitnan');
    DS_inds(ct,:) = [ind(1),data(ind(1))]; % specify index in case duplicates
end
[~,ia] = sort(DS_inds(:,2),1,'ascend'); % want DS shallow to deep
DS_inds = DS_inds(ia,:); %indices for offset calcs
tg      = DS_inds(:,2) > CpAct; % double check all P < CpAct
DS_inds = DS_inds(tg,:);

% ************************************************************************
%       NOW DO SOME SANITY CHECKS BEFORE ESTIMATING PUMP OFFSETS
% ************************************************************************
CPrange = max(data(CP_inds(:,1),1),[],1,'omitnan') - min(data(CP_inds(:,1),1),[],1,'omitnan');
DSrange = max(data(DS_inds(:,1),1),[],1,'omitnan') - min(data(DS_inds(:,1),1),[],1,'omitnan');

if size(CP_inds,1) < 3 || CPrange < 0.75 * Pwin 
    fprintf(['Not enough CP mode data detected to estimate pump offset ', ...
        '(N = %d and range = %0.1f)\n'], size(CP_inds,1), CPrange);
    return
end

if size(DS_inds,1) < 3 || DSrange < 0.75 * Pwin 
    fprintf(['Not enough DS mode data detected to estimate pump offset ', ...
        '(N = %d and range = %0.1f)\n'], size(DS_inds,1), DSrange);
    return
end

% ************************************************************************
% ************************************************************************
%                  CALCULATE PUMP OFFSET ESTIMATES
% ************************************************************************
% ************************************************************************

% ************************************************************************
% ITERATIVE POLYNOMIAL APPROACH
% Uses helper sub function at bottom in fminsearch
% ************************************************************************
poly_order = 3;
fmin_ph_offset  = NaN; % predim
wrk_inds   = [CP_inds(:,1); DS_inds(:,1)]; %indices used in iterative poly approach
wrk_P      = [CP_inds(:,2); DS_inds(:,2)];
tcp        = wrk_P < CpAct; % lower case, used for "Z like" score
P1          = data(wrk_inds,1);
Y1          = data(wrk_inds,2);
 
fmin_ph_offset     = NaN; % set default values
fmin_ph_offset_SSR = NaN;
resid_pH_std       = NaN;
pH_S2N             = NaN;
[pump_offset, FVAL, ExitFlag] = ...
    fminsearch(@(pump_offset)  get_pH_offset_ssr(pump_offset, P1, Y1, ...
    CpAct, poly_order), 0);
if ExitFlag == 1
    fmin_ph_offset     = pump_offset;
    fmin_ph_offset_SSR = FVAL;

    newY1         = Y1; % pH
    newY1(~tcp)   = newY1(~tcp) + fmin_ph_offset; % add offset to bottom
    pcoefs1       = polyfit(P1, newY1, poly_order);
    pH_resid     = newY1 - polyval(pcoefs1, P1);
    resid_pH_std = std(pH_resid); % GET STD OF newY in fit zone
    pH_S2N       = fmin_ph_offset./(resid_pH_std);
end

% ************************************************************************
% SIMPLE LINEAR 2X AVERAGE APPROACH
% Use 1st 2 points closest to CpAct & linear extrap to 1st point on other
% side of CpAct. Shoot from both directions and take the average.
% ************************************************************************

% estimate DS value from CP data
cp_p           = data(CP_inds(end-1:end,1),1); 
cp_ph          = data(CP_inds(end-1:end,1),2);
cp_ph_est      = interp1(cp_p, cp_ph, data(DS_inds(1,1),1),'linear','extrap');
cp_ph_offset   = cp_ph_est - data(DS_inds(1,1),2);

% estimate CP value from DS data
ds_p           = data(DS_inds(1:2,1),1);
ds_ph          = data(DS_inds(1:2,1),2);
ds_ph_est      = interp1(ds_p, ds_ph, data(CP_inds(end,1),1),'linear','extrap');
ds_ph_offset   = data(CP_inds(end,1),2) -ds_ph_est;

lin_ph_offset = (cp_ph_offset + ds_ph_offset)/2; % average of 2 estimates
lin_ph_diff    = cp_ph_offset - ds_ph_offset;

% signal to noise
lin_inds = [CP_inds(end-1:end,1); DS_inds(1:2,1)];
P2           = data(lin_inds,1);
Y2           = data(lin_inds,2);
tcp          = P2 < CpAct;
newY2        = Y2; % pH
newY2(~tcp)  = newY2(~tcp) + lin_ph_offset; % add offset to bottom
pcoefs2       = polyfit(P2, newY2, 1);
lin_pH_resid     = newY2 - polyval(pcoefs2, P2);
lin_resid_pH_std = std(lin_pH_resid); % GET STD OF newY in fit zone
lin_pH_S2N   = lin_ph_offset./(lin_resid_pH_std);



% ************************************************************************
% DEFINE FUNCTION OUTPUTS
out.hdr = {'poly offset' 'poly SSR' 'poly resid std' 'poly Zish', ...
    'linear offset' 'linear diff' 'linear_resid std' 'linear Zish'};
% out.hdr = {'spline offset' 'spline SSR' 'resid std' 'Zish', ...
%     'lin offset' 'lin diff' 'lin_resid std' 'lin Zish'};
out.data = [fmin_ph_offset, fmin_ph_offset_SSR, resid_pH_std, pH_S2N, ...
    lin_ph_offset, lin_ph_diff, lin_resid_pH_std, lin_pH_S2N];

out.info.CpAct = CpAct;
out.info.P  = P1;
out.info.PH = Y1;


% % *************************************************************************
% % *************************************************************************
% % BUILD SUBFUNCTION FOR ITERATIVE POLY APPROACH
% 
    function SSR = get_pH_offset_ssr(Offset, P, Y, CpAct, poly_order)
        % OPTIMIZING FUNCTION FOR PH PUMP OFFSET 3rd ORDER POLY ESTIMATE
        % Offset     = pump offset value
        % P          = Pressure subset
        % Y          = raw pH subset: includes data to above & below CpAct P
        % CpAct = pressure at which ctd goes in to cp mode from spot mode
        % poly_order = polnimal order

        %MAKE SHALLOW TO DEEP
        [~,ia] = sort(P, 'ascend');
        P      = P(ia);
        Y      = Y(ia);

        tTOP = P < CpAct;
        tBOT = P > CpAct;

        newY       = Y;
        newY(tBOT) = newY(tBOT) + Offset;
        %newY(tTOP) = newY(tTOP) - Offset;

        %GET POLY FIT & SSR
        coefs = polyfit(P(tTOP|tBOT), newY(tTOP|tBOT), poly_order);
        Yest  = polyval(coefs, P(tTOP|tBOT));
        SSR   = sum((newY(tTOP|tBOT) - Yest).^2);

    end

end