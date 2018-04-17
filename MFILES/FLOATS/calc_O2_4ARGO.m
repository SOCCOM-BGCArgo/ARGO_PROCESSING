function O2data = calc_O2_4ARGO(S, T, P, Phase, cal)
% USE: Calculate oxygen concentration from phase angle for fluorescent
% quenching optodes ( Aanderaa 3038 & 4330, SBE63). Deals with the
% different flavors of calibration based on info found in the calibration
% file. The calibration file is created with get_float_cals.m
% 
%      Last updated 12/29/15 by Josh Plant
%
% O2data = calc_O2_4ARGO(S, T,P Phase, cal)
%
% INPUTS:       arrays except for cal structure - no matrices
%     S         = Salinity (pss) 
%     T         = Temperature (C) (usually CTD T because of faster response
%                 but check CTD T vs OPT T for 1:1 relationship)
%     P         = Pressure (mbar)
%
%     Phase     = Phase angle from Optode (degrees)
%     cal       = a structure with calibration data. (Created using
%                 get_float_cals.m)
%
% OUTPUT:
%     O2data    =  [S, T, P, Phase, CalPhase, [O2], O2sol, pO2] 
%               units:
%
% EXAMPLE:
%     O2data = calc_O2_4ARGO(S, T, P, Phase, cal)

% HISTORY
%   04/04/2017 - Added code to use "ConcCoef" calibration variables to do 
%       a final correction to concentration  and pO2 in newer 4330 sensors
%   05/30/2017 - Implemented Henry Bittig's (2015) pressure correction
%       scheme. This makes  a correction in raw phase and in [O2].
%       http://dx.doi.org/10.1175/JTECH-D-15-0108.1

% ************************************************************************
% SET VARIABLES & PATHS & FORMATS
% ************************************************************************

% OXYGEN SOLUBILITY COFFICIENTS
%
% NOTE: Aanderaa uses the combined fit coefficients for the salinity terms
%       of the Co* calculation from Garcia & Gordon (1992) Equation 8.
%       Garcia & Gordon suggest using the more precise Benson & Krause
%       coefficients which is what I do here. Very minor difference!!!

% Benson & Krause (cm^3/dm^3) *** USE THESE COEFFICIENTS!!!! ***
pA = [3.88767 -0.256847 4.94457 4.05010 3.22014 2.00907]; % Temperature Coeff
pB = [-8.17083e-3 -1.03410e-2 -7.37614e-3 -6.24523e-3]; % Salinity coff
Co = -4.88682e-7;
O2_volume =  22.3916;

% ************************************************************************
% %Aanderaa default coeff (cm^3/dm^3) (Combined fit Garcia & Gordon 1992)
% pA = [1.71069 0.978188 4.80299 3.99063 3.22400 2.00856]; % Temperature Coeff
% pB = [-4.29155e-3 -6.90358e-3 -6.93498e-3 -6.24097e-3]; % Salinity coff
% Co = -3.11680e-7;
% O2_volume =  22.414; % Volume ideal gas
% ************************************************************************

TK     = T + 273.15;
Ts     = log((298.15-T)./ TK);
part1  = polyval(pB,Ts);
S_corr = S.*part1 + S.^2 * Co;
L      = polyval(pA,Ts) + S_corr;
O2sol  = (1000/O2_volume) * exp(L); % Oxygen solubility real gas, mmol/m^3 =uM/L
%O2sol  = 44.614 * exp(L); % Oxygen solubility ideal gas, mmol/m^3 =uM/L

% CALCULATE VAPOR PRESSURE H2O
pH2O = exp(52.57 -(6690.9./TK) - 4.681 .* log(TK));

% ************************************************************************
% AANDERAA 3830 OPTODE
% ************************************************************************
if strcmp(cal.type,'3830')  
    % CALCULATE OXYGEN CONCENTRATION FOR S=0
    % Calculate calibrated phase (Dphase) from Bphase
    % NOTE: this is Aanderraa's 2 pt cal correction
    
    % Bittig pressure correction (2015) Part 1, T & O2 independant
    Phase = Phase + (0.1 * P / 1000);
    
    % Aanderaa Calibrated phase (Dphase)
    CalPhase = polyval(cal.pP,Phase); % Calibrated phase (Dphase)
    
    % Calculate temperature dependent coefficients
    C0 = polyval(cal.pC0, T);
    C1 = polyval(cal.pC1, T);
    C2 = polyval(cal.pC2, T);
    C3 = polyval(cal.pC3, T);
    C4 = polyval(cal.pC4, T);
    
    pC = [C4 C3 C2 C1 C0]; % matrix w rows of coefficients, descending power
    
    [r,~] = size(Phase);
    O2 = ones(r,1)*NaN;
    for i = 1:r
        O2(i) = polyval(pC(i,:),CalPhase(i)); % umol/L, S=0
    end
    
    O2 = O2 .* exp(S_corr); % Aanderaa S correction to O2 concentration
    
   % BACK OUT pO2 from O2 CONC. 
    pO2 = O2 ./ O2sol  .* (1013.25 - pH2O) .* 0.20964 ; % partial pressure O2

% ************************************************************************
% AANDERAA 4330 OPTODE
% ************************************************************************
elseif strcmp(cal.type,'4330')
    % Bittig pressure correction (2015) Part 1, T & O2 independant
    Phase = Phase + (0.1 * P / 1000);
  
    % CALCULATE CALIBRATED PHASE
    ind = length(cal.PCoef):-1:1; % need to reverse order of coef for Matlab
    CalPhase = polyval(cal.PCoef(ind),Phase);
    
    % CALCULATE PARTIAL PRESSURE O2
    pO2   = ones(size(CalPhase))*NaN; % predim
    for i = 1:length(CalPhase)
        tmp = cal.FCoef .* (T(i).^cal.PolyDegT) .* ...
              (CalPhase(i) .^ cal.PolyDegO);
        pO2(i) = sum(tmp);
    end

    % CALCULATE O2 SATURATION (Assume a constant atm pressure)
    O2Sat = pO2 ./ ((1013.25- pH2O) * 0.20946) * 100;
    
    % CALCULATE O2 CONC FOR S = XX
    %O2 = (O2sol .* O2Sat * 44.614) / 100;
    O2 = (O2sol .* O2Sat) / 100;
    
    % IF ConcCoef variable exist need to update pO2 too
    if isfield(cal,'ConcCoef')
        O2 = cal.ConcCoef(1) + O2 * cal.ConcCoef(2);
        pO2 = (O2 ./ O2sol) .* ((1013.25- pH2O) * 0.20946);
    end
    
% ************************************************************************
% SBE 63 OPTODE
% ************************************************************************    
elseif strcmp(cal.type,'SBE63')
    disp('Need to add SBE63 code here')
    return
    
% ************************************************************************
% UNKNOWN OXYGEN SENSOR
% ************************************************************************     
else
    disp('Can not determine Oxygen optode type')
    return
end

% ************************************************************************
% ************************************************************************
% DEPTH COMPENSATION - 3.2% is an emperical value from Uchida 2008
%O2 = O2 .* ((P * 0.032/1000) + 1);

% SWITCH TO PRESSURE CORRECTION FROM BITTIG ET AL., (2015) - 05/30/2017
O2 = O2 .* ((0.00022*T + 0.0419) .* P/1000 +1);

O2data = [S, T, P, Phase, CalPhase, O2, O2sol, pO2]; % umol /L & % sat










