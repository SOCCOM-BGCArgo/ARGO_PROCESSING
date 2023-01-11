function O2data = calc_O2_4ARGO(S, T, P, Phase, O2_T, cal)
% USE: Calculate oxygen concentration from phase angle for fluorescent
% quenching optodes ( Aanderaa 3038 & 4330, SBE63). Deals with the
% different flavors of calibration based on info found in the calibration
% file. The calibration file is created with get_float_cals.m
% 
%      Last updated 04/29/19 by Tanya Maurer
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
%     O2_T      = Temperature returned from the optode
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
% 10/25/2018 Added code to calculate O2 using the SVUCoef's for 4330's that
%       have these coefficients
% 04/24/19 - Reverted back to the Garcia & Gordon (1992) combined fit coefficients for
%       calculation of oxygen solubility for 4330 optodes with polynomial
%       coefficients.  Although this doesn't conform to the SCOR WG 142
%       recommendation, it is more accurate for these cases because more
%       closely matches what Aanderaa used in derivation of cal coeffs. (And this is
%       also reflected in the Argo O2 processing doc, and verified from Henry Bittig).
% 03/10/22  - Modifications to Temperature used for Stern Volmer calculation (use Optode T)

% ************************************************************************
% SET VARIABLES & PATHS & FORMATS
% ************************************************************************

% OXYGEN SOLUBILITY COFFICIENTS
%
% NOTE: Aanderaa uses the combined fit coefficients for the salinity terms
%       of the Co* calculation from Garcia & Gordon (1992) Equation 8 for 4330 optodes.
%       Garcia & Gordon suggest using the more precise Benson & Krause
%       coefficients but need to stay be consistent with how the cal coeffs
%       were originally derived from the manufacturer.

% Benson & Krause (cm^3/dm^3) *** USE THESE COEFFICIENTS (for all floats except Aanderaa 4330 exceptions)!!!! ***
pA = [3.88767 -0.256847 4.94457 4.05010 3.22014 2.00907]; % Temperature Coeff
pB = [-8.17083e-3 -1.03410e-2 -7.37614e-3 -6.24523e-3]; % Salinity coff
Co = -4.88682e-7;
O2_volume =  22.3916;

% ************************************************************************
% %Aanderaa default coeff (cm^3/dm^3) (Combined fit Garcia & Gordon 1992)
pA_Aa = [1.71069 0.978188 4.80299 3.99063 3.22400 2.00856]; % Temperature Coeff
pB_Aa = [-4.29155e-3 -6.90358e-3 -6.93498e-3 -6.24097e-3]; % Salinity coff
Co_Aa = -3.11680e-7;
%O2_volume_Aa =  22.414; % Volume ideal gas
% ************************************************************************

TK     = T + 273.15;
Ts     = log((298.15-T)./ TK);

% Used for final O2 conversions from [O2] to pO2:
part1  = polyval(pB,Ts);
S_corr = S.*part1 + S.^2 * Co;
L      = polyval(pA,Ts) + S_corr;
O2sol  = (1000/O2_volume) * exp(L); % Oxygen solubility real gas, mmol/m^3 =uM/L

% Used for 4330s with polynomial coeffs for getting from pO2 (essentially, output from sensor after cal coeffs applied) to [O2], prior to application of conc coeffs:
part1_Aa  = polyval(pB_Aa,Ts);
S_corr_Aa = S.*part1_Aa + S.^2 * Co_Aa;
L_Aa      = polyval(pA_Aa,Ts) + S_corr_Aa;
O2sol_Aa  = 44.614 * exp(L_Aa); % Oxygen solubility ideal gas, mmol/m^3 =uM/L.  44.614 is the conversion factor listed in Argo O2 processing cookbook


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
    pO2 = O2 ./ O2sol  .* (1013.25 - pH2O) .* 0.20946 ; % partial pressure O2

% ************************************************************************
% AANDERAA 4330 OPTODE
% ************************************************************************
elseif strcmp(cal.type,'4330')
    % Bittig pressure correction (2015) Part 1, T & O2 independant
    Phase = Phase + (0.1 * P / 1000);
  
    % CALCULATE CALIBRATED PHASE
    ind = length(cal.PCoef):-1:1; % need to reverse order of coef for Matlab
    CalPhase = polyval(cal.PCoef(ind),Phase);
    
    if ~isfield(cal,'SVUFoilCoef') % 4330 polynomial returns pO2
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
        O2 = (O2sol_Aa .* O2Sat) / 100;
        
        
    elseif isfield(cal,'SVUFoilCoef') % 4330 Stern Volmer (Uchida) returns [O2]
		T = O2_T; %per cookbook!!  For 4330 polynomial - use CTD T for all calcs.
        SVU = cal.SVUFoilCoef;
        Ksv = SVU(1) + SVU(2).*T + SVU(3).*T.^2;
        P0  = SVU(4) + SVU(5).*T;
        PC  = SVU(6) + SVU(7).*CalPhase;
        O2  = (P0./PC -1) ./ Ksv;
        O2  = O2 .* exp(S_corr); % Aanderaa S correction to O2 concentration
        pO2 = O2 ./ O2sol  .* (1013.25 - pH2O) .* 0.20946 ; % partial pressure O2
    end
    
    % IF ConcCoef variable exist need to update pO2 too
	% Use B&K solubility coeffs for final conversion from [O2] back to pO2 for all models, per Henry
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










