function [pO2, pH2O, O2_uM, O2_T] = Calc_SBE63_O2(data, cal, OptTvolt)
% This function is intended to calculate SBE63 oxygen concentration for 
% NAVIS profiling float message files. given the required raw data as a
% matrix and a calibration structure. Data is corrected for salt and
% pressure
%
% INPUTS:
%   data -- an n x 5 matix [CTD P, CTD T, CTD S, O2Phase, O2 Temp volts]
%   cal  -- a structure containing calibration coefficents
%           (I create it with get_float cals.m)
%               cal.TempCoef    [TA0 TA1 TA2 TA3
%               cal.PhaseCoef   [A0 A1 A2 B0 B1 C0 C1 C2]
%	OptTvolt -- Logical input to identify proper form of optode temperature within data(:,5).  1=volts; 0=degC.
%
% OUTPUT
%   [O2_uM O2_T] -- Oxygen concentration in µmol /L, O2 sensor Temperature
%
% HISTORY
%   05/30/2017 - Implemented Henry Bittig's (2015) pressure correction
%       scheme. This makes  a correction in raw phase and in [O2].
%       http://dx.doi.org/10.1175/JTECH-D-15-0108.1
%  5/25/21  TM - added final exponential term to the conversion to PPOX_DOXY which accounts for
%                 the change in oxy-solubility with hydrostatic pressure.  
%				  Was previously missing (very small impact).  See Argo O2 processing doc, Section 1.3.2. 
%  06/02/21 TM - minor mod to O2 temperature handling for the SBE83 sensors
%  06/03/21 TM - added pO2 as fx output
%  09/29/21 TM - implemented a catch for large opt temp voltage (resulting
%                in complex number for optode temp and concentration)
%  12/15/21 JP&TM - bad data returned from sensor generates imaginary
%  (un1114 c34 for example) numbers for [O2]. This occurrs during optode
%   window T calc in volts so check valid volt range (0.001 to 3.3) & if
%   out of range set to crazy value 3.3 (-273.15)
%  3/10/22 TM - Modification to Temp used for calculation pres, sal compensation.  
%               Keep Optode T for use in calc O2 from phase delay.
%  3/29/22 TM, increased the lower RC bound on Tvolt to 0.002 (still
%  getting imaginary values from 1200.032.msg file for bad optode temp at
%  0.0015....
%  01/18/24, TM, added case for new Navis Nautilus (0063 and future) which have SBS83 optode temperature expressed as voltage
%  02/08/24, TM, Further modified this code in support of new Navis Nautilus (0063 and future) which have SBS83 optode temperature expressed as voltage in the PROFILE data only!  Optode temp is in degC for in-air....


% TESTING
% data = lr_d(:,[iP,iT,iS,iPhase,iTo]);
% Tcf = cal.O.TempCoef; % TA0-3, Temperature coefficients
% Acf = cal.O.PhaseCoef(1:3); % A0, A1, A2 Phase Coef's
% Bcf = cal.O.PhaseCoef(4:5); % B0, B1 Phase Coef's
% Ccf = cal.O.PhaseCoef(6:8); % C0, C1, C2, Ksv Coef's
% ************************************************************************
% BREAK OUT COEFFICENTS AND DATA, SET SOME CONSTANTS
Tcf = cal.TempCoef; % TA0-3, Temperature coefficients
Acf = cal.PhaseCoef(1:3); % A0, A1, A2 Phase Coef's
Bcf = cal.PhaseCoef(4:5); % B0, B1 Phase Coef's
Ccf = cal.PhaseCoef(6:8); % C0, C1, C2, Ksv Coef's

P = data(:,1); % CTD pressure
T = data(:,2); % CTD temperature
S = data(:,3); % CTD salinity
OPh = data(:,4); % SBE63 phase delay
OTV = data(:,5); % SBE63 temperature in volts

Sref = 0; %instrument usually set to 0 salinity
%E    = 0.011; % Pressure correction coefficient
E    = 0.009; % Reassesed value in Processing ARGO O2 V2.2
%clear data cal

% ************************************************************************
% CALCULATE OPTODE WINDOW TEMPERATURE (code from Dan Quittman)
% Volts to resistance, % "log" is natural log in Matlab
% Note that the SBE83 returns temperature in C, not Volts! TM 6/2/21
% THIS IS NOT ALWAYS TRUE !!! THE 3XO2 SBE83's report Volts!!!! jp 09/03/2021
tf_3X02 = strcmp(cal.type,'SBE83') && ~isempty(regexp(cal.SN,'^000[467]','once'));

% if ~strcmp(cal.type,'SBE83') || tf_3X02
if contains(cal.type,'SBE63') || tf_3X02 || (OptTvolt)
    % Check for valid voltage if out of range set to crazy value
    tbad = (OTV <0.002 | OTV > 3.3) & ~isnan(OTV); % 3.3 = -273.15K
    OTV(tbad) = 3.3; % set bad to crazy value
	L      = log(100000 * OTV  ./ (3.3 - OTV ));
	denom  = (Tcf(1) + Tcf(2).*L + Tcf(3).*L.^2 + Tcf(4).*L.^3);
	O2_T   = 1./denom - 273.15; % window temperature (degrees C)
else
	O2_T = OTV;
end

% ************************************************************************
% CALCULATE OXYGEN FROM PHASE, ml/L (modified from Dan Quittman)
CALC_T = O2_T; % Could also use CTD T, but SBE63 is pumped

% Bittig pressure correction (2015) Part 1, T & O2 independant 05/30/2017
OPh    = OPh + (0.115 * P / 1000);

V      = OPh ./ 39.4570707;
A      = Acf(1) + Acf(2).*CALC_T + Acf(3).*V.^2;
B      = Bcf(1) + Bcf(2).*V;
Ksv    = Ccf(1) + Ccf(2).*CALC_T + Ccf(3).*CALC_T.^2;
O2     = (A./B -1) ./ Ksv; % ml/L, Salinity = 0, pressure = 0

% ************************************************************************
% CALCULATE SALT CORRECTION
TK     = T + 273.15; %these should be temp from CTD
Ts    = log((298.15 - T) ./ TK); % Scaled T
% Descending powers for Matlab polyfit
SolB  = [-8.17083e-3 -1.03410e-2 -7.37614e-3 -6.24523e-3];
C0    = -4.88682e-7;

Scorr = exp((S-Sref).*polyval(SolB,Ts) + C0.*(S.^2-Sref.^2));

% ************************************************************************
% CALCULATE PRESSURE CORRECTION
%PP =(P>0).*P; % Clamp negative values to zero
%pcorr = exp(E * PP ./ (CALC_T + 273.15));

% SWITCH TO PRESSURE CORRECTION FROM BITTIG ET AL., (2015) - 05/30/2017
%pcorr = (0.00022*T + 0.0419) .* P/1000 + 1;
pcorr = (0.00022*T + 0.0419) .* P/1000 + 1; %pressure correction should always use CTD temp...!
% ************************************************************************
% FINAL OUTPUT
O2_uM = O2 .* Scorr .* pcorr .* 44.6596; % OXYGEN µM/L

% NOW COMPUTE pO2
% Benson & Krause (cm^3/dm^3) *** USE THESE COEFFICIENTS (for all floats except Aanderaa 4330 exceptions)!!!! ***
pA = [3.88767 -0.256847 4.94457 4.05010 3.22014 2.00907]; % Temperature Coeff
pB = [-8.17083e-3 -1.03410e-2 -7.37614e-3 -6.24523e-3]; % Salinity coff
Co = -4.88682e-7;
O2_volume =  22.3916;
part1  = polyval(pB,Ts);
S_corr = S.*part1 + S.^2 * Co;
L      = polyval(pA,Ts) + S_corr;
O2sol  = (1000/O2_volume) * exp(L); % Oxygen solubility real gas, mmol/m^3 =uM/L
theExp = (0.317.*P)./(8.314.*TK); % TM 5/25/21; in accordance with Argo O2 cookbook (for aircal pO2, P=0, this exp term does nothing)
% CALCULATE VAPOR PRESSURE H2O
pH2O = exp(52.57 -(6690.9./TK) - 4.681 .* log(TK));
pO2 = (O2_uM ./ O2sol) .* ((1013.25- pH2O) * 0.20946) .* exp(theExp);  

%Perform final check for imaginary #'s. If found set to crazy value
% O2_T(imag(O2_T) ~= 0) = NaN; %Navis SBE63 1114 on cycle 34 returning bunk values for phase and volt, resulting in complex values after voltage conversion.  :-/ TM, 9/29/21
% O2_uM(imag(O2_uM) ~= 0) = NaN;
% pO2(imag(pO2) ~= 0) = NaN;
