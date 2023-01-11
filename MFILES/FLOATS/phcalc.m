function [phfree,phtot]= phcalc(Vrs, Press, Temp, Salt, k0, k2, Pcoefs)
%
% THIS FUNCTION WAS COPIED AND VERIFIED FROM THE ARGO PH PROCESSING DOCUMENT ON 1/15/19 
% Processing BGC-Argo pH data at the DAC level (https://doi.org/10.13155/57195)
%
% INPUTS:
%   Vrs     = Voltage bewteen reference electrode and ISFET source
%   Press   = Pressures in decibars
%   Temp    = Temperature in degrees C
%   Salt    = Salinity (usually CTD salinity on the PSS)
%   k0      = Sensor reference potential (intercept at Temp = 0C) 
%   k2      = constant or pressure dependant temperature response (slope)
%   Pcoefs  = sensor dependent pressure coefficients

% ************************************************************************
% TESTING
% % Vrs = -0.953799;
% % 
% % Press  = 4.71;
% % Press  = 1000;
% % Temp   = 16.3902;
% % Salt   = 33.6131;
% 
% Vrs    = [-0.953799;-0.953799]; 
% Press  = [4.71; 1000];
% Salt   = [33.6131; 33.6131];
% Temp   = [16.3902; 5.0];
% 
% k0     = -1.4131;
% k2     = -0.0011416;
% Pcoefs = [0 0 0 0 0 0]';
% 
% k2 = [-0.00100947,-5.85041e-08,3.83975e-11,-1.98141e-14];
% %k2 = [-0.00100947];
% ************************************************************************

% CHANGE HISTORY:
%02/03/2022 JP - Added compatibilty for using k2f(P)

% ************************************************************************
%  SET SOME CONSTANTS
% ************************************************************************
if max(size(k2)) > 1
    disp('K2 input has pressure dependent terms!')
end
%Universal gas constant, (R) , http://physics.nist.gov/cgi-bin/cuu/Value?r
R    = 8.31446; % J/(mol K) 
F    = 96485; %Faraday constant Coulomb / mol
Tk   = 273.15 + Temp; % degrees Kelvin
ln10 = log(10); % natural log of 10
    
% ************************************************************************
% CALCULATE PHYSICAL AND THERMODYNAMIC DATA
% Dickson, A. G., Sabine, C. L., & Christian, J. R. (2007). Guide to best
% practices for ocean CO2 measurements.
% ************************************************************************

% IONIC STRENGTH OF SEAWATER (mol / kg H2O)
% Varified units by comparing to Dickson et al. 2007: Chap 5, p10 Table 2
% Dickson et al. 2007: Chap 5, p13 Eq 34
IonS = 19.924 .* Salt ./ (1000 - 1.005 * Salt);

% MEAN SEAWATER SULFATE CONCENTRATION (mol / kg solution)
% This wants to be mol/kg-seawater  as KHSO4 is on that scale
% Dickson et al. 2007: Chap 5, p10 Table 2
Stotal = (0.14 / 96.062) .* (Salt / 1.80655);

% MEAN SEAWATER CHLORIDE CONCENTRATION  (mol / kg H20)
% this wants to be mol/kg H2O as activity is on mol/kg H2O scale
% Dickson et al. 2007: Chap 5, p10 Table 2
Cltotal = 0.99889 / 35.453 .* Salt / 1.80655; %(mol / kg solution)
Cltotal = Cltotal ./(1 - 0.001005 .* Salt);  % (mol / kg H20)

% BISULFIDE DISSCIATION CONSTANT AT T,S AND IONIC STRENGTH(mol/kg solution)
% Dickson et al. 2007: Chap 5, p12 Eq 33
Khso4 = exp(-4276.1 ./ Tk + 141.328 - 23.093 .* log(Tk) + ...
        (-13856 ./ Tk + 324.57 - 47.986 .* log(Tk)) .* IonS .^ 0.5 + ...
        (35474 ./ Tk - 771.54 + 114.723 .* log(Tk)) .* IonS - ...
        2698 ./ Tk .* IonS .^ 1.5 + 1776 ./ Tk .* IonS .^ 2 + ...
        log(1 - 0.001005 .* Salt));

% Millero 1983 Chemical Oceanography vol 8
% partial molar volume and compressibility of HSO4 in seawater. 
deltaVHSO4 = -18.03 + 0.0466 .* Temp + 0.000316 .* Temp .^ 2;
KappaHSO4 = (-4.53 + 0.09 .* Temp) / 1000;

%%%%%%%  Press changed from dbar to bar here by / 10
lnKhso4fac = (-deltaVHSO4 + 0.5 .* KappaHSO4 .* (Press / 10)) .* ...
             (Press / 10) ./ (R * 10 .* Tk);

%  bisulfate association constant at T, S, P
Khso4TPS = Khso4 .* exp(lnKhso4fac);

% GAMMA +/- HCl, activity coefficient of HCl at T/S, P=1
% ADH is the Debye Huckel constant, calcualted as a polynomial 
% fit to data in Khoo et al. 1977, doi:10.1021/ac50009a016
% See Martz et al. 2010, DOI 10.4319/lom.2010.8.172, p175
% Typo in paper 2nd term should be e-4 not e-6
%  
ADH = (3.4286e-6 .* Temp .^ 2 + 6.7524e-4 .* Temp + 0.49172143); 

log10gammaHCl = -ADH .* sqrt(IonS) ./ (1 + 1.394 .* sqrt(IonS)) + ...
                (0.08885 - 0.000111 .* Temp) .* IonS;
% Millero 1983 partial molar volume of HCl in seawater
deltaVHCl = 17.85 + 0.1044 .* Temp - 0.001316 .* Temp .^ 2;

% effect of pressure on activity coefficient of HCl, divide by 2 because
% its a mean activity coefficient, divide by 10 for units in the cm3 to F
% conversion.

log10gammaHCLtP = log10gammaHCl + deltaVHCl.*(Press./10)./(R.*Tk.*ln10)./2./10;

%  Sensor reference potential

% ************************************************************************
% CHECK SIZE OF K2 VARIABLE
if max(size(k2)) == 1
    k0T = k0 + k2 * Temp; % Temp  in deg C
elseif max(size(k2)) > 1 % polynomial pressure dependance
%     k2pc = [flipud(Pcoefs);0]; TM: This was applying f(P) Pcoefs instead of k2??
    k2pc = [flipud(k2)];
    k0T  = k0 + polyval(k2pc,Press) .* Temp;
else
    disp('Max size should be >= 1 : Check k2 input!')
    return
end

% CALCULATE PRESSURE CORRECTION (POLYNOMIAL FUNCTION OF PRESSURE)
% ALL SENSORS HAVE A PRESSURE RESPONSE WHICH IS DETERMINED IN THE LAB
% AND CONTAINED IN THE POLYNOMIAL Pcoefs
pc    = [flipud(Pcoefs);0]; % Matlab wants descending powers & n+1 (add 0)
pcorr = polyval(pc,Press);
k0TP  = k0T + pcorr;
  

%  pH on free scale then corrected to get to pH total on mol/kg-seawater scale
%    pHinsituFree = (Vrs - k0TP) / (R * Tk / F * ln10) + ...
%                   log(Cltotal) / ln10 + 2 * log10gammaHCLtP
%  this will be mol kg H2O  need to convert to mol/kg-seawater
phfree = (Vrs - k0TP) ./ (R .* Tk ./ F .* ln10) + ...
    log(Cltotal) ./ ln10 + 2 * log10gammaHCLtP; %mol/kg-H2O scale

% CONVERT TO mol/kg-seawater scale - JP 2/4/16
phfree = phfree - log10(1 - 0.001005 .* Salt); %mol/kg-seawater scale

% convert to total proton scale
phtot = phfree - log10(1 + Stotal ./ Khso4TPS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
