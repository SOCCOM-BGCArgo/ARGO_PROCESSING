function [phfree,phtot]= phcalc_jp(Vrs, Press, Temp, Salt, k0, k2, Pcoefs)
%%%%%%%%%  this version is modified per Yui comments in two places
%%%%%%%%  lnKhso4fac - press/10 to get to bars
%%%%%%%%   modified by KJ 11/17/2015 to get units of mol/kg H2O and mol/kg sw correct
%%%%%%%%
%
% INPUTS:
%   Vrs     = Voltage bewteen reference electrode and ISFET source
%   Press   = Pressures in decibars
%   Temp    = Temperature in degrees C
%   Salt    = Salinity (usually CTD salinity on the PSS)
%   k0      = Sensor reference potential (intercept at Temp = 0C) paper says K??
%   k2      = linear temperature coefficient (slope)
%   Pcoefs  = sensor dependent pressure coefficients

% FOR TESTING ONLY !!!!
% load('C:\Users\jplant\Documents\MATLAB\ARGO\DATA\Cal_files\cal9313.mat');
% d = parse_APEXmsg4ARGO('\\atlas\chemwebdata\floats\f9313\9313.001.msg');
% %d = parse_APEXmsg4ARGO('\\atlas\chemwebdata\floats\f9313\9313.055.msg');
% LR = d.lr_d;
% 
% Vrs    = LR(:,8);
% Press  = LR(:,1);
% Temp   = LR(:,2);
% Salt   = LR(:,3);
% k0     = cal.pH.k0;
% k2     = cal.pH.k2;
% Pcoefs = cal.pH.pcoefs;

% ************************************************************************
%  SET SOME CONSTANTS
% ************************************************************************
%Universal gas constant, (R) , http://physics.nist.gov/cgi-bin/cuu/Value?r
%R    = 8.31451;
R    = 8.31446; % J/(mol K) jp 
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
% This wants to be mol/kg sw  as KHSO4 is on that scale
% Dickson et al. 2007: Chap 5, p10 Table 2
% Edited by KJ to correct units 11/17/2015
Stotal = (0.14 / 96.062) .* (Salt / 1.80655);

% MEAN SEAWATER CHLORIDE CONCENTRATION  (mol / kg H20)
% this wants to be mol/kg H2O as activity is on mol/kg H2O scale
% Dickson et al. 2007: Chap 5, p10 Table 2
Cltotal = 0.99889 / 35.453 .* Salt / 1.80655; %(mol / kg solution)
Cltotal = Cltotal ./(1 - 0.001005 .* Salt);  % (mol / kg H20)
% Where does the  (1 - xxx/S) come form?

% BISULFIDE DISSCIATION CONSTANT AT T,S AND IONIC STRENGTH(mol/kg solution)
% Dickson et al. 2007: Chap 5, p12 Eq 33
Khso4 = exp(-4276.1 ./ Tk + 141.328 - 23.093 .* log(Tk) + ...
        (-13856 ./ Tk + 324.57 - 47.986 .* log(Tk)) .* IonS .^ 0.5 + ...
        (35474 ./ Tk - 771.54 + 114.723 .* log(Tk)) .* IonS - ...
        2698 ./ Tk .* IonS .^ 1.5 + 1776 ./ Tk .* IonS .^ 2 + ...
        log(1 - 0.001005 .* Salt));

% WHERE DO THESE APROXIMATIONS COME FROM???
% Millero 1983 Chemical Oceanography vol 8
deltaVHSO4 = -18.03 + 0.0466 .* Temp + 0.000316 .* Temp .^ 2;
KappaHSO4 = (-4.53 + 0.09 .* Temp) / 1000;
%%%%%%%  per Yui Press changed from dbar to bar here by / 10
lnKhso4fac = (-deltaVHSO4 + 0.5 .* KappaHSO4 .* (Press / 10)) .* ...
             (Press / 10) ./ (R * 10 .* Tk);
%  bisulfate association constant at T, S, P
Khso4TPS = Khso4 .* exp(lnKhso4fac);


% GAMMA +/- HCl AT T AND S
% Polynomial fit to Khoo et al. 1977, doi:10.1021/ac50009a016
% See Matrz et al. 2010, DOI 10.4319/lom.2010.8.172, p175
% Typo in paper 2nd term should be e-4 not e-6
%  Debye Huckel constant A
ADH = (3.4286e-6 .* Temp .^ 2 + 6.7524e-4 .* Temp + 0.49172143); % jp
%ADH = (0.00000343 .* Temp .^ 2 + 0.00067524 .* Temp + 0.49172143);

log10gammaHCl = -ADH .* sqrt(IonS) ./ (1 + 1.394 .* sqrt(IonS)) + ...
                (0.08885 - 0.000111 .* Temp) .* IonS;
% Millero
deltaVHcl = 17.85 + 0.1044 .* Temp - 0.001316 .* Temp .^ 2;

ThermoPress = -deltaVHcl .* 0.0242 ./ (23061 * 1.01) .* Press ./ 10;
%%%%%%%%%%%%% per Yui comment original line modified so ThermoPress
% (in units of volts is added to E0 not to log10gammaHCL
%   log10gammaHCLtP = log10gammaHCl + ThermoPress
log10gammaHCLtP = log10gammaHCl;

%  Sensor constant

%    Debug.Print "loggamma  ", log10gammaHCl, "  ", log10gammaHCLtP
 
% ************************************************************************
E0T = k0 + k2 * Temp; % Temp  in deg C

% CALCULATE PRESSURE CORRECTION (POLYNOMIAL FUNCTION OF PRESSURE)
% ALL SENSORS HAVE A PRESSURE RESPONSE WHICH IS DETERMINED IN TH LAB
pc    = [flipud(Pcoefs);0]; % Matlab wants descending powers & n+1 (add 0)
pcorr = polyval(pc,Press);
E0TP  = E0T + pcorr;

  
%%%%%%%%%%% Per Yui Comment, ThermoPress added into Vrs -E0TP term
%%%%%%%%%%%%  corrected by Ken to get to pH total on mol/kg sw scale
%    pHinsituFree = (Vrs - E0TP) / (R * Tk / F * ln10) + ...
%                   log(Cltotal) / ln10 + 2 * log10gammaHCLtP
%  this will be mol kg H2O  need to convert to mol/kg sw
phfree = (Vrs - E0TP - ThermoPress) ./ (R .* Tk ./ F * ln10) + ...
    log(Cltotal) ./ ln10 + 2 * log10gammaHCLtP; %mol/kg-H2O scale

% CONVERT TO mol/kg-sw scale - JP 2/4/16
phfree = phfree - log10(1 - 0.001005 .* Salt); %mol/kg-sw scale

% Hfree = 10.^(-phfree);           %   hfree mol/kg H2O
% Hfree = Hfree.*(1 - 0.001005 .* Salt);     %  hfree mol/kg SW
% phfree = -log10(Hfree);

phtot = phfree - log10(1 + Stotal ./ Khso4TPS);


