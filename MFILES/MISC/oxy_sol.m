function O2sol = oxy_sol(T,S,flag)
% [ O2sol ] = oxy_sol(T,S,flag) 
% T = Temperature in deg C
% S = Salinity
% flag: 1 = mmol/m^3, 0 = umol/kg
% if no flag  then output = mmol / m^3
% O2sol =  Oxygen solubility in mmol/m^3 (= umol/L = uM)
%
% CALCULATE SATURATION CONCENTRATION OF OXYGEN IN SEAWATER AS FUNCTION OF
% S & T IN EQUILIBRIUM WITH STANDARD COMPOSITION MOIST AIR AT 1 ATM TOTAL
% PRESSURE. FROM GARCIA & GORDON (1992) Eq 8 (p1310) USING COEFFICIENTS OF
% BENSON & KRAUSE IN TABLE 1 (p1311) AS USED IN SARMIENTO & GRUBER'S "OCEAN
% BIOGEOCHEMICAL DYNAMICS" CHAPTER 3, P81, TABLE 3.2.4 (OMMITTING ERRONEOUS
% TERM A3*Ts^2 FROM ORIGINAL PAPER). 

% J. Plant 12/15/10

% Check for proper inputs
    
if nargin ~= 3 || flag > 1 || flag <0
    s1 = 'IMPROPER INPUTS - Check flag value, S & T';
    disp(s1)
    return
    
elseif flag == 1  %  mmol/m^3 
   %Benson & Krause cm^3/dm^3 coefficients
   A = [3.88767 -0.256847 4.94457 4.05010 3.22014 2.00907];
   B = [-8.17083e-3 -1.03410e-2 -7.37614e-3 -6.24523e-3];
   C = -4.88682e-7;
   Ts = log((298.15-T) ./ (273.15 + T )); % Scaled temperature
   L = polyval(A,Ts) + S.*polyval(B,Ts) + C*S.^2;

   O2sol = (1000/22.3916) * exp(L); % Oxygen solubility, mmol/m^3
   
elseif flag ==0 % umol/kg  
   %Benson & Krause umol/kg coefficients
   A = [3.80369 -9.86643e-2 5.10006 4.17887 3.20291 5.80871];
   B = [-9.51519e-3 -1.13864e-2 -7.70028e-3 -7.01577e-3];
   C = -2.75915e-7;
   Ts = log((298.15-T) ./ (273.15 + T )); % Scaled temperature
   L = polyval(A,Ts) + S.*polyval(B,Ts) + C*S.^2;

   O2sol = exp(L); % Oxygen solubility, % umol/kg  
end
   




