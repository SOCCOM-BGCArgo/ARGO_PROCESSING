function [betasw,beta90sw,bsw]= betasw_ZHH2009(lambda,Tc,theta,S,delta)
% Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
% seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710 
%
% INPUTS:
%   lambda = wavelength (nm) (700 for  FLBB)
%   Tc     = temperature in degree Celsius, must be a scalar
%   theta  = angle scattering measured at ( DAC = 124(?), FLBBAP2 = 142)
%   S      = salinity, must be scalar
%   delta  = depolarization ratio 
%            if not provided, default = 0.039 will be used.
%
% OUPUTS:
%   betasw   = volume scattering at angles defined by theta. 
%              Its size is [x y],
%              x is the number of angles (x = length(theta))
%              y is the number of wavelengths in lambda (y = length(lambda))
%   beta90sw = volume scattering at 90 degree. Its size is [1 y]
%	bw       = total scattering coefficient. Its size is [1 y]
%
%   FOR BACKSCATTERING COEFICIENTS, DIVIDE TOTAL SCATTERING BY 2

%
% Xiaodong Zhang, March 10, 2009
% 

% values of the constants
Na = 6.0221417930e23 ;   %  Avogadro's constant
Kbz = 1.3806503e-23 ;    %  Boltzmann constant
Tk = Tc+273.15 ;         %  Absolute tempearture
M0 = 18e-3;              %  Molecular weigth of water in kg/mol

error(nargchk(4, 5, nargin));
if nargin == 4
    delta = 0.039; % Farinato and Roswell (1976)
end

if ~isscalar(Tc) || ~isscalar (S)
    error('Both Tc and S need to be scalar variable');
end

lambda = lambda(:)'; % a row variable
rad = theta(:)*pi/180; % angle in radian as a colum variable

% nsw: absolute refractive index of seawater
% dnds: partial derivative of seawater refractive index w.r.t. salinity
[nsw dnds] = RInw(lambda,Tc,S);
 
% isothermal compressibility is from Lepple & Millero (1971,Deep
% Sea-Research), pages 10-11
% The error ~ +/-0.004e-6 bar^-1
IsoComp = BetaT(Tc,S);

% density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
density_sw = rhou_sw(Tc, S);
 
% water activity data of seawater is from Millero and Leung (1976,American
% Journal of Science,276,1035-1077). Table 19 was reproduced using
% Eq.(14,22,23,88,107) then were fitted to polynominal equation.
% dlnawds is partial derivative of natural logarithm of water activity
% w.r.t.salinity
dlnawds = dlnasw_ds(Tc, S);

% density derivative of refractive index from PMH model
DFRI = PMH(nsw);  %% PMH model
 
% volume scattering at 90 degree due to the density fluctuation
beta_df = pi*pi/2*((lambda*1e-9).^(-4))*Kbz*Tk*IsoComp.*DFRI.^2*(6+6*delta)/(6-7*delta);
% volume scattering at 90 degree due to the concentration fluctuation
flu_con = S*M0*dnds.^2/density_sw/(-dlnawds)/Na;
beta_cf = 2*pi*pi*((lambda*1e-9).^(-4)).*nsw.^2.*(flu_con)*(6+6*delta)/(6-7*delta);
% total volume scattering at 90 degree
beta90sw = beta_df+beta_cf;
bsw=8*pi/3*beta90sw*(2+delta)/(1+delta);
for i=1:length(lambda)
    betasw(:,i)=beta90sw(i)*(1+((cos(rad)).^2).*(1-delta)/(1+delta));
end


% ************************************************************************* 
% *************************************************************************
function [nsw dnswds]= RInw(lambda,Tc,S)
% refractive index of air is from Ciddor (1996,Applied Optics)
n_air = 1.0+(5792105.0./(238.0185-1./(lambda/1e3).^2)+167917.0./(57.362-1./(lambda/1e3).^2))/1e8;

% refractive index of seawater is from Quan and Fry (1994, Applied Optics)
n0 = 1.31405; n1 = 1.779e-4 ; n2 = -1.05e-6 ; n3 = 1.6e-8 ; n4 = -2.02e-6 ;
n5 = 15.868; n6 = 0.01155;  n7 = -0.00423;  n8 = -4382 ; n9 = 1.1455e6;

nsw = n0+(n1+n2*Tc+n3*Tc^2)*S+n4*Tc^2+(n5+n6*S+n7*Tc)./lambda+n8./lambda.^2+n9./lambda.^3; % pure seawater
nsw = nsw.*n_air;
dnswds = (n1+n2*Tc+n3*Tc^2+n6./lambda).*n_air;

function IsoComp = BetaT(Tc, S)
% pure water secant bulk Millero (1980, Deep-sea Research)
kw = 19652.21+148.4206*Tc-2.327105*Tc.^2+1.360477e-2*Tc.^3-5.155288e-5*Tc.^4;
Btw_cal = 1./kw;

% isothermal compressibility from Kell sound measurement in pure water
% Btw = (50.88630+0.717582*Tc+0.7819867e-3*Tc.^2+31.62214e-6*Tc.^3-0.1323594e-6*Tc.^4+0.634575e-9*Tc.^5)./(1+21.65928e-3*Tc)*1e-6;

% seawater secant bulk
a0 = 54.6746-0.603459*Tc+1.09987e-2*Tc.^2-6.167e-5*Tc.^3;
b0 = 7.944e-2+1.6483e-2*Tc-5.3009e-4*Tc.^2;

Ks =kw + a0*S + b0*S.^1.5;

% calculate seawater isothermal compressibility from the secant bulk
IsoComp = 1./Ks*1e-5; % unit is pa

% *************************************************************************
% *************************************************************************
function density_sw = rhou_sw(Tc, S)

% density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
a0 = 8.24493e-1;  a1 = -4.0899e-3; a2 = 7.6438e-5; a3 = -8.2467e-7; a4 = 5.3875e-9;
a5 = -5.72466e-3; a6 = 1.0227e-4;  a7 = -1.6546e-6; a8 = 4.8314e-4;
b0 = 999.842594; b1 = 6.793952e-2; b2 = -9.09529e-3; b3 = 1.001685e-4;
b4 = -1.120083e-6; b5 = 6.536332e-9;
 
% density for pure water 
density_w = b0+b1*Tc+b2*Tc^2+b3*Tc^3+b4*Tc^4+b5*Tc^5;
% density for pure seawater
density_sw = density_w +((a0+a1*Tc+a2*Tc^2+a3*Tc^3+a4*Tc^4)*S+ ...
             (a5+a6*Tc+a7*Tc^2)*S.^1.5+a8*S.^2);

% *************************************************************************
% *************************************************************************
function dlnawds = dlnasw_ds(Tc, S)
% water activity data of seawater is from Millero and Leung (1976,American
% Journal of Science,276,1035-1077). Table 19 was reproduced using
% Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
% dlnawds is partial derivative of natural logarithm of water activity
% w.r.t.salinity
% lnaw = (-1.64555e-6-1.34779e-7*Tc+1.85392e-9*Tc.^2-1.40702e-11*Tc.^3)+......
%            (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3).*S+......
%            (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^1.5+......
%            (-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S.^2;

dlnawds = (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3) +...
           1.5*(1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^0.5+......
           2*(-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S;

% *************************************************************************
% *************************************************************************
% density derivative of refractive index from PMH model
function n_density_derivative=PMH(n_wat)
n_wat2 = n_wat.^2;
n_density_derivative=(n_wat2-1).*(1+2/3*(n_wat2+2).*(n_wat/3-1/3./n_wat).^2);
