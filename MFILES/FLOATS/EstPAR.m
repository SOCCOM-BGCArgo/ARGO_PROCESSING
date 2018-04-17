function [PAR,ZPAR20] =  EstPAR(lon, lat, sdn, chl, Depth, dirs)
% This function returns  the Par profile & the depth at PAR = 20 given the
% inputs. sdn is matlab serial date number in GMT, chla is raw float chlorophyll
% depth can be Depth or pressure it doesn't matter that much for the
% sueface waters

% TESTING[PAR, ZPAR20] = EstPAR(gps(2), gps(3), gps(1), CHLAm, P);
% lon = gps(2);
% lat = gps(3);
% sdn = gps(1);
% chl = CHLAm;
% Depth =  P;

% ************************************************************************
% SET SOME DIRS / CONSTANTS
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

    dirs.cal       = [user_dir,'CAL\'];
    dirs.temp      = 'C:\temp\';
end

PAR_range  =  [400 700];
%MM01_dir   = 'C:\Users\jplant\Documents\MATLAB\Xing\';
%MM01_dir   = 'E:\Documents\MATLAB\Xing\';

planck     = 6.6260755e-34; % 10e-34 kg m^2 s^-1
lightspeed = 299792458;  % m s^-1
Avogadro   = 6.02214086*10^23;

% % ************************************************************************
% % ************************************************************************
% % TESTING
% FV_dir = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\';
% %FV_dir = 'E:\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\';
% tmp    = get_FloatViz_data([FV_dir,'9254SOOCN.TXT']);
% d      = tmp.data;
% hdr    = tmp.hdr;
% clear tmp
% 
% iLAT  = find(strcmp('Lat [°N]',hdr)  == 1);
% iLON  = find(strcmp('Lon [°E]',hdr)  == 1);
% iSDN  = find(strcmp('SDN',hdr)  == 1);
% iP    = find(strcmp('Pressure[dbar]',hdr)  == 1);
% iT    = find(strcmp('Temperature[°C]',hdr)  == 1);
% iS    = find(strcmp('Salinity[pss]',hdr)  == 1);
% iSDN  = find(strcmp('SDN',hdr)  == 1);
% iSTA  = find(strcmp('Station',hdr)  == 1);
% iCHL  = find(strcmp('Chl_a[mg/m^3]',hdr)  == 1);
% iB    = find(strcmp('b_bp700[1/m]',hdr)  == 1);
% 
% t1 = d(:,iSTA) == 1;
% dd = d(t1,:);
% 
% sdn = d(1,iSDN);
% lat = d(1,iLAT);
% lon = d(1,iLON);
% 
% Depth   = dd(:,iP);
% % CHL SHOULD BE FILTERED FOR THE NPQ CORRECTION EVENTUALLY - MAYBE??
% chl = dd(:,iCHL); % RAW CHL, NO NPQ CORRECTION
% 
% % ************************************************************************
% % ************************************************************************


P = Depth;
CHL = chl * 0.5; % GLOBAL AVERAGE CALIBRATION GAIN CORRECTION
CHL(CHL<0) = 0; %QUICK FIX SO MMO1 doesn't produce imaginary #'s


% % CREATE A SURFACE VALUE FOR INTEGRATION & SAMPLES WITH MISSING SURFACE
% % VALUES (ICE DETECTION)
% P         = [0;P];
% CHL       = [CHL(1);CHL];

IXO = (1:max(size(P)))';

% t_order = 0;
% if P(end) < P(1) % deep to shallow
%     t_order = 1;
%     P   = flipud(P);
%     CHL = flipud(CHL);
% end



[~,IX] = sort(P);
P = P(IX);
IXO = IXO(IX);
CHL = CHL(IX);

clear dd t1 d hdr


% INTERP ANY MISSING CHL VALUES
tnan      = isnan(CHL);
if sum(tnan) > 0
    CHL(tnan) = interp1(P(~tnan), CHL(~tnan), P(tnan));
end

% ************************************************************************
% CALCULATE ATTENUATION COEFICIENTS FOR PAR WAVELENGTHS FOLLOWING
% Morel, A., and S. Maritorena (2001) doi:10.1029/2000JC000319

% Load spectral values from text file
% Data from Table 2, p7168 in Morel, A., and S. Maritorena (2001)
% doi:10.1029/2000JC000319
% fid  = fopen([MM01_dir, 'MM01.txt']); % Lambda	Kw	e	chi
% tmp = textscan(fid,'%f%f%f%f','Delimiter','\t', 'HeaderLines', 1, ...
%        'CollectOutput', 1);
% data = tmp{1,1};
% 
% fclose(fid);

load([dirs.cal,'MM01.mat']);
data = MM01_data';
clear MM01_data
t_WL   = data(1,:) >= PAR_range(1) & data(1,:) <= PAR_range(2);

lambda = data(1,t_WL); % nm
Kw     = data(2,t_WL); % m^-1
e      = data(3,t_WL);
chi    = data(4,t_WL);
dWL    = mode(diff(lambda)); % wavelength step nm
clear  tmp data fid

fill_C = ones(size(CHL));
fill_R = ones(size(lambda));

% MAKE Kd matrix pressure x wavelength
Kd = (fill_C*Kw) + (fill_C*chi) .* (CHL*fill_R) .^ (fill_C*e); % m^-1

% NOW CALCULATE DEPTH INTEGRATED Kd PROFILE MATRIX
% THIS ACCOUNTS FOR ACCUMLATIVE ATTENUATION WITH DEPTH
INT_Kd = cumsum(Kd .* ([P(1); diff(P)] * fill_R),1); % m/m unit-less

% ************************************************************************
% GET SPECTRAL IRRADIANCES FROMM GREGG & CARDER ADAPTATION OF THE
% BIRD & RIORDAN SPECTRAL IRRADIANCE MODEL

% SETTING METEOROLOGICAL VALUES TO SOME AVERAGE VALUES
met.psp     = 100;     % NOT USED! pyranometer readings (W m-2) 
met.slvp    = 1013.25; % sea-level pressure (mb)
met.dryt    = 0;       % dryt = dry air temperature (deg C)
met.relhum  = 80;      % relhum = relative humidity (%)
met.windspd = 7;       % windspd = wind speed (m s-1)
met.windavg = 5;       % windavg = 24hr mean wind speed (m s-1)

[wave dectime zenang IT ID IS] = radmod_jp(sdn, lat, lon, met);

% IT ID & IS ARE IN UNITS OF (J/(m^2 s um) = W / (m^2 µm))
% CONVERT TO OCEAN OPTICS COMMUNITY PREFERED UNITS (uW m-2 nm-1) (* 100)
% ALSO SLIGHT CORRECTION FOR HOW MUCH SURFACE SUNLIGHT PENETRATES AIR/SEA
% INTERFACE (* 0.98)
% wave * 1000 for um to nm

% Downwelling energy at surface for each PAR wavelength
% Interpolated to Kd MMO1 wave lengths
Ed0 = interp1(wave* 1000, IT*0.98*100, lambda); %(W m^-2 um^-1) Downwelling
% Ed0 = interp1(wave* 1000, IT, lambda); %(W m^-2 um^-1) at each nm Downwelling
% Ed0 = Ed0 * 0.98 / 1000; %(W m^-2 nm^-1)
% Ed0 = Ed0 * 10000; %(W cm^-2 nm^-1)

% PAR profiles for each wavelength
PAR = (fill_C * Ed0) .* exp(-INT_Kd); %(W m^-2 nm^-1)



% NOW INTEGRATE ACROSS WAVE LENGTHS
WL = fill_C * lambda;
PARxWL = PAR .* WL; %(W m^-2) at each lambda 

%INT_PAR = sum((PAR(:,1:end-1) + PAR(:,2:end)) ./ 2,2) * dWL;
INT_PAR = sum((PARxWL(:,1:end-1) + PARxWL(:,2:end)) ./ 2,2) * dWL/1000;

INT_PAR = INT_PAR / planck/ lightspeed/ 10^5 /Avogadro;
%INT_PAR = INT_PAR / planck/ lightspeed /10^4/Avogadro;

ind = find(INT_PAR > 20, 1, 'last');
if isempty(ind) % night time
    PAR20 = NaN;
else
%      PAR20 = P(ind);
    tnan = isnan(INT_PAR);
    t_diff = diff(INT_PAR) == 0; % duplicate value
    t_diff = logical([0; t_diff]);
    PAR20 = interp1(INT_PAR(~tnan & ~t_diff),P(~tnan & ~t_diff), 20);
end

P      = P(IXO);
CHL    = CHL(IXO);
PAR    = INT_PAR(IXO);
ZPAR20 = PAR20;




















