function PAR =  EstPAR(lon, lat, sdn, chl, Depth, dirs)
% This function returns  the Par profile given the inputs. sdn is the 
% matlab serial date number in GMT, chla is raw float chlorophyll
% depth can be Depth or pressure it doesn't matter that much for the
% surface waters

% TESTING PAR = EstPAR(gps(2), gps(3), gps(1), CHLAm, P);
% lon = gps(2);
% lat = gps(3);
% sdn = gps(1);
% chl = CHLAm;
% Depth =  P;

% ************************************************************************
% SET SOME DIRS / CONSTANTS
dirs = []; % testing
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

    dirs.cal       = [user_dir,'CAL\'];
    dirs.temp      = 'C:\temp\';
end

PAR_range  =  [400 700];
%PAR_range  =  [450 570]; % blue-green range
%MM01_dir   = 'C:\Users\jplant\Documents\MATLAB\Xing\';
%MM01_dir   = 'E:\Documents\MATLAB\Xing\';

planck     = 6.626070040e-34; % 10e-34 kg m^2 s^-1
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
% GLOBAL AVERAGE CALIBRATION GAIN CORRECTION. Due to NPQ this will be a
% minimum in attenuation

%CHL = chl * 0.5; % GLOBAL AVERAGE CALIBRATION GAIN CORRECTION

CHL = chl; % GLOBAL AVERAGE CALIBRATION GAIN CORRECTION
CHL(CHL<0) = 0;  % QUICK FIX SO MM01 doesn't produce imaginary #'s


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

met.vis  = 15; % visibility, km
met.am   = 1;  % atmospher factor, 1 for marine
met.wsm  = 4;  % windavg = 24hr mean wind speed (m s-1)
met.ws   = 6;  % windspd = wind speed (m s-1)
met.pres = 1013.25; % sea-level pressure (mb)
met.rh   = 80; % relhum = relative humidity (%)
met.wv   = 1.5; % water vapor , cm

% GreggCarder returns W /m^2/ nm
[lam, ID, IS, IT] = GreggCarder(sdn, lon, lat, PAR_range, met);

% CONVERT FROM W/m^2/nm to umol photons / m^2 / s-1  nm-1
% Used following website as a guide:
% https://www.berthold-bio.com/service-support/support-portal/ ...
% knowledge-base/how-do-i-convert-irradiance-into-photon-flux.html
%1E-9: nm => m in E = h*c/lambda; 1E6: moles => umoles
ID = ID * 1E-9 /(planck * lightspeed) / Avogadro *1E6;
IT = IT * 1E-9 /(planck * lightspeed) / Avogadro *1E6;
IS = IS * 1E-9 /(planck * lightspeed) / Avogadro *1E6;

% ALSO SLIGHT CORRECTION FOR HOW MUCH SURFACE SUNLIGHT PENETRATES AIR/SEA
% INTERFACE (* 0.98)
% Downwelling energy at surface for each PAR wavelength
% Interpolated to Kd MMO1 wave lengths
%Ed0 = interp1(lam, IT*0.98, lambda); % umol photons / m^2 / s-1  nm-1
% Removed 0.98 because reflectance is now included
Ed0 = interp1(lam, IT, lambda); % umol photons / m^2 / s-1  nm-1

% PAR profiles for each wavelength (Attenuate Ed0)
PAR = (fill_C * Ed0) .* exp(-INT_Kd); %umol photons / m^2 /s-1  nm-1


% NOW INTEGRATE ACROSS WAVE LENGTHS TO GET AVERAGE PAR PROFILE
WL      = fill_C * lambda;
PARxWL  = PAR .* WL; %(umol photons / m^2 / s-1) at each lambda 
INT_PAR = sum((PARxWL(:,1:end-1) + PARxWL(:,2:end))./2, 2) * dWL;

%pause

P      = P(IXO);
CHL    = CHL(IXO);
PAR    = INT_PAR(IXO);
