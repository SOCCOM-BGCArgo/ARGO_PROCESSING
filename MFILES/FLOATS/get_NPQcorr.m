function CHL = get_NPQcorr(gps, data, dirs)
% data = indata; %Testing
% gps =gps;

% gps  = [INFO.sdn, nanmean(INFO.gps,1)];
% %data = [lr_d(:,[iP,iT,iS]),LR.CHLA_ADJUSTED]; % APEX
% data = [hr_d(:,[iP,iT,iS]),HR.CHLA_ADJUSTED]; % NAVIS
% dirs = [];

% data = [tmp(:,[iP,iT,iS]),CHL_EXP];
% gps = tmp(1,[iSDN,iLON,iLAT]);
% ************************************************************************
% PURPOSE: 
%   This function corrects a chlorophyll fluorescence profile for the
%   effects of non-photochemical fluorescence quenching (NPQ) using 
%   Xing et al., 2012, doi 10.4319/lom.2012.10.483 with various flavors
%   or reference depths to start looking  for the CHL max above
%   Corrections are only applied when the sun angle is greater than 0
%   degrees.
%
% INPUTS:
%   gps         1x3 profile position array [matlab SDN(GMT), longitude, lattitude]
%   data        nx4 data matrix [Pressure, Temperature, Salinity, Adj Chl]
%                   Chl should be the chl calculated from the manufacturers
%                   calibration equation / 2
%   dirs        sturcture of dir paths or empty, if empty defaults are used
%
%
% OUTPUT:
%   CHL     a structure
%           .gps  	[matlab SDN(GMT), longitude, lattitude]
%           .hdr    header for data array
%           .data 	[P CHL  Xing_MLD  Xing_Zeu  Xing_Zchl  Xing PAR2]
%
%           .Dmld 	mixed layer depth (density)
%           .Tmld 	mixed layer depth (temperature)
%           .Zeu 	Euphotic zond depth estimate
%           .Zchl   Shallow chl max (max above shallower of Zeu or MLD)
%           .ZPAR20 Depth at PAR = 20 µmol photons / m^2 / seconds
%           .sun    sun angle
%
% REQUIRED FUNCTIONS:
%   SolarAzElq.m    a function to derive the sun angle

% ************************************************************************
% DO SOME PREP WORK
% ************************************************************************
% SET DEFAULT DIRS - ONLY CARRIED ALONG SO IT CAN BE USED IN EstPAR.m
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

    dirs.cal       = [user_dir,'CAL\'];
    dirs.temp      = 'C:\temp\';
end

% PREDIMMENSION OUTPUT

% CHL.hdr   = {'Pressure', 'CHL', 'Xing_MLD', 'Sack_MLD', 'Xing_Zeu', ...
%              'Sack_Zeu', 'Xing_Zchl', 'Sack_Zchl','Xing_Zchl2', 'Sack_Zchl2'};
CHL.hdr   = {'Pressure', 'CHL', 'CHLm','CHLspike', 'Xing_MLD', ...
             'Xing_Zeu','Xing_Zchl','Xing_Zchl2'};         
%CHL.data  = ones(size(data,1),size(CHL.hdr,2))* NaN;
CHL.data  = [];
CHL.Dmld  = [];
CHL.Tmld  = [];
CHL.Zeu   = [];
CHL.Zchl  = [];
CHL.ZPAR20  = [];
CHL.sun   = [];
CHL.slope = []; % CHL:BBP slope test diagnostics
CHL.XMLDZ = [];

if max(size(gps)) ~= 3 && min(size(gps)) ~= 1
    disp('GPS data is not valid - maybe no position fix')
    CHL.gps = gps;
    return
else
    CHL.gps   = gps;
end
 
% DEFINE indata column indices
iP = 1; iT = 2; iS = 3; iC = 4; iB = 5;

if all(isnan(data(:,iC)))
    CHL.data = [];
    disp('All chl data = NaN, exiting')
    return
end

if min(data(:,iP)) > 25
    CHL.data = [];
    fprintf('No surface data found for MLD reference point: Zmin = %0.1f\r\n',...
        min(data(:,iP)))
    return
end

% ************************************************************************
% ************************************************************************
% FIGURE OUT PRESSURE DIRECTION - MAKE SHALLOW TO DEEP BUT DON'T FORGET TO
% FLIP BACK CALCULATE MIXED LAYER DEPTH
pressure_dir = 0;
if data(1,iP) > data(end,iP) % data is deep too shallow
    data = flip(data,1);
    pressure_dir = 1;
end

% FIND DENSITY MLD
mld_den_threshold = 0.03; % Dong et al 2008
den      = density(data(:,iS),data(:,iT));
tDMLD    = den - nanmean(den(1:2)) <= mld_den_threshold; % mixed layer
CHL.Dmld = max(data(tDMLD,iP)); % depth of mixed layer
clear tDMLD mld_den_threshold den

% FIND ISOTHERMAL MLD
mld_temp_threshold = 0.2; % Dong et al 2008
tTMLD    = abs(data(:,iT) - nanmean(data(1:2,iT))) <= mld_temp_threshold; 
CHL.Tmld = max(data(tTMLD,iP)); % depth of mixed layer
clear tTMLD mld_temp_threshold

% tMLD = data(:,iP) <= min([CHL.Dmld CHL.Tmld]); % shallower of 2
% iMLD = data(:,iP) == max(data(tMLD,iP)); % base of mld index

% !!! FOR NOW ONLY USING DMLD PER ADMT18 CONSENSUS - jp 12/14/17 !!!
tMLD = data(:,iP) <= CHL.Dmld; % shallower of 2
iMLD = data(:,iP) == max(data(tMLD,iP)); % base of mld index

% ************************************************************************
% ************************************************************************
% PUT CHL THROUGH A 3 or 5 POINT MEDIAN FILTER TO REMOVE SPIKES
% THIS WILL BE USED FOR CHL:BBP RATIOS & OTHER TASKS
P      = data(:,iP);
%BBP    = data(:,iB);
CHLA   = data(:,iC);
dp     = nanmean(abs(diff(P(P<100))));

if dp > 3
    win    = 3; % APEX filter window, keep odd
else
    win    = 5; % NAVIS / PROVOR filter window, keep odd
end

offset = (win-1)/2;
%tmpB   = BBP * ones(1, win);
tmpC   = CHLA * ones(1, win);
for i = 1:win
    %tmpB(offset+1:end-offset,i) = BBP(i:end-win+i);
    tmpC(offset+1:end-offset,i) = CHLA(i:end-win+i);
end
%BBPm   = median(tmpB,2);
%CHLAm  = nanmedian(tmpC,2);
% NOTE: I believe nanmedian has a bug. If more 0 vlues than non zero values
% nanmedian will return a NaN. Using the updated matlab functionm which now
% has a nan flag 12/28/17
CHLAm  = median(tmpC,2,'omitnan');

clear tmpB tmpC win i

if all(isnan(CHLAm))
    CHL.data = [];
    disp('Median filterd Chl is all NaN - check raw data')
    return
end

% ************************************************************************
% FIND EUPHOTIC DEPTH ESTIMATE
% GET QUICK & DIRTY EUPHOTIC ZONE ESTIMATE FROM CHL
% Kim et al 2015, doi:10.5194/bgd-12-3905-2015
% I/Io = 1/100 = e^-(Kd*z), Zeu = 1n(100)/Kd
% Kd = 0.0232 + 0.074*chl.^0.674;Kd(Blue), Use Max CHL for AVG
% This will be a miniumum estimate for Zeu
%
% FLBB SENSORS USING FACTORY CALIBRATIONS CAN REPORT 6X TOO HIGH IN THE
% SOUTHERN OCEAN. THE BEST GLOBAL APROXIMATION is 2X TOO HIGH

% DO A BETTER JOB ATTENUATING LIGHT
tnan  = isnan(CHLAm);
KdCHL = CHLAm(~tnan);
PX    = P(~tnan);
if isempty(PX)
    disp('NO data left after Nan removal')
    return
end
    
%KdCHL(tnan) = median(KdCHL(~tnan)); % fill nan's with median if any
KdCHL(KdCHL < 0) = 0; % DOn't want any neg values in log

Kd     = 0.0232 + 0.074*(KdCHL).^0.674; % 1% light level (Kim et al 2015)

[~,IX] = sort(PX); % SORT BEFORE INTERP TO GET Zeu at 1%
Kd = Kd(IX);
PX = PX(IX);

PdP    = [PX(1); diff(PX)];
t_good = PdP ~= 0 & ~isnan(Kd); % remove duplicates before interpolation too
INT_KD = cumsum(Kd .* PdP);
CHL.Zeu = interp1(INT_KD(t_good), PX(t_good), -log(0.05));

tZeu    = data(:,iP) <= CHL.Zeu;

if sum(tZeu) == 0 % 1st sample deeper than Zeu
    iZeu = 1;
else
    iZeu    = data(:,iP) == max(data(tZeu,iP)); % base of Zeu index
end
clear CHLmax Kd

% ************************************************************************
% NOW GET Z AT PAR = 20
[PAR, ZPAR20] = EstPAR(gps(2), gps(3), gps(1), CHLAm, P,dirs);
CHL.ZPAR20    = ZPAR20;
tZPAR         = data(:,iP) <= CHL.ZPAR20;


% FIND DEPTH OF SHALLOW CHL MAX (ABOVE SHALLOWEST Zref POSIBILITY)
% Zmin   = nanmin([CHL.Dmld, CHL.Tmld, CHL.Zeu]);
% Zmin2  = nanmin([CHL.Dmld, CHL.Tmld, CHL.ZPAR20]);

Zmin   = nanmin([CHL.Dmld, CHL.Zeu]);
Zmin2  = nanmin([CHL.Dmld, CHL.ZPAR20]);

chl_chk = ~isnan(data(:,iC)) & data(:,iP) <= Zmin2;
if all(~chl_chk) 
    disp('No surface chl data in selected range - exiting')
    CHL.data =[];
    return
end


tZchl    = data(:,iP) <= Zmin; % 1's shalow , 0's deep
chl_tmp  = CHLAm;
chl_tmp(~tZchl) = NaN; % Set values below MLD to NaN
Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % find shallow chl max index
if isempty(Cmax)
    disp('No surface chl max found in selected range - exiting')
    CHL.Zchl = NaN; % depth of base of shallow chl max
    return
else
    tZchl  = data(:,iP) <= data(Cmax,iP);
    iZchl    = data(:,iP) == max(data(tZchl,iP)); % base of Zchl index
    CHL.Zchl = max(data(iZchl,iP)); % depth of base of shallow chl max
end


tZchl2    = data(:,iP) <= Zmin2; % 1's shalow , 0's deep
chl_tmp  = CHLAm;
chl_tmp(~tZchl2) = NaN; % Set values below MLD to NaN
Cmax2 = find(chl_tmp == max(chl_tmp),1, 'last'); % find shallow chl max index
tZchl2  = data(:,iP) <= data(Cmax2,iP);
iZchl2    = data(:,iP) == max(data(tZchl2,iP)); % base of Zchl index
CHL.Zchl2 = max(data(iZchl2,iP)); % depth of base of shallow chl max


% GET SUN ANGLE FIRST - IF < 5 deg, work is done
[Az,El] = SolarAzElq(gps(:,1),gps(:,3),gps(:,2),0); % sdn lat lon alt
%El(El<0) = 0; %SET NEGATIVE OR UNDERFOOT ELEVATIONS TO ZERO
CHL.sun  = El;
%if El >= 0 % temporary - find dark samples only
if El < 0
    CHL.data =[];
    return
end

% % CHECK  CHL:Z SLOPE for significance
% if sum(tZchl) > 2
%     T = test0slope(data(tZchl,iP),CHLA(tZchl),0.05);
%     if isempty(T)
%         disp('not enough unique data, can not compute statistics')
%         CHL.slope = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
%     else
%         CHL.slope = [T.test, T.my, T.by, T.t, T.p, T.N, T.alpha, T.DF];
%     end
% else
%     disp('sample size <= 2, can not compute statistics')
%     CHL.slope = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
% end

% CHECK  CHL:BBP vs Z SLOPE for significance
% if sum(tZchl2) > 2 % use chlmax above Zpar20
%     T = test0slope(data(tZchl2,iP), CHLAm(tZchl2)./BBPm(tZchl2), 0.05);
%     if isempty(T)
%         disp('not enough unique data, can not compute statistics')
%         CHL.slope = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
%     else
%         CHL.slope = [T.test, T.my, T.by, T.t, T.p, T.N, T.alpha, T.DF];
%     end
% else
%     disp('sample size <= 2, can not compute statistics')
%     CHL.slope = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
% end


    
%   T.test  = test outcome 1 = different from zero, 0 = not
%   T.alpha = signifigance level;
%   T.N     = number of samples
%   T.DF    = Degrees of freedom (N-2)
%   T.my    = regression slope, dY/dX;
%   T.by    = regression Y intercept;
%   T.t     = t_stat for slope;
%   T.p     = p value;

% ************************************************************************
% REF DEPTHS FOUND NOW DO VARIOUS FLAVORS OF CORRECTIONS
% ************************************************************************
% XING 2012 - MLD
chl_tmp = CHLAm;
chl_tmp(~tMLD) = NaN; % Set values below MLD to NaN
Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % max chl above MLD index
Xing_MLD  = CHLA; % Get original chl profile
Xing_MLD(1:Cmax) = chl_tmp(Cmax);
CHL.XMLDZ = P(Cmax);

chl_tmp = CHLAm;
chl_tmp(~tZeu) = NaN; % Set values below Zeu to NaN
Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % max chl above Zeu index
Xing_Zeu  = CHLA; % Get original chl profile
Xing_Zeu(1:Cmax) = chl_tmp(Cmax);

chl_tmp = CHLAm;
chl_tmp(~tZchl) = NaN; % Set values below Zchl to NaN
Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % max chl above Zchl index
Xing_Zchl  = CHLA; % Get original chl profile
Xing_Zchl(1:Cmax) = chl_tmp(Cmax);

chl_tmp = CHLAm;
chl_tmp(~tZchl2) = NaN; % Set values below Zchl to NaN
Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % max chl above Zchl index
Xing_Zchl2  = CHLA; % Get original chl profile
Xing_Zchl2(1:Cmax) = chl_tmp(Cmax);



% ************************************************************************
% SACKMANN 2008
% CHL_BBP = CHLAm(iMLD) ./ BBPm(iMLD); % ratio at base of mixed layer
% Sack_MLD = CHLA;
% Sack_MLD(tMLD) = BBPm(tMLD)*CHL_BBP;
% 
% CHL_BBP = CHLAm(iZeu) ./ BBPm(iZeu); % ratio at base of euphotic depth
% Sack_Zeu = CHLA;
% Sack_Zeu(tZeu) = BBPm(tZeu)*CHL_BBP;
% 
% CHL_BBP = CHLAm(iZchl) ./ BBPm(iZchl); % ratio at base of shallow chl max
% Sack_Zchl = CHLA;
% Sack_Zchl(tZchl) = BBPm(tZchl)*CHL_BBP;
% 
% % ZING 2017
% CHL_BBP = CHLAm(iZchl2) ./ BBPm(iZchl2); % ratio at base of shallow chl max
% Sack_Zchl2 = CHLA;
% Sack_Zchl2(tZchl2) = BBPm(tZchl2)*CHL_BBP;

% ************************************************************************
% BUILD OUTPUT ARRAY
% out = [data(:,iP), data(:,iC), Xing_MLD, Sack_MLD, Xing_Zeu, Sack_Zeu, ...
%     Xing_Zchl, Sack_Zchl, Xing_Zchl2, Sack_Zchl2];

out = [data(:,iP), data(:,iC), CHLAm, data(:,iC)-CHLAm, Xing_MLD, ...
      Xing_Zeu, Xing_Zchl, Xing_Zchl2];

if pressure_dir == 1
    out = flip(out,1);
end

CHL.data = out;






