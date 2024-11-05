function CHL = get_NPQcorr(gps, data, dirs)
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
%
% CHNAGE HISTORY:
%  06/12/2019 - improved MLD calculation code to better deal with missing 
%     refrence depth samples or missing data below the MLD. -jp
%  03/02/2020 - JP - added a shallow profile check to account for 0889 surf
%  only measurements
%  05/15/24 - TM modification to smoothing and MLD specs per most recent Argo QC doc.  Argo CHLA QC doc Sept2023 : http://dx.doi.org/10.13155/35385 

% ************************************************************************
% TESTING:
% data = indata; %Testing
% gps =gps;

% gps  = [INFO.sdn, nanmean(INFO.gps,1)];
% % %data = [lr_d(:,[iP,iT,iS]),LR.CHLA_ADJUSTED]; % APEX
% data = [hr_d(:,[iP,iT,iS]),HR.CHLA_ADJUSTED]; % NAVIS
% dirs = [];

% data = [tmp(:,[iP,iT,iS]),CHL_EXP];
% gps = tmp(1,[iSDN,iLON,iLAT]);
% dirs = [];

% d = parse_NAVISmsg4ARGO('\\atlas\ChemWebData\floats\n0889\0889.031.msg');
% gps  = [d.sdn, nanmean(d.gps,1)];
% data = d.lr_d(:,[1:3,7]); % NAVIS
% dirs =[];

% ************************************************************************

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
% CHL.hdr   = {'Pressure', 'CHL', 'CHLm','CHLspike', 'Xing_MLD', ...
%              'Xing_Zeu','Xing_Zchl','Xing_Zchl2'}; 
CHL.hdr   = {'Pressure', 'CHL', 'CHLm','CHLspike', 'Xing_MLD', ...
             'Xing_Zeu','Xing_Zchl'};
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
CHL.CHLAm = [];
CHL.exit_msg = '';

if max(size(gps)) ~= 3 && min(size(gps)) ~= 1
    disp('GPS data is not valid - maybe no position fix')
    CHL.gps = gps;
    return
else
    CHL.gps   = gps;
end
 
% DEFINE indata column indices
iP = 1; iT = 2; iS = 3; iC = 4; iB = 5;

% CHECK FOR ANY NaN's IN PRESS TEMP OR PSAL & REMOVE TEMPORARILY
orig_data = data;
%P_nan     = isnan(orig_data(:,iP));
P_nan     = isnan(orig_data(:,iP)) | isnan(orig_data(:,iT)) | ...
            isnan(orig_data(:,iS));
        
if sum(P_nan > 0)
    disp(' ');
    disp(['NaNs in P, T or S for cycle date ',datestr(gps(1)), ...
        ' temporarily removed for data processing purposes']);
    data = data(~P_nan,:);
end

if all(isnan(data(:,iC)))
    CHL.data = [];
    CHL.exit_msg = 'All chl data = NaN, exiting';
    disp(CHL.exit_msg)
    return
end

if min(data(:,iP)) > 25
    CHL.data = [];
    CHL.exit_msg = sprintf(['No surface data found for MLD reference ',...
        'point: Zmin = %0.1f\r\n'], min(data(:,iP)));
    disp(CHL.exit_msg)
    return
end

if max(data(:,iP))< 75
    CHL.data = [];
    CHL.exit_msg = sprintf(['Shallow profile! MLD calculation ',...
        'questionable: Zmax = %0.1f meters\r\n'], max(data(:,iP)));
    disp(CHL.exit_msg)
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

% There will be at least one data point <=25m
% data has already been set to shallow to deep
ref_ind  = find(data(:,iP)<=25, 2, 'first'); % find ref points for mld calc
if isempty(ref_ind) % redundant but just to be safe
    CHL.data = [];
    CHL.exit_msg = sprintf(['No reference indices found for MLD reference ',...
        'point: Zmin = %0.1f\r\n'], min(data(:,iP)));
    disp(CHL.exit_msg)
    return
end

tDMLD    = den - nanmean(den(ref_ind)) <= mld_den_threshold; % mixed layer
if all(tDMLD ==1)  % no data below mld
    CHL.exit_msg = sprintf(['MLD reference exists but no data found ', ...
        'below MLD. MLD could not be determined Zmax = %0.1f\r\n'], ...
        max(data(:,iP)));
    disp(CHL.exit_msg)
    return
end


CHL.Dmld = max(data(tDMLD,iP)); % depth of mixed layer
clear tDMLD mld_den_threshold den

% FIND ISOTHERMAL MLD
mld_temp_threshold = 0.2; % Dong et al 2008
tTMLD    = abs(data(:,iT) - nanmean(data(ref_ind,iT))) <= mld_temp_threshold; 
CHL.Tmld = max(data(tTMLD,iP)); % depth of mixed layer
clear tTMLD mld_temp_threshold

% tMLD = data(:,iP) <= min([CHL.Dmld CHL.Tmld]); % shallower of 2
% iMLD = data(:,iP) == max(data(tMLD,iP)); % base of mld index

% !!! FOR NOW ONLY USING DMLD PER ADMT18 CONSENSUS - jp 12/14/17 !!!
%tMLD = data(:,iP) <= CHL.Dmld; % shallower of 2
tMLD = data(:,iP) <= CHL.Dmld.*0.9; % per Argo CHLA QC doc Sept2023 : http://dx.doi.org/10.13155/35385 pg 17.
% %iMLD not used??  And causing error in cases where 0.9*CHL.Dmld is
% %shallower than shallowest measurement... ie 3902559 cycle 87.  TM 5/16/24
% iMLD = data(:,iP) == max(data(tMLD,iP)); % base of mld index

% ************************************************************************
% ************************************************************************
% PUT CHL THROUGH A 3 or 5 POINT MEDIAN FILTER TO REMOVE SPIKES
% THIS WILL BE USED FOR CHL:BBP RATIOS & OTHER TASKS
P      = data(:,iP);
%BBP    = data(:,iB);
CHLA   = data(:,iC);
dp     = nanmean(abs(diff(P(P<100))));

% TM modified window sizes per Argo CHLA QC doc Sept2023 : http://dx.doi.org/10.13155/35385 pg 17.
if dp > 3
    win    = 5; % APEX filter window, keep odd
else
    win    = 7; % NAVIS / PROVOR filter window, keep odd
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
CHLAm     = median(tmpC,2,'omitnan');

% ADD CHLAm to output structure for dark data returns
Ctmp = CHLAm;
if pressure_dir == 1
    Ctmp = flip(CHLAm,1);
end
CHL.CHLAm = ones(size(orig_data(:,1)))*NaN;
CHL.CHLAm(~P_nan) = Ctmp;

clear tmpB tmpC win i Ctmp

if all(isnan(CHLAm))
    CHL.data = [];
    CHL.exit_msg = 'Median filtered Chl is all NaN - check raw data';
    disp(CHL.exit_msg)
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
    CHL.exit_msg = 'NO data left after NaN removal';
    disp(CHL.exit_msg)
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
% [PAR, ZPAR20] = EstPAR(gps(2), gps(3), gps(1), CHLAm, P,dirs);
% CHL.ZPAR20    = ZPAR20;
% tZPAR         = data(:,iP) <= CHL.ZPAR20;


% FIND DEPTH OF SHALLOW CHL MAX (ABOVE SHALLOWEST Zref POSIBILITY)
% Zmin   = nanmin([CHL.Dmld, CHL.Tmld, CHL.Zeu]);
% Zmin2  = nanmin([CHL.Dmld, CHL.Tmld, CHL.ZPAR20]);

Zmin   = nanmin([CHL.Dmld, CHL.Zeu]);
%Zmin2  = nanmin([CHL.Dmld, CHL.ZPAR20]);

%chl_chk = ~isnan(data(:,iC)) & data(:,iP) <= Zmin2;
chl_chk = ~isnan(data(:,iC)) & data(:,iP) <= Zmin;
if all(~chl_chk) 
    CHL.exit_msg = 'No surface chl data in selected range - exiting';
    disp(CHL.exit_msg)
    CHL.data =[];
    return
end


tZchl    = data(:,iP) <= Zmin; % 1's shalow , 0's deep
chl_tmp  = CHLAm;
chl_tmp(~tZchl) = NaN; % Set values below MLD to NaN
Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % find shallow chl max index
if isempty(Cmax)
    CHL.exit_msg = 'No surface chl max found in selected range - exiting';
    disp(CHL.exit_msg)
    CHL.Zchl = NaN; % depth of base of shallow chl max
    return
else
    tZchl  = data(:,iP) <= data(Cmax,iP);
    iZchl    = data(:,iP) == max(data(tZchl,iP)); % base of Zchl index
    CHL.Zchl = max(data(iZchl,iP)); % depth of base of shallow chl max
end


% tZchl2    = data(:,iP) <= Zmin2; % 1's shalow , 0's deep
% chl_tmp  = CHLAm;
% chl_tmp(~tZchl2) = NaN; % Set values below MLD to NaN
% Cmax2 = find(chl_tmp == max(chl_tmp),1, 'last'); % find shallow chl max index
% tZchl2  = data(:,iP) <= data(Cmax2,iP);
% iZchl2    = data(:,iP) == max(data(tZchl2,iP)); % base of Zchl index
% CHL.Zchl2 = max(data(iZchl2,iP)); % depth of base of shallow chl max


% GET SUN ANGLE FIRST - IF < 5 deg, work is done
[Az,El] = SolarAzElq(gps(:,1),gps(:,3),gps(:,2),0); % sdn lat lon alt
%El(El<0) = 0; %SET NEGATIVE OR UNDERFOOT ELEVATIONS TO ZERO
CHL.sun  = El;
%if El >= 0 % temporary - find dark samples only
if El < 0
    CHL.data =[];
    CHL.exit_msg = ['Sun elevation < 0 degress (',num2str(El), ...
        ') - No NPQ correction'];
    disp(CHL.exit_msg)
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
if isempty(Cmax)
    CHL.exit_msg = 'No surface chl max found in selected range - exiting';
    disp(CHL.exit_msg)
    CHL.data =[];
    return
end
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

% chl_tmp = CHLAm;
% chl_tmp(~tZchl2) = NaN; % Set values below Zchl to NaN
% Cmax = find(chl_tmp == max(chl_tmp),1, 'last'); % max chl above Zchl index
% Xing_Zchl2  = CHLA; % Get original chl profile
% Xing_Zchl2(1:Cmax) = chl_tmp(Cmax);



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

% out = [data(:,iP), data(:,iC), CHLAm, data(:,iC)-CHLAm, Xing_MLD, ...
%       Xing_Zeu, Xing_Zchl, Xing_Zchl2];

out = [data(:,iP), data(:,iC), CHLAm, data(:,iC)-CHLAm, Xing_MLD, ...
      Xing_Zeu, Xing_Zchl];
  
if pressure_dir == 1
    out = flip(out,1);
end

final_data = ones(size(orig_data,1), size(out,2))* NaN;
final_data(~P_nan,:) = out; % put nan's back in

    


CHL.data = final_data;







