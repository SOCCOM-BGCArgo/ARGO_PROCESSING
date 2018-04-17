function d = get_GLODAPv2_local(track, tol, depth_bnds)
% PURPOSE: 
%   Extract all GLODAPv2 nutrient data within a given distance from each
%   track position if it exists. Intended use is to find cross overs for
%   BGC ARGO float profiles for data quality control purposes
%   variable and then interpolate data along the provided track. 
%
% USAGE:
%	d = get_GLODAPv2_local(track, tol, depth_bnds)
%
% INPUTS:
%   track      = n x 4 matrix [Matlab_SDN, float cycle, Lat, Lon]
%   tol        = km distance to expand search from track position
%   depth_bnds = depth bounds [min depth  max depth]
%
% OUTPUTS:
%	d =   a stucture
%       d.hdr  = colum headers
%                [track row index, Date, cruise, station, cast, latitude,
%                 longitude, pressure, oxygen, nitrate, phts25p0'}
%       d.data = GLODAPv2 data near track positions
%                
%
%
%              
% EXAMPLES:
%  

% EDIT LOG
%   11/30/2016 6960 waas not grabbing any GLODAP data. Fixed bug in lat
%       & lon bounds max-tol should have been max+tol

%plot_it = 0; % 0 to turn off plotting

% ***********************************************************************
%   TESTING
% ***********************************************************************
% depth_bnds = [0 2000]; % depth bounds in meters
% tol = 30; % km from profile
% d = get_FloatViz_data(['C:\Users\jplant\Documents\MATLAB\ARGO\DATA\', ...
%      'FloatViz\6960ETNP.TXT']);
% [~,ia,~] = unique(d.data(:,2));
% track = d.data(ia,[1,2,4,3]); 
% clear d

% ***********************************************************************
%   SET NAMES DIRS AND PATHS
% ***********************************************************************
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
GLODAP_dir          = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\', ...
                      'DATA\GLODAP\'];
GLODAP_fn           = 'GLODAPv2 Merged Master 102716.mat';
%GLODAP_expocode_fn  = 'GLODAPv2_ExpoCodes_102716.txt';

% ***********************************************************************
%   DEFAULT OUTPUT
% ***********************************************************************
data =[];
hdr = {};

% ***********************************************************************
%   LIST OF DESIRED GLODAP VARIABLES
% ***********************************************************************
G_INFO(1,:)  = {'G2cruise'};      % REQIURED   % ID # WHICH DEFINES EXPO #
G_INFO(2,:)  = {'G2station'};     % REQIURED
G_INFO(3,:)  = {'G2cast'};        % REQIURED
G_INFO(4,:)  = {'G2year'};        % REQIURED
G_INFO(5,:)  = {'G2month'};       % REQIURED
G_INFO(6,:)  = {'G2day'};         % REQIURED
G_INFO(7,:)  = {'G2hour'};        % REQIURED
G_INFO(8,:)  = {'G2minute'};      % REQIURED
G_INFO(9,:)  = {'G2latitude'};    % REQIURED
G_INFO(10,:) = {'G2longitude'};   % REQIURED
G_INFO(11,:) = {'G2pressure'};    % REQIURED
G_INFO(12,:) = {'G2temperature'}; % REQIURED
G_INFO(13,:) = {'G2salinity'};    % REQIURED

G_INFO(14,:) = {'G2oxygen'};      % USER CAN ADJUST THESE
G_INFO(15,:) = {'G2nitrate'};
G_INFO(16,:) = {'G2silicate'};
G_INFO(17,:) = {'G2phosphate'};
G_INFO(18,:) = {'G2talk'};  
G_INFO(19,:) = {'G2tco2'};  
G_INFO(20,:) = {'G2phts25p0'}; % pH at total scale, 25c p =0 (convert)

% ***********************************************************************
% CHECK TRACK INPUT
% ***********************************************************************
s2 = ['Check track input: must have 4 rows or columns ', ...
     '[SDN float_cycle Lat Lon]'];
s3 = ' Check depth range input: should be [min_depth max_depth]';


[r,c] = size(track);
if r ~= 4 && c ~= 4
    disp(s2)
    return
end

if r == 4 && c ~= 4
    track = track'; % now  n x 4
end

if length(depth_bnds) ~= 2 || depth_bnds(1) > depth_bnds(2)
    disp(s3)
    return
end
clear s2 s3 r c t1

% ***********************************************************************
% GET / CHECK BOUNDS OF TRACK 
% ***********************************************************************
cross180 = 0; % flag for meridian crossing

% SDN LAT LON
t1 = track(:,3) == -1e10; % SET MISSING VALUES IN LAT to NaN
track(t1,3:4) = NaN;
clear d t1

% GLODAP LON -180/+180 => CONVERT LON TO -180 to +180 IF NECESSARY
t1 = track(:,4) > 180; % WOA2013 -180 + 180
track(:,4) = track(:,4) - (t1*360);

% CONVERT TOL TO DEGREES LAT AND LON
lat_tol = tol/ 100.574; % aprox degrees lat
lon_tol = abs(tol/ (111.320*cos(deg2rad(nanmean(track(:,3)))))); % aprox deg lon

% GET LAT BOUNDS - ADD 2X TOL FOR 1st SUBSET
lat_bnds = [min(track(:,3))- 2*lat_tol, max(track(:,3)+ 2*lat_tol)];
if lat_bnds(1) < -90
    lat_bnds(1) = -90;
elseif lat_bnds(2) > 90;
    lat_bnds(2) = 90;
end

%pause
% GET LON BOUNDS - ADD 2X TOL FOR 1st SUBSET
lon_bnds = [min(track(:,4))- 2*lon_tol, max(track(:,4)+ 2*lon_tol)];
if lon_bnds(1) < -180
    disp('longitude bounds - tol crosses +180 / -180 merdian')
    cross180 = 1;
    lon_bnds(1) = lon_bnds(1)+360;
    lon_bnds = sort(lon_bnds); % this will reverse order
elseif lon_bnds(2) > 180;
    disp('longitude bounds + tol crosses +180 / -180 merdian')
    cross180 = 1;
    lon_bnds(2) = lon_bnds(2)-360;
    lon_bnds = sort(lon_bnds); % this will reverse order
end
    
% CHECK TRACK CROSSES -180 / +180 MERIDIAN 
if max(abs(diff(track(:,4)))) > 340 % hard to immage a 20 degree step
    disp('Trackline crosses +180 / -180 merdian')
    cross180 = 1; % float crosses 0 / 360 longitude line
end

% ***********************************************************************
% LOAD GLODAP AND SUBSET DATA WITHIN TRACK BOUNDS
% CHECK GLODAP QUALITY FLAGS 'f' AND 'qc' & SET BAD = NaN
% ***********************************************************************lat
% BUILD VARIABLE LIST AND LOAD VARIABLES INTO A STRUCTURE. THIS WAY I
% CAN REUSE INFO LIST TO ACCESS VARIABLES & USER CAN CHANGE LIST
% FOR NUTRIENTS LOAD QUALITY FLAGS *.f AND *.qc
load_vars ={};
for i = 1: size(G_INFO,1)
    if i < 14
        load_vars = [load_vars, G_INFO{i,1}];
    elseif regexp(G_INFO{i,1},'^G2phts') % pH
        load_vars = [load_vars, G_INFO{i,1},[G_INFO{i,1},'f'],'G2phtsqc'];        
    else
        load_vars = [load_vars, G_INFO{i,1},[G_INFO{i,1},'f'], ...
            [G_INFO{i,1},'qc']];
    end
end

% LOAD DESIRED VARIABLES & CONTAIN IN A STRUCTURE
d = load([GLODAP_dir,GLODAP_fn],load_vars{:});

% NOW SUBSET AND PACK INTO A MATRIX
tP = d.G2pressure >= depth_bnds(1) & d.G2pressure <= depth_bnds(2);
tlat = d.G2latitude >= lat_bnds(1) & d.G2latitude <= lat_bnds(2);
if cross180 == 0
    tlon = d.G2longitude >= lon_bnds(1) & d.G2longitude <= lon_bnds(2);
else % cross180 = 1
    tlon = (d.G2longitude >= -180 & d.G2longitude <= lon_bnds(1)) |...
        (d.G2longitude >= lon_bnds(2) & d.G2longitude <= 180);
end
tdata = tP & tlat & tlon; 
Gdata = ones(sum(tdata),size(G_INFO,1))* NaN;
for i = 1:size(G_INFO,1)
    if i < 14
        Gdata(:,i) = d.(G_INFO{i,1})(tdata);
    elseif regexp(G_INFO{i,1},'^G2phts') % pH
        tmp = d.(G_INFO{i,1})(tdata); % GRAB PARAMETER DATA
        tnan = d.([G_INFO{i,1},'f'])(tdata) > 2 & ...
            d.('G2phtsqc')(tdata) > 1; % CHECK FLAGS, BAD = NaN
        tmp(tnan) = NaN;
        Gdata(:,i) = tmp;
    else
        tmp = d.(G_INFO{i,1})(tdata); % GRAB PARAMETER DATA
        tnan = d.([G_INFO{i,1},'f'])(tdata) > 2 & ...
            d.([G_INFO{i,1},'qc'])(tdata) > 1; % CHECK FLAGS, BAD = NaN
        tmp(tnan) = NaN;
        Gdata(:,i) = tmp; 
    end
end
sdn   = datenum(Gdata(:,4),Gdata(:,5),Gdata(:,6),Gdata(:,7),Gdata(:,8),0);
Gdata = [sdn,Gdata(:,1:3), Gdata(:,9:end)];
Ghdr  = {'Date',G_INFO{1:3}, G_INFO{9:end}};
clear d tmp tnan tdata tP tlat tlon i ia t1 d


% ***********************************************************************
% ***********************************************************************
% STEP THROUGH TRACK POINTS AND LOOK FOR CROSS OVER GLODAP DATA
% BASED ON LAT AND LON AND TOL VARIABLE
%
% NOTE: A GIVEN GLODAP PROFILE MAY BE A CROSSOVER FOR MULTIPLE POSITIONS OR
%       OR A TRACK POINT MAY HAVE MULTIPLE GLODAP PROFILES.
% ***********************************************************************
% ***********************************************************************
stations = size(track(:,1),1); % get # of points to interpolate
rr  = size(Gdata,1);

data = [];
for i = 1: stations
    if isnan(track(i,3)) % NO position (under ice, lost comms, etc)
        continue
    end
    % Haversine distance in KM from GLODAP point to track position
    [d1, ~] = lldistkm(Gdata(:,5:6),ones(rr,1)* track(i,3:4));
    t1 = d1 <=  tol;
    if sum(t1) > 0
        data = [data; ones(sum(t1),1)*track(i,2), Gdata(t1,:)];
    end
end

% ***********************************************************************
% NOW CALCULATE IN SITU PH USING SOCCOM CONSTANTS from pH at 25C
% ***********************************************************************

if ~isempty(data)
    hdr = ['float cycle', Ghdr];
    
    pHSCALEIN     = 1;  % TOTAL SCALE
    K1K2CONSTANTS = 10; % Lueker et al, 2000
    KSO4CONSTANTS = 3;  % KSO4 of Dickson & TB of Lee 2010 (USE THIS ONE !!!)
    %KSO4CONSTANTS = 1;  % KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED)
    
    % GET SHIPBOARD DATA INDICES INDICES
    iP   = find(strcmp('G2pressure', hdr)    == 1); % dbar
    iT   = find(strcmp('G2temperature', hdr) == 1);
    iS   = find(strcmp('G2salinity',hdr)     == 1);
    iDIC = find(strcmp('G2tco2',hdr)         == 1);
    iALK = find(strcmp('G2talk',hdr)         == 1);
    iPH  = find(strcmp('G2phts25p0',hdr)     == 1); % 25C & 0dbar
    iSI  = find(strcmp('G2silicate',hdr)     == 1);
    iPO4 = find(strcmp('G2phosphate',hdr)    == 1);
    
    % IN SITU CONDITIONS
    SAL     = data(:,iS);
    TEMPOUT = data(:,iT);
    PRESOUT = data(:,iP);
    SI      = data(:,iSI);
    PO4     = data(:,iPO4);
    
    % ************************************************************************
    % USE ALK & PH @25C P=0
    PAR1     = data(:,iALK);
    PAR1TYPE = 1;
    PAR2     = data(:,iPH);
    PAR2TYPE = 3;
    TEMPIN   = 25;
    PRESIN   = 0;
    
    [OUT] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
        TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
        pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);
    
    data = [data, OUT(:,18)]; % Col 18 in OUT is pH output
    hdr  = [hdr, 'ph_insitu'];
    
    % % ************************************************************************
% % USE ALK & TCO2
%     PAR1 = data(:,iALK);
%     PAR1TYPE   = 1; 
%     PAR2 = data(:,iDIC);
%     PAR2TYPE   = 2;
%     TEMPIN = data(:,iT);
%     PRESIN = data(:,iP);
% 
% [OUT] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
%     TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
%     pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);
% 
% data = [data, OUT(:,18)]; % Col 18 in OUT is pH output
% hdr  = [hdr, 'ph_insitu_alk&tco2'];
    
end







% CREATE OUTPUT AND CLEAN UP
d.hdr = hdr;
d.data = data;
clearvars -except d


    

    
    









