% Float_vs_GLODAPv2
% PURPOSE: 
%   Extract all GLODAPv2 nutrient data within a given distance from each
%   track position if it exists.

% ************************************************************************
% CHOOSE FLOAT VARIABLE TO COMPARE {Float GLODAP}
% ************************************************************************
%compare_var = {'Oxygen[µmol/kg]' 'G2oxygen' '[µmol/kg]'};
%compare_var = {'Nitrate[µmol/kg]' 'G2nitrate' '[µmol/kg]'};
compare_var = {'pHinsitu[Total]' 'G2phts25p0' ' '};





% ************************************************************************
% GLODAP SETTINGS
depth_bnds = [0 2000]; % depth bounds in meters
%depth_bnds = [0 2000]; % depth bounds in meters
tol = 20; % km from profile

% ***********************************************************************
%   SET NAMES DIRS AND PATHS
% ***********************************************************************
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
GLODAP_dir          = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\', ...
                      'DATA\GLODAP\'];
GLODAP_fn           = 'GLODAPv2 Merged Master 102716.mat';

load(['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
      'float_type_lists.mat'])
jp = regexpi(NON_SOCCOM_list(:,1),'RosSea|SoAtl|SoPac|SoOcn|drake');
t1 = cellfun(@isempty,jp);
FLOAT_LIST = [SOCCOM_list; NON_SOCCOM_list(~t1,:)];
[~, ind]   = sort(FLOAT_LIST(:,2)); % sort by UW ID
FLOAT_LIST = FLOAT_LIST(ind,:);
clear NON_SOCCOM_list SOCCOM_list

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
% LOAD GLODAP AND INITIAL SUBSET DEPTH & < 25S & PACK INTO A MATRIX
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
d = load([GLODAP_dir,GLODAP_fn],load_vars{:}); % this line forces into structure

tP = d.G2pressure >= depth_bnds(1) & d.G2pressure <= depth_bnds(2);
Gdata = ones(sum(tP),size(G_INFO,1))* NaN;

for i = 1:size(G_INFO,1)
    if i < 14
        Gdata(:,i) = d.(G_INFO{i,1})(tP);
    elseif regexp(G_INFO{i,1},'^G2phts') % pH
        tmp = d.(G_INFO{i,1})(tP); % GRAB PARAMETER DATA
        tnan = d.([G_INFO{i,1},'f'])(tP) > 2 & ...
            d.('G2phtsqc')(tP) > 1; % CHECK FLAGS, BAD = NaN
        tmp(tnan) = NaN;
        Gdata(:,i) = tmp;       
    else
        tmp = d.(G_INFO{i,1})(tP); % GRAB PARAMETER DATA
        tnan = d.([G_INFO{i,1},'f'])(tP) > 2 & ...
            d.([G_INFO{i,1},'qc'])(tP) > 1; % CHECK FLAGS, BAD = NaN
        tmp(tnan) = NaN;
        Gdata(:,i) = tmp; 
    end
end
sdn   = datenum(Gdata(:,4),Gdata(:,5),Gdata(:,6),Gdata(:,7),Gdata(:,8),0);
Gdata = [sdn,Gdata(:,1:3), Gdata(:,9:end)];
Ghdr  = {'Date',G_INFO{1:3}, G_INFO{9:end}};
clear d tmp tnan tdata tP tlat tlon i ia t1 d sdn

% GET GLODAP INDICES
iP  = find(strcmp('G2pressure',Ghdr)     == 1);
iX  = find(strcmp(compare_var{1,2},Ghdr)     == 1);

% IF pH IS CHOOSEN VAR CONVERT GLODAP @25C&0P to GLODAP INSITU
if strcmp(compare_var{1,1}, 'pHinsitu[Total]')
    % GET SHIPBOARD DATA INDICES INDICES
    iT   = find(strcmp('G2temperature', Ghdr) == 1);
    iS   = find(strcmp('G2salinity',Ghdr)     == 1);
    iDIC = find(strcmp('G2tco2',Ghdr)         == 1);
    iALK = find(strcmp('G2talk',Ghdr)         == 1);
    iSI  = find(strcmp('G2silicate',Ghdr)     == 1);
    iPO4 = find(strcmp('G2phosphate',Ghdr)    == 1);

    pHSCALEIN     = 1;  % TOTAL SCALE
    K1K2CONSTANTS = 10; % Lueker et al, 2000
    KSO4CONSTANTS = 3;  % KSO4 of Dickson & TB of Lee 2010 (USE THIS ONE !!!)
    %KSO4CONSTANTS = 1;  % KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED)

    % IN SITU CONDITIONS
    SAL     = Gdata(:,iS);
    TEMPOUT = Gdata(:,iT);
    PRESOUT = Gdata(:,iP);
    SI      = Gdata(:,iSI);
    PO4     = Gdata(:,iPO4);

    % USE ALK & PH @25C P=0
    PAR1     = Gdata(:,iALK);
    PAR1TYPE = 1;
    PAR2     = Gdata(:,iX);
    PAR2TYPE = 3;
    TEMPIN   = 25;
    PRESIN   = 0;
    
    [OUT] = CO2SYSSOCCOM(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
        TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
        pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);

%     [OUT] = CO2SYSSOCCOM_old(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, ...
%         TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4, ...
%         pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS); 
    
    Gdata = [Gdata, OUT(:,18)]; % Col 18 in OUT is pH output
    Ghdr  = [Ghdr, 'G2ph insitu'];
    iX    = find(strcmp('G2ph insitu',Ghdr)     == 1); % redo pH index
    tQC   = Gdata(:,iX) < 7.3 | Gdata(:,iX) > 8.5; % GLODAP RANGE CHECK
    Gdata(tQC,iX) = NaN;
    
    clear SAL TEMPOUT PRESOUT SI PO4 PAR1 PAR1TYPE PAR2 PAR2TYPE
    clear TEMPIN PRESIN OUT
end

 
% ************************************************************************
% NOW LOAD FLOAT DATA
FV_dir = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
          'DATA\FLOATVIZ\QC\'];
      
if strcmp(compare_var{1,1},'Oxygen[µmol/kg]')
    hdr_add_str = 'Float O2 [µmole kg^{-1}]';
elseif strcmp(compare_var{1,1},'Nitrate[µmol/kg]')
    hdr_add_str = 'Float NO3 [µmole kg^{-1}]';
elseif strcmp(compare_var{1,1},'pHinsitu[Total]')
    hdr_add_str = 'Float pH';
end
    
hdr = ['UW_ID' 'Float Cycle','Cycle Date', Ghdr, hdr_add_str]; 
data = [];

for flt_ct = 1 : size(FLOAT_LIST,1)
    fname = [FLOAT_LIST{flt_ct,1},'QC.TXT'];
    UW_str = regexp(fname, '^\d{3}\d+', 'once', 'match');
    UW_ID  = str2double(UW_str);
    disp(['Processing float ', fname, ' .....'])
    
    d = get_FloatViz_data([FV_dir,fname]);
    if isempty(d)
        disp (['No file found for ',fname, ' moving to next float'])
        continue
    end
        
    iPF  = find(strcmp('Pressure[dbar]',d.hdr)     == 1);
    iXF  = find(strcmp(compare_var{1,1},d.hdr)     == 1);
    
    if isempty(iXF)
        disp (['No variable data for ', fname, ' moving to next float'])
        continue
    end
    
    
    % GOOD FLOAT DATA ONLY
    tQC = d.data(:,iXF+1) == 0;  % good data only - checking quality flag
    d.data(~tQC,iXF) = NaN;
    tMIV = d.data(:,iXF) == -1e10; % good data only - 2nd check on value
    d.data(tMIV,iXF) = NaN;
    
    % NOW DO VARIABLE SPECIFC RANGE CHECKS TOO
    if strcmp(compare_var{1,1},'Oxygen[µmol/kg]')
        tQC = d.data(:,iXF) < -5 | d.data(:,iXF) > 450;
        d.data(tQC,iXF) = NaN;
    elseif strcmp(compare_var{1,1},'Nitrate[µmol/kg]')
        tQC = d.data(:,iXF) < -5 | d.data(:,iXF) > 55;
        d.data(tQC,iXF) = NaN;
    elseif strcmp(compare_var{1,1},'pHinsitu[Total]')
        tQC = d.data(:,iXF) < 7.3 | d.data(:,iXF) > 8.5;
        d.data(tQC,iXF) = NaN;
    end
    

    [~,ia,~] = unique(d.data(:,2)); % GET TRACK
    track = d.data(ia,[1,2,4,3]);
    clear ia tQC
    % ***********************************************************************
    % GET / CHECK BOUNDS OF TRACK
    % ***********************************************************************
    cross180 = 0; % flag for meridian crossing
    
    % SDN LAT LON
    t1 = track(:,3) == -1e10; % SET MISSING VALUES IN LAT to NaN
    track(t1,3:4) = NaN;
    clear t1
    
    % GLODAP LON -180/+180 => CONVERT LON TO -180 to +180 IF NECESSARY
    t1 = track(:,4) > 180; % WOA2013 -180 + 180
    track(:,4) = track(:,4) - (t1*360);
    
    % CONVERT TOL TO DEGREES LAT AND LON
    lat_tol = tol/ 100.574; % aprox degrees lat
    lon_tol = abs(tol/ (111.320*cos(deg2rad(nanmean(track(:,3)))))); % aprox deg lon
    
    % GET LAT BOUNDS - Do 2nd subset
    lat_bnds = [min(track(:,3)) - lat_tol, max(track(:,3) + lat_tol)];
    if lat_bnds(1) < -90
        lat_bnds(1) = -90;
    elseif lat_bnds(2) > 90;
        lat_bnds(2) = 90;
    end
    
    % GET LON BOUNDS - 2nd subset
    lon_bnds = [min(track(:,4)) - lon_tol, max(track(:,4) + lon_tol)];
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
    
    tlat = Gdata(:,5) >= lat_bnds(1) & Gdata(:,5) <= lat_bnds(2);
    if cross180 == 0
        tlon = Gdata(:,6) >= lon_bnds(1) & Gdata(:,6) <= lon_bnds(2);
    else % cross180 = 1
        tlon = (Gdata(:,6) >= -180 & Gdata(:,6) <= lon_bnds(1)) |...
            (Gdata(:,6) >= lon_bnds(2) & Gdata(:,6) <= 180);
    end
    
    GD = Gdata(tlat&tlon, :);
    
    clear lat_tol lon_tol cross180 tlat tlon lon_bnds lat_bnds
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
    rr  = size(GD,1);
    
    
    
    for cycle = 1: stations
        if isnan(track(cycle,3)) % NO position (under ice, lost comms, etc)
            continue
        end
        % Haversine distance in KM from GLODAP point to track position
        [d1, ~] = lldistkm(GD(:,5:6),ones(rr,1)* track(cycle,3:4));
        t1 = d1 <=  tol;
        
        if sum(t1) > 0
            tFD = d.data(:,2) == track(cycle,2); % get float data
            FD  = d.data(tFD,:);
            tnan = isnan(FD(:,iXF));
            FD(tnan,:) =[]; % Get rid of NaN's in Nitrate especially for NAVIS FLOATS
            if isempty(FD(:,iXF)) || size(FD,1) < 2% NO float data - move on
                continue
            end
            
            Xi  = interp1(FD(:,iPF),FD(:,iXF), GD(t1,iP)); % interp float on GLODAP P
            
%             if  sum(Xi < -5)  > 0 % sum(isnan(Ni)) > 0 ||
%                 plot(FD(:,iXF), FD(:,iPF),'bo-', 'MarkerSize', 6)
%                 hold on
%                 plot(Xi,GD(t1,iP),'ko', 'MarkerSize', 4, ...
%                     'MarkerFaceColor', 'r')
%                 plot(GD(t1,iX),GD(t1,iP),'ko', 'MarkerSize', 6, ...
%                     'MarkerFaceColor', 'g');        
%                 
%                 set(gca, 'Ydir', 'reverse')
%                 xlabel('campare_var{1,1}')
%                 ylabel('Pressure')
%                 str = [num2str(UW_ID),'  Cast: ',num2str(track(cycle,2))];
%                 title(str)
%                 legend('Float','interp float', 'GLODAP')
%                 hold off
%                 
%                 pause
%             end
                

            data = [data; ones(sum(t1),1)*UW_ID, ...
                    ones(sum(t1),1)*track(cycle,2), ...
                    ones(sum(t1),1)*track(cycle,1) GD(t1,:), ...
                    Xi];
        end
    end
    
end


iXg  = find(strcmp(compare_var{1,2},hdr)     == 1); % GLODAP
if regexp(compare_var{1,2},'^G2ph','once')
    iXg  = find(strcmp('G2ph insitu',hdr)     == 1); % GLODAP
end


iXf  = find(strcmp(hdr_add_str,hdr)     == 1); % FLOAT

range = [(min([data(:,iXg);data(:,iXf)])), ...
    (max([data(:,iXg);data(:,iXf)]))];

% GET REGRESSION
tnan = isnan(data(:,iXg)) | isnan(data(:,iXf));
[m,b,r,sm,sb] = lsqfitma(data(~tnan,iXg), data(~tnan,iXf));
reg_info   = [m,b,r,sm,sb]; % DIFF REGRESSION OUTPUT
title_str  = {sprintf('Float = %1.3f * GLODAPv2 + %3.4f',reg_info(1:2)), ...
    sprintf('R*R = %1.4f', reg_info(3))};

plot(data(:,iXg),data(:,iXf), 'bo', 'MarkerSize', 2)
xlim(range)
ylim(range)
title(title_str)
ylabel(hdr_add_str)
xlabel(hdr{iXg})

hold on
plot(xlim, ylim, 'k--', 'LineWidth', 2)
plot(xlim, xlim .* reg_info(1) + reg_info(2),'r-');
legend('data', '1:1', 'model II', 'Location','NorthWest')
hold off

out = data;
out(:,3:4) = excelsdn(out(:,3:4));
tmp = regexp(hdr_add_str,' ','split'); % use this to build file namw


%Print out file
fid = fopen([GLODAP_dir, 'FLOATvsGLODAPv2_',tmp{1,2},'.txt'], 'w');

fprintf(fid,'File created on %s by %s\r\n', datestr(now),getenv('USERNAME'));
fprintf(fid,'Created with %s\r\n', mfilename);
fprintf(fid,'GLODAP depth range %4.0f to %4.0f meters\r\n', depth_bnds);
fprintf(fid,'Max GLODAP distance from float profile position = %4.0fkm\r\n', tol);

% BUILD DATA FORMAT STR
df ='';
for i = 1:size(out,2)
    if i == size(out,2);
        df = [df, '%f\r\n'];
    else
        df = [df,'%f\t'];
    end
end

for i = 1:size(hdr,2)
    if i == size(hdr,2)
        fprintf(fid,'%s\r\n', hdr{i});
    else
        fprintf(fid,'%s\t', hdr{i});
    end
end

disp('printing data to file')
for i = 1: size(out,1)
    str = sprintf(df,out(i,:));    
    fprintf(fid, '%s', strrep(str,'NaN',''));
end

fclose(fid);

  








    

    
    









