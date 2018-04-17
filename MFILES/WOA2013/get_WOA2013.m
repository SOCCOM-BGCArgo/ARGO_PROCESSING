function data = get_WOA2013(track, depth_bnds, ocean_var)
% PURPOSE: 
%   Extract a subset from the monthly WOA2013 climatology for a given
%   variable and then interpolate data along the provided track. Access
%   different data sets by adjusting paths and entries in the WOA_info
%   cell array. My use is for profiling float data. Must Build 4D matrix
%   from 3D monthly variable files first.
%   [month x depth x lat x lon] => [time x depth x lat x lon]
%   [1 x m x n x p]             => [12 x m x n x p]
%
% USAGE:
%	data = get_WOA2013(track, depth_bnds, ocean_var)
%
% INPUTS:
%   track      = n x 3 matrix [Matlab_SDN, Lat, Lon]
%   depth_bnds = depth bounds [min depth  max depth]
%	ocean_var  = a string, ocean parameter
%                must be one of these: T S O2 O2sat NO3 Si PO4 AOU
%
% OUTPUTS:
%	data =   a data structure.
%       data.d   = data matrix subset for WOA2013 variable [depth x time]
%       data.z   = depth array subset
%              
% EXAPLES:
%   data = get_WOA2013(track, [0 1000], 'NO3')
%
%   SCRIPT REQUIRES THESE HELPER FUNCTIONS / SCRIPTS / TOOLBOXES:
%       NCTOOLBOX (from : https://code.google.com/p/nctoolbox/)

       
% ***********************************************************************
% **********************       WOA 2013 INFO       **********************
% WOA2013 file name format:  woa13_all_[v][tp][ft][gr].[form_end]
%    [v] = oceanographic variable (t s o n i p)
%           t = temperature
%           s = salinity
%           o = oxygen
%           O = % Oxygen saturation
%           n = nitrate
%           i = silicate
%           p = phosphate
%    [tp]= time averaging period
%           00 – annual statistics, all data used;
%           01 to 12 – monthly statistics (01 = Jan, 12 – Dec);
%           13 to 16 – seasonal statistics:
%               13 – North Hemisphere winter (January - March);
%               14 – North Hemisphere spring (April - June);
%               15 – North Hemisphere summer (July - September);
%               16 – North Hemisphere autumn (October - December);
%    [ft] = field type = an = Objectively analyzed climatological mean
%    [gr] = the grid size (01 – 1-degree grid resolution this use)
%
%  NOTE: Monthly nitrate, phosphate, silica only to 500m
% dimensions = time x depth x lat x lon
% lat -89.5 to 89.5
% lon -179.5 to +179.5

% ***********************************************************************
% SET SOME VARIABLES, PATHS, STRINGS, TEMPLATES
% ***********************************************************************

plot_it = 0; % 0 to turn off plotting
data_site  = 'http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA13/DATA/';
ST_path    = '/netcdf/A5B2/1.00/';  % ARGO FLOAT BASED 2005 - 2012, 1 deg
NUT_path   = '/netcdf/all/1.00/';   % All nutrients, 1 deg
fn_str     = 'woa13_XXX_VNN_01.nc'; % file name template to modify 
%O2_volume  =  22.3916; % L/ mole O2 @ STP

% VARIABLES EXTRACTION TABLE
WOA_info(1,:) = {'nitrate'    , 'n', NUT_path, 'all', 'NO3'};
WOA_info(2,:) = {'oxygen'     , 'o', NUT_path, 'all', 'O2' };
WOA_info(3,:) = {'phosphate'  , 'p', NUT_path, 'all', 'PO4'};
WOA_info(4,:) = {'silicate'   , 'i', NUT_path, 'all', 'Si' };
WOA_info(5,:) = {'temperature', 't', ST_path , 'A5B2', 'T' };
WOA_info(6,:) = {'salinity'   , 's', ST_path , 'A5B2', 'S' };
WOA_info(7,:) = {'o2sat'      , 'O', NUT_path, 'all', 'O2sat'};
WOA_info(8,:) = {'AOU'        , 'A', NUT_path, 'all', 'AOU'};

month_label = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' ...
               'Oct' 'Nov' 'Dec'};

% ***********************************************************************
% CHECK INPUTS
% ***********************************************************************
s1 = ['ocean_var must be a character string from this list: ', ...
      'T S O2 O2sat NO3 Si PO4 AOU'];  
s2 = 'Check track input: must have 3 rows or columns [SDN Lat Lon]';
s3 = ' Check depth range input: should be [min_depth max_depth]';

% if nargin ~= 3
%     disp(['Check inputs (3 required): track matrix, depth bounds', ...
%           ' and variable string']);
%     return
% end
[r,c] = size(track);
if r ~= 3 && c ~= 3
    disp(s2)
    return
end
if r == 3 && c ~= 3
    track = track'; % now  n x 3
end
if length(depth_bnds) ~= 2 || depth_bnds(1) > depth_bnds(2)
    disp(s3)
    return
end
if ~ischar(ocean_var) % not a sting
    disp(s1)
    return
end
var_ind = find(strcmp(ocean_var, WOA_info(:,5)) == 1);
if isempty(var_ind) % not a variable from the designated list
    disp(s1)
    return
end

% SDN LAT LON
t1 = track(:,2) == -1e10; % SET MISSING VALUES IN LAT to NaN
track(t1,2:3) = NaN;
clear d t1

% CONVERT LON TO -180 to +180 IF NECESSARY
t1 = track(:,3) > 180; % WOA2013 -180 + 180
track(:,3) = track(:,3) - (t1*360);

% CHECK FOR -180 / +180 MERIDIAN TRACK CROSSING, GET AREA BOUNDS
cross180 = 0;
if max(abs(diff(track(:,3)))) > 340 % hard to immage a 20 degree step
    disp('Trackline crosses +180 / -180 merdian')
    cross180 = 1; % float crosses 0 / 360 longitude line
    t1 = track(:,3) < 0; % break  - / +
    lon_bnds = [max(track(t1,3)) min(track(~t1,3))]; % ) cross bounds
else
    lon_bnds = [min(track(:,3)) max(track(:,3))]; % lon
end
lat_bnds = [min(track(:,2)) max(track(:,2))]; % lat    

clear s1 s2 r c t1

% ***********************************************************************
%   EXTRACT WOA2013 MONTHLY CLIMATOLOGY SUBSETS & BUILD ANNUAL MATRIX (4D)
%       FROM INDIVIDUAL MONTHS: [TIME x DEPTH x LAT x LON]
%             Time for each month is singlton so squeeze
% ***********************************************************************
data_path = [data_site, WOA_info{var_ind,1},WOA_info{var_ind,3}];
fn_str     = regexprep(fn_str,'V',WOA_info{var_ind,2}); % build name: add var
fn_str     = regexprep(fn_str,'XXX',WOA_info{var_ind,4});% add data set
var_name  = [WOA_info{var_ind,2},'_an']; % variable to extract

annual_flag = 0;
for i = 1:12 % STEP THROUGH MONTHLY CLIMATOLOGIES
    fname = regexprep(fn_str,'NN',num2str(i,'%02.0f')); % add month
    ds    = ncdataset([data_path,fname]); % Create netCDF object   
    
    % GET EXTRACTION INDICES
    if i == 1 %
        z    = double(ds.data('depth'));
        lat  = double(ds.data('lat'));
        lon  = double(ds.data('lon'));
        
        z_ind1 = find(z <= depth_bnds(1),1,'last');
        z_ind2 = find(z >= depth_bnds(2),1,'first');
        if isempty(z_ind2)
            z_ind2 = length(z);
            disp(['Seasonal data set max depth < requested max depth(', ...
                num2str(depth_bnds(2)),' m). Extracting seasonal data to ', ...
                num2str(z(z_ind2)),'m.','Deeper data will be extracted ',...
                ' annual data set.']);
            annual_flag = 1;
        end
        Z      = z(z_ind1:z_ind2); % lattitude subset

        lat_ind1 = find(lat <= lat_bnds(1),1,'last');
        lat_ind2 = find(lat >= lat_bnds(2),1,'first');
        LAT      = lat(lat_ind1:lat_ind2); % lattitude subset
        
        x_track = track(:,3); % for interpolating
        if cross180 == 0
            lon_ind1 = find(lon <= lon_bnds(1),1,'last');
            lon_ind2 = find(lon >= lon_bnds(2),1,'first');
            LON      = lon(lon_ind1:lon_ind2); % longitude subset
            x        = LON; % for interpolation
        else % Crosses -180 / +180
            lon_ind1 = find(lon <= lon_bnds(2),1,'last');
            lon_ind2 = find(lon >= lon_bnds(1),1,'first');
            nlon     = length(lon);
            LON      = lon([lon_ind1:nlon,1:lon_ind2]);
            t1       = LON < 0;
            x        = LON + (t1*360); % put in 0 to 360
            t1      = x_track < 0;
            x_track = x_track + t1*360;
        end
        
        % PREDIMMENSION 4D SUBSET MATRIX [TIME x DEPTH x LAT x LON]
        d = ones(12,length(Z),length(LAT),length(LON)); 
        disp(['Extracting and merging WOA2013 monthly climatologies', ...
              ' for ', WOA_info{var_ind,1},' .......takes some time.....']);
    end
  
    % EXTRACT SUBSET
    disp(['     Extracting subset for ', month_label{i}, ' .....'])
   
    if cross180 == 0
        d(i,:,:,:) = squeeze(double(ds.data(var_name,...
            [1, z_ind1, lat_ind1, lon_ind1],...
            [1, z_ind2, lat_ind2, lon_ind2])));
    else % CROSSES -180 / +180, do to data grabs and combine
        d1 = squeeze(double(ds.data(var_name,...
            [1, z_ind1, lat_ind1, lon_ind1],...
            [1, z_ind2, lat_ind2, nlon])));  
        
        d2 = squeeze(double(ds.data(var_name,...
            [1, z_ind1, lat_ind1, 1],...
            [1, z_ind2, lat_ind2, lon_ind2])));
        
        d(i,:,:,:) = cat(3, d1, d2); % join along lon dimension
    end
end

if annual_flag == 1 % GET DEEP ANNUAL NUTRIENT VALUES
    fn_str     = 'woa13_XXX_VNN_01.nc';
    fn_str     = regexprep(fn_str,'V',WOA_info{var_ind,2}); % build name: add var
    fn_str     = regexprep(fn_str,'XXX',WOA_info{var_ind,4});% add data set
    var_name  = [WOA_info{var_ind,2},'_an']; % variable to extract
    fname = regexprep(fn_str,'NN','00'); % add month
    ds    = ncdataset([data_path,fname]); % Create netCDF object
    
    z    = double(ds.data('depth'));
    lat  = double(ds.data('lat'));
    lon  = double(ds.data('lon'));
    z_ind3 = find(z > Z(end),1,'first');
    z_ind4 = find(z >= depth_bnds(2),1,'first');
    
    Z1      = z(z_ind3:z_ind4); % lattitude subset
    % CAN REUSE LAT AND LON INDICES
    
    % PREDIMMENSION 4D SUBSET MATRIX [TIME x DEPTH x LAT x LON]
    d3 = ones(12,length(Z1),length(LAT),length(LON));
    disp(['Extracting WOA2013 annual climatology', ...
        ' for ', WOA_info{var_ind,1},' .......takes some time.....']);
    
    % EXTRACT SUBSET
    if cross180 == 0
        d3(1,:,:,:) = squeeze(double(ds.data(var_name,...
            [1, z_ind3, lat_ind1, lon_ind1],...
            [1, z_ind4, lat_ind2, lon_ind2])));
    else % CROSSES -180 / +180, do to data grabs and combine
        d4 = squeeze(double(ds.data(var_name,...
            [1, z_ind3, lat_ind1, lon_ind1],...
            [1, z_ind4, lat_ind2, nlon])));  
        
        d5 = squeeze(double(ds.data(var_name,...
            [1, z_ind3, lat_ind1, 1],...
            [1, z_ind4, lat_ind2, lon_ind2])));
        
        d3(1,:,:,:) = cat(3, d4, d5); % join along lon dimension
    end
    
    % NOW REPLICATE IN 4D TO MORGE WITH SEASONAL VALUES
    for  i = 2:12
        d3(i,:,:,:) = d3(1,:,:,:);
    end
    d = cat(2, d, d3);  % Merge shallow seasonal with deep annual data
    Z = [Z;Z1];
end

       
clear z_ind1 z_ind2 lat_ind1 lat_ind2 lon_ind1 lon_ind2 t1 nlon
clear lat lon z 

% ***********************************************************************
% ***********************************************************************
% NOW INTERPOLATE WOA2013 ALONG THE FLOAT TRACK
% ***********************************************************************
% ***********************************************************************
stations = size(track(:,1),1); % get # of points to interpolate
% GET YEAR DAY FOR EACH POSITION
[Y, ~, ~, ~, ~, ~] = datevec(track(:,1)); % year for each position
yrday    = track(:,1) - datenum(Y,1,1);   % year day for each position
WOA_day  = (datenum(2014,1:12,15))' - datenum(2014,1,0); % month centered

d_interp = ones(length(Z), stations) * NaN; % predim final output matrix

for i = 1: stations
    if isnan(track(i,2)) % NO position (under ice, lost comms, etc)
        continue
    end
    % ********************************************************************
    % GET BOUNDING INDICES, INTERP WEIGHTS & SUBSET BEFORE INTERPOLATION
    % ********************************************************************
    t1 = find(WOA_day < yrday(i),1,'last'); % bound float time
    if isempty(t1) % early Jan so jday >1 & < 15
        t1 = 1;
        t2 = 12;
        wt = (yrday(i) + 15) / 30; % 2/11/15 jp 15~100%, 1~50:
    elseif t1 == 12 % late dec so jday >350 & < 365
        t2 = 1;
        wt   = (365 - yrday(i) + 15)/ 30; % 2/11/15 jp 350~100%, 365~50%
    else
        t2 = t1+1;
        dx   = WOA_day(t2) - WOA_day(t1); % time between grid points
        dx1  = yrday(i) - WOA_day(t1);
        wt   = (dx-dx1)./dx;
    end
    
    lat1 = find(LAT < track(i,2),1,'last'); % bound float lat
    lat2 = lat1 + 1;
    dx   = LAT(lat2) - LAT(lat1); % degrees lat between grid points
    dx1  = track(i,2) - LAT(lat1);
    lat_wt  = (dx-dx1)./dx; % lat weight
    
    lon1 = find(x < x_track(i), 1, 'last'); % bound float lat
    lon2 = lon1 + 1;
    dx   = x(lon2) - x(lon1); % degrees lat between grid points
    dx1  = x_track(i) - x(lon1);
    lon_wt  = (dx-dx1)./dx; % lat weight
    
    % ********************************************************************
    % INTERPOLATION PART I => GET BOUNDING PROFILES
    % 4d to 3d => [nz x 2 x 2]
    % [TIME x DEPTH x LAT x LON] => [DEPTH x LAT x LON] at time t
    D3 = squeeze(d(t1,:,lat1:lat2,lon1:lon2)*wt + ...
        d(t2,:,lat1:lat2,lon1:lon2)*(1-wt)); % 4d to 3d
    
    D_plot = squeeze(d(t1,:,:,:)*wt + ...
        d(t2,:,:,:)*(1-wt)); % FOR testing
    %*********************************************************************
    % CHECK FOR NAN PROFILES, IF EXIST DON'T INTERPOLATE - USE SIMPLE
    % AVERAGE OF GOOD PROFILES
    
    t_NaN = squeeze(sum(isnan(D3),1) > 0);% NaN's in bounding profiles? 2x2
    if sum(t_NaN(:)) ~= 0
        disp(['Boundng Climate profile(s) missing data', ...
            ' - Taking simple average of complete profiles']);
        d_interp(:,i) = mean(D3(:,~t_NaN),2); % avg non NaN profiles
        continue % done move to next position
    end
    
    % ********************************************************************
    % INTERPOLATION PART II => 3d to 2d
    % [DEPTH x LAT x LON] => [DEPTH x LON] at given lat
    D2 = squeeze(D3(:,1,:)*lat_wt + D3(:,2,:)*(1-lat_wt)); % 3d to 2d
    
    % ********************************************************************
    % PART III => 2d surface to 1d profile => depth at a given lat & lon
    % [DEPTH x LON] => Depth at a given lon
    d_interp(:,i)  = (D2(:,1)*lon_wt + D2(:,2)*(1-lon_wt)); % 2d to profile
    
    
    % ********************************************************************
    % FOR TESTING
    if plot_it == 1
        figure(1)
        
        subplot(2,4,1:2)
        p_data1 = squeeze(D_plot(1,:,:));
        levs = 40;
        cmin  = min(d(:));             % MIN MAX
        cmax  = max(d(:));
        clevs = cmin:(cmax-cmin)/levs:cmax;
        
        contourf(x,LAT,p_data1,clevs,'linecolor','none');
        xlim([min(x) max(x)])
        ylim([min(LAT) max(LAT)])
        colorbar
        hold on
        plot(x_track(1:i),track(1:i,2),'-ko','MarkerSize',6, ...
            'MarkerFaceColor','k','LineWidth',1.5)
        plot(x_track(i),track(i,2),'-ko','MarkerSize',8, ...
            'MarkerFaceColor','w','LineWidth',1.5)
        hold off
        
        subplot(2,4,5:6)
        p_data2 = squeeze(D3(1,:,:));
        levs = 40;
        cmin  = min(p_data2(:));             % MIN MAX
        cmax  = max(p_data2(:));
        clevs = cmin:(cmax-cmin)/levs:cmax;
        
        %subplot(2,1,2)
        contourf(x(lon1:lon2),LAT(lat1:lat2),p_data2,clevs, ...
            'linecolor','none');
        xlim([x(lon1) x(lon2)])
        ylim([LAT(lat1) LAT(lat2)])
        colorbar
        hold on
        scatter(x_track(i),track(i,2),64,d_interp(1,i),'filled',...
            'MarkerEdgeColor', 'k','linewidth', 2)
        hold off
        %     disp(squeeze(D3(1,:,:)))
        %     disp(d_interp(1,i))
        
        
        subplot(2,4,[3,4,7,8])
        bnd_profs = [D3(:,1,1) D3(:,1,2) D3(:,2,1) D3(:,2,2)];
        plot(d_interp(:,i),Z, 'k--')
        hold on
        plot(bnd_profs, Z, 'r-')
        set(gca,'YDir','Reverse')
        hold off
        
        pause
    end
    clear D3 D2 wt lat_wt lon_wt t_NaN t1 t2 lat1 lat2 lon1 lon2 dx dx1
end
% FUNCTION OUTPUT
data.d = d_interp;
data.Z = Z;
clearvars -except data % don't need d anymore   

























