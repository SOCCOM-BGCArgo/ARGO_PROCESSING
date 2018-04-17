function ERA = getERA_sO2Argo(inputdir,SDN,LON,LAT)
% ************************************************************************
% getERA_sO2Argo.m
% ************************************************************************
%
% ERA = getERA_sO2Argo(inputdir,SDN,lon,lat)
%   getERA extracts ERA-Interim reanalysis data along a given time track
%   (time, lat, lon).
%   ECMWF ERA INT 6 hourly data sets at:
%      http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/requests
%
%   Native ECMWF grid: 0.75 x 0.75 degree
%
%
% USE AS:
%
% INPUTS:
%    SDN = Matlab serial date number
%    LON = longitude (0<->360 or -180<->180)
%    LAT = latitude (-90<->90)
%
% OUTPUTS:
%
%
% AUTHOR: Tanya Maurer
%         MBARI
%         tmaurer@mbari.org
%
% DATE: 12/28/16
% UPDATES:
% NOTES: 4/16/18:  THE ERA REFERENCE IS NO LONGER SUPPORTED IN SAGEO2...BUT KEEP
% THIS FUNCTION IN REPOSITORY FOR USERS TO INCORPORATE UPON THEIR OWN DESIRE.
% ************************************************************************
%
% ************************************************************************

% CHECK INPUTS - Set all to columns of data
[r,c] = size(SDN); if r < c, SDN = SDN'; end
[r,c] = size(LON); if r < c, LON = LON'; end
[r,c] = size(LAT); if r < c, LAT = LAT'; end

% MAKE SURE SDN IS MONOTONIC - IF NOT FIX
t1 = diff(SDN) <=0;
% USE "ic" index to rebuild data set later
if sum(t1)> 0;                           % [C,ia,ic] = unique(A,occurrence)
    [SDNx,ia,ic] = unique(SDN,'first'); % SDNf = SDN(ia) and SDN = SDNf(ic)
    LONx = LON(ia);
    LATx = LAT(ia);
    disp([num2str(sum(t1)),...
        ' non-unique points in time removed before NCEP extraction']);
else
    SDNx = SDN;
    LONx = LON;
    LATx = LAT;
end

% ************************************************************************
% ***************  GET START AND END YEARS, LAT & LON BOUNDS  ************
% ********************** Convert lon to 0 to 360 *************************
% ************************************************************************
start_yr = datevec(SDNx(1));
start_yr = start_yr(1); % Get year for first data point

end_yr  = datevec(SDNx(end)); % ending date vec for input time
end_yr = end_yr(1);

t1 = LONx < 0;
LONx(t1) = LONx(t1) + 360; % Convert any - values to 0 to 360E

LAT_bnds = [max(LATx) min(LATx)];   % lat range of data 
LON_bnds = [min(LONx) max(LONx)]; % lon range of data

% ************************************************************************
%               ***  CHECK FOR 0 MERIDIAN CROSSING  ***
%  Assume no more than 20º of travel In 14 days  so crossing should show
%                 an abs value in the long diff > 340
% ************************************************************************
lon_diff = diff(LONx);
T_lon = find(abs(lon_diff) > 340,1,'first'); % find 1st meridian crossing
if isempty(T_lon)
    LON_flag = 1;   % get lon >= min & lon <= max
    disp(['LON_flag = ',num2str(LON_flag),' Float does not cross 0 Meridian.']);
else
    LON_flag = 0;   % get lon < min & lon > max
    t1 = (LONx >180);
    LON_bnds = [max(LONx(~t1)) min(LONx(t1))]; % re-due bounds jp 4/7/17
    
    % NORMALLY SIGN OF DIFF: E =1, W= -1, HOW EVER AT MERIDIAN CROSSING
    % IT IS THE OPPOSITE E= -1 (ex 1-359 =-358), W = 1 (ex 359-1 = 358)
    LON_dir = sign(lon_diff(T_lon)); % -1 = W, +1 =E dir of movement across meridian
    disp(['LON_flag = ',num2str(LON_flag),' Float transits across 0 Meridian.']);
end

% ************************************************************************
%               ***  Loop Through, get data  ***
% ************************************************************************

ERAdata =[];
ERAtime =[];
for yr = start_yr:end_yr           % loop through each year-files
    df = dir([inputdir,'*',num2str(yr),'*.nc']);
    fname = char(df.name);
    if size(fname,1) ~=1
        disp(['WARNING: multiple (or zero) file instances for ',num2str(yr),' ERA analysis fields. SKIPPING YEAR.']);
        continue
    else
        ERAtarget =[inputdir,fname];
        ds  = ncdataset(ERAtarget);   % create data ID for year-file
        disp(['Extracting data from: ',ERAtarget]);
    end
    
            
    % *****************************************************************
    %        DETERMINE LAT & LON INDICES FOR DATA SUBSET GRABS
    %       ONLY GET EXTRACTION INDICES 1st TIME VAR ENCOUNTERED
    %   GRID WILL BE THE SAME FOR ALL YEAR_FILES FOR A GIVEN VARIABLE
    % *****************************************************************
    if yr == start_yr
        dsv     = ds.variables{1}; % Get variable name
        lat     = double(ds.data('latitude')); % get lattitude
        lon     = double(ds.data('longitude')); % get longitude
        lat_ind = [find(lat>LAT_bnds(1),1,'last'),...
            find(lat<LAT_bnds(2),1,'first')];
        lat     = lat(lat_ind(1):lat_ind(2)); % subsetted lat variable

        if LON_flag ==1
            lon_ind = [find(lon<LON_bnds(1),1,'last'),...
                find(lon>LON_bnds(2),1,'first')];

            if length(lon_ind) == 2
                lon_1 = lon(lon_ind(1):lon_ind(2)); % subset lon
            else % float doesn't cross meridian but ERA data grab does
                LON_flag = 0; % RESET FLAG TO ZERO
                if LON_bnds(1)< min(lon)
                    LON_bnds =[LON_bnds(2) max(lon)]; %low pt = problem
                else
                    LON_bnds =[min(lon) LON_bnds(1)]; %high pt = problem
                end
                disp('ERA data grab transits 0 Meridian.');
                disp(['LON_flag reset to: ',num2str(LON_flag)]);
                disp(['LON_bnds reset to: ',num2str(LON_bnds)]);
            end
        end

        if LON_flag ==0
            lon_ind = [find(lon>LON_bnds(1),1,'first'),...
                find(lon<LON_bnds(2),1,'last')];
            lon_1 = [lon(1:lon_ind(1));lon(lon_ind(2):end)];
        end
    end

    % *****************************************************************
    % EXTRACT DATA RANGE FOR VARIABLE
    tlength = ds.size('time');         % get length of time axis
    t       = double(ds.time('time')); % get time data in Matlab SDN

    if LON_flag ==1                     % No 0 meridian crossing data
        d = double(ds.data(dsv, [1 lat_ind(1) lon_ind(1)], ...
            [tlength lat_ind(2) lon_ind(2)]));
    else
        disp('PLOT TO VERIFY MERIDIAN CROSSING CODE!!!!');
        d1 = double(ds.data(dsv, [1 lat_ind(1) 1], ...
            [tlength lat_ind(2) lon_ind(1)]));
        d2 = double(ds.data(dsv, [1 lat_ind(1) lon_ind(2)], ...
            [tlength lat_ind(2) length(lon)]));
        d = cat(3,d1,d2);          
        clear d1 d2
    end
    ERAtime = [ERAtime; t];
    ERAdata = [ERAdata; d]; %( time x lat x lon)
end


% ************************************************************************
% ***********  EXTRACT ERA CLIMATE DATA OVER FLOAT TRACK  ***************
% ************************************************************************
% Interpolate float lat & lon onto ERA time grid
% Could just interpolate on to inputs lat & lon but may want more
% detailed NCEP data time wise to play with water age, time lags and O2
% calcs
    
    if LON_flag == 0 % O MERIDIAN CROSSING - CONVERT LON to -180 to 180
        LONx(LONx >180) = LONx(LONx >180) - 360;
        lon_1(lon_1 >180) = lon_1(lon_1 >180) - 360;
        [lon_1,IX] = sort(lon_1);
        ERAdata =ERAdata(:,:,IX);
    end 
        
    LAT_i = interp1([fix(SDNx(1)); SDNx; ceil(SDNx(end))], [LATx(1); LATx; ceil(LATx(end))], ERAtime);
    LON_i = interp1([fix(SDNx(1)); SDNx; ceil(SDNx(end))], [LONx(1); LONx; ceil(LONx(end))], ERAtime);
    NO_NaNs = ~isnan(LAT_i);  
    
    LAT_i = LAT_i(NO_NaNs);
    LON_i = LON_i(NO_NaNs);
    ERA_SDN = ERAtime(NO_NaNs);
    if isempty(ERA_SDN)
        disp('WARNING: No ERA DATA exists for time period of interest.')
        ERA.PRES = [];
        return
    end

    
    % LOOP THROUGH AND EXTRACT ERA DATA
    ERAdata_pt =[];
    
    %REMNANT FROM JOSH'S getNCEP CODE:::***********************************
% % % OLD LOOP - NO INTERPOLATION. COMPARED 2 DIFFERNET NCEP PRESSURE GRIDS
% % % DIFFERENCE UP TO 15 so BETTER INTERPOLATE
% % %     for i = 1: length(NCEP_SDN)
% % %         time_dt = abs(NCEPtime - NCEP_SDN(i)); % time range - time of data pt
% % %         time_pt = find(time_dt == min(time_dt)); % time point logical index
% % %         lat_dt = abs(lat - LAT_i(i));
% % %         lat_pt = find(lat_dt == min(lat_dt)); % lat point index
% % %         lon_dt = abs(lon_1 - LON_i(i));
% % %         lon_pt = find(lon_dt == min(lon_dt)); % lon point index
% % %         
% % %         % CHOOSE the first ind if there are 2 (1/2 way bewtween)
% % %         NCEPdata_pt =[NCEPdata_pt;NCEPdata(time_pt(1),lat_pt(1),lon_pt(1))]; %extract
% % %     end
    %END, REMNANT FROM JOSH'S getNCEP CODE:::******************************

% BRUTE FORCE LINEAR INTERP - STEP DOWN THROUGH DIMMENSIONS
     for i = 1: length(ERA_SDN)
        % 3D to 2D
        time_pt = find(ERAtime >= ERA_SDN(i),1 ,'first'); % time point upper bound or =
        wt = (ERAtime(time_pt) - ERA_SDN(i)) /...
             (ERAtime(time_pt) - ERAtime(time_pt-1));
        tmp1 = squeeze(ERAdata(time_pt-1, :, :)); % get lower time surfaces
        tmp2 = squeeze(ERAdata(time_pt, :, :)); % get upper time surfaces
        tmp_s  = tmp1* wt + tmp2*(1-wt); % surface
        
%         disp([NCEPtime(time_pt-1) NCEPtime(time_pt) wt ... % for testing
%             (NCEPtime(time_pt-1)*wt + NCEPtime(time_pt)*(1-wt))  NCEP_SDN(i) ])
        
        clear tmp1 tmp2

        % 2D to 1D 90 to -90
        lat_pt = find(lat >= LAT_i(i),1 ,'last'); % 1st index , but upper bound or =
        wt = (lat(lat_pt) - LAT_i(i)) / (lat(lat_pt) - lat(lat_pt+1)); 
        tmp_l = tmp_s(lat_pt+1,:)*wt + tmp_s(lat_pt,:)*(1-wt); % line
        
%         disp([lat(lat_pt+1) lat(lat_pt) wt ... % for testing
%             (lat(lat_pt+1)*wt + lat(lat_pt)*(1-wt)) LAT_i(i)])       
        
        
        % 1D to POINT                       
        lon_pt = find(lon_1 >= LON_i(i),1 ,'first'); % lon point upper bound or =
        wt = (lon_1(lon_pt) - LON_i(i)) / (lon_1(lon_pt) - lon_1(lon_pt-1)); 
        tmp_pt = tmp_l(lon_pt-1)*wt + tmp_l(lon_pt)*(1-wt); % line
        
%         disp([lon_1(lon_pt-1) lon_1(lon_pt) wt ... % for testing
%             (lon_1(lon_pt-1)*wt + lon_1(lon_pt)*(1-wt)) LON_i(i)])    
        %pause
        % CHOOSE the first ind if there are 2 (1/2 way bewtween)
%         if isempty(tmp_pt), pause, end
        ERAdata_pt =[ERAdata_pt;tmp_pt]; % add to array
    end
       

    % NOW INTERPOLATE BACK TO INPUT SDN
    out = interp1(ERA_SDN,ERAdata_pt,SDN);
    
    %plot(flt_LONG,flt_LAT, 'bo-', jptest(:,3), jptest(:,2),'r*-')% for testing
    %pause% for testing
    
    ERA.PRES=out;

    % BUILD DATA STRUCTURE
    disp('ERA surface pressure variable created');
    disp(' ');
end

