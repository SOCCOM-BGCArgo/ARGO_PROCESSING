function ERA = getERA_NRT_sO2(inputdir,SDN,LON,LAT)
% ************************************************************************
% getERA_NRT_sO2.m
% ************************************************************************
%
% ERA = getERA_NRT_sO2(inputdir,SDN,lon,lat)
%   getERA extracts ERA-5 reanalysis data along a given time track
%   (time, lat, lon).
%
%   Native ECMWF grid: 0.25 x 0.25 degree
%   nc sp variable stored as lon x lat x time (1440 x 720 x 1460)
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
% DATE: 12/28/16 (originally created for testing ERAinterim, and modified from Josh's original getNCEP routine)
% UPDATES: 4/5/22: fixed weighting bug affecting certain floats where lat
%                  index was missing for specific ERA time points (I think
%                  due to incomplete current year's annual file)
%          6/15/23: TM, fixed another indexing bug related to NRT data
%                   processing ("expvar" variable is added for most recent year, and
%                   creates a fourth dimension in the nc file that doesn't exist in
%                   earlier years.  I thought I had fixed this already...!
%                   See 
%		   8/8/23: TM Grrr...the ERA Near Real Time (NRT) mixed variable nc file is a real pain.  My previous attempts to merge the QA with NRT had failed. 
%					Now opting for a stand-alone piece of code that solely addresses & extracts the expver=2 dimension for the NRT data.  The wrapper-code (Process_ERA_FloatRefs) 
%					will then do the merging.
% NOTES: 
% ************************************************************************
%
% ************************************************************************
NOW = datevec(now);
% THIS CODE IS ONLY FOR CURRENT YEAR!  SO, PARE DOWN THE SDN TO CURRENT
% YEAR.
TMPERA = NaN(length(SDN),1);
inputSDN = SDN;
inputLAT = LAT;
inputLON = LON;
xxx= find(inputSDN>=datenum(NOW(1),01,01));
SDN = inputSDN(xxx);
LON = inputLON(xxx);
LAT = inputLAT(xxx);

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
if yr == NOW(1) %Current year will have an extra dimension ("expvar") which is length 2, and identifies the QC'd or NRT ERA.  We want to ignore this variable; the NRT is acceptable!
    if LON_flag ==1                     % No 0 meridian crossing data
        d = double(ds.data(dsv, [1 2 lat_ind(1) lon_ind(1)], ...
            [tlength 2 lat_ind(2) lon_ind(2)]));
        d = squeeze(d);
    else
        disp('PLOT TO VERIFY MERIDIAN CROSSING CODE!!!!');
        d1 = double(ds.data(dsv, [1 2 lat_ind(1) 1], ...
            [tlength 2 lat_ind(2) lon_ind(1)]));
        d2 = double(ds.data(dsv, [1 2 lat_ind(1) lon_ind(2)], ...
            [tlength 2 lat_ind(2) length(lon)]));
        d = cat(3,d1,d2); 
        d = squeeze(d);
        clear d1 d2
    end
%     if LON_flag ==1                     % No 0 meridian crossing data
%         dd = double(ds.data(dsv, [1 2 lat_ind(1) lon_ind(1)], ...
%             [tlength 2 lat_ind(2) lon_ind(2)]));
%         dd = squeeze(dd);
%     else
%         disp('PLOT TO VERIFY MERIDIAN CROSSING CODE!!!!');
%         dd1 = double(ds.data(dsv, [1 2 lat_ind(1) 1], ...
%             [tlength 2 lat_ind(2) lon_ind(1)]));
%         dd2 = double(ds.data(dsv, [1 2 lat_ind(1) lon_ind(2)], ...
%             [tlength 2 lat_ind(2) length(lon)]));
%         dd = cat(3,dd1,dd2);
%         dd = squeeze(dd);
%         clear dd1 dd2
%     end
% %     d = [d;dd];
% else
%     if LON_flag ==1                     % No 0 meridian crossing data
%         d = double(ds.data(dsv, [1 lat_ind(1) lon_ind(1)], ...
%             [tlength lat_ind(2) lon_ind(2)]));
%     else
%         disp('PLOT TO VERIFY MERIDIAN CROSSING CODE!!!!');
%         d1 = double(ds.data(dsv, [1 lat_ind(1) 1], ...
%             [tlength lat_ind(2) lon_ind(1)]));
%         d2 = double(ds.data(dsv, [1 lat_ind(1) lon_ind(2)], ...
%             [tlength lat_ind(2) length(lon)]));
%         d = cat(3,d1,d2);          
%         clear d1 d2
%     end
    ERAtime = [ERAtime; t];
    ERAdata = [ERAdata; d]; %( time x lat x lon)
end

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
    ERAdata_pt = nan(10000,1);
    kk=1;
    bad_inds=[];

% BRUTE FORCE LINEAR INTERP - STEP DOWN THROUGH DIMMENSIONS
% keyboard 
     for i = 1: length(ERA_SDN)
        % 3D to 2D
        time_pt = find(ERAtime >= ERA_SDN(i),1 ,'first'); % time point upper bound or =
        if time_pt == 1
            wt=0;
            tmp1=0;
                    tmp2 = squeeze(ERAdata(time_pt, :, :)); % get upper time surfaces

        else
             wt = (ERAtime(time_pt) - ERA_SDN(i)) /...
             (ERAtime(time_pt) - ERAtime(time_pt-1));
                     tmp1 = squeeze(ERAdata(time_pt-1, :, :)); % get lower time surfaces
        tmp2 = squeeze(ERAdata(time_pt, :, :)); % get upper time surfaces
        end

        tmp_s  = tmp1* wt + tmp2*(1-wt); % surface
        
        clear tmp1 tmp2

        % 2D to 1D 90 to -90
        lat_pt = find(lat >= LAT_i(i),1 ,'last'); % 1st index , but upper bound or =
        wt = (lat(lat_pt) - LAT_i(i)) / (lat(lat_pt) - lat(lat_pt+1)); 
        if isempty(lat_pt)
            bad_inds = [bad_inds;i];
            continue
        end 
        tmp_l = tmp_s(lat_pt+1,:)*wt + tmp_s(lat_pt,:)*(1-wt); % line
        
        % 1D to POINT                       
        lon_pt = find(lon_1 >= LON_i(i),1 ,'first'); % lon point upper bound or =
        wt = (lon_1(lon_pt) - LON_i(i)) / (lon_1(lon_pt) - lon_1(lon_pt-1)); 
        if isempty(lon_pt)
            bad_inds = [bad_inds;i];
            continue
        end 
        tmp_pt = tmp_l(lon_pt-1)*wt + tmp_l(lon_pt)*(1-wt); % line
        ERAdata_pt(kk,1) = tmp_pt;
        kk=kk+1;
%         ERAdata_pt =[ERAdata_pt;tmp_pt]; % add to array
     end
    
     ERA_SDN(bad_inds)=[];
     ERAdata_pt(kk:end,:)=[]; %trim array
    
    % NOW INTERPOLATE BACK TO INPUT SDN
    out = interp1(ERA_SDN,ERAdata_pt,SDN);  

    % NOW, WE WANT THE OUTPUT BACK ON THE ORIGINAL TIME VECTOR
    [x,y] = intersect(inputSDN,SDN);  
    TMPERA(y) = out;
    ERA.PRES = TMPERA;

    % BUILD DATA STRUCTURE
    disp('ERA surface pressure variable created');
    disp(' ');
end

