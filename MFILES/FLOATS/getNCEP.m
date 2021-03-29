function NCEP = getNCEP(SDN,LON,LAT,dirs)
%
% NCEP = getNCP(SDN,lon,lat)
%   getNCEP extracts NCEP data along a time track given time, lon & lat
%   This is meant for float track data but "should" work for any data set
%   Inputs must be arrays and all the same size
%   Change the NCEP variables extracted by editing the cell array "NCEPname"
%   NCEP  6 hourly data sets at:
%      http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
%      README:  ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/README
% 
% INPUTS:
%    SDN = Matlab serial date number
%    lon = longitude (0<->360 or -180<->180)
%    lat = lattitude (-90<->90)
%
% OUTPUT:
%    NCEP - A structure containing NCEP SDN and climate variables
%    
% 1) nctoolbox needs to be installed
% 2) Setup_nctoolbox.m needs to be executed
%
% 03/22/2017: TM added waitbar to year loop.

% CHANGE HISTORY
% 04/07/2017 Fixed bug during 0 meridian crossing to correctly get
%   longitude data bounds. I was grabbing the min on the 0 to 180 side
%   of the crossing - now grabbing the max. ~ line 129
% 04/07/2017 Added path to extract netcdf files locally. Had to add code
%   to use matlab functions to get time stamp - nctoolbox could not find any
%   attributes (????) in the netcdf files copied directly from: 
%   https://www.esrl.noaa.gov/psd/cgi-bin/db_search/ ...
%   DBListFiles.pl?did=192&tid=58359&vid=28
% 01/10/2018 Modified code to utilize built-in MATLAB netcdf tools; ran
% into a problem with nctoolbox's ability to resolve the appropriate
% openDAP paths for NCEP surface pres data.  Went back and forth via email
% with esrl representatives who did not have an explanation.  Tested
% various installations of nctoolbox as well (per Brian Schlining) which
% was not the issue.

% *************************************************************************
% TEST
% d = get_FloatViz_data(['C:\Users\jplant\Documents\MATLAB\', ...
%     'ARGO_PROCESSING\DATA\FLOATVIZ\9094SOOCN.TXT']);
% % d = get_FloatViz_data(['C:\Users\jplant\Documents\MATLAB\', ...
% %     'ARGO_PROCESSING\DATA\FLOATVIZ\5143STNP.TXT']);
% [~,ia,~] = unique(d.data(:,2));
% track = d.data(ia,[1,4,3]); % reverse lat lon order for function
% t1 = track(:,2) == -1e10;
% track(t1,2:3) = NaN;
% track(t1,:) = [];
% load('C:\temp\test_track.mat')
% track = Wtrack;
% SDN = track(:,1);
% LAT = track(:,2);
% LON = track(:,3);
% *************************************************************************



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

%CHECK FOR  NAN's IN LAT OR LON
% t1 = isnan(LATx) | isnan(LONx);
% no_pos_sdn = [SDNx, t1];
% if sum(t1) > 0
%     disp(['NaN found in Lat or LON. No NCEP data will be returned for ',...
%         ' these time points. Consider cleaning up input data first!'])
%     SDNx = SDNx(t1);
%     LONx = LONx(t1);
%     LATx = LATx(t1);
% end
% clear t1

% ************************************************************************
% ****************************  PATHS and VARIABLES **********************
%        NCEPPath can be a URL, network dir (ie CHEM) or a local dir
% ************************************************************************
% NCEPpath1 = dirs.NCEP_TEMP;
NCEPpath1   = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/';
NCEPpath2   = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/';

          
%   (1)matlab var name, (2)NCEP file,    (3)file path 
NCEPname(1,:) ={'PRES' ,'pres.sfc.gauss.',NCEPpath1 };   % Surface  pressure, Pascals
%NCEPname(2,:) ={'RH'  ,'rhum.sig995.',    NCEPpath2};    % Relative humidity (pH2O/pH2O sat)

%NCEPname(2,:) ={'PRES2','pres.sfc.',     NCEPpath2};    % Surface  pressure, Pascals
landmask = [dirs.NCEP_TEMP,filesep,'land.sfc.gauss.nc'];


[num_vars,~] = size(NCEPname); % get number of variables to extract

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
    LON_dir = sign(lon_diff(T_lon)); % 
    disp(['LON_flag = ',num2str(LON_flag),' Float transits across 0 Meridian.']);
end

% ************************************************************************
% ************************************************************************
% BEGIN LOOPS yearfile and file name loops HERE 4/24/14 -jp
%              (1)matlab var, (2)NCEP file, (3)file path
% for fn = 1: num_vars                   % loop through variable names
    fn=1
    NCEPdata =[];
    NCEPtime =[];
    waitmsg = sprintf(['Please Wait. \nGrabbing ',num2str(end_yr-start_yr+1),...
        ' years of NCEP data from the web. ',...
        '\nSpeed depends on connectivity strength.']);
    WB = waitbar(0,waitmsg); 
    for yr = start_yr:end_yr           % loop through each year-files
        NCEPtarget =[NCEPname{fn,3},NCEPname{fn,2},num2str(yr),'.nc'];

%         %%%01/12/2018.  nctoolbox no longer works with this opendap
%         server!!!  Not sure why...went back and forth with numerous reps.
%          Matlab's netcdf tools seem to work; use those for now.%%%

        disp(['Extracting data from: ',NCEPname{fn,2},num2str(yr),'.nc']);
        waitbar((yr-start_yr+1)/(end_yr-start_yr+1),WB)
        % *****************************************************************
        %        DETERMINE LAT & LON INDICES FOR DATA SUBSET GRABS
        %       ONLY GET EXTRACTION INDICES 1st TIME VAR ENCOUNTERED
        %   GRID WILL BE THE SAME FOR ALL YEAR_FILES FOR A GIVEN VARIABLE
        % *****************************************************************
        
        if yr == start_yr
            lat = ncread(NCEPtarget,'lat');
            lon = ncread(NCEPtarget,'lon');
            lat_ind = [find(lat>LAT_bnds(1),1,'last'),...
                find(lat<LAT_bnds(2),1,'first')];
            lat     = lat(lat_ind(1):lat_ind(2)); % subsetted lat variable
            
            if LON_flag ==1 % NO FLOAT CROSSING
                lon_ind = [find(lon<LON_bnds(1),1,'last'),...
                    find(lon>LON_bnds(2),1,'first')];
                
                if length(lon_ind) == 2
                    lon_1 = lon(lon_ind(1):lon_ind(2)); % subset lon
                else % float doesn't cross meridian but NCEP data grab does
                    LON_flag = 0; % RESET FLAG TO ZERO
                    if LON_bnds(1)< min(lon)
                        LON_bnds =[LON_bnds(2) max(lon)]; %low pt = problem
                    else
                        LON_bnds =[min(lon) LON_bnds(1)]; %high pt = problem
                    end
                    disp('NCEP data grab transits 0 Meridian.');
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
        mytime = ncread(NCEPtarget,'time');
        t = mytime/24 + datenum(1800,1,1,0,0,0);
        tlength = length(t);
%         t_ind = [find(t>nanmin(SDNx),1,'first'),find(t<nanmax(SDNx),1,'last')];
        
        if LON_flag ==1                     % No 0 meridian crossing data
            d = ncread(NCEPtarget,'pres', [lon_ind(1) lat_ind(1) 1],[lon_ind(2)-lon_ind(1)+1 lat_ind(2)-lat_ind(1)+1 tlength]);
            LM = ncread(landmask,'land', [lon_ind(1) lat_ind(1) 1],[lon_ind(2)-lon_ind(1)+1 lat_ind(2)-lat_ind(1)+1 1]);  
            LMask = logical(LM);
            for ilm = 1:tlength
                Dtmp = d(:,:,ilm);
                Dtmp(LMask) = nan;
                d(:,:,ilm) = Dtmp;
            end
        else
            disp('PLOT TO VERIFY MERIDIAN CROSSING CODE!!!!');
            d1 = ncread(NCEPtarget,'pres', [1 lat_ind(1) 1],[lon_ind(2)-1 lat_ind(2)-lat_ind(1)+1 tlength]);
            LM1 = ncread(landmask,'land', [1 lat_ind(1) 1],[lon_ind(2)-1 lat_ind(2)-lat_ind(1)+1 1]);
            LMask1 = logical(LM1);
            for ilm = 1:tlength
                Dtmp = d1(:,:,ilm);
                Dtmp(LMask1) = nan;
                d1(:,:,ilm) = Dtmp;
            end
%             d2 = ncread(NCEPtarget,'pres', [lon_ind(2) lat_ind(1) 1],[length(lon)-lon_ind(2) lat_ind(2)-lat_ind(1)+1 tlength]);
            d2 = ncread(NCEPtarget,'pres', [lon_ind(2) lat_ind(1) 1],[length(lon)-lon_ind(2)+1 lat_ind(2)-lat_ind(1)+1 tlength]);
            LM2 = ncread(landmask,'land', [lon_ind(2) lat_ind(1) 1],[length(lon)-lon_ind(2)+1 lat_ind(2)-lat_ind(1)+1 1]);
            LMask2 = logical(LM2);
            for ilm = 1:tlength
                Dtmp = d2(:,:,ilm);
                Dtmp(LMask2) = nan;
                d2(:,:,ilm) = Dtmp;
            end
            d = cat(1,d1,d2);          
            clear d1 d2
        end
        NCEPtime = [NCEPtime; t];
%         NCEPdata = [NCEPdata; d']; %( time x lat x lon)
        NCEPdata = cat(3,NCEPdata,d); %( lon x lat x time)
    end
    %TRANSPOSE TO ( time x lat x lon)
    NCEPdata = permute(NCEPdata,[3,2,1]);

    % ************************************************************************
    % ***********  EXTRACT NCEP CLIMATE DATA OVER FLOAT TRACK  ***************
    % ************************************************************************
    % Interpolate float lat & lon onto NCEP time grid
    % Could just interpolate on to inputs lat & lon but may want more
    % detailed NCEP data time wise to play with water age, time lags and O2
    % calcs
    
    if LON_flag == 0 % O MERIDIAN CROSSING - CONVERT LON to -180 to 180
        LONx(LONx >180) = LONx(LONx >180) - 360;
        lon_1(lon_1 >180) = lon_1(lon_1 >180) - 360;
        [lon_1,IX] = sort(lon_1);
        NCEPdata =NCEPdata(:,:,IX);
    end 
        
%     LAT_i = interp1([fix(SDNx(1)); SDNx], [LATx(1); LATx], NCEPtime);
%     LON_i = interp1([fix(SDNx(1)); SDNx], [LONx(1); LONx], NCEPtime);
    
% Tanya fix  031417
    LAT_i = interp1([fix(SDNx(1)); SDNx; ceil(SDNx(end))], ...
            [LATx(1); LATx; LATx(end)], NCEPtime);
    LON_i = interp1([fix(SDNx(1)); SDNx; ceil(SDNx(end))], ...
            [LONx(1); LONx; LONx(end)], NCEPtime);

    NO_NaNs = ~isnan(LAT_i);  
    
    LAT_i = LAT_i(NO_NaNs);
    LON_i = LON_i(NO_NaNs);
    NCEP_SDN = NCEPtime(NO_NaNs);
   
    
    % LOOP THROUGH AND EXTRACT NCEP DATA
    NCEPdata_pt =[];
% OLD LOOP - NO INTERPOLATION. COMPARED 2 DIFFERNET NCEP PRESSURE GRIDS
% DIFFERENCE UP TO 15 so BETTER INTERPOLATE
%     for i = 1: length(NCEP_SDN)
%         time_dt = abs(NCEPtime - NCEP_SDN(i)); % time range - time of data pt
%         time_pt = find(time_dt == min(time_dt)); % time point logical index
%         lat_dt = abs(lat - LAT_i(i));
%         lat_pt = find(lat_dt == min(lat_dt)); % lat point index
%         lon_dt = abs(lon_1 - LON_i(i));
%         lon_pt = find(lon_dt == min(lon_dt)); % lon point index
%         
%         % CHOOSE the first ind if there are 2 (1/2 way bewtween)
%         NCEPdata_pt =[NCEPdata_pt;NCEPdata(time_pt(1),lat_pt(1),lon_pt(1))]; %extract
%     end

% BRUTE FORCE LINEAR INTERP - STEP DOWN THROUGH DIMMENSIONS
     for i = 1: length(NCEP_SDN)
        % 3D to 2D
        time_pt = find(NCEPtime >= NCEP_SDN(i),1 ,'first'); % time point upper bound or =
        if time_pt == 1
            wt = 0;
            tmp2 = squeeze(NCEPdata(time_pt, :, :)); % get upper time surfaces
            tmp_s  = tmp2*(1-wt); % surface
        else
            wt = (NCEPtime(time_pt) - NCEP_SDN(i)) /...
                 (NCEPtime(time_pt) - NCEPtime(time_pt-1));
            tmp1 = squeeze(NCEPdata(time_pt-1, :, :)); % get lower time surfaces
            tmp2 = squeeze(NCEPdata(time_pt, :, :)); % get upper time surfaces
            tmp_s  = tmp1* wt + tmp2*(1-wt); % surface
        end

        
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
        if isempty(tmp_pt), pause, end
%         tmp_pt
        NCEPdata_pt =[NCEPdata_pt;tmp_pt]; % add to array
%         if isempty(tmp_pt),pause, end
    end
       
    
    
    
    % NOW INTERPOLATE BACK TO INPUT SDN
    out = interp1(NCEP_SDN,NCEPdata_pt,SDN);
    
    %plot(flt_LONG,flt_LAT, 'bo-', jptest(:,3), jptest(:,2),'r*-')% for testing
    %pause% for testing
    
    s1 =['NCEP.',NCEPname{fn,1}, '=out;'];
    eval(s1);
    %whos NCEP.PRES
    
    % BUILD DATA STRUCTURE

    
    disp(['NCEP variable ',NCEPname{fn,1}, ' created']);
    disp(' ');
% end
close(WB)

