function NCEP = getNCEPtest(SDN,LON,LAT,dirs)
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
% 05/01/21 TM Clean-up of original code (interpolation, prime meridion crossing, etc).  

% *************************************************************************

NPRES = nan(length(SDN),1);

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
% ****************************  PATHS and VARIABLES **********************
%        NCEPPath can be a URL, network dir (ie CHEM) or a local dir
% ************************************************************************
%NCEPpath1   = 'C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\DATA\NCEP_TEMPORARY\';
NCEPpath1 = dirs.NCEP_TEMP;
NCEPpath1   = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss/';
% NCEPpath1   = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis2/gaussian_grid/';
%NCEPpath2   = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/';
landmask = [dirs.NCEP_TEMP,filesep,'land.sfc.gauss.nc'];
          
%   (1)matlab var name, (2)NCEP file,    (3)file path 
NCEPname(1,:) ={'PRES' ,'pres.sfc.gauss.',NCEPpath1 };   % Surface  pressure, Pascals
%NCEPname(2,:) ={'RH'  ,'rhum.sig995.',    NCEPpath2};    % Relative humidity (pH2O/pH2O sat)

%NCEPname(2,:) ={'PRES2','pres.sfc.',     NCEPpath2};    % Surface  pressure, Pascals


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

% LAT_bnds = [max(LATx) min(LATx)];   % lat range of data 
% LON_bnds = [min(LONx) max(LONx)]; % lon range of data

% % % % ************************************************************************
% % % %               ***  CHECK FOR 0 MERIDIAN CROSSING  ***
% % % %  Assume no more than 20º of travel In 14 days  so crossing should show
% % % %                 an abs value in the long diff > 340
% % % % ************************************************************************
% % % lon_diff = diff(LONx);
% % % T_lon = find(abs(lon_diff) > 340,1,'first'); % find 1st meridian crossing
% % % if isempty(T_lon)
% % %     LON_flag = 1;   % get lon >= min & lon <= max
% % %     disp(['LON_flag = ',num2str(LON_flag),' Float does not cross 0 Meridian.']);
% % % else
% % %     LON_flag = 0;   % get lon < min & lon > max
% % %     t1 = (LONx >180);
% % %     LON_bnds = [max(LONx(~t1)) min(LONx(t1))]; % re-due bounds jp 4/7/17
% % %     
% % %     % NORMALLY SIGN OF DIFF: E =1, W= -1, HOW EVER AT MERIDIAN CROSSING
% % %     % IT IS THE OPPOSITE E= -1 (ex 1-359 =-358), W = 1 (ex 359-1 = 358)
% % %     LON_dir = sign(lon_diff(T_lon)); % 
% % %     disp(['LON_flag = ',num2str(LON_flag),' Float transits across 0 Meridian.']);
% % % end

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
        
%         if yr == start_yr
            lat = ncread(NCEPtarget,'lat');
            lon = ncread(NCEPtarget,'lon');
            
           % max lon in ncep is ~358 deg.  That means float track with
           % longitude between 358-360 cannot be interpolated to NCEP.
           % So...try making a replicate datapoint of lon = 0 for lon = 360
           % (as these are equivalent), so can trick the code for the
           % interpolation.  Only do this if lon crosses merid
           if nanmax(lon)>358;
            lon = [lon;360];
           end
            time = ncread(NCEPtarget,'time');
            t = time./24 + datenum(1800,1,1,0,0,0);
            tlength = length(t);
            
            d = ncread(NCEPtarget,'pres');
%             mytmp(1,:,:) = d(1,:,:);
%             D = d;
            d(193,:,:) = d(1,:,:);
            

            LM = ncread(landmask,'land');  
            LM(193,:) = LM(1,:);
            LMask = logical(LM);
            for ilm = 1:tlength
                Dtmp = d(:,:,ilm);
                Dtmp(LMask) = nan;
                d(:,:,ilm) = Dtmp;
            end
           [X,Y,Z] = meshgrid(double(lat),double(lon),t);
           Vq = interp3(X,Y,Z,double(d),LATx,LONx,SDNx);
           xx = find(~isnan(Vq));
           ntmp = Vq(xx);
           NPRES(xx) = ntmp;
           
    end
    NCEP.PRES = NPRES;
    close(WB)

