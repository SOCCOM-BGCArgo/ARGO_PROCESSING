function Process_ERA5_float_refs(wmo, proc_yr, var)
% ************************************************************************
% Process_ERA5_float_refs.m
% ************************************************************************
% Short script to download and extract ERA5 surface pressure reference data    
% for all floats, all cycle locations.  Data is saved as mat files for  
% easier ingestion into SAGE-O2 (air-cal reference for calibrating the 
% oxygen sensor).
%
% GENERAL FUNCTION PROCEDURE:
%   - Download ERA5 reference data for specified years
%   - Use ERA5data to process specified floats
%   - Save data as mat files 
%
% USE AS: 
% Process_ERA5_float_refs('active', 'current',  {'surface_pressure', 'sp'})
% Process_ERA5_float_refs('all', 'all',  {'surface_pressure', 'sp'})
% Process_ERA5_float_refs({'5906491'} , 'current',  {'surface_pressure', 'sp'})
%
% INPUTS:
%    wmo        - 'all'         : process all flaots 
%                 'active'      : only process active floats (exlude dead)
%                 {'5906491'}   : a single or array of wmo str to process
%                 if wmo is empty, will default to 'all'
%    proc_yr    - 'current': only download current year
%                 'all'   : download all years since 2014
%                 (ie)[2021, 2022]: single or array of years to download
%                 if proc_yr is empty, will skip download and just process 
%    var        - variable name to grab on Copernicus
%                  First column: name on Copernicus server, 
%                  Second column: name of subdirectory in ncroot for the 
%                  respective parameter, and will be used in the filename 
%               - (ie) {'surface_pressure', 'sp'
%                       'mean_sea_level_pressure', 'msl'}; 
%
% OUTPUTS: 
%   The matfile saved includes two structure variables:
%       FLT (with WMO, and also SDN, CYC, LAT, LON for each cycle along the float track)
%       ERA (with PRES = surface pressure at corresponding FLT LAT, LON, SDN.
%       units are in pascals)
%
% SUPPORTING FUNCTIONS:
%    download_ERA5   

% AUTHOR:
%   Tanya Maurer
%   MBARI
%   10/31/22
% UPDATES: SB, 8/16/24   Automated and turned into function.
%--------------------------------------------------------------------------
% ----------------------- CONFIGURE EMAIL ---------------------------------
 
setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
setpref('Internet','E_mail','tmaurer@mbari.org'); % define sender
%email_list ={'tmaurer@mbari.org';'sbartoloni@mbari.org'};
email_list = {'sbartoloni@mbari.org'};
email_msg = {};

tic
kk=1;

% ************************************************************************
% Set Dirs
% user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
% ERA_repo = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\DATA\ERA5\FLOAT_REF_test\'];
% load([user_dir,'\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.mat']);

user_dir = ['\\atlas\Chem\ARGO_PROCESSING\DATA'];
ERA_repo = [user_dir,'\ERA5\FLOAT_REF\'];
load([user_dir,'\CAL\MBARI_float_list.mat']);

iDEAD = find(strcmp('tf Dead',d.hdr) == 1);
iWMO = find(strcmp('WMO',d.hdr) == 1);

% ************************************************************************
% Begin downloading ERA5
% Call fxn download_ERA5.m to get ERA5 data for 'current' or 'all' years.
if ~isempty(proc_yr)
    disp(['Begin ERA5 downloads...'])
    proc_status = download_ERA5(proc_yr, var)
    if proc_status == 0
        email_msg = [email_msg;['ERROR downloading ERA5 for year = ',proc_yr,'!']];
    end
end
% Send notification email:
if isempty(email_msg) == 0
    sendmail(email_list,'ARGOSY: ERA5 DATA DOWNLOAD', email_msg) 
end 

% ************************************************************************
% Begin float processing
% if WMO is not specified, process all wmo's
if isempty(wmo) | strcmp(wmo, 'all') == 1
    disp(['Processing all WMOs..'])
    wmo_list = d.list(:,iWMO);
elseif strcmp(wmo,'active') == 1
    for i = 1:length(d.list)
        if d.list{i,iDEAD} == 0
            wmo_list = d.list(:,iWMO);
        end
    end
else
    wmo_list = wmo;
end

for ii = 1:length(wmo_list)
    WMO = char(wmo_list(ii));

%     thedir = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\',WMO,'\'];
    thedir = [user_dir, '\FLOATS\',WMO,'\'];
    profs = dir([thedir,WMO,'*.mat']);
    thefiles = char(profs.name);
    k=1;
    SDN=[];
    LAT = [];
    LON = [];
    CYC = [];
    for i = 1:size(thefiles,1)
        FN = [thedir,thefiles(i,:)];
        load(FN)
        if isempty(INFO.gps)
            continue
        else
            SDN(k,1) = median(INFO.gps(:,1),'omitnan');
            LAT(k,1) = median(INFO.gps(:,3),'omitnan');
            LON(k,1) = median(INFO.gps(:,2),'omitnan');
            CYC(k,1) = INFO.cast;
            k=k+1;
        end
    end

%     inputdir = ['\\atlas\Chem\ARGO_PROCESSING\DATA\ERA5\',var{2},'\']
    inputdir = [user_dir,'\ERA5\',var{2},'\']
    try
        [ERA] = getERA_sO2(inputdir,SDN,LON,LAT);
        ERA.PRES
        try
            [ERA_NRT] = getERA_NRT_sO2(inputdir,SDN,LON,LAT);
        catch
            disp(['NRT ERA grab failed for float ', WMO])
        end
        tt = ERA_NRT.PRES;
        if exist('ERA_NRT') && ~isempty(tt)
            xxx = find(~isnan(ERA_NRT.PRES));
            ERA.PRES(xxx) = ERA_NRT.PRES(xxx);
%             clear ERA_NRT
        end
    catch
        ERA.PRES = [];
        disp(['ERROR!!! ERA5',var{2},'grab failed for float ',WMO,'.'])
        errorflts{kk} = WMO;
        kk = kk+1;
    end
    FLT.WMO = WMO;
    FLT.SDN = SDN;
    FLT.LAT = LAT;
    FLT.LON = LON;
    FLT.CYC = CYC;
    if strcmp(var{2},'sp') == 1
        save([ERA_repo,WMO,'_ERA5ref.mat'],'FLT','ERA');
    else
        save([user_dir,'\ERA5\FLOAT_REF_',var{2},'\',WMO,'_ERA5ref.mat'],'FLT','ERA');
    end
end

toc