function INSTALL_sage
% ************************************************************************
% INSTALL_sage.m
% ************************************************************************
%
% Saves all relevant paths to your MATLAB directory, and stores your
% working directory path for sage to reference.
%
%
% USE AS:  INSTALL_sage
%
% INPUTS:  
%
% OUTPUTS:
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 9/17/18
% UPDATES: 
% NOTES:   
% ************************************************************************
%
% ************************************************************************

fp = filesep; % File separator for current platform
workingDIR = pwd;
XX = strfind(workingDIR,'ARGO_PROCESSING');
if isempty(XX)
    disp('ERROR: ARGO_PROCESSING DIRECTORY NOT FOUND.  DID YOU MAINTAIN THE ARGO_PROCESSING TOP LEVEL DIRECTORY STRUCTURE?')
    disp('INSTALL INCOMPLETE!')
    return
else
    topdir = workingDIR(1:XX-1);
    disp('INSTALLING "ARGO_PROCESSING\MFILE" PATHS...')
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'FLOATS',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'GLODAP',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'GUIS',fp,'SAGE_O2Argo',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'MISC',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'WOA2013',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,]);

    disp('INSTALLING "ARGO_PROCESSING\DATA" PATHS...')
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'SHIPBOARD',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'ARGO_REPO',fp]);
    if ~isempty(strfind(computer,'PC'))
        tmploc     = ['C:',filesep,'temp',filesep];
    else
        tmploc = [getenv('HOME'),filesep,'temp',filesep];
    end
    if ~exist(tmploc,'dir')
        mkdir(tmploc)
    end
    addpath(tmploc)
    savepath
end
save([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'GUIS',fp,'SAGE',fp,'sage_workingDIR.mat'],'topdir');

% CHECK FOR WOA DATA.  IF DOESN'T EXIST, DOWNLOAD IT.
disp('CHECKING FOR LOCAL WOA2013 OXYGEN FILES...')
WOAdir = [topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'WOA2013',fp,'oxygen',fp];
wdir = dir([WOAdir,'woa13_all_o*_01.nc']);
woafiles = char(wdir.name);
if size(woafiles,1)<17 %not all files exist on repo
    disp('POPULATING LOCAL WOA2013 REPOSITORY FOR OXYGEN...')
    disp('May take a few minutes.')
    disp(['You can monitor download progress at ',WOAdir])
    try
        f = ftp('ftp.nodc.noaa.gov');
        cd(f);
        sf = struct(f);
        sf.jobject.enterLocalPassiveMode();
        cd(f,'pub/woa/WOA13/DATAv2/oxygen/netcdf/all/1.00/')
        disp('FTP CONNECTION TO ftp.nodc.noaa.gov/pub/woa/WOA13/DATAv2/oxygen/netcdf/all/1.00/ WAS SUCCESSFUL.')
        dir_list = dir(f);
        mget(f,'woa13_all_o*.nc',WOAdir)
        disp(['WOA2013 files saved to ',WOAdir])
        char(dir_list.name)
        close(f)
    catch
        disp('FTP CONNECTION TO ftp.nodc.noaa.gov FAILED.  ENDING INSTALL.')
        return
    end
else
    disp('WOA files were found:')
    woafiles
end

%pause(30)
disp('CHECKING FOR LOCAL WOA2013 NITRATE FILES...')
WOAdir = [topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'WOA2013',fp,'nitrate',fp];
wdir = dir([WOAdir,'woa13_all_n*_01.nc']);
woafiles = char(wdir.name);
if size(woafiles,1)<17 %not all files exist on repo
    disp('POPULATING LOCAL WOA2013 REPOSITORY FOR NITRATE...')
    disp('May take a few minutes.')
    disp(['You can monitor download progress at ',WOAdir])
    try
        f = ftp('ftp.nodc.noaa.gov');
        cd(f);
        sf = struct(f);
        sf.jobject.enterLocalPassiveMode();
        cd(f,'pub/woa/WOA13/DATAv2/nitrate/netcdf/all/1.00/')
        disp('FTP CONNECTION TO ftp.nodc.noaa.gov/pub/woa/WOA13/DATAv2/nitrate/netcdf/all/1.00/ WAS SUCCESSFUL.')
        disp('COPYING FILES....')
        dir_list = dir(f);
        mget(f,'woa13_all_n*.nc',WOAdir)
        disp(['WOA2013 files saved to ',WOAdir])
        char(dir_list.name)
        close(f)
    catch
        disp('FTP CONNECTION TO ftp.nodc.noaa.gov FAILED.  ENDING INSTALL.')
        return
    end
else
    disp('WOA files were found:')
    woafiles
end    
 


    
disp('INSTALL COMPLETE.')
