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
% UPDATES: 01/08/21: TM ftp dir no longer exists, moved to websave to
%          11/13/24:  TM update to WOA2023
% extract files from ncei.noaa.gov.
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
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'LIAR',fp]);
    addpath(genpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'ESPER',fp]));
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'GUIS',fp,'SAGE_O2Argo',fp]);
	addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'GUIS',fp,'SAGE',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'MISC',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'WOA',fp]);
	addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'CANYON',fp]);
	addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,'CANYON_B',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'MFILES',fp,]);

    disp('INSTALLING "ARGO_PROCESSING\DATA" PATHS...')
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'SHIPBOARD',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'ARGO_REPO',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'CANYON',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'GLODAP',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'FLOATS',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'FLOATVIZ',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'LIAR',fp]);
    addpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'SHIPBOARD',fp]);
    addpath(genpath([topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'WOA2023',fp]));

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
disp('CHECKING FOR LOCAL WOA2023 OXYGEN FILES...')
WOAdir = [topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'WOA2023',fp,'oxygen',fp];
wdir = dir([WOAdir,'woa23_all_o*_01.nc']);
woafiles = char(wdir.name);
if size(woafiles,1)<17 %not all files exist on repo
    disp('POPULATING LOCAL WOA2023 REPOSITORY FOR OXYGEN...')
    disp('May take a few minutes.')
    disp(['You can monitor download progress at ',WOAdir])
    try
        % Jan2021 -- TM. ftp site no longer exists.  Use websave with
        % 'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/oxygen/netcdf/all/1.00/'
        % to download all 17 files.
        WOAurl = 'https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen/netcdf/all/1.00/';
        for i = 1:17 %hard-coded for now so that error will inform of any data file removals.  Could parse html from webread to get filenames, but could be risky if html format suddenly changes.
        	woaFname = ['woa23_all_o',num2str(i,'%02d'),'_01.nc'];
            websave([WOAdir,woaFname],[WOAurl,woaFname])
            disp([woaFname,' SUCCESSFULLY DOWNLOADED FROM https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen/netcdf/all/1.00/.'])
        end
    catch
        disp('CONNECTION TO https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/oxygen/netcdf/all/1.00/ FAILED.  ENDING INSTALL.')
        return
    end
else
    disp('WOA files were found:')
    woafiles
end

%pause(30)
disp('CHECKING FOR LOCAL WOA2023 NITRATE FILES...')
WOAdir = [topdir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'WOA2023',fp,'nitrate',fp];
wdir = dir([WOAdir,'woa23_all_n*_01.nc']);
woafiles = char(wdir.name);
if size(woafiles,1)<17 %not all files exist on repo
    disp('POPULATING LOCAL WOA2023 REPOSITORY FOR NITRATE...')
    disp('May take a few minutes.')
    disp(['You can monitor download progress at ',WOAdir])
    try
        % Jan2021 -- TM. ftp site no longer exists.  Use websave with
        % 'https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/oxygen/netcdf/all/1.00/'
        % to download all 17 files.
        WOAurl = 'https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/nitrate/netcdf/all/1.00/';
        for i = 1:17 %hard-coded for now so that error will inform of any data file removals.  Could parse html from webread to get filenames, but could be risky if html format suddenly changes.
        	woaFname = ['woa23_all_n',num2str(i,'%02d'),'_01.nc'];
            websave([WOAdir,woaFname],[WOAurl,woaFname])
            disp([woaFname,' SUCCESSFULLY DOWNLOADED FROM https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/nitrate/netcdf/all/1.00/.'])
        end
    catch
        disp('CONNECTION TO https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/nitrate/netcdf/all/1.00/ FAILED.  ENDING INSTALL.')
        return
    end
else
    disp('WOA files were found:')
    woafiles
end    
 


    
disp('INSTALL COMPLETE.')
