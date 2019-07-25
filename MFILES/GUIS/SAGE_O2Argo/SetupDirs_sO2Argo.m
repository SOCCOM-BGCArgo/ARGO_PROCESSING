function dirs = SetupDirs_sO2Argo
% ************************************************************************
% SetupDirs_sO2.m
% ************************************************************************
%
% Sets up directories referenced throughout GUI.
%
%
% USE AS:  DIRS = SetupDirs_sageO2;
%
% INPUTS:  
%
% OUTPUTS: DIRS:  structure containing directory paths used throughout the
%                 sageO2 GUI.
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 12/13/16
% UPDATES: Updated 3/8/17 to allow gui user to more easily work from
%          preferred directory location.
%          Updated 10/25/2017 to improve path install process.
% NOTES:   Before launching SageO2, be sure you have run INSTALL_sageO2
% from ../ARGO_PROCESSING/MFILES/GUIS/SAGE_O2/.
% ************************************************************************
%
% ************************************************************************

fp = filesep; % File separator for current platform

try
    load sageO2Argo_workingDIR.mat %loads 'topdir', 'msgfile_dir, 'MBARInet'
catch
    warning('Problem loading sageO2_workingDIR.mat.  Did you run INSTALL_sageO2Argo from ../ARGO_PROCESSING/MFILES/GUIS/SAGE_O2Argo/ ?')
end

dirs.user_dir = topdir;
dirs.mfiles    = [dirs.user_dir,fp,'ARGO_PROCESSING',fp,'MFILES',fp];
dirs.bottle    = [dirs.user_dir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'SHIPBOARD',fp];
dirs.QCadj     = [dirs.user_dir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'CAL',fp,'QC_LISTS',fp]; % empty, but DACs may populate.  GUI checks for QC_LISTs (pre-stored adjustments to populate table)
dirs.NCEP_TEMP = [dirs.user_dir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'NCEP_TEMPORARY',fp];
dirs.Argo      = [dirs.user_dir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'ARGO_REPO',fp];
dirs.woa       =[dirs.user_dir,'\ARGO_PROCESSING\DATA\WOA2018\'];

% ERA data reference has been deprecated in this version of the software.
% Keep file path definition in comments in case of future re-incorporation.
% dirs.ERA       = [dirs.user_dir,fp,'ARGO_PROCESSING',fp,'DATA',fp,'ERA_INT',fp];
% dirs.temp      = 'C:\temp\';
%dirs.Argo = DATAdir;
end


