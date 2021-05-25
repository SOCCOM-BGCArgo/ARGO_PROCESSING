function dirs = SetupDirs_sage
% ************************************************************************
% SetupDirs_sage.m
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
    load sage_workingDIR.mat %loads 'topdir', 'msgfile_dir, 'MBARInet'
catch
    warning('Problem loading sage_workingDIR.mat.  Did you run INSTALL_sage from ../ARGO_PROCESSING/MFILES/GUIS/SAGE/ ?')
end

dirs.user_dir = topdir;
dirs.mfiles    = [dirs.user_dir,'ARGO_PROCESSING',filesep,'MFILES',filesep];
dirs.woa       = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'WOA2018',filesep];
dirs.glodap    = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'GLODAP',filesep];
dirs.mat       = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'FLOATS',filesep];
dirs.cal       = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'CAL',filesep];
dirs.FVlocal   = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'FLOATVIZ',filesep];
dirs.FV        = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'FLOATVIZ',filesep];
%dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
dirs.QCadj     = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'CAL',filesep,'QC_LISTS',filesep];
dirs.bottle    = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'SHIPBOARD',filesep];
dirs.QC_images = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'QC_images',filesep];
dirs.CANY = [dirs.user_dir,'ARGO_PROCESSING',filesep,'DATA',filesep,'CANYON',filesep];
if ~isempty(strfind(computer,'PC'))
    dirs.temp      = ['C:',filesep,'temp',filesep];
else
    dirs.temp = [getenv('HOME'),filesep,'temp',filesep];
end
% msg file directories (MBARI).  Comment out if not used.
dirs.msg   = '\\seaecho.shore.mbari.org\floats\';
% dirs.alt       = '\\atlas\ChemWebData\floats\alternate\'; 
% dirs.msg_comb  = '\\atlas\ChemWebData\floats\combined\';
end


