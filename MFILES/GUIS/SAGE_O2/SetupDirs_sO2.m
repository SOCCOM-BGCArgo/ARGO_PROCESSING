function dirs = SetupDirs_sO2
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

try
    load sageO2_workinDIR.mat %loads 'topdir', 'msgfile_dir, 'MBARInet'
catch
    warning('Problem loading sageO2_workinDIR.mat.  Did you run INSTALL_sageO2 from ../ARGO_PROCESSING/MFILES/GUIS/SAGE_O2/ ?')
end

dirs.user_dir = topdir;
dirs.mfiles    = [dirs.user_dir,'\ARGO_PROCESSING\MFILES\'];
dirs.mat       = [dirs.user_dir,'\ARGO_PROCESSING\DATA\FLOATS\'];
dirs.cal       = [dirs.user_dir,'\ARGO_PROCESSING\DATA\CAL\'];
dirs.NO3config = [dirs.user_dir,'\ARGO_PROCESSING\DATA\CAL\'];
dirs.FVlocal   = [dirs.user_dir,'\ARGO_PROCESSING\DATA\FLOATVIZ\'];
dirs.FV        = [dirs.user_dir,'\ARGO_PROCESSING\DATA\FLOATVIZ\'];
dirs.QCadj     = [dirs.user_dir,'\ARGO_PROCESSING\DATA\CAL\QC_LISTS\'];
dirs.bottle    = [dirs.user_dir,'\ARGO_PROCESSING\DATA\SHIPBOARD\'];
dirs.QC_images = [dirs.user_dir,'\ARGO_PROCESSING\DATA\QC_images\'];
% dirs.ERA       = [dirs.user_dir,'\ARGO_PROCESSING\DATA\ERA5\'];
%dirs.ERA       = ['\\atlas\Chem\ARGO_PROCESSING\DATA\ERA5\sfcpres\']; %Big files -- keep on network? Or can modify this to local repo!
dirs.ERA = '\\atlas\Chem\ARGO_PROCESSING\DATA\ERA5\FLOAT_REF\';
dirs.NCEP_TEMP       = [dirs.user_dir,'\ARGO_PROCESSING\DATA\NCEP_TEMPORARY\'];
dirs.woa       =[dirs.user_dir,'\ARGO_PROCESSING\DATA\WOA2023\'];
dirs.temp      = 'C:\temp\';

if MBARInet == 1
    dirs.msg       = '\\atlas\ChemWebData\floats\'; %main msg file directory (MBARI).  Modify if needed.
    dirs.dup       = '\\atlas\ChemWebData\floats\duplicate\'; %duplicate msg dir.
    dirs.alt       = '\\atlas\ChemWebData\floats\alternate\'; %alternate msg file directory (MBARI).  Comment out if not used.
    dirs.msg_comb  = '\\atlas\ChemWebData\floats\combined\';  %combined msg file directory (MBARI).  Comment out if not used.
    dirs.config    = '\\atlas\Chem\ISUS\Argo\'; %MBARI-specific
else
    dirs.msg = msgfile_dir;
end

