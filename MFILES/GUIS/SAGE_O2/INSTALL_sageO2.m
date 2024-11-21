function INSTALL_sageO2(MBARInet,msgfile_dir)
% ************************************************************************
% INSTALL_sageO2.m
% ************************************************************************
%
% Saves all relevant paths to your MATLAB directory, and stores your
% working directory path for sageO2 to reference.
%
%
% USE AS:  INSTALL_sageO2(1,'\\seaecho.shore.mbari.org\floats\');
%          INSTALL_sageO2(0,'C:\mymsgfilelocation\')
%
% INPUTS:  msgfile_dir = path to msg files.
%          MBARInet   = 1: working at MBARI, connected to msg files on network
%                       0: do not have access to MBARI network; have stored
%                       msg files elsewhere, as defined by 'msgfile_dir'
%
% OUTPUTS:
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 10/25/17
% UPDATES: 
% NOTES:   
% ************************************************************************
%
% ************************************************************************

workinDIR = pwd;
XX = strfind(workinDIR,'ARGO_PROCESSING');
if isempty(XX)
    disp('ERROR: ARGO_PROCESSING DIRECTORY NOT FOUND.  DID YOU MAINTAIN THE ARGO_PROCESSING TOP LEVEL DIRECTORY STRUCTURE?')
else
    topdir = workinDIR(1:XX-1);
    disp('INSTALLING "ARGO_PROCESSING\MFILE" PATHS...')
    addpath([topdir,'\ARGO_PROCESSING\MFILES\CANYON\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\FLOATS\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\GLODAP\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\GUIS\SAGE\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\GUIS\SAGE_O2\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\LIAR\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\MISC\']);
%     addpath([topdir,'\ARGO_PROCESSING\MFILES\NetCDF\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\WOA\']);
    addpath([topdir,'\ARGO_PROCESSING\MFILES\']);

    disp('INSTALLING "ARGO_PROCESSING\DATA" PATHS...')
    addpath([topdir,'\ARGO_PROCESSING\DATA\FLOATS\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\CAL\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\CAL\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\FLOATVIZ\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\FLOATVIZ\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\CAL\QC_LISTS\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\SHIPBOARD\']);
    %addpath([topdir,'\ARGO_PROCESSING\DATA\QC_images\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\NCEP_TEMPORARY\']);
    addpath([topdir,'\ARGO_PROCESSING\DATA\WOA2018\']);
    addpath('C:\temp\');
    addpath(msgfile_dir);
    savepath
end
save('sageO2_workinDIR.mat','topdir','msgfile_dir','MBARInet');
disp('INSTALL COMPLETE.')
