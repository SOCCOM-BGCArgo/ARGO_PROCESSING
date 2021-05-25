% run_SOCCOM_daily_tasks.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT TO RUN VARIOUS SOCCOM-RELATED TASKS, INCLUDING:
%  1. UPDATES TO SHIPBOARD DATABASE USED IN SOCCOM QC
%  2. UPDATES MAIN SOCCOM TRACK MAP
%  3. UPDATES TO SENSOR HEALTH COVERAGE MAPS.
% THIS SCRIPT GETS CALLED ONCE PER DAY
% BY A BAT FILE USED IN WINDOWS TASK SCHEDULER.
%
% TANYA MAURER
% MBARI
% 07/24/2019
%
% UPDATED 8/5/2019
% UPDATED 11/2/2020; TM added routine to create main SOCCOM track map (took
% job over from Ken).  

% -------------------------------------------------------------------------
% % Check for bottle file updates and regenerate shipboard data lookup table
% % (gets used in GUIs)
process_data_ref_table
pause(30)

% -------------------------------------------------------------------------
% Generate CHL sensor map (per Emmanuel's request, shows which soccom floats have bio-optics)
%currentdir = pwd;
%cd C:\Users\bgcargo\Documents\CHLsensormap\
% map_SOCCOM_chl
%pause(30) 
%cd(currentdir)


% -------------------------------------------------------------------------
cd C:\Users\bgcargo\Documents\MAProom\
map_SOCCOM_forWEBSITE

% -------------------------------------------------------------------------
%Generate sensor health coverage maps and save to sirocco
pause(120) %pause for 2 min
addpath C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\DATA\SENSOR_HEALTH_MAPS\
disp('GENERATING SENSOR HEALTH MAPS...')
savedir = 'C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\DATA\SENSOR_HEALTH_MAPS\';
Plot_sensor_coverage('N',savedir)
close all
keep savedir
disp('NITRATE MAPS COMPLETE.')
pause(60)
Plot_sensor_coverage('PH',savedir)
close all
keep savedir
disp('PH MAPS COMPLETE.')
pause(60)
Plot_sensor_coverage('O',savedir)
close all
keep savedir
disp('OXYGEN MAPS COMPLETE.')
pause(60)
Plot_sensor_coverage('CHL',savedir)
close all
keep savedir
disp('CHLOROPHYLL MAPS COMPLETE.')
pause(60)
Plot_sensor_coverage('B',savedir)
close all
keep savedir
disp('BACKSCATTER MAPS COMPLETE.')
pause(60)

MATLABfinish
% 
% copyfile C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\DATA\SENSOR_HEALTH_MAPS\*.png S:\soccom\images\SENSOR_HEALTH_MAPS\ f
