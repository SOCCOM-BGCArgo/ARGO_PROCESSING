% Loop_ARGO_float_auto.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT TO RUN Loop_ARGO_float.m IN AUTOMATIC MODE ("UPDATE").  THIS SCRIPT GETS CALLED
% BY THE BAT FILE USED IN WINDOWS TASK SCHEDULER.
%
% TANYA MAURER
% MBARI
% 06/05/2017 
%
% **** UPDATED 5/2/18 TO ALSO RUN SENSOR HEALTH COVERAGE-MAPS.  ONLY RUN THESE ONCE PER DAY ****	
% **** UPDATED 10/2/18 TO UPDATE THE SHIPBOARD DATA REFERENCE TABLE.  CHECKS FOR UPDATES ON CCHDO.  ONLY NEEDS TO RUN PERIODICALLY (WEEKLY) ****
% **** UPDATED 07/24/19 REMOVED ADDED JOBS, NOW ONLY RUNNING MAIN FLOAT
%                       UPDATE DRIVER.  ADDITIONAL TASKS (SENSOR HEALTH COVERAGE MAPS, SHIPBOARD
%                       DATA, INTERACTIVE MAP JSON FILE) ARE NOW PROCESSED THROUGH SCRIPT
%                       "run_SOCCOM_daily_tasks.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Loop_ARGO_float('update','',{' '})

MATLABfinish


