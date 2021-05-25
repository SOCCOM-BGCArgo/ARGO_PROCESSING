% run_MBARI_daily_tasks.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT TO RUN VARIOUS MBARI FLOAT-PROCESSING AND TRACKING -RELATED TASKS, INCLUDING:
%  1. UPDATES TO MBARI NON-SOCCOM FLOAT STATS TABLE, https://www3.mbari.org/chemsensor/MBARI_float_table/mbari_float_table.html
%  2. UPDATES TO MAT-FILE THAT STORES INFO ON DATE-OF-RECEIPT OF MOST RECENT MSG FILE.
%  3. UPDATES TO SHIPBOARD DATABASE USED IN SOCCOM QC
%  4. ?? OTHERS ??
%
% THIS SCRIPT GETS CALLED ONCE PER DAY
% BY A BAT FILE USED IN WINDOWS TASK SCHEDULER.
%
% TANYA MAURER
% MBARI
% 01/21/21
%
% NOTE THAT THIS SCRIPT WAS MODIFIED FROM run_SOCCOM_daily_tasks.m.  MANY 
% OF THE DAILY MAP-GENERATION TASKS WERE MERGED INTO A SINGLE SEPARATE SCRIPT, 
% C:\Users\bgcargovm\Documents\MATLAB\ARGO_PROCESSING\MFILES\MAP_ROOM\update_MBARI_floatmaps.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (1) -------------------------------------------------------------------------
% % Update NON-SOCCOM MBARI stats table
NONSOCCOM_make_sorttable_float_stats
pause(30)

% (2) -------------------------------------------------------------------------
% % Update mat file storing time-of-last-receipt of msg files for each float.  Main purpose is to know: "when did each float last talk to us?"
%get_latest_msgfile_times
%pause(30)

% (3) -------------------------------------------------------------------------
% % Check for bottle file updates and regenerate shipboard data lookup table
% % (gets used in GUIs)
process_data_ref_table
pause(30)

MATLABfinish
% 
