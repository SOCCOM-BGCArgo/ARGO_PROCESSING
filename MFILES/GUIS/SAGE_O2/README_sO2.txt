README_sO2.TXT

This is a text file with notes on how to effectively utilize the SAGE-O2 GUI.

1. This GUI was generated using the GUI Layout Toolbox (must install before running).  
	Download at: https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
2. Other toolboxes essential to GUI function: 
	nansuite (https://www.mathworks.com/matlabcentral/fileexchange/6837-nan-suite) 
	SEAWATER (https://www.mathworks.com/matlabcentral/linkexchange/links/741)
	nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox)
3. Paths: The GUI runs best if you place the "ARGO_PROCESSING" directory in C:\Users\XXXX\Documents\MATLAB\ (where XXXX is your username)
	If you place it elsewhere, be sure to edit the dirs.user_dir variable assignment in ARGO_PROCESSING\MFILES\GUIS\SAGE_O2\SetupDirs_sO2.m
	Make sure all appropriate ARGO_PROCESSING subdirectories are added to your path.
4. This GUI accesses float msg files, stored in float directories fXXXX (APEX) or nXXXX (NAVIS), during operation (XXXX = UW float ID).  
	At MBARI, these are located on a Network drive.  For other users, change the dirs.msg variable assignment in SetupDirs_sO2.m to the appropriate path.
5. This GUI is most useful for floats that have produced 5 or more profiles.
6. Comparisons to climatology are done using NCEP and ERA-interim reanalysis products.  Note that ERA-interim has an update delay of close to 4 months (longer than NCEP).  Due to the update delay, and connection speeds, ERA datasets are not downloaded realtime.  Annual datasets are stored in \ARGO_PROCESSING\DATA\ERA_INT.  To augment this repository with more recent data than the time of this writing, visit http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/ 

This GUI gets updated periodically with added functionality and enhancements.  It is provided as-is. 
If you need further assistance, please contact:

Tanya Maurer
Monterey Bay Aquarium Research Institute
tmaurer@mbari.org


