README_SAGE.TXT

This is a text file with notes on how to effectively utilize the SAGE GUI.

1. This GUI was generated using the GUI Layout Toolbox (must install before running).  
	Download at: https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
2. Other toolboxes essential to GUI function: 
	nansuite (https://www.mathworks.com/matlabcentral/fileexchange/6837-nan-suite) 
	SEAWATER (https://www.mathworks.com/matlabcentral/linkexchange/links/741)
	nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox)
3. Paths: The GUI runs best if you place the "ARGO_PROCESSING" directory in C:\Users\XXXX\Documents\MATLAB\ (where XXXX is your username)
	If you place it elsewhere, be sure to edit the dirs.user_dir variable assignment in ARGO_PROCESSING\MFILES\GUIS\SAGE\sage.m
4. Supporting Data:  WOA2013 data (for oxygen, O2sat, and nitrate) is not included in the repo due to file size.  Data should be downloaded manually from ftp://ftp.nodc.noaa.gov/pub/woa/WOA13/DATAv2/ and added to data directories in ...ARGO_PROCESSING\DATA\WOA2013\ prior to launching the GUI.  
5. This GUI accesses floatviz text files, stored in dirs.FVlocal.  To view or change this directory definition, see the first few lines of sage.m.  
6. This GUI is most useful for floats that have produced 5 or more profiles.  Do not use this GUI unless float oxygen data has been quality controlled (see SAGE-O2 GUI).  The empirical algorithms used to correct nitrate and pH within SAGE rely on accurate oxygen data as an input variable.
7. Once properly installed, navigating to ...ARGO_PROCESSING\MFILES\GUIS\SAGE\ and typing 'sage' at the command prompt will launch the GUI, adding all necessary paths.

This GUI gets updated periodically with added functionality and enhancements.  It is provided as-is. 
If you need further assistance, please contact:

Josh Plant & Tanya Maurer
Monterey Bay Aquarium Research Institute
jplant@mbari.org
tmaurer@mbari.org


