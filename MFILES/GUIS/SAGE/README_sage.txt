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
	Make sure all appropriate ARGO_PROCESSING subdirectories are added to your path.
4. This GUI accesses floatviz text files, stored in dirs.FVlocal.  To view or change this directory definition, see the first few lines of sage.m.  
5. This GUI is most useful for floats that have produced 5 or more profiles.  Do not use this GUI unless float oxygen data has been quality controlled (see SAGE-O2 GUI).  The empirical algorithms used to correct nitrate and pH within SAGE rely on accurate oxygen data as an input variable.

This GUI gets updated periodically with added functionality and enhancements.  It is provided as-is. 
If you need further assistance, please contact:

Josh Plant & Tanya Maurer
Monterey Bay Aquarium Research Institute
jplant@mbari.org
tmaurer@mbari.org


