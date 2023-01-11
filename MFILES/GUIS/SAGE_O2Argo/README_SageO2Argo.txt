README_sO2Argo.TXT

This is a text file with notes on how to effectively utilize the SAGE-O2 GUI for Argo floats.

1. This GUI's framework was structured using the GUI Layout Toolbox (must install before running).  
	Download at: https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
2. Other toolboxes essential to GUI function: 
	nansuite (https://www.mathworks.com/matlabcentral/fileexchange/6837-nan-suite) 
	SEAWATER (https://www.mathworks.com/matlabcentral/linkexchange/links/741)
	nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox)
3. Paths: The GUI runs best if you place the "ARGO_PROCESSING" directory in C:\Users\XXXX\Documents\MATLAB\ (where XXXX is your username)
4. This GUI accesses the following Argo float file types, you must have each for each float of interest (although *BRtraj.nc is not needed for floats without air-cal).  I suggest storing them in float-specific subdirectories within a master folder repository.  For example, I keep mine in, C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\DATA\ARGO_REPO\, with subdirectories ie 6900889, 6900890, etc.
	*Sprof.nc
	*meta.nc
	*BRtraj.nc
5. This GUI is most useful for floats that have produced 5 or more profiles.
6. Comparisons to climatology are done using NCEP and WOA products. NCEP data is stored C:\Users\...\Documents\MATLAB\ARGO_PROCESSING\DATA\NCEP_TEMPORARY.  To change this directory, or download NCEP data in real time within the GUI, see example NCEP path assignments in C:\Users\...\Documents\MATLAB\ARGO_PROCESSING\MFILES\FLOATS\getNCEP.m
7.  Steps to running GUI:
	PREP GUI (TO BE DONE PRIOR TO FIRST USE):
		a. install appropriate toolboxes (see 1., 2.)
		b. copy ARGO_PROCESSING repository and check paths and directory structures, as necessary (see 3., 4.)
		c. navigate to ...\ARGO_PROCESSING\MFILES\GUIS\SAGE_O2Argo and run the install (see INSTALL_sO2Argo.m)
		d. consider climatology reference preferences (if your float does not take air-cal measurements, GUI will default to World Ocean Atlas).
	PREP DATA (TO BE DONE FOR EACH NEW FLOAT YOU'D LIKE TO VISUALIZE)
		a. organize files defined in 4., may need to download from GDAC (ftp.ifremer.fr; or usgodae.org)
		b. convert *Sprof.nc to ODV*.TXT by running ARGOsprofmat2ODV.m (located in ...\ARGO_PROCESSING\MFILES\GUIS\SAGE_O2Argo\SProf_Conversion\, see inputs in file).  NOTE that this process takes some time...so you may want to run a loop for all floats you are interested in viewing in the background at your convenience, in advance of a DMQC session.  
		c. be sure TXT file produced in b. is stored in float-specific directory described in 4.
	RUN GUI
		a. type >> sageO2Argo
		b. click "Select Float" (top left) and click on directory for float of interest.
		c. you should now be able to visualize your data, toggle between reference datasets, enter specific profile or depth ranges, etc.  Note that the 'profile' view is not fully supported for all Argo floats due to the high format variability in shipboard bottle datasets supplied to various DACs.
	
This GUI gets updated periodically with added functionality and enhancements.  It is provided as-is.  See ...\ARGO_PROCESSING\MFILES\GUIS\SAGE_O2Argo\acknowledgements_sO2_Argo.txt for license information.
If you need further assistance, please contact:

Tanya Maurer
Monterey Bay Aquarium Research Institute
tmaurer@mbari.org


