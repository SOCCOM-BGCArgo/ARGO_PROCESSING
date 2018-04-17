README_sO2Argo.TXT

This is a text file with notes on how to effectively utilize the SAGE-O2 GUI for Argo floats.

1. This GUI's framework was structured using the GUI Layout Toolbox (must install before running).  
	Download at: https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
2. Other toolboxes essential to GUI function: 
	nansuite (https://www.mathworks.com/matlabcentral/fileexchange/6837-nan-suite) 
	SEAWATER (https://www.mathworks.com/matlabcentral/linkexchange/links/741)
	nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox)
3. Paths: The GUI runs best if you place the "ARGO_PROCESSING" directory in C:\Users\XXXX\Documents\MATLAB\ (where XXXX is your username)
4. This GUI accesses the following Argo float file types, you must have each for each float of interest (although *BRtraj.nc is not needed for floats without air-cal).  I suggest storing them in float-specific subdirectories within a master folder repository.  For example, I keep mine in, C:\Users\tmaurer\Documents\ARGO\EXTERNAL_DACs\DAC_files\, with subdirectories 6900889, 6900890, etc.
	*Mprof.nc
	*meta.nc
	*BRtraj.nc
5. This GUI is most useful for floats that have produced 5 or more profiles.
6. Comparisons to climatology are done using NCEP and ERA-interim reanalysis products.  Note that ERA-interim has an update delay of close to 4 months (longer than NCEP).  Due to the update delay, connection speeds, and access requirements at ECMWF, ERA datasets are no longer actively supported in this GUI.  Pre-downloaded datasets until March2017 are stored in \ARGO_PROCESSING\DATA\ERA_INT.  To augment this repository with more recent data than the time of this writing, you may try visiting http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/ NCEP data is stored C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\DATA\NCEP_TEMPORARY.  To change this directory, or download NCEP data in real time within the GUI, see lines 90-91 in C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\MFILES\FLOATS\getNCEP.m
7.  Steps to running GUI:
	PREP GUI (TO BE DONE PRIOR TO FIRST USE):
		a. install appropriate toolboxes (see 1., 2.)
		b. copy ARGO_PROCESSING repository and check paths and directory structures, as necessary (see 3., 4.)
		c. navigate to ...\ARGO_PROCESSING\MFILES\GUIS\SAGE_O2Argo and run the install (see INSTALL_sO2Argo.m)
		d. consider climatology reference preferences (if your float does not take air-cal measurements, GUI will default to World Ocean Atlas).
	PREP DATA (TO BE DONE FOR EACH NEW FLOAT YOU'D LIKE TO VISUALIZE)
		a. organize files defined in 4., may need to download from GDAC (ftp.ifremer.fr; or usgodae.org)
		b. convert *Mprof.nc to ODV*.TXT by running mprof2mat.m followed by mprofmat2ODV.m (function, see inputs in file).  
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


