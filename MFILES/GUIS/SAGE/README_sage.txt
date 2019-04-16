README_SAGE.TXT

This is a text file with notes on how to effectively utilize the SAGE GUI outside of MBARI.  Please note that this has been minimally tested on external floats!  While a pdf manual exists for SageO2Argo (oxygen QC gui), currently only this readme exists for the SAGE GUI (used for pH and nitrate QC).

1. This GUI was generated using the GUI Layout Toolbox (must install before running).  
	Download at: https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
2. Other toolboxes essential to GUI function: 
	nansuite (https://www.mathworks.com/matlabcentral/fileexchange/6837-nan-suite) 
	SEAWATER (https://www.mathworks.com/matlabcentral/linkexchange/links/741)
	nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox) - follow install instructions.
3. Paths: The GUI runs best if you place the "ARGO_PROCESSING" directory in C:\Users\XXXX\Documents\MATLAB\ (where XXXX is your username)
	Once you've downloaded the software and supporting toolboxes, navigate to ...\ARGO_PROCESSING\MFILES\GUIS\SAGE\ and run INSTALL_sage.m  The first time you run the install will take some time (~10 min?) since the software will need to download supporting reference data from the web( ie NCEP,WOA).
4. Supporting Data:  As mentioned in (3) WOA2013 data (for oxygen, O2sat, and nitrate) is not included in the repo due to file size.  Data will be downloaded automatically once the install is run.  
5. The SAGE GUI works off of ODV-compatible ascii files.  At MBARI, we generate these routinely for the SOCCOM array.  If the intent is to use this GUI for external BGC Argo floats, data from Sprof*.nc files must first be converted to ODV format.  The file converter, ARGOsprofmat2ODV.m, is provided and located in ARGO_PROCESSING\MFILES\GUIS\SAGE_O2Argo\SProf_Conversion\ (the converter is also required for SAGE_O2Argo). It is suggested that you create a working directories specific to your floats within ARGO_PROCESSING/DATA/ARGO_REPO/ for which to hold your converted Sprof files.
6. This GUI is most useful for floats that have produced 5 or more profiles.  DO NOT USE THIS GUI UNLESS OXYGEN QC HAS BEEN PERFORMED!! (IT WILL NOT WORK, THE GUI WILL BE SEARCHING FOR AN ODV*QC.TXT FILE). (see SAGE-O2 GUI).  The empirical algorithms used to correct nitrate and pH within SAGE rely on accurate oxygen data as an input variable.  If you have performed your own oxygen adjustment, you can circumvent the use of SAGE-O2, and manually make a QC directory with an ODV*QC.TXT file similar to the raw file but with your O2 corrections applied.
7. Once properly installed, navigating to ...ARGO_PROCESSING\MFILES\GUIS\SAGE\ and typing 'sage' at the command prompt will launch the GUI, adding all necessary paths.


This GUI gets updated periodically with added functionality and enhancements.  It is provided as-is! It is used extensively on APEX and NAVIS floats within the SOCCOM program but has not been tested on all BGC Argo floats!!  Modifications may be needed to suit you and your floats' needs!
If you need further assistance or would like to provide feedback, please contact:

Tanya Maurer & Josh Plant
Monterey Bay Aquarium Research Institute
tmaurer@mbari.org
jplant@mbari.org



