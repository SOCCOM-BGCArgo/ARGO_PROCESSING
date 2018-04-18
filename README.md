# ARGO_PROCESSING
Code library to support processing and quality control of BGC Argo floats.

Welcome to MBARI's ARGO_PROCESSING code repository.  The repo has two main subdirectories, one containing our MATLAB code base (MFILES), and one containing supporting reference data (DATA).  Within MFILES/FLOATS are a number of functions and data parsers that are internal to processing at MBARI, and that also support the Graphical User Interfaces (GUIs) we developed for quality control of BGC data. 

The /MFILES/GUIS/SAGE_O2Argo/ directory holds code used to run our oxygen-specific GUI, SAGE-O2Argo (SOCCOM Assessment and Graphical Evaluation for Oxygen, for use with Argo floats).  This version of SAGEO2 was developed specifically for use with Argo netCDF files and assists the user in deriving optode gain correction factors by comparing float oxygen with specific reference datasets (ie WOA2013, and NCEP reanalysis).  The software accesses files and routines within other subdirectories of the repository, so be sure to clone the entire repo before launching the GUI software.  Additionally, a SAGE-O2Argo user manual with detailed install inscructions can be found within MFILES/GUIS/SAGE_O2Argo/.
 
The SAGE GUI (used for correcting pH and NO3- data), as well as a version of SAGEO2 that utilizes .msg files (rather than Argo *.nc files) may be added to this repository in the future under /MFILES/GUIS/SAGE/ and /MFILES/GUIS/SAGE_O2/, respectively.  

Please note that this code is provided as-is and is subject to periodic updates and improvements.  If you are interested in contributing to this repository, please contact Tanya Maurer at tmaurer@mbari.org.


