# ARGO_PROCESSING
Code library to support processing and quality control of BGC Argo floats.

Welcome to MBARI's ARGO_PROCESSING code repository.  The repo has two main subdirectories, one containing our supporting MATLAB code base (MFILES), and one containing supporting reference data (DATA).  Within MFILES/FLOATS you will find many functions that are internal to processing at MBARI (ie various data parsers and subroutines).  However, many of these may be useful to external users.  The MFILES/GUIS directory holds code used to run our oxygen-specific Graphical User Interface (GUI), SAGE-O2Argo (SOCCOM Assessment and Graphical Evaluation for Oxygen, for use with Argo floats).  This version of SAGEO2 was developed specifically for use with Argo netCDF files.  It also accesses files and routines within other subdirectories of the repository, so be sure to clone the entire repo before launching the GUI software.  Additionally, a SAGE-O2Argo manual can be found within MFILES/GUIS/SAGE_O2Argo/.

If you are interested in contributing to this repository, please contact Tanya Maurer at tmaurer@mbari.org.


