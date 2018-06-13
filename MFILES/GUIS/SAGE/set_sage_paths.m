function tf = set_sage_paths(user_dir)

% This function checks to see if the required paths are set so that SAGE
% and the float processing functions can find all the appropriate files

% user_dir -    highest level dir where all the processing and GUI code is
%               stored
%               exmple: user_dir = getenv('USERPROFILE'); 
%                   returns user path,i.e. 'C:\Users\jplant'
%

% ************************************************************************
% TESTING
%user_dir = [getenv('USERPROFILE'),'\Documents\MATLAB\']; 

% ************************************************************************

 tf = 0;
% ************************************************************************
% CREATE A STRUCTURE OF DIRS TO CHECK / ADD
% ************************************************************************
% MFILE PATHS
% ADD NEW PATHS TO ThE gd STRUCTURE

gd.mCANYON  = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'CANYON'];
gd.mFLOATS  = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'FLOATS'];
gd.mGLODAP  = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'GLODAP'];
gd.mSAGE    = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'GUIS',filesep,'SAGE'];
gd.mSAGE_O2 = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'GUIS',filesep,'SAGE_O2'];
gd.mLIAR    = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'LIAR'];
gd.mMFILES  = ['ARGO_PROCESSING',filesep,'MFILES'];
gd.mMISC    = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'MISC'];
gd.mNetCDF  = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'NetCDF'];
gd.mSAGE    = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'GUIS',filesep,'SAGE'];
gd.mSAGE_O2 = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'GUIS',filesep,'SAGE_O2'];
gd.mWOA2013 = ['ARGO_PROCESSING',filesep,'MFILES',filesep,'WOA2013'];

% DATA PATHS
gd.dCANYON    = ['ARGO_PROCESSING',filesep,'DATA',filesep,'CANYON'];
gd.dERA_INT   = ['ARGO_PROCESSING',filesep,'DATA',filesep,'ERA_INT'];
gd.dGLODAP    = ['ARGO_PROCESSING',filesep,'DATA',filesep,'GLODAP'];
gd.dLIAR      = ['ARGO_PROCESSING',filesep,'DATA',filesep,'LIAR'];
gd.dSHIPBOARD = ['ARGO_PROCESSING',filesep,'DATA',filesep,'SHIPBOARD'];
gd.dWOA       = ['ARGO_PROCESSING',filesep,'DATA',filesep,'WOA2013'];

% MISC PATHS
%gd.nctoolbox  = 'nctoolbox';

% ************************************************************************
% GET EXISTING MATLAB PATHS & STORE AS CELL ARRAY
matpaths = strread(path,'%s','delimiter', pathsep);

% ************************************************************************
% CHECK PATHS IN DIRS STRUCTURE FIRST
dirs_names = fieldnames(gd);
for i = 1 : size(dirs_names)
    t1 = strcmpi([user_dir,gd.(dirs_names{i})], matpaths);
    if sum(t1) == 0; % NO MATCH, ADD PATH
        addpath([user_dir,gd.(dirs_names{i})])
        disp(['Path added to Matlab path: ', user_dir,gd.(dirs_names{i})])
        tf = 1;
    end
end

if tf == 0
    disp('All required paths already exist - no updates required')
end

        
    
    
    
    
    