% GET MATFILE LIST

WMO_dir = '5903754'; %(UW ID 7593)  
 
% SET UP DIRS
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir = [user_dir, '\Documents\MATLAB\'];
dirs.mat = [user_dir,'ARGO_PROCESSING\DATA\FLOATS\'];
dirs.msg = '\\atlas\ChemWebData\floats\';
dirs.OB  = 'C:\temp\FLOATS\'; % ODD BALL

OB_path  = [dirs.OB, WMO_dir,'\']

% COPY WMO DIR TO ODD BALL DIR
copyfile([dirs.mat,WMO_dir,],'destination')

% STEP THROUGH *.MAT PROFILES
%   GET ACCOMPANYING MESSAGE FILE
%   PARSE LR DATA INCLUDING NAN PRESSURE LINES