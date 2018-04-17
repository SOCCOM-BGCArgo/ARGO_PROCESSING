function getfile_runtOcean(float,prof,varargin)
% ************************************************************************
% getfile_runtOcean.m
% ************************************************************************
%
% Retrieves ARGO float data files from :
% http://runt.ocean.washington.edu/argo/data/
% URL maintained by Dana Swift <swift@ocean.washington.edu>
%
%
% USE AS:  **sample command-line options for float 8514, profile 15**:
%          getfile_runtOcean(8514,15)
%          getfile_runtOcean(8514,15,[],{'msg','dura','isus'})
%          getfile_runtOcean(8514,15,'\\atlas\Chem\ISUS\Tanya',[])
%          getfile_runtOcean(8514,15,'\\atlas\Chem\ISUS\Tanya',{'msg','dura','isus'})
%
% INPUTS:  float      = SOCCOM float ID
%          prof       = profile number
%        ((mydir))    = directory path to write to (optional)
%        ((getfiles)) = cell array of file-types to retrieve (optional)
%
% OUTPUTS: Files are saved to local directory.  
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 12/8/16
% UPDATES:
% NOTES: If optional arguments are
% left unspecified, only .msg file is retrieved and saved to local Matlab
% work dir.
% ************************************************************************
%
% ************************************************************************

nvarg = length(varargin); %Variable argument list

%__________________________________________________________________________
%Unless specified, output dir = working directory
if nvarg == 0 || isempty(varargin{1}); 
    mydir = pwd;
else
    mydir = varargin{1};
end

%Unless specified, retrieve only .msg file for select float.prof
if nvarg ==0 || isempty(varargin{2})
    getfiles = {'msg'};
else 
    getfiles = varargin{2};
end

%__________________________________________________________________________
%Retrieve files from URL
danaURL = 'http://runt.ocean.washington.edu/argo/data/';
for i = 1:length(getfiles);
    myfile = [sprintf('%4.4d',float),'.',sprintf('%3.3d',prof),'.',getfiles{1,i}];
    try
        websave([mydir,'\',myfile],[danaURL,sprintf('%4.4d',float),'/',myfile]);
    catch
        warning(['Error retrieving .',getfiles{1,i},' file for float ',...
            sprintf('%4.4d',float),'.  Check URL for existence of file.'])
    end
    disp([myfile,' successfully saved to ',mydir,'\'])
end

%__________________________________________________________________________
end %end getfile_runtOcean.m
    
    
