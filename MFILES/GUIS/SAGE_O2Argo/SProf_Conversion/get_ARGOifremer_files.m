function tf = get_ARGOifremer_files(WMOs,DAC,savedir)
% ************************************************************************
% get_ARGOifremer_files.m
% ************************************************************************
%
% Function to grab *BRtraj.nc, *Sprof.nc, and *meta.nc files from user-specifed
% external floats through the GDAC at ifremer.  Files are saved to
% user-specified local directory.
%
%
% USE AS:  tf = get_ARGOifremer_files({'6901582'},'coriolis','C:\Users\tmaurer\Documents\ARGO\CORIOLIS\CORIOLIS_files\')
%          tf = get_ARGOifremer_files({'6901582';'6900889';'6900890'},'coriolis','C:\Users\tmaurer\Documents\ARGO\CORIOLIS\CORIOLIS_files\')
%
% INPUTS:  
%       WMOs    = Cell array of WMO number(s) for float(s) of interest
%       DAC     = Data Assembly Center housing float of interest
%       savedir = Local directory within which to create the WMO dir to save files to 
%
%       *current DAC options at ifremer/argo/dac/ (maintain lower-case):
%           - aoml - atlantic oceanographic & meteorological laboratory
%           - bodc - british oceanographic data center
%           - coriolis - a component of the french and euro operational oceanography
%           - csio - china second institute of oceanography
%           - csiro - commonwealth scientific and industrial research organization (australia)
%           - incois - indian national center for ocean information services
%           - jma - japan meteorological agency
%           - kma - korean meteorological administration
%           - kordi - korea ocean research and development institute
%           - meds - marine environmental data service (canada)
%           - nmdis - national marine data and information service (china)
%
% OUTPUTS: tf = length(WMOs)x1 cell structure, 1 = file download success, 0 = file download failure (not found)
%               tf{i}(1) is associated with *BRtraj.nc
%               tf{i}(2) is associated with *Mprof.nc
%               tf{i}(3) is associated with *meta.nc
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 10/05/2017
% UPDATES: 7/31/19 Updated to extract Sprof file instead of Mprof file.
% NOTES: *Acquiring data from floats from separate DACs requires separate
% calls to the function.
%
% ************************************************************************
%
% ************************************************************************

if nargin < 3
  help get_ARGOifremer_files
  return
end

% Loop Through each float of interest--------------------------------------
%--------------------------------------------------------------------------
tf = cell(length(WMOs),1);
for i = 1:length(WMOs)
    
    WMO = WMOs{i};
    if isnumeric(WMO)
        WMO = num2str(WMO);
    end

    tf{i} = [1 1 1];

    % Connect to FTP server------------------------------------------------
    %----------------------------------------------------------------------
    ftp_dir = ['ifremer/argo/dac/',DAC,'/',WMO];
    ftp_target = 'ftp.ifremer.fr';
    disp(['CONNECTING TO FTP TARGET AT ',ftp_target,'...'])
    disp('-----------------------------------------------')
    f = ftp(ftp_target); 
    binary(f)

    % Enter passive mode by accessing the java object----------------------
    %----------------------------------------------------------------------
    cd(f);
    sf = struct(f);
    sf.jobject.enterLocalPassiveMode();
    ftp_path = ftp_dir;
    cd(f, ftp_path);   
    argo_dir   = dir(f);   
    argo_files = {argo_dir.name}'; % Cell array of ARGO netcdf files

    % get BRtraj file------------------------------------------------------
    %----------------------------------------------------------------------
    IndexC = strfind(argo_files, [WMO,'_BRtraj.nc']);
    Index = find(not(cellfun('isempty', IndexC)), 1);
    if ~isempty(Index)
        disp(['Grabbing ',WMO,'_BRtraj.nc from ',ftp_target,'...'])
        str = mget(f, [WMO,'_BRtraj.nc'], [savedir,WMO,'/']);
        disp([WMO,'_BRtraj.nc saved to ',savedir,WMO,'\'])
    else
        disp([WMO,'_BRtraj.nc NOT FOUND AT ',ftp_dir])
        tf{i}(1) = 0;
    end
    disp(' ')

    % get Mprof file-------------------------------------------------------
    %----------------------------------------------------------------------
    IndexC = strfind(argo_files, [WMO,'_Sprof.nc']);
    Index = find(not(cellfun('isempty', IndexC)), 1);
    if ~isempty(Index)
        disp(['Grabbing ',WMO,'_Sprof.nc from ',ftp_target,'...'])
        str = mget(f, [WMO,'_Sprof.nc'], [savedir,WMO,'/']);
        disp([WMO,'_Sprof.nc saved to ',savedir,WMO,'\'])
    else
        disp([WMO,'_Sprof.nc NOT FOUND AT ',ftp_dir])
        tf{i}(2) = 0;
    end
    disp(' ')

    % get meta file--------------------------------------------------------
    %----------------------------------------------------------------------
    IndexC = strfind(argo_files, [WMO,'_meta.nc']);
    Index = find(not(cellfun('isempty', IndexC)), 1);
    if ~isempty(Index)
        disp(['Grabbing ',WMO,'_meta.nc from ',ftp_target,'...'])
        str = mget(f, [WMO,'_meta.nc'], [savedir,WMO,'/']);
        disp([WMO,'_meta.nc saved to ',savedir,WMO,'\'])
    else
        disp([WMO,'_meta.nc NOT FOUND AT ',ftp_dir])
        tf{i}(3) = 0;
    end
    disp(' ')
    

    close(f)
end
% END ---------------------------------------------------------------------
%--------------------------------------------------------------------------


