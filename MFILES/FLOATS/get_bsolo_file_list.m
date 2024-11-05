function d = get_bsolo_file_list(cal_info, file_type)
%
% PURPOSE: 
%    This function returns a structure of float sensor file names, directories
%    & file time stamps. If a secondary float directory exists and the two 
%    directories share common files, only the most recent files are
%    retained. If unique files in the secondary directory exist, they are
%    merged with the primary directory file listing. 
%
% INPUTS:
%   cal_info  = cal.info structure from cal file which contains msg dir
%   file_type = string identifying file type (i.e 'ALK', 'CTD', 'DOX, 'ECO', 'NO3', 'OCR')
%
% OUTPUTS:
%   d.hdr  = cell array of column ID's (file name, dir path, file mod date)
%   d.list = msg file listing
%
% EXAMPLES:
% load(['C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\DATA\cal\',...
%     'calss0001.mat']);
% cal_info = cal.info;
%   d        = get_bsolo_file_list(cal_info, 'NO3')
%
% CHANGES:
% 01/25/2022 TM - Initialization of code for handling files from bsolo from SIO, modified from get_msg_list.m

% ************************************************************************
% ************************************************************************
% TESTING


% ************************************************************************
% ************************************************************************
%
% Extend path name based on file type assignment.
SOLOftype(1,:) = {'ALK alk'};
SOLOftype(2,:) = {'CTD phy'};
SOLOftype(3,:) = {'DOX dox'};
SOLOftype(4,:) = {'ECO eco'};
SOLOftype(5,:) = {'NO3 no3'};
SOLOftype(6,:) = {'OCR ocr'};

tCAL = ~cellfun(@isempty,regexp(SOLOftype, file_type,'once'));
tmp  = regexp(SOLOftype(tCAL),'\s+','split'); 
tmp = tmp{:};
dirext = tmp{1};
fileext = tmp{2};

filepath_root = cal_info.msg_dir;
filepath_full = [filepath_root,dirext,'\'];

% GET PRIMARY DIRECTORY LISTING

if isfolder(filepath_full)
    DirList1 = struct2cell(dir([filepath_full,'*.', fileext]))';
else
%     keyboard
    disp(['Could not find raw sensor file directory for: ',cal_info.name])
    d.hdr  = {'sensor file', 'directory' ,'file SDN'};
    d.list = {};
    return
end

d.hdr  = {'sensor file', 'directory' ,'file SDN'};
d.list = DirList1(:,[1,2,6]);
clearvars -except d








