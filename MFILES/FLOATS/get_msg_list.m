function d = get_msg_list(cal_info, file_type)
%
% PURPOSE: 
%    This function returns a structure of float msg file names, directories
%    & file time stamps. If a secondary float directory exists and the two 
%    directories share common files, only the most recent files are
%    retained. If unique files in the secondary directory exist, they are
%    merged with the primary directory file listing. 
%
% INPUTS:
%   cal_info  = cal.info structure from cal file which contains msg dir
%   file_type = string identifying file type (i.e 'msg', 'isus', 'dura')
%
% OUTPUTS:
%   d.hdr  = cell array of column ID's (file name, dir path, file mod date)
%   d.list = msg file listing
%
% EXAMPLES:
% load(['C:\Users\bgcargo\Documents\MATLAB\ARGO_PROCESSING\DATA\cal\',...
%     'calua12551.mat']);
% cal_info = cal.info;
%   d        = get_msg_list(cal_info, 'msg')
%
% CHANGES:
% 02/23/2021 JP - A complete overhaul in prep for GO-BGC file processing

% ************************************************************************
% ************************************************************************
% TESTING
% load(['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING_TESTDIR\DATA\cal\',...
%     'calua12396.mat']);
% cal_info = cal.info;
% file_type = 'msg';

% ************************************************************************
% ************************************************************************
msg_path = cal_info.msg_dir;

% GET PRIMARY DIRECTORY LISTING
if isfolder(msg_path)
    DirList1 = struct2cell(dir([msg_path,'*.', file_type]))';
else
    disp(['Could not find msg file directory for: ',cal_info.name])
    d.hdr  = {'msg file', 'directory' ,'file SDN'};
    d.list = {};
    return
end

% CHECK FOR SECONDARY (ALTERNATE) DIRECTORY LISTING.
inst      = regexp(msg_path,'(?<=floats\\)\w+','match','once'); % institute
msg_path2 = regexprep(msg_path, inst,[inst,'\\alternate']);
DirList2  = {};
if isfolder(msg_path2)
    DirList2 = struct2cell(dir([msg_path2,'*.', file_type]))';
end
% IF SECONDARY (ALTERNATE) DIRECTORY LISTING EXISTS, COMPARE AND MERGE
% USE FILE MODIFICATION SIZE TO DECIDE
if ~isempty(DirList2)
    DirListApnd = cell(size(DirList2)); % predim for appending list
    alt_ct = 0;
    for fct = 1:size(DirList2,1) % step through alternate dir files
        t1 = strcmp(DirList1(:,1),DirList2{fct,1});
        if sum(t1) == 1 % file exists in both directories 
            if DirList2{fct,4} > DirList1{t1,4} % alt file is larger - keep
                alt_ct = alt_ct + 1;
                DirList1(t1,:) = DirList2(fct,:); 
            end
        elseif sum(t1) == 0 % file exists only in alternate - append
            DirListApnd(fct,:) = DirList2(fct,:);
        else
            disp(['Trouble with primary & secondary directory ', ....
                'comparison for: ',DirList2{fct,1}]);
        end
    end
    % Add any file only in secondary listing to primary listing
    tg = ~cellfun(@isempty, DirListApnd(:,1));
    if sum(alt_ct) > 0
        fprintf(['%0.0f newer %s file(s) in the secondary directory ', ...
            'listing have replaced older files in the primary ', ...
            'listing\n'],sum(alt_ct), file_type);
    end
    if sum(tg) > 0
        DirList1 = [DirList1; DirListApnd(tg,:)];
        fprintf(['%0.0f %s file(s) from the secondary directory ', ...
            'listing have been appended to the primary listing\n'], ...
            sum(tg), file_type);
        [~,IX]   = sort(DirList1(:,1));
        DirList1 = DirList1(IX,:);
    end
end

d.hdr  = {'msg file', 'directory' ,'file SDN'};
d.list = DirList1(:,[1,2,6]);
clearvars -except d








