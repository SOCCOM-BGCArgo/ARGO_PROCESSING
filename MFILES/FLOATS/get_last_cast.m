function [last_cast, missing_casts] = get_last_cast(data_dir,WMO)
% HELPER FUCTION FOR ARGO DATA PROCESSING
% This function looks for *.mat cast files in a given float directory
% named by its WMO ID # to determine the most recent cast and if there are
% any missing casts

% TEST ******
% user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
% user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\'];
% data_dir = [user_dir,'DATA\FLOATS\'];
% WMO = '5904188'; % 9095
% *********

last_cast     = [];
missing_casts = [];

if exist([data_dir,WMO,'\'],'dir') == 7 % 7 means it is a dir
    tmp_fnames = ls([data_dir, WMO,'\*.mat']);
    rr         = size(tmp_fnames,1); % rows
    cast_nums  = ones(rr,1) * NaN; 
    for i = 1:rr
        cast_str = regexp(tmp_fnames(i,:), ...
            '(?<=\.)\d+(?=\.)','once','match'); % cast # from file name
        if ~isempty(cast_str)
            cast_nums(i) = sscanf(cast_str,'%f');
        end
    end
    
    cast_nums(isnan(cast_nums)) =[];
    last_cast = max(cast_nums);
    for i = 1:max(cast_nums)
         t1 = sum(cast_nums == i);
         if t1 == 0
             missing_casts = [missing_casts; i];
         end
    end
    clear t1 i cast_nums cast_str rr tmp_fnames WMO data_dir
end