function WMO  = uw2wmo(UW_num)
% Quick function to return WMO % for MBARI floats based on table built by
% get_MBARI_WMO_list.m
%
% UW_ID is UW float number as number or string


data_file = ['C:\Users\jplant\Documents\MATLAB\ARGO\DATA\Cal_files\', ...
             'float_WMO_list.mat'];
         
if exist(data_file,'file')
    load(data_file)
else
    disp(['Could not find: ',data_file]);
end

if ~ischar(UW_num)
    flt_num_str = sprintf('%04.0f',UW_num);
else
    flt_num_str = UW_num;
end

tf = strcmpi(flt_num_str,UW_ID);

WMO = WMO_ID{tf};

if ~isempty(WMO)
    disp([flt_num_str,'  ',float_names{tf}, '  ', WMO]);
else
    disp(['WMO ID # not found for ',flt_num_str]);
end



