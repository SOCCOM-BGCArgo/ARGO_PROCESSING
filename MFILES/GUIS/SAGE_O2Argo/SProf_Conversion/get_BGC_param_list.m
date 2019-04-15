function list = get_BGC_param_list(filepath)
% This function parses a txt file made the 1st sheet in the
% argo-parameters-list-core-and-b.xlsx file
% found at: http://www.argodatamgt.org/Documentation
%
% To do this (only need to do this once unless the list changes):
%   1. open argo-parameters-list-core-and-b.xlsx
%   2. Select all => copy => paste special => values into a new sheet
%   3. Edit column 3 of header line "Parameter name short version Backup"
%       In the formula bar you will see it is actually 2 lines. Use
%       backspace before "Backup" to remove linefeed and get it back to 1
%       line
%   4. Delete all columns after "core/bio/intermediate" col
%   5. Delete all rows after the last "Order" number
%   6. Save as a tab delimited text file
%
% INPUT:
%   filepath = path\filename of text file
%
% OUTPUT:
%   list = a structure
%       .hdr  =  colum names
%       .list = cell array of core and BGC variables 
%
% CREATED 09/12/2017 bu JP


% PARSE PARAM LIST AND GET CORE AND BGC VARIABLE NAMES
fid = fopen(filepath);
tline = ' ';
while ischar(tline) % find header line
    if regexpi(tline,'^Order','once') % header line starts with "Order"
        break
    end
    tline = fgetl(fid);
end
param_hdr  = regexp(tline,'\t','split');
f_str      = repmat('%s', 1,size(param_hdr,2)); % build format str
param_list = textscan(fid,f_str,'Delimiter','\t','CollectOutput', 1);
param_list = param_list{1,1}; % remove outer shell

iType = find(strcmp(param_hdr,'core/bio/intermediate') == 1);
t1    = strcmp('c',param_list(:,iType)) | strcmp('b',param_list(:,iType));
param_list  = param_list(t1,:);

% ASIGN OUTPUT
list.hdr = param_hdr;
list.list = param_list;
clearvars -except list




