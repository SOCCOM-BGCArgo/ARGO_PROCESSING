function out = parse_bad_ODVsensor_list(file_path)

% CHANGE HISTORY
% 08/28/2017 Added code to look for "NO_WMO" floats in addition to those with
%   WMO #'s -jp
% 01/10/2018 converted for use on non MBARI floats

% TEST
%file_path  = 'C:\Users\jplant\Documents\MATLAB\ARGO\ODVbad_sensor_list.txt';

% DO SOME PREP WORK
out.hdr    = []; % Output if function fails
out.list   = [];    

format_str = '%s%s%s%s';
cell_list  = cell(1000,4); % for pasing text file
list  = cell(1000,5); % for final output
max_cycles = 1000; % # way greater than the max # of cycles a float can do 

% PARSE BAD SENSOR LIST
if exist(file_path, 'file') ~= 2
    disp(['Bad sensor list not found at: ',file_path])
    return
end
fid   = fopen(file_path);
tline = ' ';
ct    = 0;

while ischar(tline)
    if regexp(tline, '^%', 'once')
        tline = fgetl(fid);
        continue
    elseif regexp(tline, '^WMO', 'once') % header line starts with WMO
        txt_hdr = regexp(tline, '\t', 'split');
    elseif regexp(tline, '^\d|^NO_WMO', 'once') % data lines start with WMO#
        ct = ct+1;
        cell_list(ct,:) = regexp(tline, '\t', 'split');
    end
    tline = fgetl(fid);
end
fclose(fid);
cell_list = cell_list(1:ct,:);
list = list(1:ct,:);

% CREATE FINAL OUTPUT
hdr  = [txt_hdr, 'CYCLE BLOCKS'];

% GET HEADER INDICES
iWMO = find(strcmp('WMO #',hdr) == 1);
iM   = find(strcmp('MBARI ID STR',hdr) == 1);
iSEN = find(strcmp('SENSOR',hdr) == 1);
iCYC = find(strcmp('CYCLES',hdr) == 1); % idividual bad profiles
iCB  = find(strcmp('CYCLE BLOCKS',hdr) == 1); % start block of bad

for i = 1:ct
    list(i,[iWMO,iM,iSEN]) = cell_list(i,1:3);
    cycle_str = regexprep(cell_list{i,iCYC}, ' ', ''); % remove any spaces
    nodes = regexp(cycle_str, ',', 'split');
    cycle_blocks =[]; cycles =[];
    for j = 1:size(nodes,2)
        if isempty(nodes{j})
            continue
        end
        if regexp(nodes{j},'\d+\-\d+')
            cycle_blocks = [cycle_blocks; (sscanf(nodes{j},'%f-%f'))'];
        elseif regexp(nodes{j},'\d+\-')
             cycle_blocks = [cycle_blocks; (sscanf(nodes{j},'%f-'))', 2000];  
        elseif regexp(nodes{j},'\-\d+')
             cycle_blocks = [cycle_blocks; 2000, (sscanf(nodes{j},'%-f'))'];                       
        else
            cycles = [cycles, sscanf(nodes{j},'%f')];
        end
        %pause
    end
    list{i,iCYC} = cycles;
    list{i,iCB} = cycle_blocks;
end

out.hdr = hdr;
out.list = list;
clearvars -except out
%         
