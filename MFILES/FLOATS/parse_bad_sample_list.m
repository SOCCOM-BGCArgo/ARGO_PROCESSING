function out = parse_bad_sample_list(file_path)

% parse_bad_sample_list.m
% Tanya Maurer
% MBARI% 3/31/22
%
% Modified from parse_bad_sensor_list.match
%
% This file parses the bad_sample_list.txt which lists a record of bad sample points for a given, profile, given sensor, given float				
% The intent is that this updated code and txt file will replace the older method of manually modifying sirocco files!
% In the bad sample list:
%------------------------
% SENSOR OPTIONS ARE:
%	S    = Salinity
%	T    = Temperature
%	O    = Oxygen
%	N    = Nitrate
%	PH   = pH
%	CHL  = Chlorophyll
%	BBP  = Particle backscatter
%	CDOM = Colored Disolved Organic Matter
% CYCLE is the cycle number holding bad data (note: one cycle number per entry!!)
% DEPTH is a character string representing depth range of bad data for a given cycle, sensor
% Valid string examples for the DEPTH column can be mixed:
%------------------------

% TEST
% fdir = 'C:\Users\bgcargovm\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL';
% fn = 'bad_sample_list.txt';
% file_path = fullfile(fdir, fn);

% DO SOME PREP WORK
out.hdr    = []; % Output if function fails
out.list   = [];    

format_str = '%s%s%s%s';
cell_list  = cell(1000,6); % for pasing text file
list  = cell(1000,6); % for final output
max_depth = 10000; % # way greater than the max depth a float can do 

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

% CREATE FINAL OUTPUT
hdr  = [txt_hdr(1:5), 'DEPTH BLOCKS', txt_hdr(6)];

% GET HEADER INDICES
iWMO = find(strcmp('WMO #',hdr) == 1);
iM   = find(strcmp('MBARI ID STR',hdr) == 1);
iSEN = find(strcmp('SENSOR',hdr) == 1);
iCYC = find(strcmp('CYCLE',hdr) == 1); % idividual bad profiles
iD = find(strcmp('DEPTH',hdr) == 1); % individual fliers
iDB = find(strcmp('DEPTH BLOCKS',hdr) == 1); % start block of bad
iFLAG = find(strcmp('FLAG',hdr) == 1); %flag (ODV style, 8 = 'bad', 4 = 'questionable')

for i = 1:ct
    list(i,[iWMO,iM,iSEN,iCYC]) = cell_list(i,1:4);
    list(i,iFLAG) = cell_list(i,end);
    depth_str = regexprep(cell_list{i,iD}, ' ', ''); % remove any spaces
    
    % Get indvidual profile list - look for a block of #'s separated by
    % commas - flexible enough to catch cases where trailing comma is
    % forgotten on last entry.

    if ~isempty(regexp(depth_str,'(\d+(?=,))|((?<=,)\d+)', 'once'))
        list{i,iD} = cellfun(@str2double,regexp(depth_str, ...
            '(\d+(?=,))|((?<=,)\d+)','match'));
    end   
    % OK NOW PROCESS BAD PROFILE BLOCKS IF ANY - LOOK for dash
    find_dash = regexp(depth_str,'-');
    dash_ct   = max(size(find_dash));
    depth_blocks =[];
    if max(size(find_dash)) > 0
        depth_blocks = ones(dash_ct,2) * NaN;
        for j = 1:dash_ct;
            if ~isempty(regexp(depth_str,'\d+\-\d+', 'once')) % X to X
                depth_blocks(j,1) = str2double(regexp(depth_str, ...
                    '\d+(?=-)','match','once'));
                depth_blocks(j,2) = str2double(regexp(depth_str, ...
                    '(?<=-)\d+','match','once'));
            elseif ~isempty(regexp(depth_str,'\-\d+', 'once')) % 1 to X
                depth_blocks(j,1) = 1;
                depth_blocks(j,2) = str2double(regexp(depth_str, ...
                    '(?<=-)\d+', 'match','once'));
            elseif ~isempty(regexp(depth_str,'\d+\-', 'once')) % X to end
                depth_blocks(j,1) = str2double(regexp(depth_str, ...
                    '\d+(?=-)','match','once'));
                depth_blocks(j,2) = max_depth; % set to unreal high cycle #
            end
            depth_str(find_dash(j)) = 'X'; % Replace dash with 'X'
        end
    end
    list{i,iDB} = depth_blocks;

        
end
out.hdr = hdr;
out.list = list(1:ct,:);
clearvars -except out
        
