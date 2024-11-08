function [tf, theflag] = isbadsensor(list, MBARI_ID_str, cycle, sensor)
% Check bad sensor list parsed with parse_bad_sensor_list.m
% and see if the current profile is on the list for the given sensor.
% Possitve test = 1, negative test = 0
%
% INPUTS:
%   list         -  a structure (hdr & cell array of bad sensor list)
%                   output or subset from parse_bad_sensor_list.m
%   MBARI_ID_str - 	mbari float name
%   cycle        - 	float cycle (profile #)
%   sensor       - 	sensor ID string
%                   valid options: (S, T, O, N, PH, CHL, BBP, CDOM)
% OUTPUT
%   tf           - 	a flag, 1 if true, 0 if flase

% TESTING
% file_path  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
%     'DATA\CAL\bad_sensor_list.txt'];
% list = parse_bad_sensor_list(file_path);
% %MBARI_ID_str = '7615CALCURRENT';
% MBARI_ID_str = '9275SOOCN';
% cycle = 1;
% sensor = 'N';
%
% 09/11/2018:  TM update:  bad sensor list now supports QF = 'questionable', as well as 'bad' (QF = 3 is quest, QF = 4 is bad, ARGO FORMAT!!!)

% ************************************************************************
tf =  0;
theflag = 0;

% no data?  - happens if list subsetted or float not on list - job done
if isempty(list.list) 
    return
end

hdr  = list.hdr;
list = list.list;


% GET HEADER INDICES
iWMO = find(strcmp('WMO #',hdr) == 1);
iM   = find(strcmp('MBARI ID STR',hdr) == 1);
iSEN = find(strcmp('SENSOR',hdr) == 1);
iCYC = find(strcmp('CYCLES',hdr) == 1); % idividual bad profiles
iCB  = find(strcmp('CYCLE BLOCKS',hdr) == 1); % start block of bad
iFLAG = find(strcmp('FLAG',hdr) == 1 ); %flag to apply, per list

% TEST TO SEE IF FLOAT & SENSOR ARE ON THE LIST
t1 = strcmpi(list(:,iM), MBARI_ID_str)  & strcmp(list(:,iSEN), sensor);
t1I = find(t1 == 1);
if sum(t1) == 0 % list exsists but float + sensor are not on it
    return
end

clear t1
for iF = 1:length(t1I)
    t1 = t1I(iF);
    % FIND FLAGS FOR ENTRIES. 
    tmpflag = str2num(list{t1,iFLAG}); %supports either 3 (quest) or 4 (bad)
    
    % BUILD CYCLE COMPARISON ARRAY
    cycle_list   = list{t1,iCYC};

    % BUILD ARRAY FOR BLOCK END POINTS
    cycle_blocks = list{t1,iCB};
    if ~isempty(cycle_blocks)
        blocks = [];
        for i = 1:size(cycle_blocks,1)
            blocks = [blocks, cycle_blocks(i,1):cycle_blocks(i,2)];
        end
        cycle_list = unique([cycle_list, blocks]);
    end

    cycle_chk  = intersect(cycle,cycle_list); % is profile in list?
    if ~isempty(cycle_chk)
        tf =1;
        if tmpflag>theflag %if on list twice, once for bad, once for questionable, then the 'bad' flag trumps.
            theflag = tmpflag;
        end
    end
end
clearvars -except tf theflag


