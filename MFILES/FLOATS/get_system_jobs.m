function d = get_system_jobs(image_name)
% THIS FUNCTION RETURNS A CELL ARRAY OF CURRENT PROCESSES RUNNING ON 
% YOUR COMPUTER BASED ON THE SUPPLIED image_name (ex 'MATLAB.EXE')
% IF image_name IS EMPTY ALL JOBS WILL BE RETURNED
%
% INPUTS:
%   image_name: can be empty or process name (ex 'MATLAB.EXE')
%
% OUTPUTS:
%   A stucture d
%       d.hdr - colum headers
%       d. list - a cell array of cell arrays (one for each column)


% TESTING
%image_name = 'MATLAB.exe';
%image_name = '';

% ************************************************************************
% GET JOB INFO FROM SYSTEM
if ~isempty(image_name)
    [~,tasks] = system(['tasklist /v /fi "IMAGENAME eq ',image_name,'"']);
else
    [~,tasks] = system('tasklist /v');
end

% ************************************************************************
% PARSE STRING TO LINES AND FIND COLUMN BREAKS
jobs = textscan(tasks, '%s', 'Delimiter', '\n'); % Parse lines to cell array
jobs = jobs{1,1};                                  % Remove outer skin
t1   = cellfun(@isempty,jobs); jobs(t1) = [];  % remove empty cells

% look for header break line, a bunch of '=' signs
% use this line to define line parsing format
t1 = strncmp('=', jobs, 1);  
breaks = (regexp(jobs{t1},'(?<=\=)\s{1}(?=\=)'))'; % '= space ='
breaks = [1; breaks; size(jobs{t1},2)];
jobs(t1) =[]; % remove line after use

% ************************************************************************
% BUILD CELL ARRAY FILLED WITH PARSED LINES FOR EACH COL
ncol  = size(breaks,1) - 1;
nrow = size(jobs,1); % ommit header and header divider rows
list = cell(1, ncol);
hdr = cell(1, ncol);
for i = 1 : nrow
    str = jobs{i};
    for j = 1 : ncol
        if i == 1
            hdr{1,j} = strtrim(str(breaks(j):breaks(j+1)));
        else
            list{1,j}{i-1,1} = strtrim(str(breaks(j):breaks(j+1)));
        end
    end
end

d.hdr = hdr;
d.list = list;

clearvars -except d
        
        
    

