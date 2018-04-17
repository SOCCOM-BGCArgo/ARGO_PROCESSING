%getfloatlist.m

% ************************************************************************
% FILTER EXPRESSIONS
% ************************************************************************
refine_expr  = ''; % REFINE FLOAT LIST FOR REGION
%refine_expr  = 'SoOcn|Drake|Ross|SoAtl'; % REFINE FLOAT LIST FOR REGION
flt_exp      = '^\d{4}\w+\.TXT'; % 4#'s,char,.txt(all float names)
exclude_expr = '(MTY)|(cor)|(Surface0)|(\d+\.txt)';

% ************************************************************************
% PATHS & NAMES
% ************************************************************************
file_name = 'FloatVIZConfig.txt';
float_dir = '\\SIROCCO\wwwroot\lobo\Data\FloatVizData\';
local_dir = ['C:\Documents and Settings\jplant\My Documents\' ...
    'MATLAB\PWP MODEL\Maps\data\'];

s1 = [float_dir,file_name]; % build target string
s2 = [local_dir,'floatlist.txt']; % build destination string
copyfile(s1,s2);

% ************************************************************************
% PARSE & FILTER & SORT
% ************************************************************************
fid  = fopen(s2);
d    = textscan(fid,'%s'); 
d    = d{1,1}; %txt file into cell array cell
fclose(fid);
delete(s2)

t1 = ~cellfun(@isempty,regexpi(d, flt_exp)); %(floats)
t2 = cellfun(@isempty,regexpi(d, exclude_expr)); % (unwanted floats)
d  = d(t1&t2);
if ~isempty(refine_expr)
   t3 = ~cellfun(@isempty,regexpi(d,refine_expr)); % (refined list)
   d  = d(t3);
end
d  = sort(d); % Output for script use
jp = d;       % will be cut and paste output

[r,c] = size(jp);
for i = 1:r
    jp(i) = regexprep(jp(i), '(^\d{1})','fname = ''$1');
    jp(i) = regexprep(jp(i), '\.TXT|\.txt|\.Txt',''';');
end

jp = char(jp); %List for getWOAdata

dd =d;
for i = 1:length(d)
%     dd(i) = {['float_info(',num2str(i),',:) ={''',d{i},'''', ...
%         '        1, 0, 5, 0, 0};']}; % for getfloat data
    
    dd(i) = {['flt_info(',num2str(i),',:) ={''',d{i},'''', ...
        '        0, 0, 0};']}; % for getfloat data
end
jp2 = char(dd);


