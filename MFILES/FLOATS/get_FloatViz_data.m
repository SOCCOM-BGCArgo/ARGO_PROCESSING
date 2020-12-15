function floatviz_data = get_FloatViz_data(floatviz_file)
% PURPOSE: 
%   This function parses the ODV compatable text files used by FloatViz
%   into Matlab variables
%
% USAGE:
%	data = get_FloatViz_data(floatviz_file)
%   data = get_FloatViz_data(path\floatviz_file)
%
% INPUTS:
%   floatviz_file = MBARI floatname ID or direct path\MBARI floatname. If a 
%   direct path is used it will override the default directory in network
%   mode
%
% OUTPUTS: a structure
%             .hdr  = cell array of column headers
%             .data = a matrix of the float data
%
% EXAMPLES:
%    get_FloatViz_data('9031SoOcnQC')
%    get_FloatViz_data('C:\temp\6091SOOCNQC.txt')
%
% 
% CODE PROVIDER: 
%    MONTEREY BAY AQUARIUM RESEARCH INSTITUE (MBARI)
%    WWW.MBARI.ORG
%    MAR2017
%
% CHANGEES:
% 04/20/2017 - added code to automatically detect HR and HRQC files and
%   point to the correct directories on sirocco to get them. lines ~57-64
% 09/19/2018 - TM, lots of mods, enhance functionality across platforms and
%              for different usage (network, internet...)
% 10/07/2019 - JP added a close(fid) before some return statement if file exist
%     but no valid data. around lines 160 & 187
% 12/13/20 - TM - Forced fopen read to UTF-8, because that is the
%    new default for Matlab 2020 writes and better cross platform sharing
%    and all ODV compatible TXT files are  now UTF-8
%    

% *************************************************************************
% TESTING
% floatviz_file = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\ADJ\ODV1901467ADJ.TXT';

% *************************************************************************
% SET PATHS & COPY FILE TO LOCAL & OPEN
% *************************************************************************
%data_source = 'internet';
data_source = 'network'; %MBARI USE ONLY 

floatviz_url  = 'http://www3.mbari.org/lobo/Data/FloatVizData/';

% create temporary folder for housing files pulled from internet or network
if ~isempty(strfind(computer,'PC'))
    temp_dir = ['C:',filesep,'temp',filesep]; % change this filepath to wherever you want your temp dir to be
else
    temp_dir = [getenv('HOME'),filesep,'temp',filesep]; % creates temp dir in home path.  Can change this if preferred
end

if ~exist(temp_dir,'dir')
    mkdir(temp_dir) 
end

    
floatviz_data =[];

switch data_source
    % GET TEXT FILE FROM DIRECT PATH OR SIROCCO AND STORE LOCALY
    case 'network'     
        floatviz_dir  = '\\sirocco\wwwroot\lobo\Data\FloatVizData\'; %MBARI USE ONLY
        if regexp(floatviz_file,filesep) % direct path (dir included)
            from_str = floatviz_file;
            to_str   = [temp_dir, regexpi(floatviz_file, ...
                        '\d{3}\d+\w+\.txt','once','match')];
        elseif regexp(floatviz_file,'HRQC','once') % QC DIR
            floatviz_dir  = [floatviz_dir,'HRQC',filesep];
            from_str = [floatviz_dir, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt'];                    
        elseif regexp(floatviz_file,'HR','once') % QC DIR
            floatviz_dir  = [floatviz_dir,'HR',filesep];
            from_str = [floatviz_dir, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt'];                                
        elseif regexp(floatviz_file,'QC','once') % QC DIR
            floatviz_dir  = [floatviz_dir,'QC',filesep];
            from_str = [floatviz_dir, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt']; 
        else
            from_str = [floatviz_dir, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt'];   
        end

        if strcmpi(from_str, to_str) % destination and source equal
            disp(['File already exists at ',to_str, ...
                  ' .... Parsing existing file'])
        elseif exist(from_str,'file') == 2
            copyfile(from_str, to_str)
        else
            disp(['Could not find: ',from_str]);
%             disp(['SOURCE: ',from_str])
%             disp(['DESTINATION: ',to_str])
            return
        end
        
    case 'internet'
         if regexp(floatviz_file,filesep) % direct path (dir included)
            from_str = floatviz_file;
            to_str   = [temp_dir, regexpi(floatviz_file, ...
                        '\w{3}\d+\w+\.txt','once','match')];
         elseif regexp(floatviz_file,'HRQC','once') % QC DIR
            floatviz_url  = [floatviz_url,'HRQC/'];
            from_str = [floatviz_url, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt'];                    
        elseif regexp(floatviz_file,'HR','once') % QC DIR
            floatviz_url  = [floatviz_url,'HR/'];
            from_str = [floatviz_url, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt'];                                
        elseif regexp(floatviz_file,'QC','once') % QC DIR
            floatviz_url  = [floatviz_url,'QC/'];
            from_str = [floatviz_url, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt']; 
        else
            from_str = [floatviz_url, floatviz_file, '.txt'];
            to_str   = [temp_dir, floatviz_file, '.txt'];   
        end
        
        [~,url_chk] = urlread(from_str); % See if file exisit on the web
        if url_chk == 1
            f = urlwrite(from_str,to_str);
            disp(' ');
            disp(['Data for float ',floatviz_file,' retrieved from:  ',...
                floatviz_url]);
            disp(['Saved as  ',f]);
        else
            disp('No file found online!')
            
            if strcmpi(from_str, to_str) % destination and source equal
                disp(['File already exists at ',to_str, ...
                    ' .... Parsing existing file'])
            elseif exist(from_str,'file') == 2
                copyfile(from_str, to_str)
            else
                disp(['Could not find: ',from_str]);
                %             disp(['SOURCE: ',from_str])
                %             disp(['DESTINATION: ',to_str])
                return
            end
        end
end

%fid = fopen(to_str);
fid = fopen(to_str,'r','n','UTF-8');
% *************************************************************************
% BUILD FORMAT STRING AND PARSE DATA
% *************************************************************************
tline = ' ';
while ischar(tline)
    if regexp(tline,'^Cruise', 'once') % stop at header line
        break
    end
    tline = fgetl(fid);
end
if ~ischar(tline)
    disp('No header line found')
    fclose(fid);
    return
end

% Adjusted CHL files from U. Maine have an extra tab in header so remove if
% encontered
tline = regexprep(tline,'\t\t', '\t'); % remove extra tab if there
hdr      = regexp(tline,'\t','split'); % CELL ARRAY OF HEADER VARIABLES 
hdr_rows = size(hdr,2);


d_format = '';
rm_cols  = [];
for i = 1: hdr_rows
    if regexp(hdr{i},'^Cruise|^Type|^Bot\.', 'once')
        d_format = [d_format,'%*s'];
        rm_cols  = [rm_cols,i];
    elseif regexp(hdr{i},'^mon|^hh', 'once')
        d_format = [d_format,'%s'];
    else
        d_format = [d_format,'%f'];
    end
end
hdr(rm_cols) =[];
d     = textscan(fid,d_format,'Delimiter','\t','CollectOutput',1);

if isempty(d{1,1}) || isempty(d{1,3}) % no data exists only hdr
    disp(['Header exists but no data found for ',floatviz_file]);
    fclose(fid); %10/07/19 need to close down 1st
    return
end

d_tmp = strcat(d{1,2}(:,1),regexprep(d{1,2}(:,2), '(\d+:\d+)',' $1'));

% some float have too many tabs (eg 5146) on the profile separation line
% this creates a row of empty values across all cells. this causes datenum
% to burb so remove the offending lines
t_NaN = isnan(d{1,1}); % cast row - NaNs mark the spots
d{1,1} = d{1,1}(~t_NaN);
d_tmp = d_tmp(~t_NaN);
d{1,3} = d{1,3}(~t_NaN,:);

% for jp = 1:length(d_tmp) % FOR TESTING
%     disp(d_tmp{jp})
%     date_str = datenum(d_tmp{jp},'mm/dd/yyyy HH:MM');
% end

try
    sdn   = datenum(d_tmp,'mm/dd/yyyy HH:MM');
catch
    disp(['Problem with date strings in FloatViz file - Getting date ',...
          'number line by line'])
    sdn = ones(size(d_tmp)) * NaN;
    for i = 1 : size(d_tmp,1)
        if ~isempty(d_tmp{i})
            sdn(i) = datenum(d_tmp{i},'mm/dd/yyyy HH:MM');
        else
            disp(['Not a valid date string for cast # ', num2str(d{1,1}(i))])
        end
    end
end
% COMBINE DATA & UPDATE HEADER
hdr(2:3) = []; % remove date and time headers
hdr  = ['SDN', hdr];
data = [sdn, d{1,1} d{1,3}];

iP   = find(strcmp('Depth[m]',hdr) == 1, 1,'first');
t1   = isnan(data(:,iP));
data(t1,:) =[]; % Remove ODV profile sparation lines

% ASSIGN TO STRUCTURE
floatviz_data.hdr  = hdr;
floatviz_data.data = data;

% CLEAN UP
fclose(fid);
delete(to_str);
clearvars -except floatviz_data





