function status = rewrite_XViz_txt(dirs)

% This function rewrites all FloatViz associated TXT file file name lists
% and returns a cell array indicating the status of each rewrite job
% Tanya Maurer & Josh Plant, MBARI
% This code is internal to the MBARI system (FloatViz web data plotting interface that is built and maintained at mbari)
% 3/24/21

%--------------------------------------------------------------------------
% SET DIRECTORIES
%--------------------------------------------------------------------------
%dirs = [];
if isempty(dirs) % TESTING SET UP
    user_dir        = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir        = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.cal        = [user_dir,'CAL\']; % save created lists here
    dirs.temp       = 'C:\temp\'; % path to temporary working dir
    dirs.sirocco    = '\\sirocco\wwwroot\';
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

%--------------------------------------------------------------------------
% CELL ARRAY OF FILES TO UPDATE (file name & subdir) & MISC PREP
%--------------------------------------------------------------------------
html_list = ...
   {'FloatVIZConfig.txt' , 'lobo\Data\FloatVizData\'; ... % FloatViz
    'SOCCOMVIZConfig.txt', 'lobo\Data\FloatVizData\'; ... % SOCCOMViz
    'GOBGCVIZConfig.txt' , 'lobo\Data\FloatVizData\'; ... % GOBGC
    'EVIZConfig.txt'     , 'lobo\Data\EVizData\'; ...     % EViz
    'NVIZConfig.txt'     , 'lobo\Data\NVizData\'};        % ADOPT-A-FLOAT

rhtml       = size(html_list,1);
status      = cell(rhtml,2);
status(:,1) = html_list(:,1); % add html file names
status(:,2) = {0};             % set to 0 & change to 1 for success

%--------------------------------------------------------------------------
% LOAD AND ORGANIZE MBARI FLOAT LIST
%--------------------------------------------------------------------------
load([dirs.cal,'MBARI_float_list.mat']);
FLOAT_LIST = d.list;

%iMSG    = strcmp('msg dir', d.hdr);
iWMO    = strcmp('WMO',d.hdr);
iMB     = strcmp('MBARI ID',d.hdr);
iINST   = strcmp('INST ID',d.hdr);
%iFLT    = strcmp('float type',d.hdr);
iPRG    = strcmp('Program',d.hdr);
iTFPH   = strcmp('tf pH',d.hdr);
iTFNO3  = strcmp('tf NO3',d.hdr);

[~,IX]     = sort(FLOAT_LIST(:,iWMO));
FLOAT_LIST = FLOAT_LIST(IX,:); % SORT BY WMO

% LIST SUBSET TESTS
tf_SOCCOM = strcmp(FLOAT_LIST(:,iPRG),'SOCCOM');
tf_GOBGC  = strcmp(FLOAT_LIST(:,iPRG),'GO-BGC');
tf_pH     = cell2mat(FLOAT_LIST(:,iTFPH)) == 1 & ...
            strncmp(FLOAT_LIST(:,iMB),'ua',2);
tf_NO3    = cell2mat(FLOAT_LIST(:,iTFNO3)) == 1;   

%--------------------------------------------------------------------------
% LOAD SHARON'S MAT STRUCTURE FILE INFO, BUT TREAT AS TEXT DATA - FASTER
%--------------------------------------------------------------------------
% sharonMasterURL = 'http://soccom.ucsd.edu/FLOAT_INFO/AllFloatInfo.m'; %SOCCOM only
sharonMasterURL = 'http://go-bgc.ucsd.edu/FLOAT_INFO/SOCCOM_GOBGC_FloatInfo.m';  % This file loc includes SOCCOM & GOBGC (Sharon cats the 2 files together for us now).

adopt_success = 0;
try
    outfile = websave([dirs.temp,'Sharon_SOCCOMMasterFloatInfo.m'],sharonMasterURL);
    adopt_success = 1;
catch
    disp(['ERROR CONNECTING TO ',sharonMasterURL,'. CANNOT UPDATE SOCCOM ADOPTAFLOATVIZ HTML.'])
end

if adopt_success == 1
    %run([dirs.temp,'Sharon_SOCCOMMasterFloatInfo.m'])
    
    % PARSE TEXT FILE LOADED AS CELL ARRAY - MUCH FASTER THAN GENERATING STRUCTURE
    % HAD to OPEN IN UTF-8 ENCODING BECAUSE OF A COUPLE CHARACTER DIFFERNECES
    fid = fopen([dirs.temp,'Sharon_SOCCOMMasterFloatInfo.m'],'r','n','UTF-8');
    tmp = textscan(fid,'%s','Delimiter','\r\n');
    fclose(fid);
    
    % LOOK FOR filled Adopt field line & subset to only those lines
    t1 = ~cellfun(@isempty, regexp(tmp{1,1},'^FLOAT\.F\d{7}\.Adopt(?=\s+\=\s*{'')','once'));
    sharon_info = tmp{1,1}(t1);
    
    % FIND WMO
    adopt_WMO = regexp(sharon_info,'\d{7}','once','match');
    
    % GET SCHOOL USING LAZY EXPRESSION (GET'S FIRST OCCURANCE)
    % char(39) is a single quote, looking for " ',' "
    school_pat   = ['(?<={',char(39),').*?(?=',char(39),',',char(39),')'];
    adopt_school = regexp(sharon_info, school_pat,'match','once');
    
    % GET ADOPTED FLOAT NAME bounded  by " ',' " and " ','http "
    name_pat = ['(?<=',char(39),'\s*,',char(39),').*?(?=',char(39),'\s*,',char(39),'http)'];
    adopt_name = regexp(sharon_info, name_pat,'match','once');
    % REPLACE ALL CHARS UP TO LAST " ',' ", LEAVES NAME REMAINING
    adopt_name   = regexprep(adopt_name,['.+',char(39),',',char(39)],'');
    adopt_info = [adopt_WMO, adopt_name, adopt_school];
    % SORT BY SHARON'S LIST BY WMO FIRST (MASTER LIST AKREADY IS)
    [~,IX]       = sort(adopt_info(:,1)); % sort by adopted name
    adopt_info   = adopt_info(IX,:);

    tf_adopt = ismember(FLOAT_LIST(:,iWMO),adopt_WMO); % logical aray size of FLOAT_LIST
    adopt_list = FLOAT_LIST(tf_adopt,:);
    
    [~,IX]  = sort(adopt_info(:,2)); % sort by adopted name
    adopt_list = adopt_list(IX,:); % main adopt list now sorted by name
    
    clear fid tmp ans t1 sharon_info school_pat name_pat
    clear tf_adopt adopt_school adopt_WMO adopt_name IX
end

%--------------------------------------------------------------------------
% LOOP THROUGH HTML FILES & REWRITE
%--------------------------------------------------------------------------
for file_ct = 1:rhtml
    %fprintf('Rewriting %s ....\n',html_list{file_ct,1});
    % Define files and move to local
    FVfile   = [dirs.sirocco, html_list{file_ct,2}, html_list{file_ct,1}];
    FVlocal  = [dirs.temp, html_list{file_ct,1}];
    TMPfile  = [dirs.temp,'temp_', html_list{file_ct,1}];
    
    % SUBSET LISTS TO PROGRAM IF NEEDED & DEFINE FILE SUFFUX
    LIST = FLOAT_LIST;
    file_end = '.TXT'; % file appender for EViz / NViz
    if strcmp(html_list{file_ct,1},'SOCCOMVIZConfig.txt') % SOCCOM
        LIST = LIST(tf_SOCCOM,:);
    elseif strcmp(html_list{file_ct,1},'GOBGCVIZConfig.txt') % GOBGC
        LIST = LIST(tf_GOBGC,:);
    elseif strcmp(html_list{file_ct,1},'EVIZConfig.txt') % EViz
        LIST = LIST(tf_pH,:);
        file_end = '_DURA.TXT';
    elseif strcmp(html_list{file_ct,1},'NVIZConfig.txt') % NViz
        LIST = LIST(tf_NO3,:);
        file_end = '_NO3.TXT';
    end
    
    % COPY HTML FILE TO LOCAL & REWRITE FLOAT INFO
    tf       = copyfile(FVfile, FVlocal,'f');
    if tf ~=1
        fprintf(['%s COPY ERROR FROM SIROCCO TO LOCAL!! SKIPPING OVER ', ...
            'REBUILD OF FLOATVIZ.HTM\n'], html_list{1})
        continue
    end
    
    fid1 = fopen(FVlocal,'r','n','windows-1252');
    fid  = fopen(TMPfile,'w','n','windows-1252');
    k = 0;
    while ~feof(fid1)
        line = fgetl(fid1);
        if regexp(line,'\#VARIABLES|\#END','once')
        %if contains(line,'ENDLIST')
            k = 2;
        end
        if k == 1
            continue
        end
        fprintf(fid, '%s\r\n', line) ;
        if contains(line,'#Files')
            k=1;
            for fct = 1:size(LIST,1)
                WMOID     = LIST{fct,iWMO};

                % BUILD FILE NAME STR
                file_name = [WMOID, file_end];
                
                % PRINT LINE TO FILE
                fprintf(fid, '%s\r\n', file_name);
            end
        end
    end
    fclose(fid);
    fclose(fid1);
    status{file_ct,2} = 1;
    if status{file_ct,2} == 1 % Successful HTML file re-write
        fprintf('%s successfully updated at: %s\n', ...
            html_list{file_ct,1}, FVfile);
        copyfile(TMPfile, FVfile); % THIS UPDATES FILE ON SIROCCO
    else
        fprintf('%s Could not be updated at: %s\n', ...
            html_list{file_ct,1}, FVfile);
    end
end 

    
     
     