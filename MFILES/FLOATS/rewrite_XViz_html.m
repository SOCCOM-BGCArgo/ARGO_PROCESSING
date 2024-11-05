function status = rewrite_XViz_html(dirs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function rewrites all FloatViz associated HTML files and returns
% a cell array indicating the status of each rewrite job
% (This script compiles all float information to update the master list)
%
%
%
% Tanya Maurer & Josh Plant, MBARI
% This code is internal to the MBARI system (FloatViz web data plotting
% interface that is built and maintained at mbari)
% 3/24/21
% UPDATES:
% 10/27/21, TM, improved indexing on mapping of MBARI float WMO to Sharon's
% adopt-a-float list.
%
% UPDATES:
% 7/7/21 TM -- added float location to drop-down wmo-select in floatviz
% interface html
% 9/8/21 EC -- added lines that remove ua19061(WMO 5906341) and un1248
% (WMO 5906026) from float list and sharon_info
% 10/18/21 TM -- reverted exception for un1248 as this float was heard from
% last week for the first time!! :0)
% 11/04/2021 EC -- added OCR floats to the exception string
% 8/1/2022 TM - modified handling of Sharon's Adopt-a-float master list.
%   Also moved from including school & name on html, to name & wmo (per
%   George request)
% 09/12/2024 JP - updated excusion list to remove test OCR's and and a few
%                 more DOA float
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% SET DIRECTORIES
%--------------------------------------------------------------------------

if isempty(dirs) % TESTING SET UP
    user_dir        = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\bgcargovm'
    user_dir        = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.cal        = [user_dir,'CAL\']; % save created lists here
    dirs.temp       = 'C:\temp\'; % path to temporary working dir
    dirs.sirocco    = '\\sirocco\wwwroot\';
    % dirs.sirocco    = dirs.temp
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

%--------------------------------------------------------------------------
% CELL ARRAY OF FILES TO UPDATE (file name & subdir) & MISC PREP
%--------------------------------------------------------------------------
html_list = ...
    {'FloatViz.htm'       , 'chemsensor\'; ... % FloatViz
    'default.htm'        , 'soccom\'; ...     % SOCCOMViz
    'GOBGCViz.htm'       , 'gobgc\'; ...      % GOBGC
    'EViz.htm'           , 'chemsensor\'; ... % EViz
    'NViz.htm'           , 'chemsensor\'; ... % NViz
    'AdoptAFloatViz.htm' , 'soccom\'};        % ADOPT-A-FLOAT

rhtml       = size(html_list,1); % save dimensions of html list
status      = cell(rhtml,2); % make another cell array of the same size
status(:,1) = html_list(:,1); % add html file names
status(:,2) = {0};             % set to 0 & change to 1 for success

% \t PUTS TABS IN THE TEMPLATE - TEMPLATE USED IN FPRINTF
template = [repmat('\t',1,9),'<option value="%s" >%s</option>\r\n'];
% creates a structure

%--------------------------------------------------------------------------
% LOAD AND ORGANIZE MBARI FLOAT LIST
%--------------------------------------------------------------------------
load([dirs.cal,'MBARI_float_list.mat']); % call float list builder
FLOAT_LIST = d.list; % populate a float list cell array

% create T/F values for location of column headers in float list
%iMSG    = strcmp('msg dir', d.hdr);
iWMO    = strcmp('WMO',d.hdr);
iMB     = strcmp('MBARI ID',d.hdr);
iINST   = strcmp('INST ID',d.hdr);
%iFLT    = strcmp('float type',d.hdr);
iPRG    = strcmp('Program',d.hdr);
iTFPH   = strcmp('tf pH',d.hdr);
iTFNO3  = strcmp('tf NO3',d.hdr);
iREG    = strcmp('Region',d.hdr);

[~,IX]     = sort(FLOAT_LIST(:,iWMO));
FLOAT_LIST = FLOAT_LIST(IX,:); % SORT BY WMO


% *************************************************************************
% subset float list to exclude certain floats: 1248, 19061
% deployed but doa

% LIST SUBSET TESTS
% Dead on deployment floats: ua19061|ua18340|ua21960|wn1353|un1526 JP 09/2024
exclude_str = '5906341|5906543|2903871|4903460|5907057'; % 

%exclude_str = '5906341|5906320|5906446'; % ua19061, WMO 5906341
%exclude_str = '5906026|5906341'; % ua19061, WMO 5906341 % un1248, WMO 5906026
tf_exclude = cellfun(@isempty, regexp(FLOAT_LIST(:,iWMO), exclude_str,'once'));
FLOAT_LIST = FLOAT_LIST(tf_exclude,:);
% *************************************************************************

tf_SOCCOM = strcmp(FLOAT_LIST(:,iPRG),'SOCCOM');
tf_GOBGC  = strcmp(FLOAT_LIST(:,iPRG),'GO-BGC');
tf_pH     = cell2mat(FLOAT_LIST(:,iTFPH)) == 1 & ...
    strncmp(FLOAT_LIST(:,iMB),'ua',2);
tf_NO3    = cell2mat(FLOAT_LIST(:,iTFNO3)) == 1;

%--------------------------------------------------------------------------
% LOAD SHARON'S MAT STRUCTURE FILE INFO, BUT TREAT AS TEXT DATA - FASTER
%--------------------------------------------------------------------------
% sharonMasterURL = 'http://soccom.ucsd.edu/FLOAT_INFO/AllFloatInfo.m'; % This file loc just includes SOCCOM
%sharonMasterURL = 'http://go-bgc.ucsd.edu/FLOAT_INFO/SOCCOM_GOBGC_FloatInfo.m';  % This file loc includes SOCCOM & GOBGC (Sharon cats the 2 files together for us now).
sharonMasterURL = 'https://www3.mbari.org/gobgc/tables/FLOAT_INFO/SOCCOM_GOBGC_FloatInfo.m';
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
    % TM 8/1/22, I think loading array is more straightforward; less prone
    % to regexp error when Sharon changes things.  The indexing was all
    % mucked up.  
    Sharon_SOCCOMMasterFloatInfo;
    Fnm = fieldnames(FLOAT);
    for idn = 1:length(Fnm)
        tmpflt = FLOAT.(Fnm{idn});
%         try
        fwmo = tmpflt.WMO;
%         catch
%             keyboard
%         end
        if isempty(tmpflt.Adopt) %not adopted!
            fschl = ' ';
            fnam = ' ';
        else
        fschl = tmpflt.Adopt{1};
        fnam = tmpflt.Adopt{4};
        end
        SHARON{idn,1} = fwmo;
        SHARON{idn,2} = fnam;
        SHARON{idn,3} = fschl;
    end
    
    
% % %     fid = fopen([dirs.temp,'Sharon_SOCCOMMasterFloatInfo.m'],'r','n','UTF-8');
% % %     tmp = textscan(fid,'%s','Delimiter','\r\n');
% % %     fclose(fid);
% % %     
% % %     % LOOK FOR filled Adopt field line & subset to only those lines
% % %     t1 = ~cellfun(@isempty, regexp(tmp{1,1},'^FLOAT\.F\d{7}\.Adopt(?=\s+\=\s*{'')','once'));
% % %     
% % %     sharon_info = tmp{1,1}(t1);
% % % 
% % %     tf_exclude = cellfun(@isempty, regexp(sharon_info, exclude_str,'once')); % JP 09/09/21
% % %     sharon_info = sharon_info(tf_exclude, :);
% % %     
% % %     % FIND WMO
% % %     adopt_WMO = regexp(sharon_info,'\d{7}','once','match');
% % %     
% % %     % GET SCHOOL USING LAZY EXPRESSION (GET'S FIRST OCCURANCE)
% % %     % char(39) is a single quote, looking for " ',' "
% % %     school_pat   = ['(?<={',char(39),').*?(?=',char(39),',',char(39),')'];
% % %     
% % %     adopt_school = regexp(sharon_info, school_pat,'match','once');
% % %     
% % % % % 
% % % % %     % GET ADOPTED FLOAT NAME bounded  by " ',' " and " ','http "
% % % % %     name_pat = ['(?<=',char(39),'\s*,',char(39),').*?(?=',char(39),'\s*,',char(39),'http)'];
% % % % %     adopt_name = regexp(sharon_info, name_pat,'match','once');
% % % % %     TM 8/1/2022 - this regexp is too prone to error if Sharon does a
% % % % quick - modification (ie if no html!  she will enter something else.
% % % % Using "split" function might be more robust?
% % %     for ian = 1:size(sharon_info,1)
% % %         tmpss = sharon_info{ian,:};
% % %         adopt_name_tmp = split(tmpss,',');
% % %         ANT = adopt_name_tmp{4};
% % %         ANT2 = ANT(2:end-1);
% % %         Aname{ian,:} = ANT2;
% % %     end
% % % 
% % %     % REPLACE ALL CHARS UP TO LAST " ',' ", LEAVES NAME REMAINING
% % % %     adopt_name   = regexprep(adopt_name,['.+',char(39),'',char(39)],'');
% % %     adopt_name = Aname;
% % %     adopt_info = [adopt_WMO, adopt_name, adopt_school];
% % %     adopt_info
% % %     pause
    % SORT BY SHARON'S LIST BY WMO FIRST (MASTER LIST ALREADY IS)
    [~,IX]       = sort(SHARON(:,1)); % 
    adopt_info   = SHARON(IX,:);
    % REMOVE empties!  Floats with no adoption.
    [nrows,ncols] = cellfun(@size,adopt_info(:,2));
    xxx = find(ncols==1);
    adopt_info(xxx,:) = [];
%     [tf_adopt] = ismember(FLOAT_LIST(:,iWMO),adopt_WMO); % logical aray size of FLOAT_LIST
    [~,IA,IB] = intersect(FLOAT_LIST(:,iWMO),adopt_info(:,1)); %TM 10/27/21 use intersect in case of duplicate float nums!
    adopt_list = FLOAT_LIST(IA,:);
    TMP_adopt_info = adopt_info(IB,:);
    [~,IX]  = sort(TMP_adopt_info(:,2)); % sort by adopted name
    
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
    if strcmp(html_list{file_ct,1},'default.htm') % SOCCOM
        LIST = LIST(tf_SOCCOM,:);
    elseif strcmp(html_list{file_ct,1},'GOBGCViz.htm') % GOBGC
        LIST = LIST(tf_GOBGC,:);
    elseif strcmp(html_list{file_ct,1},'EViz.htm') % EViz
        LIST = LIST(tf_pH,:);
        file_end = '_DURA.TXT';
    elseif strcmp(html_list{file_ct,1},'NViz.htm') % NViz
        LIST = LIST(tf_NO3,:);
        file_end = '_NO3.TXT';
    elseif strcmp(html_list{file_ct,1},'AdoptAFloatViz.htm') && ...
            adopt_success == 1% SOCCOM Adopt
        LIST = adopt_list;
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
        if contains(line,'ENDLIST')
            k = 2;
        end
        if k == 1
            continue
        end
        fprintf(fid, '%s\r\n', line) ;
        if contains(line,'BEGINLIST')
            k=1;
            for fct = 1:size(LIST,1)
                WMOID     = LIST{fct,iWMO};
                MBARI_ID  = LIST{fct,iMB};
                REGION = LIST{fct,iREG};
                html_ID   = repmat('.',1,30); %  predim str with periods
                ind       = size(html_ID,2) - size(REGION,2)+1;
                
                % BUILD FILE NAME STR
                file_name = [WMOID, file_end];
                
                % BUILD HTML LIST ID STR
                if strcmp(html_list{file_ct,1},'AdoptAFloatViz.htm') && ...
                        adopt_success == 1% SOCCOM Adopt
                    t1 = strcmp(adopt_info(:,1),WMOID); % a match on both lists
                    if sum(t1) == 1
%                         html_ID = [adopt_info{t1,2},'....',adopt_info{t1,1}];
                        htsize = length(adopt_info{t1,2})+7;
                        numdot = 33-htsize;
                        html_ID = [adopt_info{t1,2},repmat('.',1,numdot),adopt_info{t1,1}];
                    end
                elseif regexp(WMOID,'^NO','once')
                    html_ID(1:6) = 'NO WMO';
                    html_ID(14:13+size(MBARI_ID,2)) = MBARI_ID;
                    html_ID(ind:end) = REGION;
                else
                    html_ID(1:7) = WMOID;
                    html_ID(15:14+size(MBARI_ID,2)) = MBARI_ID;
                    html_ID(ind:end) = REGION;
                end
                
                
                % PRINT LINE TO FILE
                fprintf(fid, template, file_name, html_ID);
            end
        end
    end
    fclose(fid);
    fclose(fid1);
    status{file_ct,2} = 1;
    if status{file_ct,2} == 1 % Successful HTML file re-write
        fprintf('%s successfully updated at: %s\n', ...
            html_list{file_ct,1}, FVfile);
            copyfile(TMPfile, FVfile); % THIS UPDATES FILE ON SIROCCO %TEST
    else
        fprintf('%s Could not be updated at: %s\n', ...
            html_list{file_ct,1}, FVfile);
    end
end



