function out = process_data_ref_table
% Parse Sharon's SOCCOM data ref table % get data files

% HISTORY:
% 04/10/2019 Fixed CCHDO file time stamp accounting to properly choose
%   the most recent file in unmerged section if two valid chioices exist
% 05/14/2020 Multiple fixes to deal with non stadard "station number"
%   (CUSTARD ex = CTD030SS), cruise name with a period in name (i.e.
%   A13.5). -9999 vs - vs none as entries. BROKE CRUISE NOW CCAMLR in
%   special case assignments. Removv fill 0's from numeric strings
% 12/29/2020 - JP - Sharon's table changed so I updated some indices so
%   cruise, station & cycle values where being extracted properly again
% 03/31/21 - TM - tweaks due to changes in Sharon's table structure; plus
%    the addition of the GOBGC data reference table.
% 08/22/22 - TM - Temporary (?) fix for Marion Isl station entries in
%    Sharon's table (currently being listed as "AM012**" and that is breaking
%    the code.
% 01/10/23 - JP - updates to catch INST ID which was no longer extracted
%    properly because of some html formating changes / variants. Also some
%    general clean up of returned data
% 07/17/23 - TM - update links to data ref tables now that MBARI is managing.

out.ref_table        = 0;
out.mbari_names      = 0;
out.cchdo_file_urls  = 0;

% ************************************************************************
% SET UP & DEFINITIONS
% ************************************************************************
dirs = [];
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.cal       = [user_dir,'CAL\'];
    dirs.temp      = 'C:\temp\';
    
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end
save_dir = [user_dir,'SHIPBOARD\'];

non_soccom_fn = 'NON_SOCCOM_BottleData_lookup_table.txt';
data_url      = 'https://cchdo.ucsd.edu';
% ************************************************************************

% ************************************************************************
% GET SHARON'S DATA REFERENCE TABLE AS A STRING
% BOTH SOCCOM AND GOBGC DATA REF TABLES ARE REQUIRED!
% ************************************************************************
%data_ref_target1  = 'http://soccom.ucsd.edu/floats/SOCCOM_data_ref.html'; % This link is pre- Sharon-retirement
data_ref_target1  = 'https://www3.mbari.org/soccom/tables/SOCCOM_data_reference.html';
big_str1 = webread(data_ref_target1); % Sharon's HTML page as big string
%data_ref_target2  = 'http://go-bgc.ucsd.edu/GOBGC_data_reference.html'; % This link is pre- Sharon-retirement
data_ref_target2  = 'https://www3.mbari.org/gobgc/tables/GOBGC_data_reference.html';
big_str2 = webread(data_ref_target2); % Sharon's HTML page as big string
big_str = [big_str1 big_str2]; %Cat them together!

ind     = regexp(big_str,'\n'); % find newline characters
ind     = [0,ind]; % zero because it will be i+1 when subsetting
Nind    = size(ind,2);

hdr  = {'MBARI ID' 'INST ID' 'WMO' 'CRUISE' 'STATION' 'CAST' 'Data file' ...
    'CCHDO URL' 'Bottle file URL' 'CCHDO upload flag' 'CTD file URL'};

% DEFINE COLUMN INDEXES FROM HEADER
iM   = find(strcmp('MBARI ID',hdr)          == 1);
iINST  = find(strcmp('INST ID',hdr)             == 1);
iWMO = find(strcmp('WMO',hdr)               == 1);
iCRU = find(strcmp('CRUISE',hdr)            == 1);
iSTA = find(strcmp('STATION',hdr)           == 1);
iCST = find(strcmp('CAST',hdr)              == 1);
iDF  = find(strcmp('Data file',hdr)         == 1);
iCCH = find(strcmp('CCHDO URL',hdr)         == 1);
iBF  = find(strcmp('Bottle file URL',hdr)   == 1);
iCCF = find(strcmp('CCHDO upload flag',hdr) == 1);
iCTD = find(strcmp('CTD file URL',hdr)     == 1);

% ************************************************************************
% STEP THROUGH STRING AND FIND FLOATS
try
    disp(['Parsing data reference table from: ',data_ref_target1,' and ',data_ref_target2]);
    list = cell(1000,size(hdr,2)); % predim
    ct   = 0;
    for i = 1:Nind-1
        test_str = big_str(ind(i)+1: ind(i+1)-1);
        %         pause
        if regexp(test_str,'^\s+<tr|^\s+</tr', 'once') % format lines/no info/skip
            continue % Don't really need this line
            % Useful info lines start with "<td"
        elseif regexp(test_str,'^\s+<td', 'once')
            ct = ct+1;
            col_inds = regexp(test_str,'<td');
            % EXTRACT USEFUL INFO
            
            WMO = regexp(test_str,'(?<=>)\d{7}(?=<)', 'once', 'match');
            if isempty(WMO)
                WMO = 'NO WMO';
            end
            INST_str = test_str(col_inds(6):col_inds(7));
            %INST_ID = regexp(INST_str,'(?<=html" >).+(?=</a)', 'once', 'match');
            INST_ID = regexp(INST_str,'(?<=s*html\s*" >).+(?=</a)', 'once', 'match'); % jp 01/10/23
            if isempty(INST_ID) % maybe no shtml link because float DOA
                INST_ID = regexp(INST_str,'(?<=>\s*)\d+(?=<)', 'once', 'match');
            end

            cchdo = regexp(test_str,'https://cchdo\.\w+\.\w+/\w+/\w+', ...
                'once', 'match');
            if isempty(cchdo)
                cchdo = 'NO URL YET';
            end
            
            cruise_str = test_str(col_inds(2):col_inds(3));
            %cruise = regexp(cruise_str,'(?<=>)\w+(?=<)', 'once', 'match');
            cruise = regexp(cruise_str,'(?<=>).+(?=</td)', 'once', 'match');
            
            sta_str = test_str(col_inds(10):col_inds(11));
            sta_str = regexprep(sta_str,'none|uwDrop|-999|-','-999');
            station = regexp(sta_str,'(?<=>)\w+|(?<=>)\-\w+', 'once', 'match');
            
            % EXCEPTIIONS
            % CUSTARD CRUISE
            %             if regexp(station,'S+','once') % non standard "station number"
            %                 station = regexp(station,'[1-9]\d*(?=S+)','match','once');
            %             end
            
            % REMOVE FILL ZEROS IN NUMBER IF THEY EXIST & DEAL WITH
            % ALPHA-NUMERIC STATION NUMBERS FOR CUSTARD CRUISE
            if regexp(station,'^0|S+','once')
                station = regexp(station,'[1-9]\d*','match','once');
            elseif regexp(station,'AM','once')  % TM 8/22/22; temporary (?) fix for Marion Isl entries in Sharon's table. 
                station = '-999';
            elseif regexp(station,'P','once')  % TM 8/22/22; temporary (?) fix for Sikuliaq entry in Sharon's table. These are hacks!
                station = '-999';
            end
            
            cast_str = test_str(col_inds(11):col_inds(12));
            cast_str = regexprep(cast_str,'none|-999|-','-999');
            cast = regexp(cast_str,'(?<=>)\w+|(?<=>)\-\w+', 'once', 'match');
            
            % REMOVE FILL ZEROS IN NUMBER IF THEY EXIST
            if regexp(cast,'^0','once')
                cast = regexp(cast,'[1-9]\d*','match','once');
            end
            
            % ASSIGN OUTPUTS
            list(ct,[iINST, iWMO]) = {INST_ID, WMO};
            list(ct,[iCRU, iSTA, iCST]) = {cruise, station, cast};
            list(ct,iCCH) = {cchdo};
            
        end
        clear  WMO UW_ID cchdo cruise cast station test_str
    end
    list = list(1:ct,:);
    rlist = size(list,1);
    clist = size(list,2);

    out.ref_table        = 1;
catch
    disp('Reference tables could not be parsed - check link & try again.');
    disp(data_ref_target1);
    disp(data_ref_target2);
    return
end
clear big_str


% ***********************************************************************
% ADD MBARI FLOAT NAMES FROM MASTER LIST
% ***********************************************************************
try
    load([dirs.cal,'MBARI_float_list.mat']);
    MBARI = d;
    for i = 1:rlist
        tf = strcmp(list{i, iWMO}, MBARI.list(:,3)); % WMO
        tf2 = strcmp(list{i, iINST}, MBARI.list(:,2)); % INST ID
        if sum(tf) > 0
            list{i, iM} = MBARI.list{tf,1};
        elseif sum(tf2) > 0
            list{i, iM} = MBARI.list{tf2,1};
        else
            list{i, iM} = 'NO MBARI ID';
        end
    end
    out.mbari_names      = 1;
catch
    disp(['Could not load: ', dirs.cal,'MBARI_float_list.mat']);
    return
end

tg    = ~strcmp(list(:,1),'NO MBARI ID');
list  = list(tg,:);
rlist = size(list,1);

clear MBARI tf i

% ***********************************************************************
% TABLE MADE - NOW TRY AND GET LINKS TO DATA FILES AT CCHDO
% ***********************************************************************
try
    fprintf('Finding  CCHDO URLs to bottle & ctd files\n')
    loop_ct = 0;
    for i = 1:rlist
        loop_ct = loop_ct +1;
        if loop_ct > 60
            loop_ct = 1;
            fprintf('.\n')
        else
            fprintf('.')
        end
        bottle_flag = 1; % Unmerged Data as Received = 0, 1 toggles bopttle search
        %         disp(list{i,iM})
        if strcmp(list{i,iCCH},'NO URL YET') % No URL  so no data path or file
            list{i,iBF} = 'NO DATA PATH YET';
            list{i,iCCF} = '0';
            list{i,iCTD} = 'NO DATA PATH YET';
            %list{i,iDF} = 'NO DATA FILE';
            list{i,iDF} = '';
            continue
        end
        
        data_ref_target  = list{i,iCCH}; % CCHDO url
        big_str = webread(data_ref_target); % CCHDO page as big string
        
        % TRY TO FIND SOME RELAVENT INDICES BASED ON SECTIONS IN THE HTML CODE
        % bottle , ctd, unmerged data, documentation, history
        ind1 = regexp(big_str,'(?<=<li><h4>)bottle(?=</h4>\n)', 'once');
        ind2 = regexp(big_str,'(?<=<li><h4>)ctd(?=</h4>\n)', 'once');
        ind3 = regexp(big_str,'(?<=>)Unmerged Data as Received(?=</h3>\n)', 'once');
        ind4 = regexp(big_str,'(?<=<li><h4>)documentation(?=</h4>\n)', 'once');
        ind5 = regexp(big_str,'(?<=>)Data History(?=</h3>\n)', 'once');
        
        % ****************************************************************
        % LOOK FOR UNMERGED BOTTLE DATA FILES FIRST. IF THEY EXIST THESE ARE
        % THE ONES TO USE. IF NONE ARE FOUND TRY LOOKING FOR FILES IN THE
        % "BOTTLE" SECTION. hy1.csv (CCHDO file) normally takes priority
        % over exc.csv (Bob Key file) but this is not always the case
        % ****************************************************************
        % TRY UNMERGED DATA FILE SECTION FIRST
        if ~isempty(ind3) && ~isempty(ind5)
            cchdo_upload = 0;
            sub_str = big_str(ind3:ind5);
            ind     = regexp(sub_str,'\n'); % find newline characters
            ind     = [0,ind]; % zero because it will be i+1 when subsetting
            Nind    = size(ind,2);
            
            tmp_cell = cell(10,2); %predim
            ct = 0;
            time_flag = 0;
            for j = 1:Nind-1 % collect file names & dates
                test_str = sub_str(ind(j)+1: ind(j+1)-1);
                t1 = regexp(test_str,'hy1\.csv', 'once'); %priority
                t2 = regexp(test_str,'exc\.csv', 'once');
                
                if ~isempty(t1)
                    ct = ct+1;
                    quote_ind = regexp(test_str,'"');
                    file_path = test_str(quote_ind(1)+1:quote_ind(2)-1);
                    tmp_cell{ct,1} = file_path;
                    time_flag = 1;
                elseif ~isempty(t2)
                    ct = ct+1;
                    quote_ind = regexp(test_str,'"');
                    file_path = test_str(quote_ind(1)+1:quote_ind(2)-1);
                    tmp_cell{ct,1} = file_path;
                    time_flag = 1;
                end
                
                % It looks to me like the file path string could also be
                % used to determine file priority (if date doesn't always
                % work), for example: "/data/12717/320620161224_hy1.csv"
                % The dir between data and file name appears to be
                % incremental so the highest value should take priority
                % I will check with Sharon at some point. 4/10/19 - jp
                if time_flag == 1 % Get file time stamp
                    time_str = regexp(test_str,'\d{4}\-\d{2}\-\d{2}', ...
                        'once', 'match');
                    if ~isempty(time_str)
                        sdn = datenum(time_str,'yyyy-mm-dd');
                        tmp_cell{ct,2} = sdn;
                        time_flag = 0;
                    end
                end
                clear t1 t2 file_path time_str sdn
            end
            
            % OK IF UNMERGED FILE(S) EXIST CHOOSE THE RIGHT ONE, OTHERWISE
            % LOOK IN THE BOTTLE SECTION
            if ct == 1
                list{i,iBF} = [data_url, tmp_cell{ct,1}];
                bottle_flag = 0; % no need to look in bottle section
            elseif ct > 1 % multiple files - choose most recent
                tmp_cell = tmp_cell(1:ct,:);
                [~,IX] = sort(cell2mat(tmp_cell(:,2)),'descend');
                list{i,iBF} = [data_url, tmp_cell{max(IX),1}];
                bottle_flag = 0; % no need to look in bottle section
            end
        end
        
        %         if regexp(data_ref_target,'320620161224', 'once') % testing
        %             pause
        %         end
        
        clear sub_str ind Nind tmp_cell
        
        % ****************************************************************
        % TRY LOOKING IN BOTTLE SECTION IF NOT FOUND IN UNMERGED SECTION
        % ****************************************************************
        if bottle_flag == 1 && ~isempty(ind1) % BOTTLE DATA HEADER LINE
            cchdo_upload = 1;
            end_ind = min([ind2 ind3 ind4 ind5]); % Make string as small as possible
            sub_str = big_str(ind1:end_ind); % constrain better?
            ind1A = regexp(sub_str,'"');
            if isempty(ind1A)
                list{i,iBF} = 'NO DATA PATH YET';
                list{i,iCCF} = '0';
            else
                file_path = sub_str(ind1A(1)+1:ind1A(2)-1); % 1st one should be file
                if regexp(file_path,'hy1\.csv', 'once') % but double check
                    list{i,iBF} = [data_url, file_path];
                else
                    file_path = sub_str(ind1A(3)+1:ind1A(4)-1); % if not, try second one (ie cruise 74JC20151217 it's the second file)
                    if regexp(file_path,'hy1\.csv', 'once') % but double check
                        list{i,iBF} = [data_url, file_path];
                    else
                        list{i,iBF} = 'NO DATA PATH YET';
                        list{i,iCCF} = '0';
                    end
                end
                clear file_path sub_str ind1A
            end
        elseif isempty(list{i,iBF}) % still no
            list{i,iBF} = 'NO DATA PATH YET';
            list{i,iCCF} = '0';
        end
        list{i,iCCF} = num2str(cchdo_upload);
        
        % GET BOTTLE DATA FILE NAME FROM FILE PATH
        bottle_path = list{i,iBF};
        if ~strcmp(bottle_path,'NO DATA PATH YET') % BOTTLE
            file_name   = regexp(bottle_path,'\w+\.csv|\w+\.\w+\.csv', ...
                'once','match');
            list{i,iDF} = file_name;
        else
            %list{i,iDF} = 'NO DATA FILE';
            list{i,iDF} = '';
        end
        
        % ****************************************************************
        % TRY GETTING CTD FILE PATH NEXT - ONLY LOOK IN CTD SECTION
        if ~isempty(ind2) % CTD file section exists
            end_ind = min([ind3 ind5 ind4]); % Make string as small as possible
            %sub_str = big_str(ind1:end_ind); % constrain better?
            sub_str = big_str(ind2:end_ind); % better? jp 5/31/19
            ind2A = regexp(sub_str,'"' ,'once');
            if isempty(ind2A)
                list{i,iCTD} = 'NO DATA PATH YET';
            elseif regexp(sub_str,'ct1.zip','once') % ctd zipped file exists
                file_path = regexp(sub_str,'/\w+/\w+/\w+ct1.zip','once', ...
                    'match');
                list{i,iCTD} = [data_url,file_path];
            else
                list{i,iCTD} = 'NO DATA PATH YET';
            end
        else
            list{i,iCTD} = 'NO DATA PATH YET';
        end
        %if i == 77, pause,end % TESTING
        clear file_path sub_str ind4A end_ind
    end
    fprintf('\n')
    out.cchdo_file_urls = 1;
catch
    %catch ME
    disp('CCHDO URL search failed!.');
    %disp(ME.message)
    return
end

% ***********************************************************************
% ADD NON SOCCOM BOTTLE DATA LINKS TO THE LIST
fid = fopen([user_dir,'SHIPBOARD\',non_soccom_fn],'r','n','UTF-8');
tline = ' ';
while ischar(tline)
    if regexp(tline,'^MBARI_ID', 'once')
        ns_hdr = regexp(tline,  '\t', 'split');
        f_str = repmat('%s', 1, size(ns_hdr,2));
        break
    end
    tline = fgetl(fid);
end
d = textscan(fid, f_str,'Delimiter', '\t', 'CollectOutput',1);
d = d{1,1};
fclose(fid);

for i = 1:size(d,1)
    list(rlist+i,[iM,iINST,iWMO,iCRU,iSTA,iCST,iDF]) = ...
        d(i,[iM,iINST,iWMO,iCRU,iSTA,iCST,iDF]);
    list(rlist+i,[iCCH,iBF,iCCF,iCTD]) = {'NO URL YET' ...
        'NO DATA PATH YET','-1' 'NO DATA PATH YET'};
end
clear fid d i tline ns_hdr f_str

% ***********************************************************************
% ***********************************************************************
% !!! EXCEPTIONS SECTION !!!
% Add data file names to the list that exist but are not available for
% public use and should not be destributed outside of SOCCOM
%
% Also add float file name if the normal file selection precidence does not
% apply (i.e HAZMAT)
% ***********************************************************************
% ***********************************************************************
% keyboard
% ((TM 2/14/24)) SOCCOM 17S; final is on cchdo, but code defaults to hy1 file (??) which (in this case only?) has strange entries for bottle (mixed char-numeric).  The exc. file is the most recent and best formatted so hardcode for this!
% t1 = strcmp(list(:,iCRU),'I7S') & cellfun(@isempty,(list(:,iDF)));
t1 = strcmp(list(:,iCRU),'I7S');
if sum(t1) > 0
    list(t1,iDF) = {'49NZ20191229_usethis.exc.csv'};
    list(t1,iBF) = {'NO DATA PATH YET'};
end

% ((TM 2/15/24)) And similar for P17!!!
% t1 = strcmp(list(:,iCRU),'P17S') & cellfun(@isempty,(list(:,iDF)));
% keyboard
t1 = strcmp(list(:,iCRU),'P17S');
if sum(t1) > 0
    list(t1,iDF) = {'49NZ20170208_usethis.exc.csv'};
    list(t1,iBF) = {'NO DATA PATH YET'};
end

% Kaxis (Deployment crusie aboard AURORA AUSTRALIS - South central indian ocean
t1 = strcmp(list(:,iCRU),'Kaxis') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'09AR20160111.exc.csv'};
end

% Investigator (IN_2018v01) bottle data for the Austr. occupation of SR03
% DATA appear to be at CCHDO now - jp 04/10/19
t1 = strcmp(list(:,iCRU),'SR03') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'096U20180110.exc.csv'};
end

%ACE cruise.
t1 = strcmp(list(:,iCRU),'ACE') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'RUB320161220.exc.csv'};
end

%NCAOR cruise.
t1 = strcmp(list(:,iCRU),'NCPOR') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'91AA20171209.exc.csv'};
end

%BROKE cruise is now CCAMLR. jp 5/14/20
t1 = strcmp(list(:,iCRU),'CCAMLR') & cellfun(@isempty,(list(:,iDF)));
%t1 = strcmp(list(:,iCRU),'BROKE') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'490S20181205.exc.csv'};
end

%PS117 cruise.
t1 = strcmp(list(:,iCRU),'PS117') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'06AQ20181215.exc.csv'};
end

%ANDREXII cruise.  4/30/19.  Bottle file came direct from Bob Key on 4/17/19.  Not online yet.  May remove exception once available/updated online.
t1 = strcmp(list(:,iCRU),'ANDREXII') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'74JC20190221.exc.csv'};
end

%GOBGC OOI Irminger Sea Cruise
t1 = strcmp(list(:,iCRU),'OOI') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'33VB20210803_hy1.csv'};
end

%A16N 7/17/23 - Preliminary file available for the southern leg of cruise (south of 28N).  This includes stations for 4 of the 7 floats deployed.  
t1 = strcmp(list(:,iCRU),'A16N') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
    list(t1,iDF) = {'33RO20230306.csv'};
end

%Aaron, Solomon cruise, see email from Sharon on Aug 15, 2023; linking here https://soccompu.princeton.edu/DeploymentCruises/Pacific/Araon/2023/24HU20230110/Master_Solomon.html
t1 = strcmp(list(:,iCRU),'Solomon') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'24HU20230110_preliminary_hy1.csv'};
end

%Thwaites cruise, see email from Sharon on September 6, 2023; linking here https://soccompu.princeton.edu/DeploymentCruises/Pacific/Palmer/2020/320620200125/Master_Thwaites.html
t1 = strcmp(list(:,iCRU),'Thwaites') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'320620200125_hy1.csv'};
end

%PIPERS cruise; see Sharon's email from 12/4/2023 linking here: https://soccompu.princeton.edu/DeploymentCruises/Pacific/Palmer/2017/320620170410/Master_PIPERS.html
t1 = strcmp(list(:,iCRU),'PIPERS') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'320620170410_hy1.csv'};
end

%GoaFrem, RR2307 cruise; see Sharon's email from 6/5/2024 linking here: https://soccompu.princeton.edu/DeploymentCruises/Indian/Revelle/2023/33RR20230629/Master_GoaFrem_RR2307.html
t1 = strcmp(list(:,iCRU),'GoaFrem') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'33RR20230629_hy1.csv'};
end

%Nitrite_2023 cruise; see Sharon's email from 01/23/2024 linking here: https://soccompu.princeton.edu/DeploymentCruises/Pacific/Revelle/2023/33RR20231118/Master_Nitrite_2023.html
t1 = strcmp(list(:,iCRU),'Nitrite_2023') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'33RR20231118_prelim_hy1.csv'};
end

%A16N_Leg2  cruise; see Sharon's email from 01/25/2024 linking here:https://soccompu.princeton.edu/DeploymentCruises/Atlantic/Brown/2023/33RO20230413/Master_A16N_Leg2.html
t1 = strcmp(list(:,iCRU),'A16N') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'33RO20230413_A16N_L2_prelim_v1_hy1.csv'};
end

%Resing cruise, see email from Sharon on Dec 19, 2023; 
t1 = strcmp(list(:,iCRU),'Resing') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'33RR20210918_preliminary_hy1.csv'};
end

%Resing cruise, see email from Sharon on Feb 14, 2024; 
t1 = strcmp(list(:,iCRU),'A12(A13.5)') & cellfun(@isempty,(list(:,iDF)));
if sum(t1) > 0
   list(t1,iDF) = {'A13-A12_stations_checked_by_Leti_hy1.csv'};
end

%P02E; error in process_data_ref after 4/11/24 when TRD file added to 'unmerged data' section of https://cchdo.ucsd.edu/cruise/33RR20220613 ...?
t1 = strcmp(list(:,iCRU),'P02E');
if sum(t1) > 0
   list(t1,iDF) = {'33RR20220613_hy1.csv'};
end

%%SOCCOM A13.5 cruise, 2020; listed as 33RO20200321 but odd file name?
%t1 = strcmp(list(:,iCRU),'A13.5') & cellfun(@isempty,(list(:,iDF)));
%if sum(t1) > 0
%    list(t1,iDF) = {'A13-A12_2020_RB2002_2020_stationsfinal2021.csv'};
%end


%HAZMAT cruise.
% Presently there is a hy1 & exc files in the unmerged section
% hy1 usually takes priority but Bob's file (exc is more recent in this
% case). NOT NEEDED CORRECT TIME STAMP ACCOUNTING TAKES CARE OF IT NOW !!
% t1 = strcmp(list(:,iCRU),'HAZMAT')
% if sum(t1) > 0
%     list(t1,iDF) = {'320620161224.exc.csv'};
% end

clear nM nUW nWMO nCRU nSTA nCST nDF nCCH nBF nCCF nCTD

out.list_hdr = hdr;
out.list     = list;

% ***********************************************************************
% PRINT FULL & SHORT LISTS TO TEXT FILE
fid = fopen([save_dir,'BottleData_lookup_table_full.txt'], 'w','n','UTF-8');
if fid == -1
    disp('File may be open in anothe program:');
    disp([save_dir,'BottleData_lookup_table_full.txt']);
    disp('Close other instance and try program again');
    return
end
f_str = repmat('%s\t', 1, clist);
f_str = regexprep(f_str,'\\t$','\\r\\n');
fprintf(fid, f_str, hdr{:});
print_list = list'; % transpose so no loop needed for printing cell array
fprintf(fid, f_str, print_list{:});
fclose(fid);

short_list = list(:,1:iDF);
print_short_list = short_list';
fid = fopen([save_dir,'BottleData_lookup_table.txt'], 'w','n','UTF-8');
f_str = repmat('%s\t', 1, iDF);
f_str = regexprep(f_str,'\\t$','\\r\\n');
fprintf(fid, f_str, hdr{1:iDF});
fprintf(fid, f_str, print_short_list{:});
fclose(fid);

% ***********************************************************************
% OK FILE PATHS FOUND NOW COPY FILES LOCALLY
if ~isdir(save_dir)
    mkdir(save_dir);
end

if ~isdir([save_dir,'ctd\'])
    mkdir([save_dir,'ctd\']);
end

% COPY BOTTLE FILES TO LOCAL
% Make a Unique list by file name
fprintf('\nCopying files locally\n')
[~,ia,~] = unique(list(:,iDF));
blist    = list(ia,:);
loop_ct  = 0;
for i = 1:size(blist,1)
    loop_ct = loop_ct +1;
    if loop_ct > 60
        loop_ct = 1;
        fprintf('.\n')
    else
        fprintf('.')
    end
    
    bottle_path = blist{i,iBF};
    if ~strcmp(bottle_path,'NO DATA PATH YET') % BOTTLE
        try
            outfilename = websave([save_dir, blist{i,iDF}],bottle_path);
        catch
            disp(['Websave failed on: ',bottle_path])
            disp('Trying next file.')
        end
    end
end
fprintf('\n');
clear blist ia


% COPY CTD FILES LOCALLY
[~,ia,~] = unique(list(:,iCTD));
clist    = list(ia,:);
loop_ct  = 0;
for i = 1:size(clist,1)
    loop_ct = loop_ct +1;
    if loop_ct > 60
        loop_ct = 1;
        fprintf('.\n')
    else
        fprintf('.')
    end
    
    ctd_path = clist{i,iCTD};
    if ~strcmp(ctd_path,'NO DATA PATH YET') % CTD
        ctd_file_name   = regexp(ctd_path,'\w+\.zip', 'once','match');
        try
            outfilename = websave([save_dir,'ctd\', ctd_file_name], ctd_path);
        catch
            disp(['Websave failed on: ',ctd_path])
            disp('Trying next file.')
        end
    end
end
fprintf('\n');
clear clist ia






