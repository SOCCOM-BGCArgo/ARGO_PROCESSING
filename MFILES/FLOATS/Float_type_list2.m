% Float_type_list2.m
% Build a float list based on type for Annie Wong
%"f dirs"  only No NAVIS YET

% *************************************************************************
% SET DIRS AND PATHS
% *************************************************************************
%dirs.msg = '\\atlas\chemwebdata\floats\';
dirs = [];
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
save_dir = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    
    dirs.mat       = [user_dir,'FLOATS\'];
    dirs.cal       = [user_dir,'CAL\'];
    dirs.NO3config = [user_dir,'CAL\'];
    dirs.FVlocal   = [user_dir,'FLOATVIZ\'];
    dirs.FV        = [user_dir,'FLOATVIZ\'];
    %dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
    
    dirs.QCadj     = [user_dir,'CAL\QC_LISTS\'];
    dirs.temp      = 'C:\temp\';
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    dirs.config    = '\\atlas\Chem\ISUS\Argo\';
    
    dirs.log = [user_dir,'Processing_logs\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

% *************************************************************************
% DEFINE FLOAT TYPES FOR ANNIE
% *************************************************************************
% hdr_list = {'p';'t';'s';'Topt';'BPhase';'TPhase';'RPhase';'no3'; ...
%             'pH(V)'; 'FSig';'BbSig';'TSig'; 'O2ph'; 'O2tV';'phV';...
%             'phT';'Fl'; 'Bb'; 'Cd'; 'Cdm'}; % suna
           
% APEX FLOATS 

ftype{1,1} = {'p';'t';'s'; ...  % TYPE 1 or 2
              'Topt';'TPhase';'RPhase';'no3';'pH(V)';'FSig';'BbSig';'TSig'}; 
ftype{1,2} = {'p';'t';'s'; ...  % TYPE 3
              'Topt';'TPhase';'RPhase';'no3';'FSig';'BbSig';'TSig'};       
ftype{1,3} = {'p';'t';'s'; ...  % TYPE 4
              'Topt';'TPhase';'RPhase';'pH(V)'}; 
ftype{1,4} = {'p';'t';'s'; ...  % TYPE 5
              'BPhase';'Topt';'no3';'FSig';'BbSig';'TSig'}; 
ftype{1,5} = {'p';'t';'s'; ...  % TYPE 6
              'TPhase';'Topt';'no3'};
ftype{1,6} = {'p';'t';'s'; ...  % TYPE 7
              'BPhase';'Topt';'no3'};
ftype{1,7} = {'p';'t';'s'; ...  % TYPE 8
              'TPhase';'Topt';'no3';'FSig';'BbSig';'TSig'};  
ftype{1,8} = {'p';'t';'s'; ...  % TYPE 9
              'Topt';'TPhase';'RPhase';'no3'};   
ftype{1,9} = {'p';'t';'s'; ...  % TYPE 10
              'TPhase';'Topt';'pH(V)'};           


% NAVIS FLOATS
ftype{1,10} = {'p';'t';'s'; ... % TYPE 20 (ALL)
              'no3';'O2ph';'O2tV';'Fl';'Bb';'Cdm';'phV';'phT'};        
ftype{1,11} = {'p';'t';'s'; ... % TYPE 21 (NO MCOMS)
              'no3';'O2ph';'O2tV';'phV';'phT'};     
ftype{1,12} = {'p';'t';'s'; ... % TYPE 22 (NO pH) & Cd not Cdm in hdr
              'no3';'O2ph';'O2tV';'Fl';'Bb';'Cd'};
ftype{1,13} = {'p';'t';'s'; ... % TYPE 23 (NO NITRATE)
              'O2ph';'O2tV';'Fl';'Bb';'Cdm';'phV';'phT'}; 
  
% *************************************************************************
% UPDATE FLOAT LIST
% *************************************************************************
[float_names, UW_ID, WMO_ID] = get_MBARI_WMO_list(dirs);

% *************************************************************************
% GET SOCCOM FLOAT LIST SUBSET USING HTML FILE
% *************************************************************************
% GET SOCCOM FLOAT HTML FILE
urlwrite('http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html', ...
 [save_dir, 'sensor_stats.txt']);
     
fid = fopen([save_dir,'sensor_stats.txt']);

% *************************************************************************
% FIND UW ID's FOR SOCCOM FLOATS
tline = ' ';
ct    = 0;
SOCCOM_flts = cell(1000,3);  % over predim and clean up later
while ischar(tline) 
    UW_ID_str = regexp(tline,'(?<=>)\d+(?=<)', 'once', 'match');
    if  ~isempty(UW_ID_str)
        ct = ct+1;
        SOCCOM_flts{ct,1} = UW_ID_str;
    end
    tline = fgetl(fid);
end
SOCCOM_flts = SOCCOM_flts(1:ct,1); % cell array of SOCCOM UW_ID's
clear ct tline UW_ID_str 

% *************************************************************************
% BUILD SOCCOM NAME AND ID LIST
SOCCOM_list = cell(1000,5);
ct = 0;% predim
for i = 1:size(SOCCOM_flts,1)
    t1 = strcmp(SOCCOM_flts{i},UW_ID);
    if sum(t1) == 1
        ct = ct+1;
        SOCCOM_list(ct,1:3) = [float_names(t1) UW_ID(t1) WMO_ID(t1)];
    elseif sum(t1) > 1 % two float with same UW_ID
        ind = find(t1 == 1);
        for j = length(ind)
            if regexpi(float_names{ind(j)},'SoOCN','once')
                ct = ct+1;
                SOCCOM_list(ct,1:3) = [float_names(ind(j)) UW_ID(ind(j)) WMO_ID(ind(j))];
            end
        end
    end
end
SOCCOM_list = SOCCOM_list(1:ct,:);

% *************************************************************************
% BUILD NON SOCCOM NAME AND ID LIST
NON_SOCCOM_list = cell(1000,5);
ct = 0;% predim
for i = 1:size(UW_ID,1)
    t1 = strcmp(UW_ID{i}, SOCCOM_list(:,2));
    if sum(t1) == 0 % NOT A SOCCOM FLOAT
        ct = ct+1;
        NON_SOCCOM_list(ct,1:3) = [float_names(i) UW_ID(i) WMO_ID(i)];
    end
end
NON_SOCCOM_list = NON_SOCCOM_list(1:ct,:);
  
% *************************************************************************
% *************************************************************************
% SOCCOM FLOATS
% DETERMINE FLOAT TYPE AND ADD TO LIST
% *************************************************************************
% *************************************************************************
float_type = ones(size(SOCCOM_list(:,1),1),1)*0; % predim
%disp('CHECKING SOCCOM FLOATS FOR FLOAT TYPE ..................')
for i = 1 : size(SOCCOM_list,1)
    %flt_name = float_names{i};
    flt_name = SOCCOM_list{i,1};
    SOCCOM_list{i,4} = 0;
    
    % CHECK FOR FLOAT DIRECTORY
    if exist([dirs.msg,'f',SOCCOM_list{i,2},'\'],'dir');
        msg_dir = [dirs.msg,'f',SOCCOM_list{i,2},'\'];
        %disp('APEX float directory detected')
        FT = 'APEX';
        SOCCOM_list{i,5} = FT;
    elseif exist([dirs.msg,'n',SOCCOM_list{i,2},'\'],'dir');
        msg_dir = [dirs.msg,'n',SOCCOM_list{i,2},'\'];
        FT = 'NAVIS';
        SOCCOM_list{i,5} = FT;
    else
        disp(['could not find float directory for : ',SOCCOM_list{i,1}]);
        continue
    end
    
    % NOW CHECK FOR MESSAGE FILE WITH a HEADER ROW
    hdr_chk = 0; % 0 = empty
    
    for j = 1:3 % check 1st 3 msg files to look for a header
        target = [msg_dir,'*',sprintf('%03.0f',j),'.msg'];
        msg_file = ls(target);
        if ~isempty(msg_file)
            if strcmp(FT,'APEX')
                d = parse_APEXmsg4ARGO([msg_dir,strtrim(msg_file)]);
                hdr_chk = ~isempty(d.lr_hdr) ;
            else % has to be a NAVIS float
                d = parse_NAVISmsg4ARGO([msg_dir,strtrim(msg_file)]);
                hdr_chk = ~isempty(d.lr_hdr) ;
            end
            if hdr_chk == 1 % header found
                break
            end
        end
        %disp(['No header or file found for ',target]);
    end
    
    if hdr_chk == 0 % tried 3 times & still no header so move on
        disp(['NO HEADER FOUND TO DETERMINE TYPE FOR ' flt_name])
        float_type(i) = -1;
    end
    
    if hdr_chk == 1
        hdr = d.lr_hdr;
        if size(hdr,1) == 1 % for navis
            pause
            hdr = hdr';
        end
        for j = 1: size(ftype,2)
            nrows = size(ftype{j},1);
            if size(hdr,1) == nrows %  size of equal -- maybe the same
                tf = sum(strcmp(hdr, ftype{j})); % compare hdr vars
                if tf == nrows % A match almost
                    if     j == 1 && d.FlbbMode == 1 % TYPE 1 APEX
                        float_type(i) = 1;
                    elseif j == 1 && d.FlbbMode == 0 % TYPE 2
                        float_type(i) = 2;
                    elseif j == 2                    % TYPE 3
                        float_type(i) = 3;
                    elseif j == 3                    % TYPE 4
                        float_type(i) = 4;
                    elseif j == 4                    % TYPE 5
                        float_type(i) = 5;
                    elseif j == 5                    % TYPE 6
                        float_type(i) = 6;
                    elseif j == 6                    % TYPE 7
                        float_type(i) = 7;
                    elseif j == 7                    % TYPE 8
                        float_type(i) = 8;                       
                    elseif j == 8                    % TYPE 9
                        float_type(i) = 9; 
                    elseif j == 9                    % TYPE 10
                        float_type(i) = 10; 
                        
                    elseif j == 10                    % TYPE 20 NAVIS
                        float_type(i) = 20;
                    elseif j == 11                    % TYPE 21
                        float_type(i) = 21;
                    elseif j == 12                    % TYPE 22
                        float_type(i) = 22;
                    elseif j == 13                    % TYPE 23
                        float_type(i) = 23;
                    else                             % UNDEFINED TYPE
                        float_type(i) = 0;
                    end
                end
            end
        end
        SOCCOM_list{i,4} = float_type(i);
        %disp([flt_name, '  TYPE = ', num2str(float_type(i))]);
    end
    
    clear flt_name msg_files file_path d type_test tf nrows j
    
end

% *************************************************************************
% *************************************************************************
% NON SOCCOM FLOATS
% DETERMINE FLOAT TYPE AND ADD TO LIST
% *************************************************************************
% *************************************************************************
float_type = ones(size(NON_SOCCOM_list(:,1),1),1)*0; % predim
%disp('CHECKING FLOAT TYPE FOR NON SOCCOM FLOATS ..................')
for i = 1 : size(NON_SOCCOM_list,1)
    %flt_name = float_names{i};
    flt_name = NON_SOCCOM_list{i,1};
    NON_SOCCOM_list{i,4} = 0;
    
    % CHECK FOR FLOAT DIRECTORY
    if exist([dirs.msg,'f',NON_SOCCOM_list{i,2},'\'],'dir');
        msg_dir = [dirs.msg,'f',NON_SOCCOM_list{i,2},'\'];
        %disp('APEX float directory detected')
        FT = 'APEX';
        NON_SOCCOM_list{i,5} = FT;
    elseif exist([dirs.msg,'n',NON_SOCCOM_list{i,2},'\'],'dir');
        msg_dir = [dirs.msg,'n',NON_SOCCOM_list{i,2},'\'];
        FT = 'NAVIS';
        NON_SOCCOM_list{i,5} = FT;
    else
        %disp(['could not find float directory for : ',NON_SOCCOM_list{i,1}]);
        continue
    end
    
    % NOW CHECK FOR MESSAGE FILE WITH a HEADER ROW
    hdr_chk = 0; % 0 = empty
    
    for j = 1:3 % check 1st 3 msg files to look for a header
        target = [msg_dir,'*',sprintf('%03.0f',j),'.msg'];
        msg_file = ls(target);
        if ~isempty(msg_file)
            if strcmp(FT,'APEX')
                d = parse_APEXmsg4ARGO([msg_dir,strtrim(msg_file)]);
                hdr_chk = ~isempty(d.lr_hdr) ;
            else % has to be a NAVIS float
                d = parse_NAVISmsg4ARGO([msg_dir,strtrim(msg_file)]);
                hdr_chk = ~isempty(d.lr_hdr) ;
            end
            if hdr_chk == 1 % header found
                break
            end
        end
        disp(['No header or file found for ',target]);
    end
    
    if hdr_chk == 0 % tried 3 times & still no header so move on
        disp(['NO HEADER FOUND TO DETERMINE TYPE FOR ' flt_name])
        float_type(i) = -1;
    end
    
    if hdr_chk == 1
        hdr = d.lr_hdr;
        for j = 1: size(ftype,2)
            nrows = size(ftype{j},1);
            if size(hdr,1) == nrows %  size of equal -- maybe the same
                tf = sum(strcmp(hdr, ftype{j})); % compare hdr vars
                if tf == nrows % A match almost
                    if     j == 1 && d.FlbbMode == 1 % TYPE 1 APEX
                        float_type(i) = 1;
                    elseif j == 1 && d.FlbbMode == 0 % TYPE 2
                        float_type(i) = 2;
                    elseif j == 2                    % TYPE 3
                        float_type(i) = 3;
                    elseif j == 3                    % TYPE 4
                        float_type(i) = 4;
                    elseif j == 4                    % TYPE 5
                        float_type(i) = 5;
                    elseif j == 5                    % TYPE 6
                        float_type(i) = 6;
                    elseif j == 6                    % TYPE 7
                        float_type(i) = 7;
                    elseif j == 7                    % TYPE 8
                        float_type(i) = 8;                       
                    elseif j == 8                    % TYPE 9
                        float_type(i) = 9; 
                    elseif j == 9                    % TYPE 10
                        float_type(i) = 10; 
                        
                    elseif j == 10                    % TYPE 20 NAVIS
                        float_type(i) = 1;
                    elseif j == 11                    % TYPE 21
                        float_type(i) = 2;
                    elseif j == 12                    % TYPE 22
                        float_type(i) = 3;
                    elseif j == 13                    % TYPE 23
                        float_type(i) = 4;
                    else                             % UNDEFINED TYPE
                        float_type(i) = 0;
                    end
                end
            end
        end
        NON_SOCCOM_list{i,4} = float_type(i);
        %disp([flt_name, '  TYPE = ', num2str(float_type(i))]);
    end
    
    clear flt_name msg_files file_path d type_test tf nrows j
    
end


% *************************************************************************
% SAVE SOME FILES
% *************************************************************************
save([save_dir,'float_type_lists'],'SOCCOM_list','NON_SOCCOM_list')
% MAKE A TEXT FILES
%r = size(float_names,1);
r = size(SOCCOM_list,1);
fid = fopen([save_dir,'SOCCOM_float_type_list.txt'], 'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\r\n','MBARI float name', 'UW ID#', ...
    'WMO ID #','FLOAT TYPE','ANNIE TYPE');
for i = 1:r
    %     fprintf(fid,'%s\t%s\t%s\t%s\r\n',float_names{i}, UW_ID{i}, WMO_ID{i}, ...
    %         num2str(float_type(i)));
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\r\n',SOCCOM_list{i,1}, SOCCOM_list{i,2}, ...
        SOCCOM_list{i,3},SOCCOM_list{i,5}, num2str(SOCCOM_list{i,4}));
end
fclose(fid);

% MAKE A TEXT FILES
%r = size(float_names,1);
r = size(NON_SOCCOM_list,1);
fid = fopen([save_dir,'NON_SOCCOM_float_type_list.txt'], 'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\r\n','MBARI float name', 'UW ID#', ...
    'WMO ID #','FLOAT TYPE','ANNIE TYPE');
for i = 1:r
    %     fprintf(fid,'%s\t%s\t%s\t%s\r\n',float_names{i}, UW_ID{i}, WMO_ID{i}, ...
    %         num2str(float_type(i)));
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\r\n',NON_SOCCOM_list{i,1}, ...
        NON_SOCCOM_list{i,2}, NON_SOCCOM_list{i,3}, ...
        NON_SOCCOM_list{i,5}, num2str(NON_SOCCOM_list{i,4}));
end
fclose(fid);
