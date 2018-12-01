%make_Sprof2ODV_raw.m
% A wrapper script to convert BGC ARGO Sprof files to ODV COMPATIBLE TEXT
% FILES

% 1) argo_floats_bysensor.m
%       generates a verified list of WMO #s and dac paths for floats
%       that carry a given sensor determined by the choosen sensor_list
%       indice (VERIFICATION TAKES a WHILE - Don't run every time
%
% 2) Sprof2mat.m
%       takes the Sprof file from NetCDf to a matrix. The data is binned
%       according to the bins matrix. See header comments for more info
%
% 3) sprofmat2ODV.m 
%       Creates an ODV compatible text file from the Sprof matrix created
%       by Sprof2mat.
site_flag = 2;
dirs = [];
fp = filesep;
Sprof_dir = '\\atlas\Chem\ARGO\DATA\Sprof\';

% NORMALLY 1, row of float list to start processing at
% SET to different value if process crashes usually due to a ftp failure
start_row = 1; 

% *************************************************************************
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir    = getenv('USERPROFILE');
    dirs.data   = [user_dir,fp,'Documents',fp,'MATLAB',fp,'ARGO',fp];
    dirs.temp   = [getenv('HOMEDRIVE'),fp,'temp',fp]; % for my computer homedrive = C:
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end
diary([dirs.data,'non_mbari_float_processing_log.txt']);

% *************************************************************************
% GENERATE FLOAT LIST FOR GIVEN SENSOR
% "DON'T RUN EVERY TIME - TAKES A WHILE
% *************************************************************************
%floatbysensor('O2';

% BUILD FORMAT STRING AND PARSE DATA
fid = fopen([dirs.data,'FLBBorMCOMSorFLNTU.txt']);
%fid = fopen([dirs.data,'ISUSorSUNA.txt']);
%fid = fopen([dirs.data,'AANDERAAorSBE63.txt']);
%fid = fopen([dirs.data,'aircal_floats.txt']);
%fid = fopen([dirs.data,'DURAorSEAFET.txt']);
tline = ' ';
while ischar(tline)
    if regexp(tline,'^ID', 'once') % stop at header line
        break
    end
    tline = fgetl(fid);
end
if ~ischar(tline)
    disp('No header line found')
    return
end

list_hdr = regexp(tline,'\t','split'); % CELL ARRAY OF HEADER VARIABLES 
hdr_rows = size(list_hdr,2);

d_format = '';
for i = 1: hdr_rows
        d_format = [d_format,'%s'];
end
list = textscan(fid,d_format,'Delimiter','\t','CollectOutput',1);
list = list{1,1}; % get rrid of outer shell
fclose(fid);

%pause

% REMOVE ANY MBARI FLOATS FROM THE LIST
t1 = strcmp(list(:,6),'"Argo UW-MBARI"');
t2 = strcmp(list(:,6),'"Argo UW-SOCCOM"');
%t3 = strcmp(list(:,4),'OPERATIONAL'); % STATUS

list(t1|t2,:) =[];

%pause

%t3 = strcmp('1901329',list(:,2));
%t3 = strcmp('5905167',list(:,2));
%t3 = strcmp('6901529',list(:,2));
%t3 = strcmp('6901494',list(:,2));
%t3 = strcmp('6901495',list(:,2));
%t3 = strcmp('3901496',list(:,2));
%t3 = strcmp('5904923',list(:,2));
%t3 = strcmp('6901473',list(:,2)); 

% t3 = strcmp('5903630',list(:,2));
% t3 = strcmp('1901329',list(:,2));
% t3 = strcmp('6901524',list(:,2));
% list = list(t3,:);

% REFINE LIST - COMMENT OUT IF NOT USING
% t3 = ~cellfun(@isempty,regexp(list(:,2),'6901647|6901180|1901339|6901654'));
% list = list(t3,:);


% *************************************************************************
% MAKE A LIST BASED ON BAD SENSOR LIST -- COMMENT OUT WHEN NOT USED
% file_path  = ['C:\Users\jplant\Documents\MATLAB\ARGO\', ...
%               'ODVbad_sensor_list.txt'];
% bad  = parse_bad_ODVsensor_list(file_path);
% bad_hdr  = bad.hdr;
% bad_list = bad.list;
% t3 = ones(size(list(:,1)))* 0;
% for i = 1:size(t3,1)
%     tf = strcmp(list{i,2},bad_list(:,1));
%     if sum(tf > 0)
%         t3(i) = 1;
%     end
% end
% t3 = logical(t3);
% list = list(t3,:);

%pause
% *************************************************************************
[rrr,ccc] = size(list);
clear t1 t2
iWMO = find(strcmp('REF',list_hdr) ==1);
iDAC = find(strcmp('DAC PATH',list_hdr) ==1);

% GET SOME LIST INDICES
% *************************************************************************
% NOW STEP THROUGH LIST AND CONVERT SPROF FILES TO ODV FILE

for fn_ct = start_row: rrr
    info.WMO        = list{fn_ct,iWMO};
    info.fn         = [info.WMO, '_Sprof.nc'];
    info.dac_path   = list{fn_ct,iDAC};
    info.dac        = regexp(info.dac_path,'(?<=dac/)\w+','once','match');
    info.local_path = [Sprof_dir,info.dac,'\'];
    
    % TEST FOR FILE
    if isempty(ls([info.local_path, info.fn]))
        disp(['FILE NOT FOUND: ',info.local_path, info.fn])
        continue
    end


    
    sprof_out = Sprof2mat(info);
    tf_odv = sprofmat2ODV(sprof_out, []);
end
diary off





