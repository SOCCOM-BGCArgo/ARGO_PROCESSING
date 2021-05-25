% Float_type_list.m
% Build a float list based on type for Annie Wong
%"f dirs"  only No NAVIS YET

data_dir = '\\atlas\chemwebdata\floats\';
save_dir = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\';

type_01 = {'p';'t';'s';'Topt';'TPhase';'RPhase';'no3';'pH(V)'; ...
           'FSig';'BbSig';'TSig'};

% Update float list
[float_names UW_ID WMO_ID] = get_MBARI_WMO_list;

% FOR SOCCOM SUBSETTING GET SOCCOM FLOAT HTML FILE
urlwrite('http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html', ...
    'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\sensor_stats.txt');
fid = fopen('C:\Users\jplant\Documents\MATLAB\ARGO\DATA\sensor_stats.txt');

tline = ' ';
ct = 0;
SOCCOM_flts = cell(1000,3);  % predim
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

SOCCOM_list = cell(1000,4);
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
                SOCCOM_list(ct,1:3) = [float_names(j) UW_ID(j) WMO_ID(j)];
            end
        end
    end
end
SOCCOM_list = SOCCOM_list(1:ct,:);
        


%float_type = ones(size(UW_ID,1),1)*0; % predim
float_type = ones(size(SOCCOM_list(:,1),1),1)*0; % predim

%for i = 1 : size(UW_ID,1)
for i = 1 : size(SOCCOM_list,1)
    %flt_name = float_names{i};
    flt_name = SOCCOM_list{i,1};
    SOCCOM_list{i,4} = 0;
    
%     if exist([data_dir,'f',UW_ID{i},'\'],'dir');
%         msg_file = ls([data_dir,'f',UW_ID{i},'\*003.msg']); % 3 is random
%     else
%         disp([flt_name,' is probably a NAVIS float skipping for now'])
%         continue
%     end
    
    if exist([data_dir,'f',SOCCOM_list{i,2},'\'],'dir');
        msg_file = ls([data_dir,'f',SOCCOM_list{i,2},'\*003.msg']); % 3 is random
    else
        disp([flt_name,' is probably a NAVIS float skipping for now'])
        continue
    end
    
    % Choose 3rd msg file on list (random choice) and look at header
    file_path = [data_dir,'f',SOCCOM_list{i,2},'\',strtrim(msg_file)];
    if exist(file_path, 'file') == 2
        d = parse_APEXmsg4ARGO(file_path);
        
        % Check Header size first - they must be equal
        if size(type_01,1) ~= size(d.lr_hdr,1)
            disp(['No TYPE 1 match for ', flt_name, ' Headers ',...
                'different sizes'])
            continue
        end
        
        if d.FlbbMode == 0; % NO FLBB SENSOR but data columns exist
                disp(['No TYPE 1 match for ', flt_name, ' FlbbMode ',...
                'flag = 0'])
            continue
        end
        
        for j = 1 : size(d.lr_hdr,1)
            if  isempty(regexp(type_01{j,1}, d.lr_hdr{j,1}, 'once'))
                disp(['No TYPE 1 match for ', flt_name, ' No ',...
                    type_01{j,1}, ' variable in header'])
                break
            end
        end
        
        float_type(i) = 1;
        SOCCOM_list{i,4} = 1;
        disp([flt_name, '  TYPE = ', num2str(float_type(i))]);
    end
    clear flt_name msg_files file_path d type_test
end

% MAKE A TEXT FILE
%r = size(float_names,1);
r = size(SOCCOM_list,1);
fid = fopen([save_dir,'float_type_list.txt'], 'w');
fprintf(fid,'%s\t%s\t%s\t%s\r\n','MBARI float name', 'UW ID#', ...
    'WMO ID #','ANNIE TYPE');
for i = 1:r
%     fprintf(fid,'%s\t%s\t%s\t%s\r\n',float_names{i}, UW_ID{i}, WMO_ID{i}, ...
%         num2str(float_type(i)));
    fprintf(fid,'%s\t%s\t%s\t%s\r\n',SOCCOM_list{i,1}, SOCCOM_list{i,2}, ...
        SOCCOM_list{i,3}, num2str(SOCCOM_list{i,4}));
end
fclose(fid);

