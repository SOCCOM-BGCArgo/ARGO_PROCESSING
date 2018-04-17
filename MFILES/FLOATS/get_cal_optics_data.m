function floatviz_data = get_cal_optics_data(floatviz_file)

%floatviz_file = '9095SoOcn';

% *************************************************************************
% SET PATHS & COPY FILE TO LOCAL & OPEN
% *************************************************************************
data_source = 'internet';
%data_source = 'network';

network_dir   = ['C:\Users\jplant\Documents\MATLAB\ARGO\DATA\', ...
                 'Calibrated_optical_data\'];
floatviz_url  = 'ftp://misclab.umeoce.maine.edu/floats/SOCCOM/';
temp_dir      = 'C:\temp\';
floatviz_data =[];

switch data_source
    case 'network' % GET TEXT FILE FROM SIROCCO AND STORE LOCALY
        from_str = [network_dir, floatviz_file, '.txt'];
        to_str   = [temp_dir, floatviz_file, '.txt'];
        
        if exist(from_str,'file') == 2
            copyfile(from_str, to_str)
        else
            disp(['Could not find: ',from_str]);
            return
        end
        
    case 'internet'
        if regexp(floatviz_file,'QC','once')
            floatviz_url  = [floatviz_url,'QC/'];
        end
        from_str = [floatviz_url,floatviz_file,'.txt']; % build target string
        to_str   = [temp_dir,floatviz_file,'.txt'];    % build destination string
        
        
        [~,url_chk] = urlread(from_str); % See if file exisit on the web
        if url_chk == 1
            f = urlwrite(from_str,to_str);
            disp(' ');
            disp(['Data for float ',floatviz_file,' retrieved from:  ',...
                floatviz_url]);
            disp(['Saved as  ',f]);
        else
            disp('No file found!')
            return
        end
end

fid = fopen(to_str);

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
    return
end

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
d_tmp = strcat(d{1,2}(:,1),regexprep(d{1,2}(:,2), '(\d+:\d+)',' $1'));
sdn   = datenum(d_tmp,'mm/dd/yyyy HH:MM');

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





