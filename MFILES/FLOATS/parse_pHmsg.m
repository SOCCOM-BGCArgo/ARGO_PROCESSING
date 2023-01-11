function d = parse_pHmsg(pH_file)
% PURPOSE:
%   This function parses an APEX pH file (*.dura)
%   file.
%
% USAGE:
%	d   = parse_pHmsg(pH_file)
%
% INPUTS:
%	pH_file = Path\calibration file name as a string
%
% OUTPUTS:
%   d   = a structure
%           d.hdr  = cell array of column headers
%        	d.data = a matrix of pH data
%           d.EOT  = complete data file flag
%
% EXAMPLE:
%   d  = parse_pHmsg('C:\temp\9095.009.dura')
%
% Created 01/22/2016 by jp
%
% REVISIONS:
% 08/11/2020 - JP - added code to determine MSC firmware versions and save
%                   it as a field in the ouput structure. Also delt with
%                   bogus time stamp issue (18110). If time stamp bogus
%                   subsitute EOP date
% 11/28/22 - JP - add code to deal with partial file case where one partial
%                 data line exists but no data & no comma sepration in line
%                 a brute force fix but should be Ok

% TESTING:
%pH_file ='C:\temp\8514.034.dura'; % TEST
%pH_file ='C:\temp\19644.081.dura'; % TEST

% ************************************************************************
% SET FORMATS AND VARIABLES
% ************************************************************************
d.hdr        = ' '; % Returns if no data
d.data       = [];
d.EOT        = 0;
d.MSC_FW_ver = 0;
d.EOP_sdn    = NaN;
d.bogus_time = 0;

% check for valid target
if exist(pH_file,'file')
    fid   = fopen(pH_file);
else
    disp(['Can not find *.dura file at: ',pH_file])
    return
end

%All *.dura ph headers
% 'CRC' 'ID' 'Date/Time' 'Epoch Sec' 'CTD Pres' 'CTD Temp' 'CTD Sal', ...
% 'Record Ctr' 'Power Cycle Ctr' 'Error Ctr', 'Internal Temp', ...
% 'Internal Humidity' 'System Power, V' 'System Power, I', ...
% 'Computed pH foobar' 'Backup Batt, V' , ...
% 'Vrs' 'SD Vrs' 'Vk' 'SD Vk' 'Ik' 'Ib'

pH_format = '%*s%*s%s%*f%f%f%f%*f%*f%*f%f%f%f%f%*f%f%f%f%f%f%f%f';

pH_hdr = {'Date/Time' 'CTD Pres' 'CTD Temp' 'CTD Sal' 'Internal Temp', ...
    'Internal Humidity' 'System Power, V' 'System Power, I', ...
    'Backup Batt, V' 'Vrs' 'std Vrs' 'Vk' 'std Vk' 'Ik' 'Ib'};

% ************************************************************************
% FIND # OF SAMPLE SPECTRA DATA ROWS
% COUNT COMMAS FOR DATA COMPETENESS
% ************************************************************************
tline    = ' ';
data_ct  = 0;
comma_ct = [];
MSC_FW_ver = 0;

% MSC firmware version App build dates as of  08/12/2020
msc_fw_versions(1,:) = {'Jan 27 2012', 0};
msc_fw_versions(2,:) = {'Apr 27 2012', 0};
msc_fw_versions(3,:) = {'Apr  7 2015', 1};
msc_fw_versions(4,:) = {'Sep 30 2015', 2};
msc_fw_versions(5,:) = {'Nov  1 2018', 3};
msc_fw_versions(6,:) = {'Apr 30 2019', 4};

while ischar(tline)
    if regexp(tline,'App Build','once') % msc build version here
        dstr = regexp(tline,'\w{3}\s+\d+\s+\d{4}','match', 'once');
        t1   = strcmp(msc_fw_versions(:,1), dstr);
        if sum(t1) == 1
            MSC_FW_ver = msc_fw_versions{t1,2};
        end
    end
    
    % get end profile time incase dura file has bad time stamp, i.e. 18110
    if regexp(tline,'^<EOP>', 'once')
        dstr = regexp(tline,'\w{3}\s+\d+\s+.+\d{4}','match', 'once');
        if ~isempty(dstr)
            d.EOP_sdn = datenum(dstr,'mmm dd HH:MM:SS yyyy');
        end
    end
    
    if  regexp(tline,'^0x','once') % DATA LINES
        data_ct = data_ct + 1;
        comma_ct = [comma_ct;length(regexp(tline,','))] ;
    end
    tline = fgetl(fid); % get 1st line
end

if isnan(d.EOP_sdn)
    disp(['Could not find EOP line in ', pH_file,'. file is ', ...
        'probably incomplete'])
end

% check for data lines
if data_ct < 1
    disp(['File exists but no data lines found for : ',pH_file]);
    fclose(fid);
    return
end

% Check total comma count. Error thrown by partial file 19644.081.dura
% where only one partial data line & comma count  = 0
if size(comma_ct,1) == 1 && comma_ct(1,1) == 0
    fprintf(['File exists but only 1 partial data line with no data ', ...
        'inside for %s\n'], pH_file)
    fclose(fid);
    return
end

% ************************************************************************
% GO BACK TO TOP OF FILE & THEN START PARSING THE DATA
% ************************************************************************
frewind(fid) % GO BACK TO TOP

d_cols     = size(pH_hdr,2);            % # data columns
comma_ct   = mode(comma_ct);            % most frequent value
data       = ones(data_ct, d_cols)*NaN; % predim
tline      = ' ';                       % predim
ct         = 0;   % line counter
bogus_time = 0;
EOT        = 0;

while ischar(tline)                              % step through lines again
    if  regexp(tline,'^0x','once')            	 % data line?
        if length(regexp(tline,',')) == comma_ct % full line?
            ct = ct+1;
            dline = textscan(tline,pH_format,'Delimiter',',',...
                'Collectoutput',1);
            % GET TIME STAMP
            if regexp(dline{1,1}{1}, '01/01/2000','once')
                bogus_time = 1;
                data(ct,1) = d.EOP_sdn;
            else
                data(ct,1) = datenum(dline{1,1}{1},'mm/dd/yyyy HH:MM:SS');
            end
            if ~isempty(dline{1,2})
                data(ct,2:d_cols) = dline{1,2};
            else
                disp('SKIPPING INCOMPLETE pH DATA LINE');
                disp(tline);
            end
        else
            disp('SKIPPING INCOMPLETE pH DATA LINE');
            disp(tline);
        end
        
    end
    
    if regexp(tline,'^<EOT>','once')
        EOT = EOT +1;
    end
    tline = fgetl(fid); % get 1st line
end
fclose(fid);

if bogus_time == 1
    disp([pH_file,' has bogus time stamps - subsituting EOP ',...
        'date from APF11'])
end

% ************************************************************************
% ASSIGN DATA TO STRUCTURE AND CLEAN UP
% ************************************************************************
d.data       = data;
d.hdr        = pH_hdr;
d.EOT        = EOT;
d.MSC_FW_ver = MSC_FW_ver;
%d.EOP_sdn    = EOP_sdn;
d.bogus_time = bogus_time;
clearvars -except d


