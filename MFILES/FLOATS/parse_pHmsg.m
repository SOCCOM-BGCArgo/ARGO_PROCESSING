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

%pH_file ='C:\temp\8514.034.dura'; % TEST


% ************************************************************************
% SET FORMATS AND VARIABLES
% ************************************************************************
d.hdr = ' '; % Returns if no data
d.data= [];
EOT = 0;
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
data_ct = 0;
comma_ct = [];

while ischar(tline) 
    if  regexp(tline,'^0x','once') % DATA LINES
        data_ct = data_ct + 1;
        comma_ct = [comma_ct;length(regexp(tline,','))] ; 
    end
    tline = fgetl(fid); % get 1st line
end

% check for data lines
if data_ct < 1
    disp(['File exists but no data lines found for : ',pH_file]);
    fclose(fid)
    return
end

% ************************************************************************
% GO BACK TO TOP OF FILE & THEN START PARSING THE DATA
% ************************************************************************
frewind(fid) % GO BACK TO TOP

d_cols   = size(pH_hdr,2);            % # data columns
comma_ct = mode(comma_ct);            % most frequent value
data     = ones(data_ct, d_cols)*NaN; % predim
tline    = ' ';                       % predim
ct = 0;                               % line counter   

while ischar(tline)                              % step through lines again
    if  regexp(tline,'^0x','once')            	 % data line?
        if length(regexp(tline,',')) == comma_ct % full line?
            ct = ct+1;
            dline = textscan(tline,pH_format,'Delimiter',',',...
                             'Collectoutput',1);
            data(ct,1) = datenum(dline{1,1}{1},'mm/dd/yyyy HH:MM:SS');
            data(ct,2:d_cols) = dline{1,2};
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

% ************************************************************************
% ASSIGN DATA TO STRUCTURE AND CLEAN UP
% ************************************************************************
d.data = data;
d.hdr  = pH_hdr;
d.EOT = EOT;
clearvars -except d


