function data = parse_NAVISmsg4ARGO(msg_file)
% PURPOSE: 
%   This function parses an NAVIS biochemical float *.msg file and
%   returns a structure of raw data. The # of colums are determined from
%   the park, Aircal, low res & high res header lines. The file is
%   insprctrf twice: onec to determine header size & to get meta data &
%   a second time to extract profile data
%
%   Parsing steps: 
%       1) scan file for header lines (low res, park, Air sequence)
%          and meta data then go back to top of file
%       2) build cp mode header & determine cp hex format
%       3) scan again - extract data (park, spot, cp, air sequence, SurfObs)
%          If a *.cp file exists in the same path as the msg file it will
%          automatically get processed.
%       4) convert PTS values to usefull numbers
%       5) Sanity checks
%
% USAGE:
%	data = parse_NAVISmsgV2(file_name)
%
% INPUTS:
%	msg_file  = string of file name or  path\file name to */msg file
%
% OUTPUTS:
%   data = a structure of data
%
% EXAMPLES:
%
% CHANGES:
%    01/16/23 JP - re-vamp of original code. Got rid of CSAE / SWITCH blocks.
%                  Now processes OCR data & deals with cp data in seperate
%                  cp file, SorfaceOBS & in-air sequence flavors extracted


% TESTING
%msg_file = 'C:\temp\0037.003.msg';
%msg_file = 'C:\temp\0566.001.msg';
%msg_file = 'C:\temp\0063.002.msg';
%msg_file = 'C:\temp\0063.001.msg';
%msg_file = 'C:\temp\0571.002.msg';
%msg_file = 'C:\temp\1512.000.msg';
%msg_file = 'C:\temp\0062.003.msg';
%msg_file = 'C:\temp\1341.002.msg'; % 3XO2 NOT WORKING YET



% ************************************************************************
% PREDIMENSION OUTPUT STRUCTURE (If no data this is the output)
% ************************************************************************
data.cast           = NaN; % profile #
data.floatID        = '';  % float ID as string
data.sdn            = NaN; % profile termination time
data.gps            = [];  % gps location fix
data.FwRev_str      = '';
data.FwRev_date_str = '';
data.FwRev          = [];
data.CpActivationP  = [];
data.FlbbMode       = NaN; % not in Navis files
data.CTDtype        = '';
data.CTDsn          = '';
data.ice_flag       = 0;                 % NEW
data.ice_hex        = []; 
data.tf_cp_infile   = 0;
data.EOT            = 0; % can be greater than 1 if multiple gps fixes

data.pk_hdr = {}; % park data header cell array
data.pk_d   = []; % park data matrix
data.lr_hdr = {}; % low res header cell array
data.lr_d   = []; % low res data matrix
data.hr_hdr = {}; % high res header
data.hr_d   = []; % high res data
data.aircal_hdr = {}; % O2 surf & just below surf sequence data
data.aircal     = []; % O2 surf & just below surf sequence data
data.air        = []; % Surf O2 obs after air bladder has been inflated

if ~isfile(msg_file)
    fprintf('No message file found for %s!\n',msg_file);
    return
end

% CAST # FROM MESSAGE FILE PATH
str          = regexp(msg_file,'\d{3}(?=\.msg$)','match','once');
data.cast    = str2double(str); % cast # from file name
data.floatID = regexp(msg_file,'\d+(?=\.\d{3}\.msg$)','match','once');
clear str

% ************************************************************************
%                                PART 1
%           FIND ALL AVAILABLE HEADER LINES & GRAB META DATA
%                header lines will determine data formats
% ************************************************************************

% navis nautilus: lower case "p" & upper case 'Time' in Unix Time in
% park header & upper case "P" & lower case "time" O2AirCal header
pk_hdr_exp1    = '^\$\s+Date\s+p\s+t\s+s\s+'; 
pk_hdr_exp2    = '^\s+UNIX\s+Time\s+p'; % 1/14/24 JP Navis Nautilus
aircal_hdr_exp = '^\s+UNIX\s+time\s+pnm'; % 1/14/24 JP Navis Nautilus 0063
lr_hdr_exp     = '^\$\s+\s+p\s+t\s+s\s+';
hr_hdr_exp     = '^\$\s+\s+p\s+t\s+s\s+';
% can add to surf_obs_exp if OCR is added
surf_obs_exp   = '[-\.\d]+(?=dbar|us|V)'; % V appears to be T (0061,0062,0063)

fid   = fopen(msg_file);
tline = fgetl(fid);% initialize with first line
% FILE EXISTS BUT EMPTY INSIDE
if tline == -1 % Go to next i
    fprintf('File exist but empty inside - moving to next message file.\n')
    fclose(fid);
    data = []; % function will return empty value if no msg data
    return;
end

tf_o2_air = 0;

while ischar(tline) % loop through & find header lines & meta data

    % *********************************************************************
    % ******  META DATA FROM TOP OF FILE ( Lines start with "$ \w")  ******
    % header lines will have several spaces after $ vs only one for meta
    if regexp(tline, '^\$ \w','once') 
        if regexp(tline,'^\$ Mission config', 'once')
            data.FwRev_str = regexp(tline,'(?<=FwRev.+)\w+\s\w+', ...
                'once','match');
            data.FwRev_date_str = regexp(tline,'(?<=FwRev.+)\d+', ...
                'once','match');
            data.FwRev = str2double(data.FwRev_date_str);
        end
    
        % GET CpActivation pressure
        if regexp(tline,'^\$ CpAct', 'once')
            data.CpActivationP = str2double(regexp(tline,'\d+', ...
                'once','match'));
        end

% NO FLBBMode line in NAVIS
%         if regexp(tline,'\$ FlbbMode', 'once') 
%             data.FlbbMode = str2double(regexp(tline,'\d+','once','match'));
%         end
        
        % GET profile termination time
        if regexp(tline,'^\$ Profile', 'once')
            d_str = regexp(tline,'(?<=terminated\:\s+\w{3}).+', ...
                'once','match');
            if isempty(d_str)
                fprintf('No termination time found for %s! No data!\n',...
                    file_name);
                fclose(fid);
                return
            else
                data.sdn = datenum(d_str,'mmm dd HH:MM:SS yyyy'); 
            end
        end

    % ********************************************************************
    % **************   META DATA LINES AT BOTTOM OF FILE  ****************
    % pertinent meta lines at bottom of file start with a letter
    elseif regexp(tline, '^[a-zA-z]','once') 

        % Under ice?
        if regexp(tline,'^IceEvasionRecord','once')
            ihex = regexp(tline,'(?<=\=)\w+','once','match'); % extract hex
            ice_hex = hex2bin(ihex, 8);             % convert hex to binary
            data.ice_hex = ice_hex;
            if ice_hex(8) == 1 % look last digit of hex 8
                data.ice_flag = 1; % change ice flag to 1 if under ice now
            end
        end

        % Get GPS fix(es) - a bunch of data checking here for complete data
        % but I can't remember the history (JP 01/14/24)
        % re-aranged so sdn lon lat & always a sdn (gps or profile)
        if regexp(tline,'^Fix:','once')
            gps = sscanf(tline,'%*s %f %f',2); % lon lat
            if ~isempty(gps)
                gps_sdn = char(sscanf(tline,'%*s %*f %*f %17c',1))'; %[mm/dd/yyyy; hhmmss]
                if size(gps_sdn,2) == 17
                    Gsdn = datenum(gps_sdn,'mm/dd/yyyy HHMMSS');
                    data.gps = [data.gps; [Gsdn gps']];
                else
                    data.gps = [data.gps; [data.sdn NaN NaN]];
                end
            else
                data.gps = [data.gps; [data.sdn NaN NaN]];
            end
        end

        if regexp(tline, '^O2AirCal\:','once') % i.e. 0061, 0062, 0063
            % No header provided for some air cal sequences so use a
            % counter as a flag (ie 0061, 0062
            tf_o2_air = tf_o2_air +1;
        end

     % *****************  FIND PARK HEADER LINE  *********************
    elseif regexp(tline, [pk_hdr_exp1,'|',pk_hdr_exp2],'once') % park header
        hdr_str = regexprep(tline, '\$', ' ');
        hdr_str = strtrim(regexprep(hdr_str, 'UNIX\s+Time', 'UNIXTime'));
        data.pk_hdr = regexp(hdr_str,'\s+', 'split');

        if strcmp('Date', data.pk_hdr{1})
            pk_format = ['%*s%s%s%s%s',repmat('%f',1,size(data.pk_hdr,2)-1)];
        else
            pk_format = repmat('%f',1,size(data.pk_hdr,2));
        end
   
    % ****************  FIND SPOT SAMPLE HEADER LINE  ******************
    elseif regexp(tline, lr_hdr_exp,'once') % low res data header found
        hdr_str = regexprep(tline, '\$', ' ');
        data.lr_hdr  = regexp(strtrim(hdr_str),'\s+', 'split');
        % below doesn't need to happen but not certain about downstream affects
        data.lr_hdr = data.lr_hdr'; 
        lr_format = repmat('%f',1, size(data.lr_hdr,1));

    % ***********   FIND CTD SERIAL # & STATS LINE   ******************
    % if it exists! 0063 has seperate cp file & no cp serial # line in
    % msg file, 0061 & 0062 & most older navis msg files have it
    elseif regexpi(tline,'^#.+sbe','once') % cp serial #  & stats line
         data.tf_cp_infile = 1;
         % if not found here try & get from cp file later
         data.CTDtype = regexpi(tline,'sbe\w+(?=serno)','match', 'once');
         data.CTDsn   = regexpi(tline,'(?<=serno[)\d+','match', 'once');

    % ***********   FIND IN-AIR SEQUENCE HEADER LINE   ******************
    % if it exists! 0063 has one (seperate cp file) 0061 & 0062 do not
    elseif regexp(tline, aircal_hdr_exp,'once') % O2 Air cal header found
        hdr_str         = strtrim(regexprep(tline, 'UNIX\s+time', 'UNIXTime'));
        data.aircal_hdr = regexp(hdr_str,'\s+', 'split');
        air_format   = repmat('%f',1,size(data.aircal_hdr,2));

    elseif regexp(tline,'^<EOT>','once') % END OF FILE COUNTER
        data.EOT = data.EOT +1;
    end

    tline = fgetl(fid);
end

if isempty(data.lr_hdr)
    fprintf('Incomplete message file! No profile data for %s\n',msg_file);
    fclose(fid);
    return
end

if isempty(data.gps) % NO GPS FIX FOUND SET DUMMY ARRAY WITH PROFILE SDN
    data.gps = [data.sdn NaN NaN];
end

% CREATE AIR CAL SEQUENCE HEADER IF NEED BE
if isempty(data.aircal_hdr) && sum(tf_o2_air) > 0 % i.e. 0061, 0062
    data.aircal_hdr = {'Date' 'pnm' 'P' 'O2ph' 'O2T'}; % no hdr provided
    air_format = ['%*s%s%s%s%s',repmat('%f',1,size(data.aircal_hdr,2)-1)];
end

frewind(fid); % go back to top of file (prep for data extraction)

% ************************************************************************
%                                  PART 2
% BUILD CP DATA HEX FORMAT LINE, HEX CONVERSION MATRIX FOR PTS & CP HEADER
% cp header built off of spot sample header so need to add bin count info
% for all cp sensors. I chose to group bin counts at the end of the matrix
% ************************************************************************

% OK NOW DETERMINE HEX FORMAT FOR CP DATA
% [> is neg number, scale factor]
pH_flag  = 1; % default
OCR_flag = 0; % default
hex_ind  = [1:3, 5:6, 8:10]; % col indices for PTS, O2 Phase&T, MCOMS 
nbin_ind = [4,7,11]; %nbins for pts, o2, mcoms
nbin_hdr = {'nbin ctd';'nbin oxygen';'nbin MCOMS'}; % putting bins at end
conv_ind = 1:3; % col indices for PTS
hex_conv = [32768 10; ...     % p  these are used to convert to pts
            61440 1000; ...   % t
            61440 1000];      % s

if size(intersect({'no3' 'O2ph' 'Mch1' 'OCRch1' 'tilt' 'phVrs'},data.lr_hdr),1) == 6
% 6 senssor Navis Nautilus
    OCR_flag = 1;
    hex_format  =['%04x%04x%04x%02x%06x%06x%02x%06x%06x%06x%02x', ...
        '%06x%06x%06x%06x%02x%02x%02x%06x%02x'];
    hex_ind  = [hex_ind, 12:15,17:19]; % PTS, O2 Phase&T, MCOMS, OCR, tilt&std, Vrs
    nbin_ind = [nbin_ind, 16,20];
    nbin_hdr = [nbin_hdr; 'nbin OCR'; 'nbin pH'];

elseif size(intersect({'phV' 'phT'},data.lr_hdr),1) == 2
% O2, MCOMS, pH (OLD)    
    hex_format  = '%04x%04x%04x%02x%06x%06x%02x%06x%06x%06x%02x%06x%04x%02x';
    hex_ind     = [hex_ind, 12:13]; % PTS, O2 Phase&T, MCOMS, pH V&T 
    nbin_ind    = [nbin_ind, 14];
    nbin_hdr    = [nbin_hdr; 'nbin pH'];
    hex_conv    = [hex_conv; 61440 1000]; % Add for pH T
    conv_ind    = [conv_ind, 10]; % add col indice for pH T 

elseif size(intersect({'phVrs' 'phVk'},data.lr_hdr),1) == 2
% O2, MCOMS, pH  (NEW 10/27/20 NO pHT & different hdr names)    
    hex_format  = '%04x%04x%04x%02x%06x%06x%02x%06x%06x%06x%02x%06x%02x';
    hex_ind     = [hex_ind, 12]; % PTS, O2 Phase&T, MCOMS, pH V
    nbin_ind    = [nbin_ind, 13];
    nbin_hdr    = [nbin_hdr; 'nbin pH'];

else % O2, MCOMS
    pH_flag     = 0;
    hex_format  = '%04x%04x%04x%02x%06x%06x%02x%06x%06x%06x%02x';
end

hex_char_ct = cumsum(str2double(regexp(hex_format,'\d+','match')));
hex_length  = max(hex_char_ct);
hex_cols    = size(regexp(hex_format,'x'),2);

% BUILD HIGH RES HEADER (CP DATA HEADER) - all but NO3 & pH diagnostics
tKEEP = cellfun(@isempty, regexp(data.lr_hdr,'no3|phVk|phIb|pHIk','once'));
data.hr_hdr = [data.lr_hdr(tKEEP); nbin_hdr]; % add bin counts to end

% No LR header ID for tilt std so add it in to HR header if needed
iT = find(strcmp(data.hr_hdr,'tilt') == 1);
if size(iT,1) == 1
    data.hr_hdr = [data.hr_hdr(1:iT,1);'STDtilt';data.hr_hdr(iT+1:end,1)];
end

% ************************************************************************
%                                  PART 3
%        PARSE DATA LINES (PARK, SPOT, CP, SURF AIR  SEQ, SURF OBS)
% ************************************************************************
% initialize counters
pk_ct = 0; spot_ct = 0; cp_ct = 0; air_ct = 0; surfobs_ct = 0;

pk_data   = ones(1000, size(data.pk_hdr,2)) * NaN;
spot_data = ones(1000, size(data.lr_hdr,1)) * NaN;
cp_data   = ones(1000, size(data.hr_hdr,1)) * NaN;
cp_mode   = 0; % 1 means in cp data mode block
air_data  = ones(100, size(data.aircal_hdr,2)) * NaN;


tline = fgetl(fid);% initialize with first line
while ischar(tline)
    tline = strtrim(tline);

    % PARK SAMPLE DATA LINES
    if regexp(tline,'^ParkObs\:','once')
        pk_tmp = textscan(tline, pk_format, 'CollectOutput',1);
        dstr   = sprintf('%s',pk_tmp{1,1}{:}); % merged str of date parts
        sdn    = datenum(dstr,'mmmddyyyyHH:MM:SS');
        pk_ct  = pk_ct+1;
        pk_data(pk_ct,:) = [sdn pk_tmp{2}];
    elseif regexp(tline,'^PPtM\:','once')
        pk_tmp = cell2mat(textscan(tline, ['%*s',pk_format], ...
            'CollectOutput',1)); % all numeric return
        % Unix time = number of seconds since Jan 1, 1970 00:00:00.
        pk_tmp(:,1)      = datenum(1970,1,1,0,0,0) + pk_tmp(:,1)/86400;
        pk_ct            = pk_ct+1;
        pk_data(pk_ct,:) = pk_tmp;

    % ***************  SPOT SAMPLE DATA LINES ********************
    % starts with press value & could be negative at the surface
    % (souldn't be though)
    elseif regexp(tline,'^-*\d+\.\d+','once') & ...
            isempty(regexp(tline,'\(Park Sample\)$','once'))
        tmp     = textscan(tline,lr_format,1,'CollectOutput',1);
        spot_ct = spot_ct+1;
        spot_data(spot_ct,:) = tmp{1}; % tmp1: p, t, s, etc

    elseif regexp(tline,'^ser1\:','once') % ENTERING CP MODE BLOCK
        cp_mode = 1;
    elseif regexp(tline,'^Resm','once')   % LEAVING CP MODE BLOCK
        cp_mode = 0;

    % ************* IN-FILE CP SAMPLE DATA LINES *********************
    elseif data.tf_cp_infile == 1 && cp_mode == 1
        % JP 12/12/19, MODIFIED TO BETTER CATCH PARTIAL HEX LINES W DATA
        %HEX CHARS BEGINING OF STR, EXPECTED LENGTH & ALL CHARS HEX CHARS?
        if size(tline,2) == hex_length && ...
                size(regexp(tline,'^[A-F0-9]+','match','once'),2) == hex_length
            tmp = sscanf(tline,hex_format);
            if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
                cp_ct = cp_ct+1;
                cp_data(cp_ct,:) = tmp'; % cp
            end

        %SIZE ~= HEX_LENGTH & ALL HEX CHARS, BUT PARTIAL OR SHORT LINE
        elseif size(tline,2) > 11 && ... % 12 chars means p t & s (14 for bin count too)
                size(tline,2) == size(regexp(tline,'^[A-F0-9]+','match','once'),2)
            % figure out # of valid measurements
            ind = find(hex_char_ct <= size(tline,2),1,'last'); % jp 7/2020

            [tmp,N] = sscanf(tline,hex_format);
            if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
                dtmp = ones(hex_cols,1)* NaN; % predim to excepted variable size input jp 7/2020
                dtmp(1:ind) = tmp(1:ind); % jp 7/2020
                cp_ct = cp_ct+1;
                cp_data(cp_ct,:) = dtmp'; % cp
            end
        end

    % ***************  IN-AIR SEQUENCE DATA LINES *******************
    elseif regexp(tline,'^O2AirCal:','once')
        air_ct  = air_ct+1;
        if strcmp(data.aircal_hdr{1},'Date') % 0061, 0062
            air_tmp = textscan(tline, air_format, 'CollectOutput',1);
            dstr   = sprintf('%s',air_tmp{1,1}{:}); % merged str of date parts
            sdn    = datenum(dstr,'mmmddyyyyHH:MM:SS');
            air_data(air_ct,:) = [sdn air_tmp{2}];
        else % 0063
            air_tmp = textscan(tline, ['%*s',air_format], 'CollectOutput',1);
            air_data(air_ct,:) = air_tmp{1};
        end

    % ******************   SINGLE SURF OBS   *********************
    % CTD is off but could be more than one line of data. There could also
    % be multiple O2 sensors or possibly OCR in the future. Size is not
    % predetermined so build array on the fly vs predim. Add if blocks for
    % more cases
    elseif regexp(tline,'^SurfaceObs','once') % return [temp, phase]
        tmp_cell = regexp(tline, surf_obs_exp,'match'); % P phase T
        tmp_data = str2double(tmp_cell);
        if size(tmp_data,2) == 3 % P, phase T
            data.air = [data.air; tmp_data([1,3,2])]; % P, T, Phase
        end
    end

    tline = fgetl(fid);
end
fclose(fid);
clear fid tline

data.pk_d   = pk_data(1:pk_ct,:);
data.lr_d   = spot_data(1:spot_ct,:);
data.hr_d   = cp_data(1:cp_ct,:);
% re-arrange cp data: move bin counts to end
data.hr_d   = data.hr_d(:,[hex_ind, nbin_ind]);
data.aircal = air_data(1:air_ct,:);

% ***********************************************************************
% *********   CHECK IF CP DATA ARE IN SEPERATE *.CP FILE ****************
% some redundant coding for hex parsing but I can't think of a better way

cp_fp = regexprep(msg_file,'msg','cp'); % create path to *.cp file
if data.tf_cp_infile == 0 && isfile(cp_fp) % No data & *.cp file exist
    cp_ct   = 0; % counter
    cp_data = ones(1000, size(data.hr_hdr,1)) * NaN; % predim

    fid2  = fopen(cp_fp);
    tline = fgetl(fid2); % init with 1st line
    while ischar(tline)
        if regexpi(tline,'^#.+sbe','once') % cp serial #  & stats line
            data.CTDtype = regexpi(tline,'sbe\w+(?=serno)','match', 'once');
            data.CTDsn   = regexpi(tline,'(?<=serno[)\d+','match', 'once');

        elseif size(tline,2) == hex_length && ...
                size(regexp(tline,'^[A-F0-9]+','match','once'),2) == hex_length
            tmp = sscanf(tline,hex_format);
            if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
                cp_ct = cp_ct+1;
                cp_data(cp_ct,:) = tmp'; % cp
            end

        %SIZE ~= HEX_LENGTH & ALL HEX CHARS, BUT PARTIAL OR SHORT LINE
        elseif size(tline,2) > 11 && ... % 12 chars means p t & s (14 for bin count too)
                size(tline,2) == size(regexp(tline,'^[A-F0-9]+','match','once'),2)
            % figure out # of valid measurements
            ind = find(hex_char_ct <= size(tline,2),1,'last'); % jp 7/2020

            [tmp,N] = sscanf(tline,hex_format); % recover partial data
            if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
                dtmp = ones(hex_cols,1)* NaN; % predim to excepted variable size input jp 7/2020
                dtmp(1:ind) = tmp(1:ind); % jp 7/2020
                cp_ct = cp_ct+1;
                cp_data(cp_ct,:) = dtmp'; % cp
            end
        end

        tline = fgetl(fid2);
    end
    fclose(fid2);
    data.hr_d   = cp_data(1:cp_ct,:);
    data.hr_d   = data.hr_d(:,[hex_ind, nbin_ind]); % move bin counts to end
end

% ************************************************************************
%                                  PART 4
%             CONVERT PTS VALUES to USEFULL NUMBERS AND/OR COUNTS
%
% THIS IS A MODIFICATION AND CONDENSING OF Dan Quittman's CODE FOR:
% hextop.m, hextot.m and hextos.m functions
% for converting 16 bit 2's complement format to decimal
% It appears this is also a way to deal with a 12 bit A/D board in a
% 16 bit world as well as signed integers in a hex world
% ************************************************************************
if ~isempty(data.hr_d)

    tmp     = data.hr_d(:,conv_ind); % p, t, s OR p, t, s, pH t OR p, t, s, pH
    hex_var1 = ones(size(tmp(:,1))) * hex_conv(:,1)'; % matrix
    hex_var2 = ones(size(tmp(:,1))) * hex_conv(:,2)'; % matrix
    tNaN = tmp*NaN; % set up NaN flags
    tNaN(tmp - hex_var1 ~= 0) = 1; % No NaN set to 1
    tHi  = tmp - hex_var1 > 0;
    tLo  = tmp - hex_var1 < 0;

    data.hr_d(:,conv_ind) = (tHi .*(tmp-65536)./ hex_var2 + ...
        tLo .* tmp ./ hex_var2).* tNaN;

    % NOW DO BIO-SENSORS
    rail_vals = data.hr_d == 2^24-1; % Check out of bounds A/D
    data.hr_d(rail_vals) = NaN;
    data.hr_d(:,4)   = (data.hr_d(:,4)/100000) - 10; % O2 Phase
    data.hr_d(:,5)   = (data.hr_d(:,5)/1000000) - 1; % O2 temperature volts
    data.hr_d(:,6:8) = data.hr_d(:,6:8) - 500; % MCOMS

    if pH_flag == 1
        % find ph column
        pH_ind = ~cellfun(@isempty, regexp(data.hr_hdr,'pHVrs|phV','once'));
        data.hr_d(:,pH_ind')   = (data.hr_d(:,pH_ind')/1000000) - 2.5; % pH volts
    end

    if OCR_flag == 1
        OCR_ind     = ~cellfun(@isempty, regexp(data.hr_hdr,'^OCR','once'));
        tilt_ind    = ~cellfun(@isempty, regexp(data.hr_hdr,'^tilt','once'));
        stdtilt_ind = ~cellfun(@isempty, regexp(data.hr_hdr,'^STDtilt','once'));
        data.hr_d(:,OCR_ind') = data.hr_d(:,OCR_ind') .* 1024 + 2013265920;
        if sum(tilt_ind) == 1
            data.hr_d(:,tilt_ind') = data.hr_d(:,tilt_ind') ./ 10;
        end
        if sum(stdtilt_ind) == 1
            data.hr_d(:,stdtilt_ind') = data.hr_d(:,stdtilt_ind') ./ 100;
        end
    end
end

% CHECK LR OCR DATA FOR '0' VALUES. THESE SHOULD BE NaN'S AS NO
% MEASUREMENTS WERE TAKEN. ALSO DO COVERSION TO GET TRUE OCR COUNTS
if ~isempty(data.lr_d)

    if OCR_flag == 1
        OCR_ind = ~cellfun(@isempty, regexp(data.lr_hdr,'^OCR','once'));
        tmp = data.lr_d(:,OCR_ind');
        tmp(tmp == 0) = NaN;
        data.lr_d(:,OCR_ind') = tmp;
        data.lr_d(:,OCR_ind') = data.lr_d(:,OCR_ind') .* 1024 + 2013265920;
    end

end

% ************************************************************************
%                                  PART 5
% **********************   CLEAN UP & SANITY CHECKS   ********************
% ************************************************************************

tnan = isnan(data.lr_d(:,1));
if isempty(data.lr_d) && isempty(data.hr_d)
    fprintf('No data found in %s\n',msg_file);
elseif isempty(data.hr_d)
    fprintf('No cp data found in %s\n',msg_file);
elseif isempty(data.lr_d)
    fprintf('No spot sample data found in %s\n',msg_file);
elseif sum(tnan) > 0
    fprintf(['%0.0f missing pressure values found in % and ', ...
        'will be discarded\n'], sum(tnan),msg_file);
    data.lr_d(tnan,:) = [];
end

if isempty(data.pk_d)
    fprintf('No park data found for %s\n',msg_file);
end

clearvars -except data




           






