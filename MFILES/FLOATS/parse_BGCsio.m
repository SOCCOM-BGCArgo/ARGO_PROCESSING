function d = parse_BGCsio(filename, verbosity, msgtype)

% PURPOSE:
%   This function parses an SIO BGC msg file.
%   A structure is returned containing the data
%
% USAGE:
%	data = parse_BGCsio(filename,merge)
%
% INPUTS:
%	filename  = calibration file name or path\file name as a string
%   verbosity = 0: data in = data out
%               1: only non fill Presure levels, shortened header ID's, &
%                  only parameters used for pH calc
%   msgtype      = type of msg: 'alk' 'dox' 'eco' 'ocr' 'phy'
%
% OUTPUTS:
%   d = a structure of data and calibration info
%       d.hdr          = data column ID's     (cell array of strings)
%       d.data         = sensor data          (numeric matrix)
%       d.aomlID       = float ID for aoml    (string)
%       d.floatID      = SIO float ID         (number as string 0 padded)
%       d.cast         = cast number          (scalar)
%       d.file_name    = source file name     (string)
%       d.MTIME_hdr    = header for extracted ascent profile timings
%       d.MTIME        = extracted ascent profile timings & MC codes
%       d.BOP_sdn      = start of profile time from MC = 500 line
%       d.EOP_sdn      = end of profile time from MC = 600 line
%       d.pres_axes_ct = number of unique pressure axes retuned by float
%       d.gps          = gps fix if available
%       d.EOP_sdn_str  = end of profile time as a date string for viewing
%       d.pres_levels  = sample levels count  (scalar)
%       d.sensors      = cell array or cell array of cell array listing
%                        sensors associated with pressure axes (sequential)
%       d.max_BGC_pres = deepest BGC sample   (scalar, dbar)
%       d.verbosity    = 0: full; 1:condensed (scalar 1 or 0)
%
% EXAMPLE:
%   jp = parse_ALKsio('09999_000997_0002.alk')
%
% Modified by jg; Created 01/13/2016 by jp
%
% REVISONS:
%      2/16/2022 TM: Modifications to reflect InAir handling of dox files
%      (to follow jg format changes).
%      07/18/2022  TM & JP: Modifications to "sensors_str" regular
%      expression to reflect jg format changes.
%      8/11/22     TM: If 'verbosity' is set to 1: ONLY remove NaN pressure lines
%       that are beyond Max Press!  If there are missing pressures along the
%       profile, preference is to record them as such).
%      01/25/2024, TM; added extraction of MC700 START OF TRANSMISSION as
%      an alternate option for EOP (end of profile...they should be only a
%      couple minutes off and allows for a timestamp to associate with our
%      ODV profiles).
%
% TESTING
%filename  = 'C:\temp\08666_000001_0042.alk'; msgtype   = 'alk';
%filename  = 'C:\temp\08666_000001_0042.dox'; msgtype   = 'dox';
%filename  = 'C:\temp\08666_000001_0042.phy'; msgtype   = 'phy';
%filename  = 'C:\temp\08666_000001_0042.eco'; msgtype    = 'eco';
%filename  = 'C:\temp\08666_000001_0042.ocr'; msgtype    = 'ocr';
%filename  = 'C:\temp\08666_000001_0042.no3'; msgtype    = 'no3';
% filename  = 'C:\temp\08666_000001_0005.phy'; msgtype   = 'phy'; %empty file
%filename  = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\MFILES\SIO_beta_code\ss0001\CTD\09040_000001_0013.phy';
% filename  = 'C:\temp\09040_000001_0000.phy';
% msgtype   = 'phy';
% verbosity = 1;

%
% ************************************************************************
% DO SOME PREP
% ************************************************************************
ftype_filter   = 'alk|dox|eco|ocr|phy';

%sensors_str    = '(?<=\[)\w+(?=\s+pressure axis\])|SBE\w+';
%sensors_str    = '(?<=PRESSURE AXIS\:)\w+|SBE\w+';
%sensors_str    = '(?<=PRESSURE AXIS\:).+(?=\])|SBE\w+';
%sensors_str    = '(?<=BGC[\w\s]*:)\w{3}.+(?=\])|SBE\w+';
%sensors_str    = '(?<=BGC[\w\s]*:)\w{3}.*(?=\])|SBE\w+'; % jp 06/15/22 .* instead of .+
% sensors_str    = '(?<=BGC[\w\s]*:)\w{3}.*(?=\s\()|SBE\w+'; % 07/18/22 tm attempt :-)
sensors_str    = '(?<=BGC:)\w{3}|SBE\w+'; %07/18/22 jp streamlining

dstart_str     = '^=========='; % data lines start shortly after this line
hdr_short_filt = '^\w+.+\)(?=\s+Offset)'; % regexp str to shorten hdr names
%Predim output if no data present
d.hdr          = {}; % cell array of strings
d.data         = []; % numeric matrix
d.aomlID       = ''; % string array
d.floatID      = ''; % string array
d.cast         = ''; % string array
d.file_name    = ''; % string array
d.pres_levels  = []; % numeric scalar
d.sensors      = {};
d.max_BGC_pres = []; % numeric scalar
d.verbosity    =  0;

if strcmp(msgtype, 'dox')
    d.InAir     = [];
elseif strcmp(msgtype, 'phy')
    d.MTIME_hdr    = {};
    d.MTIME        = [];
    d.pres_axes_ct = [];
    d.gps      = [];
    d.BOP_sdn      = [];
    d.EOP_sdn      = [];
    d.EOP_sdn_str  = '';
    dstart_str     = '^=================';
    %     hdr_short_filt = '^\w+.+\)';
    hdr_short_filt = '^\w+.+';
end

% CHECK FOR VALID FILE TYPE
ftype        = regexp(filename,'\w+$','once','match'); % get file suffix
if isempty(regexp(ftype, ftype_filter, 'once'))
    fprintf('Input file is not a valid file type : %s\n', filename);
    fprintf('Valid file extensions: %s\n', regexprep(ftype_filter,'|',' '));
    return
end

% CHECK FILE EXTENSION FOR PROPER FILE TYPE
if ~strcmp(ftype, msgtype)
    fprintf('Input msg file does not match expected file type (%s)\n', ...
        msgtype);
    fprintf('Input file: %s\n', filename);
    return
end

% CHECK FILE EXISTANCE
if ~isfile(filename)
    fprintf('Input file does not appear to exist at: %s\n',filename);
    return
end

% GET INFO FROM FILE NAME
expr_str = ['\w+(?=\.', msgtype,')'];
fn = regexp(filename, expr_str,'match','once'); % fn w\o extension
flt_info = regexp(fn,'\_','split');

if size(flt_info, 2) == 3
    d.aomlID    = flt_info{1};
    d.floatID   = flt_info{2};
    d.cast      = str2double(flt_info{3});
    d.file_name = filename;
else
    fprintf(['Could not resolve aomlID, floatID, & cast from file name:',...
        ' %s\n'],filename);
end

fid   = fopen(filename); % OPEN FILE FOR READING

% ************************************************************************
% ** IF CTD PHY FILE: GET ASCENDING TIMES, PRESSURE & MEASUREMENT CODES ***
% ************************************************************************
if strcmp(msgtype,'phy')
    MTIME_hdr = {'SDN' 'PRES' 'MC CODE'};
    date_fmt  = 'yyyymmddHHMMSS';
    MTIME     = ones(1000,size(MTIME_hdr,2))* NaN;
    m_ct      = 0;
    tline     = ' ';
    mt_toggle = 0;
    while ischar(tline)
        if regexp(tline,'^MC[567]','once') % measurement codes 5XX or 6XX
            %disp(tline)
            mc_code = str2double(regexp(tline,'(?<=MC)\d+','once','match'));
            if regexp(tline,'ASC MEAS PRES','once')
                mt_toggle = 1;
                m_ct = m_ct+1; %Increment 1st line of paired lines
                MTIME(m_ct,2) = str2double(regexp(tline,'[\d\.-]+$', ...
                    'match','once'));
                MTIME(m_ct,3) = mc_code;
            elseif mt_toggle == 1 & regexp(tline,'TIME STATUS\=','once')
                mt_toggle = 0;
                date_str = regexp(tline,'\d+$','match','once');
                MTIME(m_ct,1) = datenum(date_str, date_fmt);
            elseif regexp(tline, 'START OF TRANS','once')
                m_ct = m_ct+1;
                date_str = regexp(tline,'\d+$','match','once');
                MTIME(m_ct,1) = datenum(date_str, date_fmt);
                MTIME(m_ct,3) = mc_code;
            end
        end
        
        % if it existed assent timing is now over
        if regexp(tline,'^NUMBER OF PROFILES','once')
            d.pres_axes_ct = str2double(regexp(tline,'\d+$','match','once'));
            break
        end
        tline = fgetl(fid);
    end
    d.MTIME_hdr = MTIME_hdr;
    d.MTIME     = MTIME(1:m_ct,:);
    
    t1 = MTIME(:,3) == 600;
    t2 = MTIME(:,3) == 500; % assuming this is profile start meas code ???
    t1B = MTIME(:,3) == 700; % START OF TRANSMISSION.  Use this if MC600 EOP is 99999999999999!  (case cycle 110, 111 for solo0002).  Very close to end of profile time...
    if sum(t1) == 1
        d.EOP_sdn = MTIME(t1,1);
        d.EOP_sdn_str = datestr(MTIME(t1,1),'mm/dd/yyyy HH:MM:SS');
        if d.EOP_sdn > now % Not possible...this can occur when the date is fill value (99999999999999)
            d.EOP_sdn = MTIME(t1B,1);
            d.EOP_sdn_str = datestr(MTIME(t1B,1),'mm/dd/yyyy HH:MM:SS');
        end
    end
    if sum(t2) == 1
        d.BOP_sdn = MTIME(t2,1);
    end
    
    clear  mt_toggle mc_code m_ct MTIME MTIME_hdr date_str t1 t2
end

% ************************************************************************
%               **** GET HEADER FROM "COLUMN LINES" ****
% ************************************************************************
tline   = ' ';
hdr = cell(20,1); % columns because header ID's are long - easier to read
hdr_ct = 0;
sensor_ct = 0;
while isempty(regexp(tline, dstart_str,'once')) && ischar(tline)
    if regexp(tline,'PRESSURE LEVELS','once')
        d.pres_levels = str2double(regexp(tline,'\d+','match','once'));
    end
    
    % GET SENSOR TYPES FROM PRESSURE AXES ID
    if regexp(tline,'^\d+\.\s+VERTICAL','once') % data col ID line
        sensor_ct = sensor_ct + 1;
        str = regexp(tline, sensors_str,'match', 'once');
        d.sensors{sensor_ct} = str;
    end
    
    if regexp(tline,'^\d+\.\s+COLUMN','once') % data col ID line
        hdr_ct = hdr_ct + 1;
        str = regexp(tline,'(?<=COLUMN\s+)\w+.+','match','once');
        hdr{hdr_ct} = strtrim(str);
    end
    tline = fgetl(fid);
end
hdr  = hdr(1:hdr_ct);
rhdr = size(hdr,1);
data = ones(1000,rhdr) * NaN; % predim
airdata = ones(1000,rhdr) * NaN;
% ************************************************************************
%                     **** PARSE DATA LINES ****
% ************************************************************************
line_ct = 0;
gps_ct = 0;
air_ct = 0;
airkey = 0;
while ischar(tline)
    %For dox data ==> parse InAir:
    if airkey == 1
        tmpA = regexp(strtrim(tline),'\s+','split'); % delim = 1 or more spaces
        air_ct = air_ct+1;
        if size(tmpA,2) == rhdr % All good
            airdata(air_ct,:) = str2double(tmpA);
            % partial data line, fill in what you can
        elseif size(tmpA,2) < rhdr && size(tmpA,2) > 0
            fprintf('Partial InAir data line  detected...');
            airdata(air_ct,1:size(tmpA,2)) = str2double(tmpA);
        else
            disp('Oversized InAir data line - something is wrong!!')
        end
%     end
    % data line? zero or more spaces at begining followed by 1 or more #'s
    elseif regexp(tline,'^\s*-*\d+','once')
        tmp = regexp(strtrim(tline),'\s+','split'); % delim = 1 or more spaces
        line_ct = line_ct+1;
        if size(tmp,2) == rhdr % All good
            data(line_ct,:) = str2double(tmp);
            % partial data line, fill in what you can
        elseif size(tmp,2) < rhdr && size(tmp,2) > 0
            fprintf('Partial data line at %0.1f meter detected', ...
                str2double(tmp{1}));
            data(line_ct,1:size(tmp,2)) = str2double(tmp);
        else
            disp('Oversized data line - something is wrong!!')
        end
        
        % if phy file get gps fix make = [sdn lon lat]
    elseif strcmp(msgtype,'phy') && strncmp(tline,'MC703',5)
        tmp = regexp(tline,'\s+','split');
        if size(tmp,2) == 10 % all good
            gps_ct = gps_ct+1;
            dstr = [tmp{4},' ',tmp{5}];
            dnum = datenum(dstr,'yyyy/mm/dd HH:MM:SS');
			if str2double(tmp{3})>900 & str2double(tmp{2})>99 %JG gps fill value
            d.gps(gps_ct,1:3) = [dnum, nan, ...
                nan]; % sdn, lon, lat
				else
            d.gps(gps_ct,1:3) = [dnum, str2double(tmp{3}), ...
                str2double(tmp{2})]; % sdn, lon, lat
				end
        end
    elseif regexp(tline,'InAir','once')
        airkey = 1;
        air_ct = 0;
    end
    tline = fgetl(fid);
end

fclose(fid);
data           = data(1:line_ct,:);
tFILL          = data == -999 | data == -99;
data(tFILL)    = NaN; % replace fill values with NaN
iP             = strncmp(hdr,'PRESSURE (dbar)',15);
[d.max_BGC_pres,mpI] = max(data(:,iP),[],1,'omitnan');

% If oxygen data file, look for in air measurements & seperate
if strcmp(msgtype, 'dox')
    %     tAir    = data(:,1) == -9;
    %     d.InAir = data(tAir,:);
    %     data    = data(~tAir,:);
    tnanA = isnan(airdata(:,1));
    d.InAir.data = airdata(~tnanA,:);
    if size(airdata,2)==9
        d.InAir.hdr = {'pumpIndicator','ElapsedSec','NA','Pres','NA','PhaseDelay','NA','TempOpt','NA'};
    else
        d.InAir.hdr = {'pumpIndicator','ElapsedSec','NA','Pres','NA','PhaseDelay','NA','CalcDOXY','NA','TempOpt','NA'};
    end
end

% ************************************************************************
%                      ****  TIDY UP DATA  ****
% ************************************************************************
if verbosity == 1
    d.verbosity = 1;
    % REMOVE NaN PRESSURE LINES 
    % (TM NOTE Aug22: ONLY remove NaN pressure lines
    % that are beyond Max Press!  If there are missing pressures along the
    % profile, preference is to record them as such).
    data = data(mpI:end,:);
%     tnan = isnan(data(:,iP'));
%     data = data(~tnan,:);
    
    % REMOVE RESOLUTION,INDEX COLUMNS
    % May need to divide data by gain 1st if not always 1
    fmt  = '^RESOLUTION|^BGC|^SAMPLING|^CHECKSUM';
    tf   = cellfun(@isempty, regexp(hdr,fmt,'once'));
    hdr  = hdr(tf);
    data = data(:,tf');
    
    % SHORTEN HEADER ID NAMES
    hdr = regexp(hdr, hdr_short_filt,'once','match');
    
    % TM 2/16/22; after John's reformatting of InAir, columns headers are
    % slightly different than profile data...
    % % %     if strcmp(msgtype, 'dox')
    % % %         d.InAir = d.InAir(:,tf');
    % % %     end
end
d.data = data;
d.hdr  = hdr;

clear tFILL iP tnan data hdr