function d = parse_NO3sio(filename, verbosity)

% PURPOSE: 
%   This function parses an SIO NO3 msg file.
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
%       d.cast         = cast number          (number as string 0 padded)
%       d.file_name    = source file name     (string)
%       d.pres_levels  = sample levels count  (scalar)
%       d.max_BGC_pres = deepest BGC sample   (scalar, dbar)
%       d.verbosity    = 0: full; 1:condensed (scalar 1 or 0)
%
%       d.spectra_pix_range = pixel range returned by sensor (1x2 array)
%       d.BISTb_hdr         = msg file wl & ecoef data hdr (cell array} 
%       d.BISTb_data        = msg file wl & ecoef data (matrix)
%       d.DC                = sample DC intensity, counts (array)
%       d.UV_INTEN          =  UV sample intensities, counts (matrix)
%
% EXAMPLE:
%   jp = parse_NO3sio_jp('08666_000001_0042.no3')
%
% Modified by jg; Created 01/13/2016 by jp

%      8/11/22     TM: If 'verbosity' is set to 1: ONLY remove NaN pressure lines
%       that are beyond Max Press!  If there are missing pressures along the
%       profile, preference is to record them as such).
%
% REVISONS:
%
% TESTING
% filename  = 'C:\temp\08666_000001_0042.no3'; msgtype    = 'no3';
% verbosity = 1;

%
% ************************************************************************
% DO SOME PREP
% ************************************************************************
% SENSOR OFFSET JUST A GUESS, dbar (& I guess this would change with press
% a bit too!
d.tf_sensor_offset = 0; %TM 2/15/22 don't apply this direct to the pressure axis.  Only apply the offset in the T/S interpolation
d.sensor_offset    = -0.8; 

dstart_str     = '^=========='; % data lines start shortly after this line
%Predim output if no data present
d.hdr          = {}; % cell array of strings
d.data         = []; % numeric matrix
d.aomlID       = ''; % string array
d.floatID      = ''; % string array
d.cast         = ''; % string array
d.file_name    = ''; % string array
d.pres_levels  = []; % numeric scalar
d.max_BGC_pres = []; % numeric scalar
d.verbosity    =  0; % logical 1 or 0

d.spectra_pix_range = [NaN NaN]; % 1x2 numeric scalar
d.BISTb_hdr         = {}; % cell array of strings
d.BISTb_data        = []; % numeric matrix
d.DC                = []; % numeric array
d.UV_INTEN          = []; % numeric matrix
d.SDN               = []; 
% CHECK FOR VALID FILE TYPE
ftype        = regexp(filename,'\w+$','once','match'); % get file suffix
if isempty(regexp(ftype, 'no3', 'once'))
    fprintf('Input file is not valid no3 file : %s\n', filename);
    fprintf('Valid file extension: no3\n');
    return
end

% CHECK FILE EXISTANCE
if ~isfile(filename)
    fprintf('Input file does not appear to exist at: %s\n',filename);
    return
end

% GET INFO FROM FILE NAME
expr_str = ['\w+(?=\.no3)'];
fn = regexp(filename, expr_str,'match','once'); % fn w\o extension
flt_info = regexp(fn,'\_','split');
if size(flt_info, 2) == 3
    d.aomlID    = flt_info{1};
    d.floatID   = flt_info{2};
    d.cast      = flt_info{3};
    d.file_name = filename;
else
    fprintf(['Could not resolve aomlID, floatID, & cast from file name:',...
        ' %s\n'],filename);
end

fid   = fopen(filename); % OPEN FILE FOR READING

% ************************************************************************
% ****                IF NO3, GET PIXEL BLOCK INFO                     ****
% ************************************************************************
tline    = ' ';
tf_BISTb = 0;
pix_ct   = 0;
while ischar(tline)
    if regexp(tline,'^BIST FIT PIXEL BEGIN', 'once')
        d.spectra_pix_range(1) = ...
            str2double(regexp(tline,'\d+$','match','once'));
    elseif regexp(tline,'^BIST FIT PIXEL END', 'once')
                d.spectra_pix_range(2) = ...
            str2double(regexp(tline,'\d+$','match','once'));
    elseif tf_BISTb == 0 & regexp(tline,'^BISTb CALIBRATION HEADER', 'once')
        tf_BISTb   = 1; % toggle on for in data block
        % delim = 0 or more commas followed by 1 or more spaces
        tmp        = regexp(tline, ',*\s+','split'); 
        BISTb_hdr  = ['Pixel', 'Wavelength', tmp(8:10)];
        BISTb_data = ones(100,size(BISTb_hdr,2))* NaN;
    elseif tf_BISTb == 1 && strncmp(tline,'BISTb SPEC',10)
        pix_ct               = pix_ct+1;
        tmp                  = regexp(tline,'\s+','split');
        BISTb_data(pix_ct,:) = str2double(tmp([4,6:9]));
    elseif tf_BISTb == 1 && ~strncmp(tline,'BISTb SPEC',10)
        break % end of pixel calibration block
    end
 tline = fgetl(fid);   
end
if pix_ct > 0
    d.BISTb_hdr = BISTb_hdr;
    d.BISTb_data = BISTb_data(1:pix_ct,:);
end
clear tmp pix_ct BISTb_data BISTb_hdr tf_BISTb

% ************************************************************************
%               **** GET HEADER FROM "COLUMN LINES" ****
% ************************************************************************   
tline   = ' ';
hdr = cell(20,1); % columns because header ID's are long - easier to read
hdr_ct = 0;
while isempty(regexp(tline, dstart_str,'once')) && ischar(tline)
    if regexp(tline,'PRESSURE LEVELS','once')
        d.pres_levels = str2double(regexp(tline,'\d+','match','once'));
    end
    
    if regexp(tline,'^\d+\.\s+COLUMN','once') % data col ID line
        hdr_ct = hdr_ct + 1;
        % regexp filt: find 0 or more "<" followed by 1 or more letters or
        %              numbers followed by 1 or more of anything, but all
        %              this must be preceded by "COLUMN" and a bunch of
        %              spaces
        str = regexp(tline,'(?<=COLUMN\s+)<*\w+.+','match','once');
        hdr{hdr_ct} = strtrim(str);
    end
    tline = fgetl(fid);
end
hdr  = hdr(1:hdr_ct);

% ADD PIXEL COUNTS TO HDR
if ~isempty(d.BISTb_data)
    hdr = [hdr; cellstr(num2str(d.BISTb_data(:,1)))];
else
    warning('no NO3 data found!')
    d.data = [];
    return;
end

rhdr = size(hdr,1);
data = ones(1000,rhdr) * NaN; % predim

% ************************************************************************
%                     **** PARSE DATA LINES ****
% ************************************************************************
line_ct = 0;
gps_ct = 0;
while ischar(tline)
    % data line? zero or more spaces at begining followed by 1 or more #'s
    if regexp(tline,'^\s*-*\d+','once') 
        tmp = regexp(strtrim(tline),'\s+','split'); % delim = 1 or more spaces
        line_ct = line_ct+1;
        if size(tmp,2) == rhdr % All good
            data(line_ct,:) = str2double(tmp);
        % partial data line, fill in what you can
        elseif size(tmp,2) < rhdr && size(tmp,2) > 0
            fprintf('Partial data line at %0.1f meter detected\', ...
                str2double(tmp{1}));
            data(line_ct,1:size(tmp,2)) = str2double(tmp);
        else
            disp('Oversized data line - something is wrong!!')
        end
    end
    tline = fgetl(fid);
end

fclose(fid);
data           = data(1:line_ct,:);
tFILL          = data == -999 | data == -99;
data(tFILL)    = NaN; % replace fill values with NaN
iP             = strncmp(hdr,'PRESSURE (dbar)',15);
[d.max_BGC_pres,mpI] = max(data(:,iP),[],1,'omitnan');

% APPLY PRESSURE OFFSET TO SENSOR PRESSURE (SUNA AT BOTTOM OF FLOAT)
if d.tf_sensor_offset
    tnan = isnan(data(:,iP'));
    data(~tnan,iP') = data(~tnan,iP') + d.sensor_offset;
    fprintf(['A %0.2f dbar offset has been added to the nitrate ', ...
        'pressure axis values\n'],d.sensor_offset);
end
    

% ************************************************************************
%                      ****  TIDY UP DATA  ****
% ************************************************************************
if verbosity == 1
    d.verbosity = 1;
    % REMOVE NaN PRESSURE LINES
        data = data(mpI:end,:);
%     tnan = isnan(data(:,iP'));
%     data = data(~tnan,:);
    
    % REMOVE RESOLUTION,INDEX COLUMNS
    % May need to divide data by gain 1st if not always 1
    fmt  = '^RESOLUTION|^BGC|^SAMPLING|FIT_ERROR';
    tf   = cellfun(@isempty, regexp(hdr,fmt,'once'));
    hdr  = hdr(tf);
    data = data(:,tf');
    
    % SHORTEN HEADER ID NAMES
    hdr_short_filt = '^\w+.+\)(?=\s+Offset)|^\d+'; % regexp str to shorten hdr names
    hdr = regexp(hdr, hdr_short_filt,'once','match');
end
d.data = data;
d.hdr  = hdr;

% LAST - MAKE SOME NO3 SPECIFIC VARIABLES
iDRK = strncmp(d.hdr,'UV_INT_DARK',11)';
d.DC = d.data(:,iDRK);
d.SDN = d.DC*NaN; % This is a junk variable to play nice with calc_FLOAT_NO3.m
iUV = ~cellfun(@isempty,regexp(d.hdr,'^\d+','once'))';
d.UV_INTEN = d.data(:,iUV);

% TM: Add missing fields just for consistency across float types (although not needed in code for the
% SOLO SUNA!!)
% d.SWDC = [];
d.WL_fit_win = [217 240];
d.pix_fit_win = NaN;
d.CalTemp = NaN;
d.CalDate = NaN;
d.CalDateStr = '';

%clear tFILL iP tnan data hdr iDRK iUV
clearvars -except d