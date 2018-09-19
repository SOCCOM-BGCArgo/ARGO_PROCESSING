function data = parse_NAVISmsg4ARGO(msg_file)
% PURPOSE: 
%   This function parses an NAVIS biochemical float *.msg file and
%   returns a structure of raw data. The # of colums are determined from
%   the low resolution header line
%
% USAGE:
%	data = parse_NAVISmsg4ARGO(file_name)
%
% INPUTS:
%	file_name  = string of file name or  path\file name
%
% OUTPUTS:
%       data =  a structure of the data and headers. Size varies depending
%               on data flag.
%           data.FwRev  = Firmware revision number
%           data.CpActivationP =  Activation pressure (dbar) for cp?
%           data.FlbbMode = Flag in newer firmware 1=Yes,0=NO,[] no exist
%           data.lr_hdr = low res header cell array
%           data.lr_d   = low res data matrix
%           data.hr_hdr = high res header cell array
%           data.hr_d   = high res data matrix
%           data.pk_hdr = park sample header cell array
%           data.pk_d   = park sample data matrix
%           data.cast   = profile number
%           data.sdn    = profile termination date [matlab sdn]
%           data.gps    = gps location fix [lon, lat]
%           data.air    = Air oxygen measurements [Temp phase (Rphase)]
% EXAMPLES:
%   parse_NAVISmsg4ARGO('c:\temp\7601.003.msg')
%   parse_NAVISmsg4ARGO('0508.003.msg')
%
% CHANGES
%   01/18/2017 - cp header line search more generic & extracted CTD
%               type and serial number ultimately for ODV meta info
%   06/09/2017 - add code to look for <EOT> and record # of instances

% ************************************************************************
% FORMATS & VARIABLES
% ************************************************************************
% FOR TESTING
%msg_file = 'c:\temp\0506.003.msg'; % TESTING
%msg_file = 'c:\temp\0509.018.msg'; % TESTING
%msg_file = 'c:\temp\0276.002.msg'; % TESTING
%msg_file = 'c:\temp\0037.005.msg'; % TESTING
%msg_file = 'c:\temp\0511.008.msg';
 
% PREDIMENSION OUTPUT STRUCTURE (If no data this is the output)
data.pk_hdr = {}; % park data header cell array
data.pk_d   = []; % park data matrix
data.lr_hdr = {}; % low res header cell array
data.lr_d   = []; % low res data matrix
data.hr_hdr = {}; % high res header
data.hr_d   = []; % high res data
data.cast   = NaN; % profile #
data.sdn    = NaN; % profile termination time
data.gps    = []; % gps location fix
data.air    = []; % Air oxygen measurements
data.FwRev  = [];
data.CpActivationP = [];
data.FlbbMode = NaN;
data.CTDtype = '';
data.CTDsn   = '';
data.EOT    = 0;

% ************************************************************************
% GET PARK AND DATA HEADERS, THEN BUILD FORMAT STRINGS (VARIES W/ FLOAT)
% park header line starts with '$   Date        p       t      s'
% data header line starts with '$       p        t        s '
%
% HEX: O2,MCOMS, pH (I THINK?)
%     1       2        3         4         5        6      7      8   
% Pressure Temp(C) Salinity  nbins_pts  O2Phase  O2Temp nbins_O2  Fl  
%  9  10       11          12    13      14
% Ntu Cd  nbins_MCOMS(?)  phV   phT  nBins_pH{?)
%
% HEX: O2,MCOMS,(I THINK?)
%     1       2        3         4         5        6      7      8  
% Pressure Temp(C) Salinity  nbins_pts  O2Phase  O2Temp nbins_O2  Fl  
%  9  10     11
% Ntu Cd  nbins_MCOMS
% ************************************************************************               
fid = fopen(msg_file);
tline = ' ';
while ischar(tline) % # find header lines   
    if regexp(tline, '\$\s+Date\s+p\s+t\s+s\s+','once') % park header found
        pk_hdr = regexp(tline,'\s+', 'split');
        data.pk_hdr = pk_hdr(~strcmp('$', pk_hdr));
        pk_format  = '';
        for i =1:size(data.pk_hdr,2)
            if strcmp('Date', data.pk_hdr{i})
                pk_format = [pk_format,'%*s%s%s%s%s'];
            else
                pk_format = [pk_format,'%f'];
            end
        end
    end

    if regexp(tline, '\$\s+\s+p\s+t\s+s\s+','once') % data header found
        lr_hdr = regexp(tline,'\s+', 'split');
        data.lr_hdr = lr_hdr(~strcmp('$', lr_hdr));
        data.lr_hdr = data.lr_hdr';
        lr_format  = '';
        for i = 1:size(data.lr_hdr,1) % build format string
            lr_format = [lr_format,'%f'];
        end
        break  % data header found - leave loop
    end
    tline = fgetl(fid);
end
fclose(fid);

if isnumeric(tline) || isempty(data.lr_hdr)
    disp(['Incomplete message file! No profile data for ',msg_file])
    return
end

% OK NOW DETERMINE HEX FORMAT FOR CP DATA
% [> is neg number, scale factor]
hex_ind  = [1:3, 5:6, 8:10]; % col indices for PTS, O2 Phase&T, MCOMS 
nbin_ind = [4,7,11]; %nbins for pts, o2, mcoms
nbin_hdr = {'nbin ctd';'nbin oxygen';'nbin MCOMS'}; % putting bins at end
conv_ind = 1:3; % col indices for PTS
hex_conv = [32768 10; ...     % p  these are used to convert to pst
            61440 1000; ...   % t
            61440 1000];      % s
        
if sum(strcmp('phV', data.lr_hdr)) == 1 % O2, MCOMS, pH
    pH_flag = 1;
    hex_format = '%04x%04x%04x%02x%06x%06x%02x%06x%06x%06x%02x%06x%04x%02x';
    hex_length = 60;
    hex_cols   = size(regexp(hex_format,'x'),2);
    hex_ind    = [hex_ind, 12:13]; % PTS, O2 Phase&T, MCOMS, pH V&T 
    nbin_ind   = [nbin_ind, 14];
    nbin_hdr   = [nbin_hdr; 'nbin pH'];
    hex_conv   = [hex_conv; 61440 1000]; % Add for pH T
    conv_ind   = [conv_ind, 10]; % add col indice for pH T
else % O2, MCOMS
    pH_flag    = 0;
    hex_format = '%04x%04x%04x%02x%06x%06x%02x%06x%06x%06x%02x';
    hex_length = 48;
    hex_cols   = size(regexp(hex_format,'x'),2);
end

% BUILD HIGH RES HEADER (CP DATA HEADER)
tNO3 = strcmpi('no3', data.lr_hdr); 
data.hr_hdr = data.lr_hdr(~tNO3); % all but nitrate
data.hr_hdr = [data.hr_hdr;nbin_hdr]; % add bin counts to end
clear lr_hdr pk_hdr i tline fid ans tNO3

% ************************************************************************
% ************************************************************************
%                         PARSE MESSAGE FILE
% ************************************************************************
% ORDER = header, parkPt, termination time, cp data, gps fix, surf obs
% ************************************************************************
fid   = fopen(msg_file);
tline = fgetl(fid);% initialize with first line

% FILE EXISTS BUT NO DATA OR NON STANDARD FILE FORMAT & FILE ENDED
if tline == -1 % Go to next i
    disp('File exist but empty inside - moving to next message file.')
    fclose(fid);
    data = []; % function will return empty value if no msg data
    return;
end

% CAST # FROM MESSAGE FILE NAME
str      = regexp(msg_file,'\d{3}(?=\.msg)','match');
data.cast = str2double(str{1}); % cast # from file name
clear str

% ************************************************************************
% ************************************************************************
% PARSING MSG FILE
% ************************************************************************
% ************************************************************************
data_chk   = 0;
msg_task   = 'profile time';
low_res    = [];
high_res   = [];
CpActP     = [];
FwRev      = [];
FwRev_str  = '';
FlbbMode   = NaN;
pk_data    = [];

% STEP THROUGH PARK SAMPLES TO GET TO PROFILE TERMINATION TIME
% AND START OF LOW RES SAMPLES
while ischar(tline)
    tline     = strtrim(tline); % at top so while catches tline = -1
    
    % GET SOME INFO FOR ANNIE
    FwRev_ind = regexp(tline,'FwRev', 'once');
    if ~isempty(FwRev_ind) % Firmware version
        data.FwRev_str = regexp(tline,'(?<=FwRev.+)\w+\s\w+', ...
                        'once','match');
        data.FwRev     = str2double(regexp(tline,'(?<=FwRev.+)\d+', ...
                        'once','match'));
        clear FwRev_ind
    end
    
    CpAct_ind = regexp(tline,'\$ CpAct', 'once');
    if ~isempty(CpAct_ind) %contant profiling activation P ?
        data.CpActivationP = str2double(regexp(tline,'\d+','once','match'));
        clear CpAct_ind
    end
    
    FlbbMode_ind = regexp(tline,'\$ FlbbMode', 'once');
    if ~isempty(FlbbMode_ind) % flag may exist stating whether Flbb on
        data.FlbbMode = str2double(regexp(tline,'\d+','once','match'));
        clear FlbbMode_ind
    end
    
    % GET PARK DATA
    if regexp(tline,'^ParkObs\:', 'once')
        pk_tmp = textscan(tline, pk_format, 'CollectOutput',1);
        d_str = [pk_tmp{1}{1},' ',pk_tmp{1}{2},' ',pk_tmp{1}{3},' ', ...
                 pk_tmp{1}{4}];
        sdn = datenum(d_str,'mmm dd yyyy HH:MM:SS');
        pk_data = [pk_data; sdn pk_tmp{2}];
    end
    
    % GET PROFILE TERMINATION TIME
    if data_chk == 0 && ~isempty(regexp(tline,'^\$ Profile', 'once'))
        d_str = regexp(tline,'(?<=terminated\:\s+\w{3}).+','once','match');
        data.sdn = datenum(d_str,'mmm dd HH:MM:SS yyyy'); 
        data_chk = 1;
        clear d_str
    % LOOK FOR PRESSURE VALUE AT BEGINING OF LOW RES DATA LINE ^(\d+\.\d+)
    elseif data_chk == 1
        ind1 = regexp(tline,'^(\d+\.\d+)|^(-\d+\.\d+)','once');
        if ~isempty(ind1)
            msg_task   = 'profile data';
            break
        end
    elseif isnumeric(tline) % tline = -1, end of file w/o termination time
        disp(['No termination time found for ',file_name, ' NO DATA!'])
        data = [];
        return
    end
    tline = fgetl(fid);
end

% EXTRACT PROFILE DATA
while ischar(tline)
    tline = strtrim(tline);
    
    % GET SOME INFO FOR ANNIE - Older float, found in footer
    FwRev_ind = regexp(tline,'FwRev', 'once');
    if ~isempty(FwRev_ind) % Firmware version
        FwRev = str2double(regexp(tline,'(?<=FwRev.+)\d+','once','match'));
        clear FwRev_ind
    end
    switch msg_task % do different tasks based on position in msg file
        
        case 'profile data' % EXTRACT DISCRETE SAMPLE DATA
            if isempty(regexp(tline,'^#','once')) % # means end of low res
                ind1 = regexp(tline,'\(Park Sample\)$','once');
                if isempty(ind1) && ~isempty(tline)
                    tmp = textscan(tline,lr_format,1,...
                        'CollectOutput',1);                 
                    low_res = [low_res; tmp{1}]; % tmp1: p, t, s, etc
                end
            %elseif ~isempty(regexp(tline,'Sbe41cpSerNo','once')) %cp hdr
            elseif ~isempty(regexpi(tline,'^#.+sbe','once')) %cp hdr
                data.CTDtype = regexpi(tline,'sbe\w+(?=serno)', ...
                               'match', 'once');
                data.CTDsn = regexpi(tline,'(?<=serno[)\d+', ...
                               'match', 'once');                 
                msg_task = 'cp_data';  % finished low res, start high res
            else % No high resolution data found - look for a fix
                msg_task = 'GPS data';
            end
            
        case 'cp_data'     % EXTRACT CP SAMPLE DATA
            % find cp data line & make sure hex characters only
            if length(tline) == hex_length % check size
                if regexp(tline,['[A-F0-9]{',num2str(hex_length),'}'],...
                        'once'); % hex chars?
                    tmp = sscanf(tline,hex_format);
                    if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
                        high_res = [high_res; tmp']; % cp
                    end
                    if data_chk == 1
                        data_chk = 2; % In the cp data lines now
                    end
                end
                
            % DEAL WITH CP DATA BUG (2nd to last line)
            % BUILD UP SHORT DATA LINE WITH 0's and F's
            elseif hex_length == 60 && length(tline) == 52 % short data line
                if regexp(tline,'[A-F0-9]{52}','once'); % hex chars?
                    tmp = sscanf([tline,'FFFFFF00'],hex_format); %build up
                    high_res = [high_res; tmp']; 
                    
                    if data_chk == 1
                        data_chk = 2; % In the cp data lines now
                    end
                end
                    
            elseif data_chk == 2 && length(tline) ~= hex_length % end of cp
                msg_task = 'GPS data';
            end
            
        case 'GPS data'    % GET GPS POSITION & AIR OXYGEN VALUES
            if ~isempty(regexp(tline,'^Fix:','once'))
                gps = sscanf(tline,'%*s %f %f',2); %[lon; lat]
                if ~isempty(gps)
                    data.gps = [data.gps; gps'];
                end
            end
            if (regexp(tline,'^(SurfaceObs)','once'))
                ind = regexp(tline,'/'); % 2nd & 3rd mark O2 data region
                if size(ind,2) >= 3 % make sure data line complete
                    str = strtrim(tline(ind(2)+1:ind(3)-1)); %O2 data string
                    nct = length(regexp(str,'\d+\.')); % how many #'s in str
                    T_pos = ~isempty(regexp(str,'^\d+\.\d+C','once')); % T 1st
                    if nct == 3 % 4330 optode
                        if T_pos == 1 % 3 #'s, T at begining
                            tmp = sscanf(str,'%fC %f %f')'; % T TPh RPh
                        else % 3 #'s, T at end
                            tmp = sscanf(str,'%f %f %fC')'; % TPh RPh T
                            tmp = tmp([3,1,2]);             % T TPh RPh
                        end
                        data.air =[data.air; tmp];
                    elseif nct == 2 % 4330 optode or older optode
                        if T_pos == 1 % 2 #'s, T at begining
                            tmp = sscanf(str,'%fC %f')'; % temp, phase
                        else          % 2 #'s, T at end
                            tmp = sscanf(str,'%f %fC')'; % phase, temp
                            tmp = tmp([2,1]);            % temp, phase
                        end
                        data.air =[data.air; tmp];
                    else
                        disp('Could not parse surf obs')
                    end
                else
                    disp('Partial Surf Obs line - could not parse')
                end
            end
            
            if regexp(tline,'^<EOT>','once')
                data.EOT = data.EOT +1;
            end
    end
    tline = fgetl(fid); % Grab next text line
end
fclose(fid);

if isempty(data.gps)
    data.gps = [NaN NaN];
end

% ************************************************************************
%             CONVERT CP TO USEFULL NUMBERS AND/OR COUNTS
%    THIS IS A MODIFICATION AND CONDENSING OF Dan Quittman's CODE FOR:
% hextop.m, hextot.m and hextos.m functions
% for converting 16 bit 2's complement format to decimal
% It appears this is also a way to deal with a 12 bit A/D board in a
% 16 bit world as well as signed integers in a hex world
% ************************************************************************
if ~isempty(high_res)
    high_res = high_res(:,[hex_ind, nbin_ind]); % add bin ct cols to end
         
%     % BUILD TEST MATRIX: > CUTOFF, NEG #'s, = 1; < CUTOFF, POS #'s, = 0;
%     % = CUTOFF, BAD VALUE, = NAN
%     t_hi = high_res(:,1:3) - (t0 * hex_conv(:,1)') > 0; % neg numbers 
%     tmp(t_hi)  = 1; % neg values, test matrix = 1
%     t_low = high_res(:,1:3) - (t0 * hex_conv(:,1)') < 0; % neg numbers
%     tmp(t_low) = 0; % pos values, test matrix = 0
%     % VALUES NOT CAUGHT BY FLAGS REMAIN AS NAN's
    
    tmp     = high_res(:,conv_ind); % p, t, s OR p, t, s, pH t
    hex_var1 = ones(size(tmp(:,1))) * hex_conv(:,1)'; % matrix 
    hex_var2 = ones(size(tmp(:,1))) * hex_conv(:,2)'; % matrix 
    tNaN = tmp*NaN; % set up NaN flags
    tNaN(tmp - hex_var1 ~= 0) = 1; % No NaN set to 1
    tHi  = tmp - hex_var1 > 0;
    tLo  = tmp - hex_var1 < 0;
    
    high_res(:,conv_ind) = (tHi .*(tmp-65536)./ hex_var2 + ...
        tLo .* tmp ./ hex_var2).* tNaN;

    % NOW DO BIO-SENSORS
    rail_vals = high_res == 2^24-1; % Check out of bounds A/D
    high_res(rail_vals) = NaN;
    high_res(:,4)   = (high_res(:,4)/100000) - 10; % O2 Phase
    high_res(:,5)   = (high_res(:,5)/1000000) - 1; % O2 temperature volts
    high_res(:,6:8) = high_res(:,6:8) - 500; % MCOMS
    if pH_flag == 1
        high_res(:,9)   = (high_res(:,9)/1000000) - 2.5; % pH volts 
    end
end


clear tmp t0 t_hi t_low
% ************************************************************************
%  FILL OUTPUT STRUCTURE
% ************************************************************************
if isempty(low_res) && isempty(high_res)
    disp(['No data found in ',msg_file])
else
    if ~isempty(low_res)
        data.lr_d = low_res;
    else
        disp(['No low resolution data found for ',msg_file])
    end
    
    if ~isempty(high_res)
        data.hr_d  = high_res; % p t s + bio-sensors
    else
        disp(['No high resolution data found for ',msg_file])
    end
    
    if ~isempty(pk_data)
        data.pk_d = pk_data;
    else
        disp(['No park data found for ',msg_file])
    end
    
end


clearvars -except data



