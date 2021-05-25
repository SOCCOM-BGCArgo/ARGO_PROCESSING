function data = parse_APEXmsg4ARGO(msg_file)
% PURPOSE:
%   This function parses an UW/MBARI APEX biochemical float *.msg file and
%   returns a structure of raw data. The # of colums are determined from
%   the low resolution header line
%
% USAGE:
%	data = parse_APEXmsg(file_name)
%
% INPUTS:
%	msg_name  = string of file name or  path\file name
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
%           data.cast   = profile number
%           data.sdn    = profile termination date [matlab sdn]
%           data.gps    = gps location fix [lon, lat]
%           data.air    = Air oxygen measurements [Temp phase (Rphase)]
%           data.EOT    = flag count of EOT occurances
% EXAMPLES:
%   parse_APEXmsg4ARGO('c:\temp\7601.003.msg')
%   parse_APEXmsg4ARGO('7601.003.msg')
%
% CHANGES
%   01/18/2017 - cp header line search more generic & extracted CTD
%               type and serial number ultimately for ODV meta info
%   01/23/2017 - if no PTS remove line
%   04/06/2017 - add code to catch bad cp data lines, "00000000000000[2]",
%       in middle of cp section of msg file ~line 230
%   05/16/2017 - add code to look for <EOT> and record # of instances
%
%   09/28/2020 - Code was skipping over first line of optode air-cal
%   sequence.  Added a fix at the end of cp-mode hex extraction. (may be
%   revised at a later date).
%
%   10/26/20 - TM, Modified gps fix extraction to carry over datetime. Was
%              not included in original code.
%   10/28/20 - JP, Modified Surface Obs parseing code to include PAL sensor
%              float formats
%   04/28/21 - JP fixed 12712 issue processing cyles 61 & 111 with no
%              profile terminted time line in msg file. Get time from
%              SBE cp header line
%  05/13/21 - TM added check for hi-res pressure level of 0; exists in msg
%             file 12892.064.msg but bad data, should be removed (AOML screening for this) 

% ************************************************************************
% FORMATS & VARIABLES
% ************************************************************************
% FOR TESTING
% msg_file = '\\atlas\chemwebdata\floats\alternate\f9095\9095.011.msg';

%msg_file = 'c:\temp\9031.003.msg'; % TESTING
%msg_file = 'c:\temp\9662.003.msg'; % TESTING
%msg_file = 'c:\temp\7622.032.msg'; % TESTING
%msg_file = 'c:\temp\9752.002.msg';
%msg_file = 'c:\temp\6966.004.msg'; % TESTING
%msg_file = 'c:\temp\9031.006.msg';

% 10/28/2020 JP testing
%msg_file = 'c:\temp\12652.036.msg';
%msg_file = 'c:\temp\11090.003.msg';
%msg_file = 'c:\temp\7614.003.msg';
%msg_file = 'c:\temp\12878.006.msg'; % under ice surf obs values = 0
%msg_file = 'c:\temp\5143.002.msg'; % no surf obs line data.air = [];
%msg_file = 'c:\temp\6968.002.msg'; %T & Phase onle , no RPhase
%msg_file = 'C:\temp\12712.061.msg';
%msg_file = 'C:\temp\12712.016.msg';


% HIGH RESOLUTION SAMPLES FOR UW / MBARI APEX FLOATS
% The first 4-bytes of the encoded sample represents the pressure in
% centibars.  The second 4-bytes represents the temperature in
% millidegrees.  The third 4-bytes represent the salinity in parts per
% million.  The final 2-bytes represent the number of samples collected in
% the 2dbar pressure bin. 2's compliment


% PREDIMENSION OUTPUT STRUCTURE (If no data this is the output)
data.lr_hdr = {}; % low res header cell array
data.lr_d   = []; % low res data matrix
data.hr_hdr = {}; % high res header
data.hr_d   = []; % high res data
data.cast   = NaN; % profile #
data.sdn    = NaN; % profile termination time
data.gps    = []; % gps location fix
data.air    = []; % Air oxygen measurements
data.aircal = []; % air calibration measurements{ w/ & w/o air inflation
data.FwRev  = [];
data.CpActivationP = [];
data.FlbbMode = NaN;
data.CTDtype = '';
data.CTDsn   = '';
data.EOT     = 0; % complete message file flag

f_aircal = '%*s%s%s%s%s%*f%f%f%f%f%f';

% ************************************************************************
% GET HEADER ROW & COLUMN INDICES FOR FLOAT VARS
% ************************************************************************
fid = fopen(msg_file);

tline = '';
% header line starts with '$       p        t        s '
while isempty(regexp(tline,'p\s+t\s+s\s+','once')) % # find header line
    tline = fgetl(fid);
    if isnumeric(tline)
        disp(['Incomplete message file! No profile data for ',msg_file])
        fclose(fid);
        return
    end
end
fclose(fid);

ind1       = regexp(tline,'\$\s+p\s+t\s+s','once');
msg_hdr    = textscan(tline(ind1+1:end),'%s');
float_vars = msg_hdr{1,1}; % r x 1 cell array with header variables

% BUILD LOW RES DATA FORMAT STRING - VARIES W/ FLOAT
low_res_format  = '';
for i =1:length(float_vars)
    low_res_format = [low_res_format,'%f '];
end
clear ind1 i tline fid msg_hdr

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
alt_sdn_chk = 0;
msg_task   = 'profile time';
low_res    = [];
high_res   = [];
CpActP     = [];
FwRev      = [];
FlbbMode   = NaN;

% STEP THROUGH PARK SAMPLES TO GET TO PROFILE TERMINATION TIME
% AND START OF LOW RES SAMPLES
while ischar(tline)
    tline     = strtrim(tline); % at top so while catches tline = -1

    % GET SOME INFOR FOR ANNIE
    FwRev_ind = regexp(tline,'FwRev', 'once');
    if ~isempty(FwRev_ind) % Firmware version
        FwRev = str2double(regexp(tline,'(?<=FwRev.+)\d+','once','match'));
        clear FwRev_ind
    end
    CpAct_ind = regexp(tline,'\$ CpAct', 'once');
    if ~isempty(CpAct_ind) %contant profiling activation P ?
        CpActP = str2double(regexp(tline,'\d+','once','match'));
        clear CpAct_ind
    end
    
    FlbbMode_ind = regexp(tline,'\$ FlbbMode', 'once');
    if ~isempty(FlbbMode_ind) % flag may exist stating whether flbb on
        FlbbMode = str2double(regexp(tline,'\d+','once','match'));
        clear FlbbMode_ind
    end
    
    if data_chk == 0 && ~isempty(regexp(tline,'^\$ Profile', 'once'))
        msg_time_format = '%*s %*s %*s %*s %*s %s %s %s %s';
        d_str = textscan(tline,msg_time_format,'CollectOutput',1);
        d_str = d_str{1}; % cells: mmm dd hh:mm:ss 2008
        s1 = [d_str{2},'-',d_str{1},'-',d_str{4},' ',d_str{3}]; % = Malab format
        data.sdn = datenum(s1,0); % Profile end time as Matlab SDN
        %if ~isreal(SDN), pause, end % WHY DID I PUT THIS HERE?? -jp
        clear d_str s1
        data_chk = 1;
    
    % 04/28/2021 12712 61 & 111 strange msg files no "profile terminated line"
    % to get time from - try and  get time from SBE CP header line
    elseif data_chk == 0 && ~isempty(regexp(tline,'^\$ Discrete samples', 'once'))
        disp(['WARNING: Profile terminated line not found for ' ,msg_file]);
        disp('Attempt to estimate termination time from SBE CP header.')
        data_chk = 1;
        alt_sdn_chk = 1;
        
    % LOOK FOR PRESSURE VALUE AT BEGINING OF LOW RES DATA LINE ^(\d+\.\d+)
    elseif data_chk == 1
        ind1 = regexp(tline,'^(\d+\.\d+)|^(-\d+\.\d+)','once');
        if ~isempty(ind1)
            msg_task   = 'profile data';
            break
        end
    % 04/28/2021 JP - this test is never reached becuase the while statement
    % catches the -1 before the last elseif block
%     elseif isnumeric(tline) % tline = -1, end of file w/o termination time
%         disp(['No termination time found for ',file_name, ' NO DATA!'])
%         data = [];
%         return
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
                    tmp = textscan(tline,low_res_format,1,...
                        'CollectOutput',1);
                    low_res = [low_res; tmp{1}]; % tmp1: p, t, s, etc
                end
                %elseif ~isempty(regexp(tline,'Sbe41cpSerNo','once')) %cp hdr
            elseif ~isempty(regexpi(tline,'^#.+sbe','once')) %cp hdr
                if alt_sdn_chk == 1 % 04/28/2021 JP no teminated time, get from SBE hdr
                    str = regexp(tline,'(?<=\#\s+)[\w\s\:]+(?=\s+Sbe)', ...
                        'once','match');
                    data.sdn = datenum(str,'mmm dd yyyy HH:MM:SS');
                end
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
            if length(tline) == 14 % right % of chars
                high_res_format = '%04x%04x%04x%02x';
                ind1 = regexp(tline,'[A-Z0-9]{14}','once'); % hex chars?
                if ~isempty(ind1)
                    tmp = sscanf(tline,high_res_format);
                    high_res = [high_res; tmp']; % cp = p t s
                    
                    if data_chk == 1
                        data_chk = 2; % In the cp data lines now
                    end
                end
            elseif data_chk == 2 & regexp(tline,'0{14}\[\d+\]','once')
                disp(['No data hex line in middle of cp data: ',tline, ...
                    ' for ', msg_file])
            elseif data_chk == 2 && length(tline) ~= 14 % end of cp
                % TM 9/28/20 Code was skipping first line of aircal data.
                % Add OptodeAirCal block here after end of cpmode
                % extraction.  The remainder of the aircal array will get
                % appended further down.  (Perhaps a better way to do
                % this!)
                if (regexp(tline,'^OptodeAirCal','once'))
                    ac_tmp = textscan(tline,f_aircal,1,'collectoutput',1);
                    d_str  = [ac_tmp{1,1}{1},' ',ac_tmp{1,1}{2},' ', ...
                        ac_tmp{1,1}{3},' ',ac_tmp{1,1}{4}];
                    if ~isempty(ac_tmp{1,2})
                        data.aircal =[data.aircal; ...
                            datenum(d_str,'mmm dd yyyy HH:MM:SS'), ...
                            ac_tmp{1,2}];
                    end
                    clear ac_tmp d_str
                end
                msg_task = 'GPS data';
            end
            
        case 'GPS data'    % GET GPS POSITION & AIR OXYGEN VALUES
            if ~isempty(regexp(tline,'^Fix:','once'))
                gps = sscanf(tline,'%*s %f %f',2); %[lon; lat]
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
            
            % complete surf line is bounded by curly braces, count to check
            if ~isempty(regexp(tline,'^(SurfaceObs)','once')) && ...
                    size(regexp(tline,'{|}'),2) == 2
                
                % use slash as delimiter to break up string
                stmp = regexp(tline,'/','split');
                if size(stmp,2) == 4 % count substrings (in cell array now)
                    O2_str = strtrim(stmp{3});
                elseif size(stmp,2) == 2 % PAL SENSOR FLOATS!
                    O2_str = strtrim(stmp{2});
                    O2_str = regexprep(O2_str,'}',''); % get rid of "}"
                else
                    disp(['SurfaceObs line format not recognized. ', ...
                        'Could not parse line'])
                    O2_str = '';
                end
                
                %space charcter divides numbers, number ct = space count +1
                %Is "C" at end of string?
                num_ct  = sum(O2_str == ' ',2) + 1;
                tf_endC = ~isempty(regexp(O2_str,'C$', 'once'));
                if num_ct == 3 && tf_endC % 3 #'s & "C" at end of string
                    tmp = sscanf(O2_str,'%f %f %fC')'; % TPh RPh T
                    tmp = tmp([3,1,2]);             % T TPh RPh
                elseif num_ct == 3 % 3 #'s & "C" must be at begining
                    tmp = sscanf(O2_str,'%fC %f %f')'; % T TPh RPh
                elseif num_ct == 2 && tf_endC % 2 #'s & "C" at end of string
                    tmp = sscanf(O2_str,'%f %fC')'; % phase, temp
                    tmp = tmp([2,1]);            % temp, phase
                elseif num_ct == 2
                    tmp = sscanf(O2_str,'%fC %f')'; % temp, phase
                else
                    disp('Could not parse surf obs')
                end
                data.air =[data.air; tmp];
            end
            
            
            
            
            %                 if size(ind,2) >= 3 % make sure data line complete
            %                     str = strtrim(tline(ind(2)+1:ind(3)-1)); %O2 data string
            %                     nct = length(regexp(str,'\d+\.')); % how many #'s in str
            %                     T_pos = ~isempty(regexp(str,'^\d+\.\d+C|^-\d+\.\d+C','once')); % T 1st
            %                     if nct == 3 % 4330 optode
            %                         if T_pos == 1 % 3 #'s, T at begining
            %                             tmp = sscanf(str,'%fC %f %f')'; % T TPh RPh
            %                         else % 3 #'s, T at end
            %                             tmp = sscanf(str,'%f %f %fC')'; % TPh RPh T
            %                             tmp = tmp([3,1,2]);             % T TPh RPh
            %                         end
            %                         data.air =[data.air; tmp];
            %                     elseif nct == 2 % 4330 optode or older optode
            %                         if T_pos == 1 % 2 #'s, T at begining
            %                             tmp = sscanf(str,'%fC %f')'; % temp, phase
            %                         else          % 2 #'s, T at end
            %                             tmp = sscanf(str,'%f %fC')'; % phase, temp
            %                             tmp = tmp([2,1]);            % temp, phase
            %                         end
            %                         data.air =[data.air; tmp];
            %                     else
            %                         disp('Could not parse surf obs')
            %                     end
            %                 else
            %                     disp('Partial Surf Obs line - could not parse')
            %                 end
            %                 pause
            %             else
            %                 disp('Partial Surf Obs line was not parses')
            %             end
            
            if (regexp(tline,'^OptodeAirCal','once'))
                ac_tmp = textscan(tline,f_aircal,1,'collectoutput',1);
                d_str  = [ac_tmp{1,1}{1},' ',ac_tmp{1,1}{2},' ', ...
                    ac_tmp{1,1}{3},' ',ac_tmp{1,1}{4}];
                if ~isempty(ac_tmp{1,2})
                    data.aircal =[data.aircal; ...
                        datenum(d_str,'mmm dd yyyy HH:MM:SS'), ...
                        ac_tmp{1,2}];
                end
                clear ac_tmp d_str
            end
            
            % check for file termination, one per satelite comms/GPS cycle
            if regexp(tline,'^<EOT>','once')
                data.EOT = data.EOT +1;
            end
            
            
    end
    tline = fgetl(fid); % Grab next text line
end
fclose(fid);

if isempty(data.gps)
    data.gps = [data.sdn NaN NaN];
end

% ************************************************************************
% CONVERT CP TO USEFULL NUMBERS  --  THIS IS A MODIFICATION AND CONDENSING
% OF Dan Quittman's hextop.m, hextot.m and hextos.m functions
% for converting 16 bit 2's complement format to decimal
% It appears this is also a way to deal with a 12 bit A/D board in a
% 16 bit world as well as signed integers in a hex world
% ************************************************************************
if ~isempty(high_res)
    % [> is neg number, scale factor]
    hex_conv    = [32768 10; ...     % p  these are used to convert to pst
        61440 1000; ...   % t
        61440 1000];      % s
    
    tmp  = high_res(:,1:3) * NaN; % predim test array w NaN's
    t0   = ones(size(tmp(:,1))); % helper array for matrix building
    
    % BUILD TEST MATRIX: > CUTOFF, NEG #'s, = 1; < CUTOFF, POS #'s, = 0;
    % = CUTOFF, BAD VALUE, = NAN
    t_hi = high_res(:,1:3) - (t0 * hex_conv(:,1)') > 0; % neg numbers
    tmp(t_hi)  = 1; % neg values, test matrix = 1
    t_low = high_res(:,1:3) - (t0 * hex_conv(:,1)') < 0; % neg numbers
    tmp(t_low) = 0; % pos values, test matrix = 0
    % VALUES NOT CAUGHT BY FLAGS REMAIN AS NAN's
    
    high_res(:,1:3) = high_res(:,1:3) - (t0*[65536 65536 65536]) .* tmp;
    high_res(:,1:3) = high_res(:,1:3) ./ (t0 * hex_conv(:,2)');
end

clear tmp t0 t_hi t_low
% ************************************************************************
%  FILL OUTPUT STRUCTURE
% ************************************************************************
if isempty(low_res) && isempty(high_res)
    disp(['No data found in ',msg_file])
else
    data.FwRev         = FwRev;
    data.CpActivationP = CpActP;
    data.FlbbMode      = FlbbMode;
    
    if ~isempty(low_res)
        tnan = all(isnan(low_res(:,1:3)),2); % MISSING PT&S?
        low_res(tnan,:) = []; % If no PTS remove line
        data.lr_hdr = float_vars;
        data.lr_d = low_res;
    else
        disp(['No low resolution data found for ',msg_file])
    end
    
    if ~isempty(high_res)
        %check for any zero pressure level - remove (ie 12892.064.msg; AOML
        %screens for this)
        xx = find(high_res(:,1)==0);
        high_res(xx,:) = [];
        data.hr_d  = high_res(:,1:4); % p t s nbin
        data.hr_hdr = [float_vars(1:3);'nbin ctd'];
    else
        disp(['No high resolution data found for ',msg_file])
    end
end
clearvars -except data



