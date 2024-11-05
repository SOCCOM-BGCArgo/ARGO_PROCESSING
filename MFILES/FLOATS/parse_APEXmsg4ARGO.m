function data = parse_APEXmsg4ARGO(msg_file) 

% PURPOSE:
%   This function parses an UW/MBARI APEX biochemical float *.msg file and
%   returns a structure of raw data. The # of columns are determined from
%   the low resolution header line
%
% USAGE:
%	data = parse_APEXmsg(file_name)
%
% INPUTS:
%	msg_name  = string of file name or  path\file name
%
% OUTPUTS:
%    data =  a structure of the data and headers. Size varies depending
%           on data flag.
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
%           data.ice_flag = UnderIceEvasion y/n
%           data.ice_hex = UnderIceEvasion status for the previous 8 cycles
%           data.gps    = gps location fix [lon, lat]
%           data.irid   = iridium location fix (newer floats)
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
%   05/13/21 - TM added check for hi-res pressure level of 0; exists in msg
%             file 12892.064.msg but bad data, should be removed (AOML screening for this)
%   06/03/21 - TM added code in support of SBE83 optode on APEX test float (air-cal format spec)
%
%   02/07/22 - EC added parser section for ParkPt and ParkPtFLBB
%   02/19/22 - EC added IceEvasionRecord in 8-digit and current status flag
%   04/27/22 - TM; incorporated the iridium fix positioning into parser (newer floats).
%   05/16/23 - JP; FL2BB parsing added (park & profile), 2 fluor channels ch0 = 470nm, ch1 = 435nm
%              also added code to capture park pressure. A little clean up
%              to processing to reduce overhead & spped up a bit
%              (hopefully)
%   05/16/23 - TM; modifications for 6-sensor APEX with OCR.  OCR data
%                  exists in low-res data an optode-air-cal samples.
%   06/08/23 - JP - Minor fixes to clear bugs incurred from JP's previuos
%              updates
%   06/14/23 - TM; small bug fix to the single air-measurement extraction
%   section (formatting irregularities on the SBS83 test floats!!)
%
%   04/23/24 - TM; small update to P=0 exclusion line at end of function,
%               based on correspondence with AOML.  They are now excluding only lines
%               where ALL data returns are zero (not just P=0...).  Case in
%               point is 20329.010.msg HR first sample.
%   08/19/2024 TM; Modification to aircal parsing to include flbb data.

%----------------
% FOR TESTING
% msg_file = '\\atlas\chemwebdata\floats\alternate\f9095\9095.011.msg';

% Full PkPt
% msg_file = 'C:\Users\eclark\Documents\Matlab\ARGO_PROCESSING\MFILES\EKC commented tests\12878.054.msg'

% Full ParkPtFLBB
% wmo 5906471
% msg_file = 'C:\Users\eclark\Documents\Matlab\ARGO_PROCESSING\MFILES\EKC commented tests\20532.003.msg'
%msg_file = 'c:\temp\9031.003.msg'; % TESTING

% 10/28/2020 JP testing
% msg_file = 'c:\temp\12652.036.msg';

% 07/05/2022 9095 cycle 125 gps fix bug testing
% msg_file = 'c:\temp\9095.125.msg';

% 05/16/2023 20974 FL2BB msg file parser testing
%msg_file = 'c:\temp\20974.006.msg'; % FL2BB test msg
%msg_file = 'c:\temp\9630.004.msg';
%msg_file = 'c:\temp\9630.000.msg';
%msg_file = 'c:\temp\20974.010.msg'; % 6 sensor test msg w OCR
%msg_file = 'c:\temp\19806.040.msg';
%msg_file = 'c:\temp\12688.074.msg';
% msg_file = 'c:\temp\18169.001.msg';

%----------------




% ************************************************************************
% FORMATS & VARIABLES
% ************************************************************************
% HIGH RESOLUTION SAMPLES FOR UW / MBARI APEX FLOATS
% The first 4-bytes of the encoded sample represents the pressure in
% centibars.  The second 4-bytes represents the temperature in
% millidegrees.  The third 4-bytes represent the salinity in parts per
% million.  The final 2-bytes represent the number of samples collected in
% the 2dbar pressure bin. 2's compliment


% PREDIMENSION OUTPUT STRUCTURE
data.pk_hdr = {}; % park data header cell array
data.pk_d   = []; % park data matrix
data.lr_hdr = {}; % low res header cell array
data.lr_d   = []; % low res data matrix
data.hr_hdr = {}; % high res header
data.hr_d   = []; % high res data
data.cast   = NaN; % profile #
data.sdn    = NaN; % profile termination time
data.sdni   = NaN; %iridium timestamp
data.gps    = []; % gps location fix
data.irid   = []; % iridium location fix (only some floats)
data.air    = []; % Air oxygen measurements
data.aircal = []; % air calibration measurements{ w/ & w/o air inflation
data.FwRev  = [];
data.FwMod  = []; %Apf9 or Apf11?
data.CpActivationP = [];
data.FlbbMode = NaN;
data.CTDtype  = '';
data.CTDsn    = '';
data.EOT      = 0; % complete message file flag
data.ice_flag = [];
data.ice_hex  = []; 
data.ParkPressure  = [];

f_aircal         = '%*s%s%s%s%s%*f%f%f%f%f%f';
f83_aircal       = '%*s%s%s%s%s%*f%f%f%f%f';
% sixsensor_aircal = '%*s%s%s%s%s%*f%f%f%f%f%f%8x%8x%8x%8x';
sixsensor_aircal = '%*s%s%s%s%s%*f%f%f%f%f%f%s%s%s%8x%8x%8x%8x'; %need to include flbb!?
high_res_format  = '%4x%4x%4x%2x';

pk_ct = 0; % park sample counter
lr_ct = 0; % low res sample counter

% ************************************************************************
% GET HEADER ROW & COLUMN INDICES FOR FLOAT VARS
% ************************************************************************
fid   = fopen(msg_file);
tline = '';
% header line starts with '$       p        t        s '
% find out if the file is complete/ whether or not to proceed
while isempty(regexp(tline,'p\s+t\s+s\s+','once')) % # find header line
    tline = fgetl(fid);
    if isnumeric(tline) % end of file & fgetl returns -1
        disp(['Incomplete message file! No profile data for ',msg_file])
        fclose(fid);
        return
    end
end

% split by white spece blocks
float_vars = regexp(tline,'\s+','split');
float_vars = float_vars(2:end)'; % loose the dollar sign capture

tf_OCR     = strncmp(float_vars,'Ocr',3); % logical size of float_vars
OCR_mode   = sum(tf_OCR) >0;

% BUILD LOW RES DATA FORMAT STRING - VARIES W/ FLOAT
low_res_format = [repmat('%f', 1, size(float_vars(~tf_OCR),1)), ...
    repmat('%8x', 1, size(float_vars(tf_OCR),1))];

% ************************************************************************
% ************************************************************************
%                         PARSE MESSAGE FILE
% ************************************************************************
% ORDER = header, parkPt, termination time, cp data, gps fix, surf obs
% ************************************************************************
frewind(fid); % go back to top of file
tline = ' '; % reinitialize

% FILE EXISTS BUT NO DATA OR NON STANDARD FILE FORMAT & FILE ENDED
% *** I don't think this ever happens now jp 05/17/23 ***
% if tline == -1 % Go to next i
%     disp('File exists but empty inside - moving to next message file.')
%     fclose(fid);
%     data = []; % function will return empty value if no msg data
%     return;
% end

% CAST # FROM MESSAGE FILE NAME
data.cast = str2double(regexp(msg_file,'\d{3}(?=\.msg)','match','once'));

% ************************************************************************
% ************************************************************************
% PARSING MSG FILE
% ************************************************************************
% ************************************************************************
data_chk    = 0;
alt_sdn_chk = 0;
msg_task    = 'profile time';
low_res     = [];
high_res    = [];
CpActP      = [];
FwRev       = [];
FwMod      = [];
FlbbMode    = NaN;
pk_data     = [];
ice_flag    = 0;
ice_hex     = [];
pk_pres     = [];% $ ParkPressure(1000) [dbar]

pk_ct   = 0;
pk_chk  = 0;

% col ct from format line & APEX should never have more than 70 LR samples
low_res = nan(100, size(regexp(low_res_format,'\%'),2)); 

% *************************************************************
% STEP THROUGH PARK SAMPLES TO GET TO PROFILE TERMINATION TIME
% AND START OF LOW RES SAMPLES

while ischar(tline)                                 % while each line starts w/ a character
    tline     = strtrim(tline);                      % at top so while catches tline = -1
 
    % GET SOME INFO FOR ANNIE
    FwRev_ind = regexp(tline,'FwRev', 'once');
    if ~isempty(FwRev_ind) % Firmware version
        FwRev = str2double(regexp(tline,'(?<=FwRev.+)\d+','once','match'));
        clear FwRev_ind
    end
    CpAct_ind = regexp(tline,'\$ CpAct', 'once');
    if ~isempty(CpAct_ind) %constant profiling activation P ?
        CpActP = str2double(regexp(tline,'\d+','once','match'));
        clear CpAct_ind
    end
    % TM, added 6/20/23
    FwMod_ind = regexp(tline,'Apf', 'once');
    if ~isempty(FwMod_ind) % Firmware Model config
        FwMod = regexp(tline,'(?<=Apf)\w+','once','match');
        FwMod = ['Apf',FwMod];
        clear FwMod_ind
    end

    FlbbMode_ind = regexp(tline,'\$ FlbbMode', 'once');
    if ~isempty(FlbbMode_ind) % flag may exist stating whether flbb on
        FlbbMode = str2double(regexp(tline,'\d+','once','match'));
        clear FlbbMode_ind
    end

    if ~isempty(regexp(tline,'\$\sParkPressure','once')) % find pkPressure
        data.ParkPressure = str2double(regexp(tline,'\d+','once','match'));
    end

    % ********************************************************************
    % ********************************************************************
    
    % GET PARK DATA
    % PARKPT
    if regexp(tline,'^ParkPt', 'once') % a park data line of some flavor
        pk_cell = regexp(tline,'\S+', 'match'); % chunk up ParkPt line by char groups (Big \S)

        % Do some one time tests. There will always be a LR profile header
        % row if you get here. Earlier floats had FLBB profile data but not park data
        % Figure out park data flavor
        if pk_chk == 0
            pk_chk = 1;
            if regexp(tline,'^ParkPt\:', 'once')
                pk_hdr        = {'Date' 'p' 't'}; % T & S only
                valid_pk_cols = 9;
                keep_idx      = 8:9;
            elseif regexp(tline,'^ParkPtFlbb\:', 'once') & sum(strcmp(float_vars,'FSig')) == 1
                pk_hdr        = {'Date' 'p' 't', 'Fsig' 'Bbsig' 'Tsig'};
                valid_pk_cols = 12;
                keep_idx      = 8:12;
            elseif regexp(tline,'^ParkPtFlbb\:', 'once') & sum(strcmp(float_vars,'FSig[0]')) == 1
                pk_hdr        = {'Date' 'p' 't' 'FSig[0]' 'FSig[1]' 'BbSig'  'TSig'};
                valid_pk_cols = 13;
                keep_idx      = 8:13;
            end
            pk_data = ones(300,size(pk_hdr,2))*NaN; % predim output
        end

        if size(pk_cell,2) ~= valid_pk_cols % valid col lengths defined above loop
            disp(['Uncharacteristic number of ParkPt columns check ',...
                'for incomplete msg file'])
            continue % If not a valid park line you can go to the next line
        end

        d_str = strtrim(sprintf('%s ',pk_cell{2:5}));
        try
            sdn = datenum(d_str,'mmm dd yyyy HH:MM:SS');       % format date
        catch
            sdn = nan;
        end

        pk_ct = pk_ct+1;
        pk_data(pk_ct,:) = [sdn, str2double(pk_cell(keep_idx))];
    end
    
    % ********************************************************************
	% Below block turns on a discrete profile block switch looking for a couple flavors
    % ********************************************************************
    
    % this set of if statements is CONFUSING ME
    if data_chk == 0 && ~isempty(regexp(tline,'^\$ Profile', 'once'))
        
        %msg_time_format = Profile_#_terminated:_Day_Month_#_HH:MM:SS_year
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
    end
    tline = fgetl(fid);
end

% EXTRACT PROFILE DATA
while ischar(tline)
    tline = strtrim(tline);  
    
      % ICE INFO        
    IceEvas_ind = regexp(tline,'IceEvasionRecord','once');
    if ~isempty(IceEvas_ind) 
        
        ihex = regexp(tline,'(?<=IceEvasionRecord.+)\w+','once','match'); % extract hex
        ice_hex = hex2bin(ihex, 8);                     % convert hex to binary
        lastnum = str2double(ice_hex(8));          % extract last digit of hex 8
        if lastnum == 1
            ice_flag = 1               ;                          % change flag to 1 if under ice now
        end
                clear IceEvas_ind
    end
    
    
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
                    %lr_ct = lr_ct +1;
                    tmp = sscanf(tline,low_res_format)';
                    % jp 06/01/23 must have data & data line must have
                    % PRES
                    if ~isempty(tmp) && ~isnan(tmp(1))
                        lr_ct = lr_ct +1;
                        low_res(lr_ct,1:length(tmp)) = tmp;
                    elseif ~isempty(tmp) && isnan(tmp(1)) 
                        fprintf(['Discrete data line has no PRES value ',...
                            'for %s and will be excluded\n'],msg_file);
                    end
                    %low_res(lr_ct,1:length(tmp)) = tmp;
                end
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
                ind1 = regexp(tline,'[A-Z0-9]{14}','once'); % hex chars?
                if ~isempty(ind1)
                    tmp = sscanf(tline,high_res_format);
                    high_res = [high_res; tmp']; % cp = p t s
                    
                    if data_chk == 1
                        data_chk = 2; % In the cp data lines now
                    end
                end
            elseif data_chk == 2 & regexp(tline,'0{14}\[\d+\]','once')
%                  comment out this elseif/display when running through fleet to save time/power  
                disp(['No data hex line in middle of cp data: ',tline, ...
                    ' for ', msg_file])
            elseif data_chk == 2 && length(tline) ~= 14 % end of cp
                % TM 9/28/20 Code was skipping first line of aircal data.
                % Add OptodeAirCal block here after end of cpmode
                % extraction.  The remainder of the aircal array will get
                % appended further down.  (Perhaps a better way to do
                % this!)
                if (regexp(tline,'^OptodeAirCal','once'))
                    if OCR_mode
                        f_aircal = sixsensor_aircal; %sixsensor APEX with OCR!
                    end
                    ac_tmp = textscan(tline,f_aircal,1,'collectoutput',1);
                    d_str  = [ac_tmp{1,1}{1},' ',ac_tmp{1,1}{2},' ', ...
                        ac_tmp{1,1}{3},' ',ac_tmp{1,1}{4}];
                    if ~isempty(ac_tmp{1,2})
%                         if OCR_mode
%                             data.aircal =[data.aircal; ...
%                                 datenum(d_str,'mmm dd yyyy HH:MM:SS'), ...
%                                 ac_tmp{1,2}, double(ac_tmp{1,3})];
%                         else
                            data.aircal =[data.aircal; ...
                                datenum(d_str,'mmm dd yyyy HH:MM:SS'), ...
                                ac_tmp{1,2}];
%                         end
                    end
                    clear ac_tmp d_str
                end
                if (regexp(tline,'^Sbe83AirCal','once'))
                    ac_tmp = textscan(tline,f83_aircal,1,'collectoutput',1);
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

            if ~isempty(regexp(tline,'^IridiumFix:','once'))
                irid = sscanf(tline,'%*s %*s %f %f',2); %[lon; lat]
                if ~isempty(irid)
                    gps_sdni = char(sscanf(tline,'%*s %*s %*f %*f %*s %*s %19c',1))'; %[mm/dd/yyyy; hh:mm:ss]
                    if size(gps_sdni,2) == 19
                        Isdn = datenum(gps_sdni,'mm/dd/yyyy HH:MM:SS');
                        data.irid = [data.irid; [Isdn irid']];
                    else
                        data.irid = [data.irid; [data.sdni NaN NaN]];
                    end
                else
                    data.irid = [data.irid; [data.sdni NaN NaN]];
                end
            end			
            
            % complete surf line is bounded by curly braces, count to check
            if ~isempty(regexp(tline,'^(SurfaceObs)','once')) && ...
                    size(regexp(tline,'{|}'),2) >= 2 % why is this hear?

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
                    
                % This is a patch for surface obs typo in msg file, 19746
                % find numberC in middle
                
                elseif num_ct == 3 & regexp(O2_str,'\s[\d\.]+C\s', 'once')
                    tmp = sscanf(O2_str,'%f %fC %f')'; % T TPh RPh
                    
                elseif num_ct == 3 % 3 #'s & "C" must be at begining
                    tmp = sscanf(O2_str,'%fC %f %f')'; % T TPh RPh
                elseif num_ct == 2 && tf_endC % 2 #'s & "C" at end of string; This is for the 2 APEX test floats (19065, 19727) with SBE83 -- they had the 'C' placed at the end next to the phase (erroneously).  These are anomalies. TM 6/14/23
                    tmp = sscanf(O2_str,'%f %fC')'; % phase, temp
                    if ~strcmp(FwMod,'Apf11')
                        tmp = tmp([2,1]);            % temp, phase % TM 6/20/23, for older 3830 with surfaceobs!  These are Apf9i.  ie 7620.
                    end 
                elseif num_ct == 2
                    tmp = sscanf(O2_str,'%fC %f')'; % temp, phase %This is the configuration for the 5-sensor APEX with SBS83 (17465)
                else
                    disp('Could not parse surf obs')
                end
                data.air =[data.air; tmp];
            end
            
            
            if (regexp(tline,'^OptodeAirCal','once'))
                if OCR_mode
                    f_aircal = sixsensor_aircal; %sixsensor APEX with OCR!
                end
                ac_tmp = textscan(tline,f_aircal,1,'collectoutput',1);
                d_str  = [ac_tmp{1,1}{1},' ',ac_tmp{1,1}{2},' ', ...
                    ac_tmp{1,1}{3},' ',ac_tmp{1,1}{4}];
                if ~isempty(ac_tmp{1,2}) % TM 8/19/24, I'm not sure why we need to carry along the OCR data in the aircal??  For float 20126 data was getting parsed incorrectly (flbb data missing?!).  Do all APEX with OCR have flbb?  Need to double check...
%                     if OCR_mode
%                         data.aircal =[data.aircal; ...
%                             datenum(d_str,'mmm dd yyyy HH:MM:SS'), ...
%                             ac_tmp{1,2}, double(ac_tmp{1,3})];
%                     else
                        data.aircal =[data.aircal; ...
                            datenum(d_str,'mmm dd yyyy HH:MM:SS'), ...
                            ac_tmp{1,2}];
%                     end
                end
                clear ac_tmp d_str
            end
            
            if (regexp(tline,'^Sbe83AirCal','once'))
                ac_tmp = textscan(tline,f83_aircal,1,'collectoutput',1);
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
% CONVERT CP TO USEFUL NUMBERS  --  THIS IS A MODIFICATION AND CONDENSING
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
if pk_ct > 0 % Park data exist
    data.pk_hdr  = pk_hdr;
    data.pk_d    = pk_data(1:pk_ct,:);
else
    fprintf('No park data were detected for this file: %s\n', msg_file);
    data.pk_d =[];
end


if isempty(low_res) && isempty(high_res)
    disp(['No data found in ',msg_file])
else
    data.FwRev         = FwRev;
    data.FwMod         = FwMod;
    data.CpActivationP = CpActP;
    data.FlbbMode      = FlbbMode;
    data.OCRMode       = OCR_mode;
    data.ice_flag      = ice_flag;                 % NEW
    data.ice_hex       = ice_hex;  
    data.ice_hex       = ice_hex; % NEW
    
    if ~all(isnan(low_res(:)))
    %if ~isempty(low_res)
        tnan            = all(isnan(low_res(:,1:3)),2); % MISSING PT&S?
        low_res(tnan,:) = []; % If no PTS remove line
        data.lr_hdr     = float_vars;
        data.lr_d       = low_res(1:lr_ct,:);
    else
        disp(['No low resolution data found for ',msg_file])
    end
    
    if ~isempty(high_res)
        %check for any zero pressure level - remove (ie 12892.064.msg; AOML
        %screens for this)
        if strfind(msg_file,'12892')
            xx = find(high_res(:,1)==0);
            high_res(xx,:) = [];
        end
%       % update 4/23/24; modify logic to only exclude lines where all PTS
%       are zero (per AOML).  This won't be backwards compatible with
%       12892...so hard-code exception for now...
        xx = find(all(high_res == 0,2));
        high_res(xx,:) = [];
        data.hr_d  = high_res(:,1:4); % p t s nbin
        data.hr_hdr = [float_vars(1:3);'nbin ctd'];
    else
        disp(['No high resolution data found for ',msg_file])
    end
end


clearvars -except data 

