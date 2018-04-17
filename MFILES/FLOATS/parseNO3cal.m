function cal = parseNO3cal(cal_file)

% PURPOSE: 
%   This function parses an APEX or NAVIS biochemical float nitrate
%   calibration files and returns a structure  of calibration data. 
%   Ultimately to be used in calculating the actual nitrate concentration.
%
% USAGE:
%	data = parseNO3cal(cal_file)
%
% INPUTS:
%	cal_file  = Path\calibration file name as a string
%
% OUTPUTS:
%       cal =   a structure of calibration info. Variable names change for
%               ISUS or SUNA.
%           cal.type       = SUNA or ISUS
%           cal.SN         = Serial number(s)

%           cal.CalTemp    = temperature the instrument was calibrated
%                            in the lab
%           cal.CalDateStr = lab calibration date string
%           cal.CalSDN     = serial date number of the lab calibration date
%           cal.WL         = Wavelength array
%           cal.ESW        = Extinction coefficients for seawater at WL's
%           cal.ENO3       = Extinction coefficients for nitrate at WL's
%           cal.TSWA       = Satlantic proprietary - not used
%           cal.EHS        = Extinction coefficients for bisulfide at WL's
%           cal.Ref        = Reference intensity through pure water at WL's
%
%           cal.pixel_base = Default is 1 (1-256), 0 (0-255)
%           cal.depth_lag  = sensor depth lag, m (early floats) 0 = default
%           cal.WL_offset  = Adjustable Br wavelength offset (default = 210)
%           cal.min_fit_WL = superseeds msg file if exists
%           cal.max_fit_WL = superseeds msg file if exists
%           cal.DC_flag    = Default is 1 (can change later)
%                            1 use DC in NO3 calc, 0 use SWDC in NO3 calc 
%           cal.pres_coef  = Bromide extinction coefficient ~ 0.02
%
% EXAMPLE:
%   jp = parseNO3cal('\\atlas\Chem\ISUS\Argo\5143StnP\5143.cal')
%
% Created 12/29/2015 by jp

%TEST
%cal_file ='C:\Users\jplant\Documents\MATLAB\ARGO\DATA\Cal_files\0508.cal'; % for testing

%cal_file ='C:\Users\jplant\Documents\MATLAB\ARGO\DATA\Cal_files\7663.cal'; % for testing

% ************************************************************************
% SET FORMATS, PREDIMENSION STRUCTURE, BOOKKEEPING
% ************************************************************************

% initialize output structures with NaN's. If parsing fails this will be
% the output
cal.type = 'ISUS';
cal.SN   = '';

cal.CalTemp    = [];
cal.CalSDN     = [];
cal.CalDateStr = '';
cal.DC_flag    = 1; % 1 DC, 0 SWDC

% SET OPTICAL WAVELENGTH OFFSET - MAY COME FROM CAL FILE IN LATER VERSIONS
% OFF CALIBRATION FIL
%cal.WL_offset = 210;

% SET SOME FLAGS, Default = 1, change manually later where needed

% cal.pixel_base = 1; % 1 (1-256), 0 (0-255)

% check for valid target
if exist(cal_file,'file')
    fid   = fopen(cal_file);
else
    disp(['Can not find calibration file at: ',cal_file])
    return
end

% FIND # OF CALIBRATION DATA ROWS
tline = ' ';
data_ct = 0;
while ischar(tline) % find number of wavelengths in cal file
    if  regexp(tline,'^E','once') % DATA LINES
        data_ct = data_ct + 1;
    end
    tline = fgetl(fid); % get 1st line
end

% check for data lines
if data_ct == 0
    disp(['File exists but no data lines found for : ',cal_file]);
    return
end

% ************************************************************************
% PARSE CALIBRATION FILE
% ************************************************************************
% GO TO BEGINING OF FILE - DEAL WITH HEADER INFO FIRST
frewind(fid);
tline = fgetl(fid); % prime the engine - get 1st line
ct = 0;
while ischar(tline)
    tline = strtrim(tline); % just in case
    % GET HEADER LINE INFO
    if strncmp(tline,'H',1) % header line so check for header info
        
        % GET INSTRUMENT TYPE AND SN(s)
        if  strncmp(tline,'H,SUNA',6) % check for Satlantic instrument
            cal.type = 'SUNA';
            % get #'s after you see SUNA
            cal.SN   = regexp(tline,'(?<=SUNA)\s+\d+', 'once','match');
        end
        if  regexp(tline,'Lamp#','once') % check for ISUS instrument
            tline = regexprep(tline,',','');
            cal.SN   = tline(2:end);          
        end
        
        % GET CALIBRATION DATE
        if  regexp(tline,'^H,\d+/\d+/\d{4}','once') % look for ISUS date
            d_str = regexp(tline,'\d+/\d+/\d{4}', 'match', 'once');
            cal.CalDateStr = d_str;
            cal.CalSDN = datenum(d_str,'mm/dd/yyyy');
        end
        if  regexp(tline,'creation time','once') % look for SUNA date
            d_str = regexp(tline,'\d{2}-\w{3}-\d{4}','match', 'once');
            cal.CalDateStr = d_str;
            cal.CalSDN = datenum(d_str,'dd-mmm-yyyy');
        end
        
        % GET CALIBRATION TEMPERATURE
        if  regexp(tline,'CalTemp','once') % ISUS or SUNA cal temp
            cal.CalTemp = sscanf(tline,'H,CalTemp,%f',1);
        end

        % GET PIXEL BASE, 1 (1-256), 0 (0-255)
        if  regexp(tline,'Pixel base','once') 
            cal.pixel_base = sscanf(tline,'H,Pixel base,%f',1);
        end
        
        % GET SENSOR DEPTH OFFSET (LAG) % may be more than 1 offset value
        if  regexp(tline,'Sensor Depth','once') 
            cal.depth_lag = sscanf(tline,'H,Sensor Depth offset,%f,%f');
        end
        
        % GET BROMIDE WAVELENGTH OFFSET % default is 210
        if  regexp(tline,'Br wavelength','once')
            cal.WL_offset = sscanf(tline,'H,Br wavelength offset,%f');
        end

        % CHECK FOR UPDATED FIT WINDOW
        if  regexp(tline,'H,Min fit','once') 
            cal.min_fit_WL = sscanf(tline,'H,Min fit wavelength,%f');
        end
        if  regexp(tline,'H,Max fit','once') 
            cal.max_fit_WL = sscanf(tline,'H,Max fit wavelength,%f');
        end 
        
        % CHECK IF SW DARK CURRENT IS USED
        if  regexp(tline,'Use seawater dark','once') 
            if regexpi(tline,'yes','once')
                cal.DC_flag = 0; % use sw DC
            end
        end  
        
        % GET BROMIDE PRESSURE COEFFICIENT ~ 0.02
        if  regexp(tline,'Pressure coef','once')
            cal.pres_coef = sscanf(tline,'H,Pressure coef,%f');
        end
        
        % PARSE HEADER LINE
        if  regexp(tline,'Wave','once') & regexp(tline,'NO3', 'once')
            % shorten hdr names and/or rename to make cals easier later
            tline = regexprep(tline,'New ',''); % shorten hdr names a bit
            tline = regexprep(tline,'DI DC Corr','Ref'); % ISUS
            tline = regexprep(tline,'Reference','Ref'); % SUNA
            tline = regexprep(tline,'Wavelength','WL'); % SUNA
            tline = regexprep(tline,'WaveLen','WL'); % ISUS
            tline = regexprep(tline,',NO3',',ENO3'); % SUNA
            tline = regexprep(tline,',SWA',',ESW'); % SUNA
            tline = regexprep(tline,',ASW,',',ESW,'); % 7622 APEX w SUNA
            tline = regexprep(tline,',T\*ASW',',TSWA'); %Satlantic ISUS 7663
            hdr = textscan(tline, '%s', 'Delimiter',',');
            hdr = hdr{1,1}; 
            hdr(1) = []; % Remove 1st cell 'H'
            cal_format = '%*s';
            for i = 1:length(hdr) % Build format string
                cal_format = [cal_format,'%f'];
            end
            data = ones(data_ct, length(hdr))* NaN; % predim
        end
    end
     
    % PARSE CALIBRATION DATA LINES   
    if  regexp(tline,'^E','once') % data lines begin with "E"
        ct = ct+1;
        d = textscan(tline, cal_format, 'Delimiter',',',...
            'CollectOutput',1);
        data(ct,:) = d{1,1};
    end
    tline = fgetl(fid);
end

% ASSIGN VARIABLE NAMES BASED ON HEADER
for i = 1:length(hdr)
    s1 = ['cal.',hdr{i,1},'=data(:,i);'];
    eval(s1)
end

fclose(fid);


clearvars -except cal
%clear  fid tline hdr d CT CT_format hdr_format cal_format