function spec = parse_NO3msg(isus_file)

% PURPOSE: 
%   This function parses an APEX or NAVIS biochemical float nitrate *.isus
%   file. A structure is returned containing the seawater UV intensity
%   spectrum and other data used to determine the nitrate concentration
%
% USAGE:
%	data = parse_NO3msg(isus_file)
%
% INPUTS:
%	isus_file  = Path\calibration file name as a string
%
% OUTPUTS:
%   spec =   a structure of data and calibration info
%       spec.SDN        = sample time, Matlab sdn        (array)
%      	spec.P          = sample pressure, dbar          (array)
%     	spec.T          = sample temperature, C          (array)
%     	spec.S          = sample salinity, pss           (array)
%      	spec.DC         = sample DC intensity, counts    (scalar)
%     	spec.UV_INTEN   = Measured UV intensities        (matrix)
%      	spec.SWDC       = sea water DC intensity, counts (scalar)
%                         NaN if doesn't exist (earlier floats)
%
%      	spec.spectra_pix_range = pixel registration for the wavelengths in
%                                the calibration file
%       spec.WL_Fit_Range      = wavelength bounds of fit window - used in
%                                SUNA floats to determine spectra_pix_range
%    	spec.pix_fit_win       = pixel range for Multiple linear regression.
%                                Used to subset sample spectra to fit window
%
%     	spec.CalTemp    = calibration temperature this is the temperature
%                         the instrument was calibrated at in the lab.
%    	spec.CalDateStr = lab calibration date string
%     	spec.CalDate    = serial date number for the lab cal date
%
%
% EXAMPLE:
%   jp = parse_NO3msg('\\atlas\chemwebdata\floats\f5143\5143.005.isus')
%
% Created 01/13/2016 by jp
%
% REVISONS:
% 02/28/17 Added code to look for null data values in Salinity and remove
% all rows from the pertinent data fields
% 08/11/2020 - JP - added code to determine MSC firmware versions and save
%                   it as a field in the ouput structure. Also delt with
%                   bogus time stamp issue (18110). If time stamp bogus
%                   subsitute EOP date

% isus_file ='\\atlas\chemwebdata\floats\f5143\5143.005.isus'; % for testing
% isus_file ='\\atlas\Chem\ISUS\Argo\8497Hawaii\8497.042.isus';
% isus_file ='C:\temp\9095.068.isus';
% isus_file ='C:\temp\6403.113.isus';
% isus_file ='C:\temp\6967.225.isus';
% isus_file ='C:\temp\6976.152.isus';
% isus_file = 'C:\temp\9666.003.isus';
% ************************************************************************
% PARSE *.isus FILE
% ************************************************************************
% FORMATS TO PARSE BASE 10 PART OF *.ISUS MESSAGE FILE
% HEX WILL BE DONE ON THE FLY BELOW AS HEX BLOCK CAN CHANGE SIZE
% FORMATS TO PARSE BASE 10 PART OF *.ISUS MESSAGE FILE
% HEX WILL BE DONE ON THE FLY BELOW AS HEX BLOCK CAN CHANGE SIZE
d_format = '%s %s %s %f %f %f %f %f %f %f'; % six cols, 1st is string
d_format = [d_format,'%f %f %f %f %f %f %f %f %f %f'];
d_format = [d_format,'%f %f %f %f %s %f'];  

%Predim output if no spectra present
spec.SDN        = [];
spec.P          = [];
spec.T          = [];
spec.S          = [];
spec.DC         = [];
spec.UV_INTEN   = [];
spec.SWDC       = [];

spec.CalTemp           = NaN;
spec.CalDate           = NaN;
spec.CalDateStr        = '';
spec.pix_fit_win       = NaN;
spec.spectra_pix_range = NaN;
spec.file_name         = '';
spec.float             = '';
spec.cast              = '';
spec.EOT               = 0; % complete data file flag

spec.MSC_FW_ver = 0;
spec.EOP_sdn    = NaN; 
spec.bogus_time = 0;

% check for valid target
if exist(isus_file,'file')
    fid   = fopen(isus_file);
    spec.file_name = regexp(isus_file,'\d+\.\d+\.\D+','once','match');
    if ~isempty(spec.file_name)
        flt_info = regexp(spec.file_name,'\.','split');
        spec.float = flt_info{1};
        spec.cast  = flt_info{2};
        clear flt_info
    end
else
    disp(['Can not find *.isus file at: ',isus_file])
    return
end

% ************************************************************************
% FIND # OF SAMPLE SPECTRA DATA ROWS & HEX BLOCK LENGTHS
% FOR SAMPLE COMPLETENESS
% ************************************************************************
tline   = ' ';
data_ct = 0;
hdr_ct  = 0;
hdr_chk  = 0;
hex_len = [];
MSC_FW_ver = 0;

% % MSC firmware version App build dates as of  08/12/2020
msc_fw_versions(1,:)  = {'Jun 11 2007', 0};
msc_fw_versions(2,:)  = {'May 28 2008', 0};
msc_fw_versions(3,:)  = {'Feb  6 2009', 0};
msc_fw_versions(4,:)  = {'May  9 2011', 0};
msc_fw_versions(5,:)  = {'Aug 24 2011', 0};
msc_fw_versions(6,:)  = {'Feb 24 2012', 0};
msc_fw_versions(7,:)  = {'Aug 23 2012', 0};
msc_fw_versions(8,:)  = {'Apr  7 2015', 1};
msc_fw_versions(9,:)  = {'Sep 30 2015', 2};
msc_fw_versions(10,:) = {'Apr 30 2019', 4};


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
            spec.EOP_sdn = datenum(dstr,'mmm dd HH:MM:SS yyyy');
        end
    end
    
    if  regexp(tline,'^0x','once') % DATA LINES
        hdr_chk =1;
        data_ct = data_ct + 1;
        % GET HEX STRING FROM TLINE
        hex_str = regexp(tline,'\w{150,}','once','match'); % >= 150 char
        if ~isempty(hex_str)
            hex_len(data_ct) = length(hex_str);
        else
            hex_len(data_ct) = 0;
        end
    end
    if regexp(tline,'^H,','once') & hdr_chk == 0% DATA LINES
        hdr_ct = hdr_ct+1;
    end
    tline = fgetl(fid); % get 1st line
end

if isnan(spec.EOP_sdn)
    disp(['Could not find EOP line in ', isus_file,'. file is ', ...
            'probably incomplete'])
end


% Check for header lines
if hdr_ct ==  0
    disp(['File exists but no header lines found in file for: ', ...
        isus_file]);
    fclose(fid);
    return
end
        
% check for data lines
if data_ct == 0
    disp(['File exists but no data lines found for : ',isus_file]);
    fclose(fid);
    return
else
    hex_length      = mode(hex_len); % most common value
    spec_length     = hex_length /4;
    line_test       = hex_len ~= hex_length; % 1's are bad lines
    data_ct         = sum(~line_test); % Sum good lines
    
    spec.SDN        = ones(data_ct,1)*NaN; % predim
    spec.P          = spec.SDN;
    spec.T          = spec.SDN;
    spec.S          = spec.SDN;
    spec.DC         = spec.SDN;
    spec.SWDC       = spec.SDN;
    spec.UV_INTEN   = ones(data_ct,spec_length)*NaN; % predim

    %hex_length   = length(dline{1,25}{1})/4;
    hex_template = '%04x';
    hex_format   = hex_template;
    % BUILD HEX FORMAT STRING % all hex #'s = 4 char
    for i = 1:spec_length  - 1 % build hex format string
        hex_format = [hex_format,hex_template];
    end
end

% ************************************************************************
% GO TO BEGINING OF FILE AND START PARSING THE DATA FILE
% ************************************************************************
frewind(fid); 
tline = fgetl(fid); % prime the engine - get 1st line 
line_ct = 0;
hdr_chk = 0;
bogus_time = 0;

while ischar(tline)
    % GET HEADER LINE INFO
    if  hdr_chk == 0 && ~isempty(regexp(tline,'^H','once'))
        hdr_ct = hdr_ct +1;
        if regexp(tline,'Calibration\s+Date','once') % calibration date
            % use to match to cal file calibration date
            d_str = textscan(tline, '%*s %*s %s', 'Delimiter',',',...
                'CollectOutput',1); 
            spec.CalDateStr = d_str{1,1}{1};
            spec.CalDate = datenum(d_str{1,1},'mm/dd/yyyy');
        end        
        if regexp(tline,'Sw\s+Calibration\s+Temp','once') % cal temp
            % use to match to cal file calibration temp
            spec.CalTemp = str2double(regexp(tline,'\d+\.\d+','match')); 
        end
        if regexp(tline,'Wavelength Fit','once') % Get WL MLR (SUNA)
            % This can be used to create pix_fit_win for SUNA files
            spec.WL_fit_win = str2double(regexp(tline,'\d+\.*\d*','match')); 
        end
        if regexp(tline,'Pixel Fit','once') % Get pixel range for MLR 
            % This will be used to subset sample spectra to fit window
            % for MLR - eventually have user control to adjust
            spec.pix_fit_win = str2double(regexp(tline,'\d+','match')); 
        end
    end
    


    % PARSE DATA LINES
    if  ~isempty(regexp(tline,'^0x','once')) % data line
        dline = textscan(tline,d_format,'Delimiter',',');
        if ~isempty(dline{1,25})
            % CHECK SIZE OF HEX STRING, IF PARTIAL SPECTRA MOVE TO NEXT LINE
            % ALSO CHECK FOR NULL VALUE IN SALINITY, SKIP IF NULL
            if length(dline{1,25}{1}) == hex_length %&& dline{1,7} ~= -1
                line_ct = line_ct +1;
                
                %spec.SDN(line_ct)  = datenum(dline{1,3},'mm/dd/yyyy HH:MM:SS');
                
                % GET TIME STAMP
                if regexp(dline{1,3}{1}, '01/01/2000','once')
                    bogus_time = 1;
                    spec.SDN(line_ct) = spec.EOP_sdn;
                else
                    spec.SDN(line_ct) = datenum(dline{1,3},'mm/dd/yyyy HH:MM:SS');
                end
                
                
                spec.P(line_ct)    = dline{1,5};
                spec.T(line_ct)    = dline{1,6};
                spec.S(line_ct)    = dline{1,7};
                spec.DC(line_ct)   = dline{1,18}; % dark current
                
                % EARLIER FLOATS DID NOT RETURN SW DC AFTER HEX BLOCK so SET
                % SW DC TO DC IF NEED BE
                if isempty(dline{1,26})
                    spec.SWDC(line_ct) = spec.DC(line_ct); % Older model
                else % used if shutter stuck open
                    spec.SWDC(line_ct) = dline{1,26};
                end
                
                hex_UV_INTEN       = dline{1,25}{1}; % hex string still
                UV_INTEN = (sscanf(hex_UV_INTEN, hex_format))'; %col to row too
                spec.UV_INTEN(line_ct,1:spec_length) = UV_INTEN;
                if line_ct == 1 % get some one time info
                    hdr_chk = 1; % header before and after only need once
                    
                    % this range is the pixel registration for WL's in cal file
                    spec.spectra_pix_range = [dline{23},dline{24}];
                end
            end
        end
    end
    
    if regexp(tline,'^<EOT>','once')
        spec.EOT = spec.EOT +1;
    end
    
    tline = fgetl(fid);
end
fclose(fid);

if bogus_time == 1
    disp([isus_file,' has bogus time stamps - subsituting EOP ',...
        'date from APF11'])
end

% CHECK FOR NULL SALINITY VALUES (-1) - MEANS BAD CTD DATA
% IF THERE, REMOVE LINES FROM ALL PERTINENT FILEDS
t1 = spec.S == -1; % could return 0's, 1's, or empty
if sum(t1) > 0
    rr = size(t1,1); % # data rows
    fns = fieldnames(spec);
    for i = 1:size(fns,1)
        if size(spec.(fns{i}),1) == rr
            spec.(fns{i})(t1,:) =[];
        end
    end
end

spec.MSC_FW_ver = MSC_FW_ver;
%spec.EOP_sdn    = EOP_sdn;
spec.bogus_time = bogus_time;
            
clearvars -except spec









