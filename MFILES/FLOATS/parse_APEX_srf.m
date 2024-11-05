function d = parse_APEX_srf(srf_file)
% PURPOSE:
%   This function parses an UW/MBARI APEX biochemical float "*.srf" file
%   and returns a structure of raw data. The "*.srf" files are the in air
%   O2 compliment file to the msg file for some APEX floats with two O2
%   sensors on board (ie ua19298,ua19843)
%
% USAGE:
%	d = parse_APEX_srf(srf_file)
%
% INPUTS:
%	srf_file  = file name or  path\file name as string
%
% OUTPUTS:
%   d =  a structure of the data and headers. 
%       d.hdr      = hdr - column ID's for d.data
%       d.data     = in air data
%       d.cast     = profile number
%       d.EOT      = flag count of EOT occurances
%       d.msg_file = parsed file name
%
% EXAMPLES:
%    d = parse_APEX_srf('C:\temp\19843.004.srf')
%    d = parse_APEX_srf('19843.004.srf')
%
% CHANGES:
%
% TESTING
% srf_file = 'C:\temp\19843.004.srf'; % 2XO2 float

% ************************************************************************
% FORMATS & VARIABLES
% ************************************************************************
% PREDIMENSION OUTPUT STRUCTURE (If no data this is the output)
d.hdr  = {}; % low res header cell array
d.data = []; % low res data matrix
d.cast = '';
d.EOT  = 0; % complete message file flag

% ************************************************************************
% GET HEADER ROW & COLUMN INDICES FOR FLOAT VARS
% ************************************************************************
if isempty(regexp(srf_file,'srf$','once'))
    disp(['WARNING: ', srf_file,' does not appear to be a "*.srf" file']);
    return
else
    fid        = fopen(srf_file);
    d.srf_file = srf_file;
    d.cast     = regexp(srf_file,'\d{3}(?=\.srf)','match','once');
end

tline   = '';
tf_hdr  = 0;
line_ct = 0;
while ischar(tline)
    % Find & parse header line - % header line includes 'Date'
    if tf_hdr == 0 && ~isempty(regexp(tline,'Date','once')) 
        d.hdr = regexp(tline,'[\w+\(\)]+','match');
        chdr  = size(d.hdr,2);
        data  = ones(30, chdr)* NaN;
        tf_hdr = 1;
    end
    % Find and parse data lines
    if tf_hdr == 1 && ~isempty(regexp(tline,'^SurfO2','once')) % data line found
        dstr  = regexp(tline,'(?<=SurfO2:\s+)[a-zA-Z]{3}.+\:\d{2}',...
                'match','once'); % pull date string
        dnum  = datenum(dstr,'mmm dd yyyy HH:MM:SS');
        ind   = regexp(tline,'\d{8}\d+','once'); % start of Unix epoch seconds
        nums  = regexp(tline(ind:end),'[\w+-\.]+', 'match'); % pull numbers or NaN str
        dline = [dnum, str2double(nums)];
        if size(dline,2) == chdr % All good
            line_ct = line_ct + 1;
            data(line_ct,:) = dline;
        else
            disp(['WARNING: Partial srf data line not parsed for ', ...
                srf_file]);
        end
    end
    
    if regexp(tline,'<EOT>','once')
        d.EOT = d.EOT+1;
    end
    
    tline = fgetl(fid);
end
fclose(fid);
d.data = data(1:line_ct,:);
clearvars -except d



