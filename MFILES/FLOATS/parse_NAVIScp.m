function cp_data = parse_NAVIScp(cp_file, msg_data)
% PURPOSE: 
%   This function parses an NAVIS biochemical float *.cp file and
%   returns a structure of raw data. The # of colums are determined from
%   the low resolution header line
%
% USAGE:
%	cp_data = parse_NAVIScp(cp_file, msg_data)
%
% INPUTS:
%	cp_file  = string of file name or path\file name
%
% OUTPUTS:
%       cp data =  a matrix HR data derived from hex formated data
%       msg_data = output structure from parse_NAVISmsg4ARGO.m
%
% EXAMPLES:
%   d  = parse_NAVISmsg4ARGO('c:\temp\1516.001') % get msg file info 1st
%   d2 = parse_NAVIScp('c:\temp\1516.001.cp', d) %

%
% CHANGES:
% 12/13/23 JP initial code built

% TESTING
% msg_data  = parse_NAVISmsg4ARGO('C:\temp\1516.001.msg'); % get msg file info 1st
% cp_file   = 'c:\temp\1516.001.cp';

cp_data = [];
if ~isfile(cp_file)
    fprintf(['cp file (%s) not found - no HR data will be extracted ', ...
        'for this cycle\n'], cp_file);
    return
end

% bring in hex parsing info from msg file data structrue
pH_flag     = msg_data.hr_info.pH_flag;
hex_format  = msg_data.hr_info.hex_format;
hex_char_ct = msg_data.hr_info.hex_char_ct;
hex_length  = msg_data.hr_info.hex_length;
hex_cols    = msg_data.hr_info.hex_cols;
hex_ind     = msg_data.hr_info.hex_ind;
nbin_ind    = msg_data.hr_info.nbin_ind;
nbin_hdr    = msg_data.hr_info.nbin_hdr;
conv_ind    = msg_data.hr_info.conv_ind;
hex_conv    = msg_data.hr_info.hex_conv;

% PARSE HEX DATA & CONVERT TO DECIMAL -LINE BY LINE BECAUSE LINES CAN BE
% VARIABLE SIZE
fid = fopen(cp_file);
tmp = textscan(fid,'%s','Delimiter','\r\n','CollectOutput',1);
fclose(fid);
hex_data = tmp{1,1};
tf       = cellfun(@isempty,regexp(hex_data,'^\#|^ser','once'));
hex_data = hex_data(tf);
rhex     = size(hex_data,1);
high_res = ones(rhex, hex_cols) * NaN; % predim
clear tmp fid ans

for row_ct = 1:rhex
    hex_line = hex_data{row_ct};

    if size(hex_line,2) == hex_length && ...
            size(regexp(hex_line,'^[A-F0-9]+','match','once'),2) == hex_length
        tmp = sscanf(hex_line, hex_format);
        if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
            high_res(row_ct,:) = tmp'; % cp
        end

        %SIZE ~= HEX_LENGTH & ALL HEX CHARS, BUT PARTIAL OR SHORT LINE
    elseif size(hex_line,2) > 11 && ... % 12 chars means p t & s (14 for bin count too)
            size(hex_line,2) == size(regexp(hex_line,'^[A-F0-9]+','match','once'),2)
        % figure out # of valid measurements
        ind = find(hex_char_ct <= size(hex_line,2),1,'last'); % jp 7/2020

        tmp = sscanf(hex_line, hex_format);
        if sum(tmp(1:3)) ~= 0 % IS THERE P,S,T or EMPTY LINE?
            dtmp = ones(hex_cols,1)* NaN; % predim to accept variable size input jp 7/2020
            dtmp(1:ind) = tmp(1:ind); % jp 7/2020
            high_res(row_ct,:) = dtmp'; % cp
        end
    end
end
tg       = ~isnan(high_res(:,1));
high_res = high_res(tg,:);

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
    
    tmp     = high_res(:,conv_ind); % p, t, s OR p, t, s, pH t OR p, t, s, pH
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

    cp_data = high_res;
end

clearvars -except cp_data


