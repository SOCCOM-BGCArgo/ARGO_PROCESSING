function [LON,LAT] = get000msg_position(msg_file)

% Very simple function to retrieve the *.000.msg lat/lon position
% information.  Can be used as an estimate for a float's first profile
% position.
%
% Tanya Maurer
% MBARI
% 01/10/19
%
%
% ************************************************************************
% GET GPS FIX IN *000.MSG FILE
% ************************************************************************               

LON = [];
LAT = [];
GPS = [];

is000 = strfind(msg_file,'000.msg');
if isempty(is000)
    disp([msg_file,' IS NOT CYCLE 000.  USE FUNCTION parse_APEXmsg4ARGO.m or parse_NAVISmsg4ARGO.m'])
    return
end

EX = exist(msg_file,'file');
if EX ~= 2 %file does not exist
    disp([msg_file,' DOES NOT EXIST.  NOT ABLE TO RETRIEVE POSITION INFO.'])
    return
end

fid = fopen(msg_file);
tline = '';

while isempty(regexp(tline,'^Fix:','once')) % # find header line
    tline = fgetl(fid);
    if isnumeric(tline)
        disp(['Incomplete message file! No profile data for ',msg_file])
        fclose(fid);
        return
    end
end

gps = sscanf(tline,'%*s %f %f',2); %[lon; lat]
if ~isempty(gps)
    GPS = [GPS;gps];
else
    disp('NO POSITION DATA IN FILE.')
end

if GPS(1) < 0
    GPS(1) = GPS(1)+360;
end
    
LON = GPS(1);
LAT = GPS(2);
