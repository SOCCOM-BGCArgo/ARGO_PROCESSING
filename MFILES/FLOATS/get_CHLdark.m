function DrkStats = get_CHLdark(MBARI_ID_str, dirs, range)
%
% get_CHLdark.m returns insitu minimum CHL counts for the defined profile
% range for profiles deeper than 900 meters
%
% INPUTS:
%   MBARI_ID_str  = MBARI Float ID as a string
%   dirs          = Either an empty array or a structure with directory
%                   strings where files are located. If dirs is empty
%                   default paths will be used.
%   range         = if empty all profiles are inspected, if scalar all
%                   profiles up to number, if size = 2 a profile range
%
% OUTPUTS:
%   Darkstats   = a structure of insitu dark information
%       .DC     median DC counts
%       .N      # of profiles used
%       .std    standard deviation
%       .data   [Cycle, Pres, DC]

% ************************************************************************
% TESTING

% MBARI_ID_str = '9762SOOCN';
% MBARI_ID_str = '0507SOOCN';
% range = [5];
% dirs =[];
% ************************************************************************

DrkStats =[];
max_depth = 900;

if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\'];
    
    dirs.cal       = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
    dirs.FV        = [user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ\'];
end

% CHECK DATA FILE
d = get_FloatViz_data([dirs.FV, MBARI_ID_str, '.TXT']);
if isempty(d)
    disp(['FILE NOT FOUND: ',dirs.FV, MBARI_ID_str, '.TXT'])
    return
end
d.data(d.data == -1e10) = NaN; % MVI to NaNs

% CHECK DATA
% GET FLOATVIZ DATA INDICES
iSTA = find(strcmp('Station',d.hdr)         == 1);
iP   = find(strcmp('Pressure[dbar]',d.hdr)  == 1);
iT   = find(strcmp('Temperature[°C]',d.hdr) == 1);
iS   = find(strcmp('Salinity[pss]',d.hdr)   == 1);
iC   = find(strcmp('Chl_a[mg/m^3]',d.hdr)   == 1);

t_good = ~isnan(d.data(:,iC)) & d.data(:,iC+1) ~= 8;
d.data = d.data(t_good,:);
clear t_good

% CHECK RANGE VARIABLE DIMENSIONS
rsz = max(size(range));
if rsz == 2 % block of profiles 
    t1 = d.data(:,iSTA) >= range(1) & d.data(:,iSTA) <= range(2);
    d.data = d.data(t1,:);
end

if isempty(d.data)
    disp('File exists but no good CHL data found')
    return
end

% CHECK CAL FILE
float_cal_path = [dirs.cal,'cal',MBARI_ID_str,'.mat'];
if exist(float_cal_path,'file')
    load(float_cal_path);
    if isfield(cal,'CHL') % Convert back to counts
        d.data(:,iC) = d.data(:,iC) ./ cal.CHL.ChlScale + cal.CHL.ChlDC;
    else
        disp(['No calibration coefficents found for: ',float_cal_path])
        return
    end  
else
    disp(['Cal file not found: ',dirs.cal,'cal',MBARI_ID_str,'.mat']);
    return
end

cycles = unique(d.data(:,iSTA));
rr     = size(cycles,1);
data   = ones(rr,3)*NaN;
ct     = 0;

% FIND SHALLOWEST MINIMUM FOR EACH PROFILE
for i = 1:rr
    t1  = d.data(:,iSTA) == cycles(i);
    tmp = d.data(t1,:);
    
    [~,IX] = sort(tmp(:,iP)); % SORT
    tmp    = tmp(IX,:);
    
    if max(tmp(:,iP)) >= max_depth % full depth profile
        ct = ct+1;
        ind = find(tmp(:,iC) == min(tmp(:,iC)),1,'first');
        data(ct,:) = [cycles(i), tmp(ind,iP), tmp(ind,iC)];
    end
end
data(isnan(data(:,3)),:) =[];

% ************************************************************************
% NOW FILL OUTPUT
if rsz == 1 && size(data,1) >= range(1)
    data = data(1:range(1),:); % [cycle, shallowest P, min ct]
elseif rsz == 1
    disp('Not enough profiles to get seawater dark counts');
    return
end

DrkStats.DC   = median(data(:,3));
DrkStats.N    = size(data,1);
DrkStats.std  = std(data(:,3));
DrkStats.data = data;

%clearvars -except DrkStats

    







