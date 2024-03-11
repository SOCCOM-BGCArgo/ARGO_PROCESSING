function DrkStats = get_CHLdark(cal, dirs, range)
%
% get_CHLdark.m returns insitu minimum CHL counts for the defined profile
% range for profiles deeper than 900 meters
%
% INPUTS:
%   cal       = calibration data structure containing MBARI_ID & WMO, &
%               cal data
%   dirs      = Either an empty array or a structure with directory
%                  strings where files are located. If dirs is empty
%                  default paths will be used.
%   range     = if empty all profiles are inspected, if scalar all
%                   profiles up to number, if size = 2 a profile range
%
% OUTPUTS:
%   Darkstats   = a structure of insitu dark information for each
%       .DC     median DC counts
%       .N      # of profiles used
%       .std    standard deviation
%       .data   [Cycle, Pres, DC]
%
% CHANGES
%    07/12/23 - JP - Updates to include CHL435 - kind of a rewrite of the
%       code using loops and dynamic structures. If additional CHL channels
%       there will be a suite of additional fields. The entire cal
%       structure is now an input vs cal.info since the cal was being
%       loaded inside anyway.
% ************************************************************************
% TESTING
% cal_fn = sprintf('cal%s.mat','ua21910');
% load(['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\',cal_fn])
% range        = 5;
% dirs         = [];

% ************************************************************************
MBARI_ID_str = cal.info.name;
WMO          = cal.info.WMO_ID;
DrkStats     = [];
max_depth    = 900;
chl_cals     = {'CHL' 'CHL435'}; % All possible chl channels

if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\'];
    
    dirs.cal       = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
    dirs.FV        = [user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ\'];
end

% CHECK DATA FILE
d = get_FloatViz_data([dirs.FV, WMO, '.TXT']);
if isempty(d)
    disp(['FILE NOT FOUND: ',dirs.FV, WMO, '.TXT'])
    return
end
d.data(d.data == -1e10) = NaN; % MVI to NaNs

% CHECK DATA
% GET FLOATVIZ DATA INDICES
I.STA    = find(strcmp('Station',d.hdr)         == 1);
I.P      = find(strcmp('Pressure[dbar]',d.hdr)  == 1);
I.T      = find(strcmp('Temperature[°C]',d.hdr) == 1);
I.S      = find(strcmp('Salinity[pss]',d.hdr)   == 1);
I.CHL    = find(strcmp('Chl_a[mg/m^3]',d.hdr)   == 1);
I.CHL435 = find(strcmp('Chl_a435[mg/m^3]',d.hdr)   == 1);

% CHECK RANGE VARIABLE DIMENSIONS
rsz = max(size(range));
if rsz == 2 % block of profiles - non standard operatrions
    t1 = d.data(:,I.STA) >= range(1) & d.data(:,I.STA) <= range(2);
    d.data = d.data(t1,:);

    if isempty(d.data)
        disp('File exists but no cycles in selected range!')
        return
    end
end

for sct = 1:size(chl_cals,2)
    Cstr      = chl_cals{sct};
    Cstr_flag = sprintf('%s_flag', lower(Cstr)); % flag in info struct
    if ~isfield(cal,Cstr)
        continue
    elseif isfield(cal,Cstr) && cal.info.(Cstr_flag) == 0
        fprintf(['WARNING: Cal structure exists for %s but sensor ',...
            'info flag contradicts - check calibration!!!\n'],Cstr);
        return
    end

    t_good = ~isnan(d.data(:,I.(Cstr))) & d.data(:,I.(Cstr)+1) ~= 8;
    d.data = d.data(t_good,:);
    clear t_good

    if isempty(d.data)
        disp('File exists but no good CHL data found')
        continue
    end
    
    % Convert back to counts
    d.data(:,I.(Cstr)) = d.data(:,I.(Cstr)) ./ cal.(Cstr).ChlScale ...
        + cal.(Cstr).ChlDC;

    
    % FIND SHALLOWEST MINIMUM FOR EACH PROFILE
    cycles = unique(d.data(:,I.STA)); % get unique cycle count
    rr     = size(cycles,1);
    data   = ones(rr,3)*NaN;
    ct     = 0;

    for i = 1:rr
        t1  = d.data(:,I.STA) == cycles(i);
        tmp = d.data(t1,:);

        [~,IX] = sort(tmp(:,I.P)); % Sort shallow to deep
        tmp    = tmp(IX,:);

        if max(tmp(:,I.P)) >= max_depth % full depth profile?
            ct = ct+1;
            [~,ind] = min(tmp(:,I.(Cstr)));
            %ind = find(tmp(:,I.(Cstr)) == min(tmp(:,I.(Cstr))),1,'first');
            data(ct,:) = [cycles(i), tmp(ind,I.P), tmp(ind,I.(Cstr))];
        end
    end
    data(isnan(data(:,3)),:) =[]; % remove cycle that didn't meet max P

    % Subset to profile count limit
    if rsz == 1 && size(data,1) >= range(1)
        data = data(1:range(1),:); % [cycle, shallowest P, min ct]
    elseif rsz == 1
        disp('Not enough profiles to get seawater dark counts');
        return
    end
 
    % Fill output structure
    nstr = regexp(Cstr,'\d+','match','once'); % strip #, empty for CHL 
    DrkStats.(['DC',nstr])   = median(data(:,3));
    DrkStats.(['N',nstr])    = size(data,1);
    DrkStats.(['std',nstr])  = std(data(:,3));
    DrkStats.(['data',nstr]) = data;
end

clearvars -except DrkStats

    







