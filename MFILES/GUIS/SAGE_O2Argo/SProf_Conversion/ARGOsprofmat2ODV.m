function tf_odv = ARGOsprofmat2ODV(f_info, outputDIR)

% ************************************************************************
% PURPOSE: 
%    This function creates ODV compatitble text files used in SAGE
%    using the *Sprof NetCDF files. If pH data exists the
%    LIAR approach will be used to estimate alkalinity.  Sprof2mat.m is
%    called as as subroutine.
%
%
% USAGE EXAMPLE:
%           f_info.WMO        = '6901580';
%           f_info.fn         = '6901580_Sprof.nc';
%           f_info.dac_path   = '/ifremer/argo/etc/argo-synthetic-profile/';
%           f_info.local_path = 'C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\DATA\ARGO_REPO\6901580\';
%           f_info.dac = 'coriolis';
%           tf_odv = ARGOsprofmat2ODV(f_info, f_info.local_path)
%
% INPUTS:
%   f_info: a structure containing the following:
%       f_info.WMO = WMO number (string) of float of interest 
%       f_info.fn = Sprof file name of interest (to convert)
%       f_info.dac_path = path to GDAC where Sprof file lives
%       f_info.local_path = Local path where Sprof file lives
%       f_info.dac = dac that 'owns' the float.
%
%   outputDIR =  A string defining where the
%               ODV ascii file will be written. 

%
% CHANGE HISTORY
% 09/27/2017 - created by JP
% 10/16/2018 - modified for Sprof files
% 04/13/19  - TM, modified fx routines to make running easier for
% users. Also removed calls to MBARI network locations etc...
% 04/15/2019 - fixed QF flagging bugs, added comment line describing Argo
%              QF to ODV QF conversion, added missing position interpolation,  Also renamed to "ARGOsprofmat2ODV.m" to distinguish from MBARI internal version.
% 3/25/21 - TM - Forced all fopen writes to UTF-8, because that is the
%    new default for Matlab 2020 and better cross platform sharing
% TESTING
%dirs =[]

%
d = ARGOSprof2mat(f_info);
info = d.INFO;
rdata = d.data; % row data from Sprof netcdf
rhdr  = d.hdr;
clear d

[raw_r,raw_c] = size(rdata);


% ************************************************************************
% DO SOME PREP WORK
% ************************************************************************
% SET UP LOCAL DATA PATH
%site_flag = 0;
fp = filesep;
mytempfolder = [getenv('HOMEDRIVE'),fp,'temp',fp]; % for my computer homedrive = C:
if ~exist(mytempfolder, 'dir')
  mkdir(mytempfolder);
  disp([mytempfolder,' does not exist. Creating directory...'])
end
param_list_file = 'argo-parameters-list-core-and-b.txt';

% **********************************************************

% ************************************************************************
% SET UP DIRECTORIES AND PATHS
% ************************************************************************
% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
user_dir = getenv('USERPROFILE');
dirs.temp = mytempfolder;

% **********************************************************


tf_odv = 0;
fill_0 = ones(raw_r,1) * 0; % for adding QF arrays


% ************************************************************************
% ADD LAT QF COL 
% QC FLAGS STILL ARGO AT THIS POINT
% INTERP QF FLAG = 3, POSITIONS STILL MISSING = NaN, QF = 99, GOOD = 1
% LOOK FOR MISSING PROFILES TOO
% ************************************************************************
iLAT     = find(strcmp('Lat',rhdr) ==1);
tnan_lat = isnan(rdata(:,iLAT));
QC_lat   = fill_0 +1; % position fixes now  = 1
QC_lat(tnan_lat) = 99;

rhdr    = [rhdr(1:iLAT),'Lat_QC', rhdr(iLAT+1:raw_c)]; % raw data
rdata   = [rdata(:,1:iLAT),QC_lat , rdata(:,iLAT+1:raw_c)]; % ARGO QF's,
[raw_r,raw_c] = size(rdata);

[~,ia,~] = unique(rdata(:,2)); % unique casts
pos_fix  = rdata(ia,1:6);   % WMO cast, profile, SDN, LON, LAT (position subset)


% ************************************************************************
% SCAN FOR MISSING FLOAT PROFILES
missing_profile_str = 'No missing float profiles';
missing_profiles =[];
for i = 1 : max(pos_fix(:,2))
    t1 = sum(pos_fix(:,2) == i,1);
    if t1 == 0
        missing_profiles = [missing_profiles, i];
    end
end
if ~isempty(missing_profiles)
    missing_profile_str = ['Missing Float profile(s) for station(s): ', ...
        sprintf('%0.0f ',missing_profiles)];
    disp(missing_profile_str);
end
clear i t1 missing_profiles

% ************************************************************************
% SCAN FOR MISSING POSITION FIXES & INTERPOLATE IF POSSIBLE
missing_pos_str = '';
interp_pos_str  = '';
t_nan = isnan(pos_fix(:,6)); % any nan's in LAT? 

%  if not all nans, try to interpolate
% (ie if never a pos fix then no interpolation, move on)
if sum(t_nan,1) > 0 && sum(t_nan,1)~=length(pos_fix(:,6)) 
    diff_nan   = [0;diff(t_nan)]; % 1 = start of NaN, -1 = start of good
    nan_start  = find(diff_nan == 1);
    nan_end    = find(diff_nan == -1) - 1;
    
    % some interp may be possible - resize to complete bounds only
    % # of start indices should equal # of end indices
    if size(nan_start,1) > size(nan_end,1) % more starts than ends, no last profile
        nan_start = nan_start(1:size(nan_end,1));     
    elseif size(nan_start,1) < size(nan_end,1) % more ends than starts
        nan_end = nan_end(2:end); % 1st profile doesn't have a position
    elseif nan_end(1)<nan_start(1) && nan_start(end)>nan_end(end) % same number of starts, ends, but don't line up cuz missing positions at both beginning and end
        nan_end = nan_end(2:end);
        nan_start = nan_start(1:end-1);
    end
    
    for i = 1 : size(nan_start,1)
        if isempty(nan_end) % no upper bound, can't interpolate, move on
            break
        end
        
        % all good now, nan's bounded at this point, procede to interp  
        if size(nan_start,1) == size(nan_end,1) 
            bnds = pos_fix([nan_start(i)-1,nan_end(i)+1],:); % bounds
            bait = pos_fix(nan_start(i):nan_end(i),:); % data to interp
            % Meridian crossing - float unlikely to move more than 180 deg
            if abs(bnds(2,5) - bnds(1,5)) > 180 
                t1 = bnds(:,5) < 180;
                bnds(t1,5) = bnds(t1,5) + 360; % temp add 360 to small side
            end
            missing_pos = [bait(:,1:4), ...
                           interp1(bnds(:,4),bnds(:,5:6),bait(:,4))];
            % bring meridian crossing back to reality if need be           
            t1 = missing_pos(:,5) > 360;
            missing_pos(t1,5) = missing_pos(t1,5) - 360;
            pos_fix(nan_start(i):nan_end(i),5:6) = missing_pos(:,5:6);
            for j = 1:size(missing_pos,1) % step through casts
                t1 = rdata(:,2) == missing_pos(j,2);
                rdata(t1,4:6) = ones(sum(t1),1) * missing_pos(j,4:6); % matrix
                rdata(t1,7)   = 3; %set QF to 3, ODV coversion takes it to 4
            end          
        end
    end
    clear diff_nan nan_start nan_end bnds bait i j t1
   
    % NOW ASSES WHAT HAS BEEN INTERPOLATED AND WHAT HAS NOT
    t_nan2 = isnan(pos_fix(:,6)); % Still any nan's in LAT? NO INTERP YET
    t_nan3 = t_nan & ~t_nan2; % used to be nan's but now interpolated
    if sum(t_nan3 > 0) % Interpolated values
        interp_pos_str = ['Missing Float position interpolated for ', ...
            'station(s): ', sprintf('%0.0f ',pos_fix(t_nan3,2))];
        interp_pos_str = [interp_pos_str,'\r\n//Latitude quality flag = 4',...
            ' for interpolated float positions'];
        disp(interp_pos_str)      
    end
    if sum(t_nan2 > 0) % Could not interpolate these casts
        missing_pos_str = ['No position for station(s): ', ...
            sprintf('%0.0f ',pos_fix(t_nan2, 2))];
        disp(missing_pos_str)
        rdata(isnan(rdata(:,4)),7) = 99;  % BIO ARGO MISSING DATA QF VALUE

    end  
    clear t_nan t_nan2 t_nan3 
end
% ************************************************************************
% **********************************************************************																		  
% NOW ADD SOME DATA COLUMNS AND HEADERS NEEDED FOR THE ODV FILE
% PRES QC, TEMP QC, PSAL QC, density,  O2sat, pH25 
% ************************************************************************
% ************************************************************************
% GET SOME INDICES (will be the same for raw and adjusted data)
iLAT  = find(strncmp('Lat',rhdr,5) == 1); % set order =lat, lat QF, p, t, s  
iP  = find(strcmp('PRES',rhdr)     == 1); 
iT  = find(strcmp('TEMP',rhdr)     == 1); 
iS  = find(strcmp('PSAL',rhdr)     == 1);

PRES_QF = fill_0 + 1;             % ARGO QF = GOOD
PRES_QF(isnan(rdata(:,iP))) = 99; % BIO ARGO QF MISSING VALUE

potT   = theta(rdata(:,iP),rdata(:,iT),rdata(:,iS),0);
den    = density(rdata(:,iS),potT) -1000; % pot den anom (sigma-theta)
den_QF = fill_0 + 1;
den_QF(isnan(den)) = 99;
% den_QF(rdata(:,iT) == 4;% out of range T or S so den bad too
% den_QF(rdata(:,iS)) = 4;% out of range T or S so den bad too

% CALC DEPTH
float_z  = ones(size(rdata(:,iP)))*NaN;
nan_lat  = isnan(rdata(:,iLAT));
mean_lat = nanmean(rdata(:,iLAT));
% If pos fix available
float_z(~nan_lat) = sw_dpth(rdata(~nan_lat,iP),rdata(~nan_lat,iLAT)); 
% otherwise use average lat
float_z(nan_lat) = sw_dpth(rdata(nan_lat,iP),rdata(nan_lat,iP)*0+mean_lat);

float_z_QF = fill_0 + 1; 
float_z_QF(isnan(float_z)) = 99;
clear nan_lat mean_lat


% ADD P,T,S QF's DENSITY AND DEPTH TO THE DATA SETS
if iT > iS
    iNEXT = iT+2;
else
    iNEXT = iS+2;
end

raw_hdr = [rhdr(1:iLAT+1),rhdr(iP),'PRES_QC',rhdr(iT:iT+1), rhdr(iS:iS+1)...
    'SIGMA_THETA','SIGMA_THETA_QC','DEPTH','DEPTH_QC',rhdr(iNEXT:raw_c)];

raw_data = [rdata(:,1:iLAT+1), rdata(:,iP), PRES_QF, rdata(:,iT:iT+1), ...
    rdata(:,iS:iS+1), den, den_QF, float_z, float_z_QF, ...
    rdata(:,iNEXT:raw_c)];

[raw_r,raw_c] = size(raw_data); % get new size

clear t_MVI rhdr rdata ahdr adata den float_z float_z_QF
clear PRES_QF TEMP_QF PSAL_QF den_QF

% REDO INDICES
iP  = find(strcmp('PRES',raw_hdr) == 1); 
iT  = find(strcmp('TEMP',raw_hdr) == 1); 
iS  = find(strcmp('PSAL',raw_hdr) == 1);
iO  = find(strcmp('DOXY',raw_hdr) == 1);
i380 = find(strcmp('DOWN_IRRADIANCE380',raw_hdr) == 1);
i412 = find(strcmp('DOWN_IRRADIANCE412',raw_hdr) == 1);
i490 = find(strcmp('DOWN_IRRADIANCE490',raw_hdr) == 1);
iPAR = find(strcmp('DOWNWELLING_PAR',raw_hdr) == 1);
iHS  = find(strcmp('BISULFIDE',raw_hdr) == 1);


% ************************************************************************
% ADD O2 % SAT IF O2 exists
if ~isempty(iO)
    O2sat     = raw_data(:,iO) ./ oxy_sol(raw_data(:,iT), ...
                raw_data(:,iS),0) *100;
    O2sat_QF  = raw_data(:,iO+1);  %use oxygen conc QF
    
    raw_hdr = [raw_hdr(1:iO+1),'DOXY_%SAT', 'DOXY_%SAT_QC', ...
               raw_hdr(iO+2:raw_c)];
    raw_data = [raw_data(:,1:iO+1), O2sat, O2sat_QF, raw_data(:,iO+2:raw_c)];
    
    clear O2sat O2sat_QF O2sat_adj O2sat_adj_QF
end

[raw_r,raw_c] = size(raw_data); % get new size

clear iP iT iS %iO iN

% FILL ALL QC NAN's with ARGO FILL VALUES
% NEXT STEP WILL CATCH THESE
qc_chk  = regexp(raw_hdr,'_QC$','once');
qc_chk  = ~cellfun(@isempty, qc_chk); % logical 1 for QC col
tmp = raw_data(:,qc_chk); % QC cols
tmp(isnan(tmp)) = 99;
raw_data(:,qc_chk) = tmp; 

% data_chk = logical([qc_chk(2:end),0]); % NOW GET data cols (QC -1)
% tmp = raw_data(:,data_chk); 
% tmp(isnan(tmp)) = 99999;
% raw_data(:,data_chk) = tmp;
clear qc_chk data_chk tmp

% ************************************************************************
% CONVERT ARGO DATA QUALITY FLAGS TO ODV DATA QUALITY FLAGS
% ************************************************************************
    for i = 1:raw_c
        if regexp(raw_hdr{i},'\_QC', 'once')   
            tmp1 = raw_data(:,i);
            
            tmp1(tmp1 == 0) = 10; % temporary (Argo not inspected)
            tmp1(tmp1 == 1) = 0; % Argo good               => ODV good
            tmp1(tmp1 == 2) = 0; % Argo probably good      => ODV good
            tmp1(tmp1 == 8) = 0; % Argo interpolated value => ODV good
            tmp1(tmp1 == 4) = 8; % Argo bad                => ODV bad
            tmp1(tmp1 == 3) = 4; % Argo probably bad       => ODV questionable
            tmp1(tmp1 == 5) = 4; % Argo changed value      => ODV questionable
            tmp1(tmp1 == 10 | tmp1 == 99) = 1; % NO QC or Missing value
            
            raw_data(:,i) = tmp1;
        end
    end
    clear i tmp1 tmp2 

clear chl_tmp i  c  casts iP c Boss_data Boss_hdr Boss_tmp

% ************************************************************************
% CHECK FOR MISSING VALUES - SHOULD ALL BE SET TO NaN's NOW
% ************************************************************************
raw_missing_data_str = 'No missing data found';
nan_sum =(sum(isnan(raw_data), 2))>0; % cols so I can get casts later
if sum(nan_sum,1) > 0
    tmp = unique(raw_data(nan_sum,2));
    raw_missing_data_str = ['Missing Float data detected for raw data', ...
        ' station(s): ', sprintf('%0.0f ',tmp)]; 
    disp(raw_missing_data_str);
end

clear nan_sum tmp
   
% NOW SET ALL QF NaN's to 1
for i = 1:raw_c % SET QF NaN's to 1
    if regexp(raw_hdr{i},'\_QC', 'once')
        tmp2 = raw_data(:,i);
        tmp2(isnan(tmp2)) = 1; % NO QC or Missing value
        raw_data(:,i) = tmp2;
    end
end

% ************************************************************************
% ************************************************************************
% PRINT DATA TO FILE 
% http://blogs.mathworks.com/loren/2006/04/19/high-performance-file-io/
% Change 'w' to 'W' for the file feopens to save time
% ************************************************************************
% ************************************************************************
MVI_str = '-1e10'; % MISSING VALUE INDICATOR FOR ODV

% ************************************************************************
% BUILD LOOK UP CELL ARRAY to match variables and set format string
% ************************************************************************
%RAW ODV FILE
ODV_raw(1,:)  = {'Pressure[dbar]'        '%0.2f' 'PRES' '' '' ''}; % ?
ODV_raw(2,:)  = {'Temperature[°C]'       '%0.4f' 'TEMP' '' '' ''};   
ODV_raw(3,:)  = {'Salinity[pss]'         '%0.4f' 'PSAL' '' '' ''};   
ODV_raw(4,:)  = {'Sigma_theta[kg/m^3]'   '%0.3f' 'SIGMA_THETA' '' '' ''};
ODV_raw(5,:)  = {'Depth[m]'              '%0.3f' 'DEPTH' '' '' ''};
ODV_raw(6,:)  = {'Oxygen[µmol/kg]'       '%0.1f' 'DOXY' '' '' ''};   
ODV_raw(7,:)  = {'OxygenSat[%]'          '%0.1f' 'DOXY_%SAT' '' '' ''};
ODV_raw(8,:)  = {'Nitrate[µmol/kg]'      '%0.2f' 'NITRATE' '' '' ''}; 
ODV_raw(9,:)  = {'Chl_a[mg/m^3]'         '%0.4f' 'CHLA' '' '' ''};   
ODV_raw(10,:) = {'b_bp700[1/m]'          '%0.6f' 'BBP700' '' '' ''};
ODV_raw(11,:) = {'CDOM[ppb]'             '%0.2f' 'CDOM' '' '' ''};
ODV_raw(12,:) = {'pHinsitu[Total]'       '%0.4f' 'PH_IN_SITU_TOTAL' '' '' ''};   
ODV_raw(13,:) = {'CP660[1/m]'            '%0.4f' 'CP660' '' '' ''}; 

% ADD THESE FOR ODV FLAVOR #2 -PROVOR
ODV_raw(14,:) = {'D_IRRAD380[W/m^2/nm]'  '%4.4f' 'DOWN_IRRADIANCE380' '' '' ''}; 
ODV_raw(15,:) = {'D_IRRAD412[W/m^2/nm]'  '%4.4f' 'DOWN_IRRADIANCE412' '' '' ''}; 
ODV_raw(16,:) = {'D_IRRAD490[W/m^2/nm]'  '%4.4f' 'DOWN_IRRADIANCE490' '' '' ''}; 
ODV_raw(17,:) = {'D_PAR[W/m^2/nm]'       '%4.4f' 'DOWNWELLING_PAR' '' '' ''}; 
ODV_raw(18,:) = {'Bisulfide[µmol/kg]'    '%4.4f' 'BISULFIDE' '' '' ''}; 

% ************************************************************************
% FIGURE OUT ODV FILE FORMAT TYPE: [NO pH] ]pH] [NAVIS]
% REMOVE SPECIFC VARIABLES FROM THE LOOKUP CELL ARRAY
% % if isempty(iPH)
% %     ind1 = find(strcmp('pHinsitu[Total]',ODV_raw(:,1))   == 1); 
% %     ODV_raw(ind1,:) = [];
% % end

if isempty(i380) && isempty(i412) && isempty(i490) && isempty(iPAR)
    ind1 = find(strcmp('D_IRRAD380[W/m^2/nm]',ODV_raw(:,1))   == 1);
    ind2 = find(strcmp('D_IRRAD412[W/m^2/nm]',ODV_raw(:,1))   == 1);
    ind3 = find(strcmp('D_IRRAD490[W/m^2/nm]',ODV_raw(:,1))   == 1);
    ind4 = find(strcmp('D_PAR[W/m^2/nm]',ODV_raw(:,1))   == 1);
    ODV_raw([ind1,ind2,ind3,ind4],:) = [];
end

if isempty(iHS)
    ind1 = find(strcmp('Bisulfide[µmol/kg]',ODV_raw(:,1))   == 1);
    ODV_raw(ind1,:) = [];
end

raw_var_ct = size(ODV_raw,1);

% ************************************************************************
% ************************************************************************
% CREATE RAW ACII FILE - LOW RES
% ************************************************************************
% ************************************************************************

% PRINT META DATA HEADER LINES FIRST
disp(['Printing raw data to: ',outputDIR, 'ODV',strtrim(info.WMO(1,:)),'.TXT']);
fid_raw  = fopen([outputDIR, 'ODV',strtrim(info.WMO(1,:)),'.TXT'],'W','n','UTF-8');  

fprintf(fid_raw,'//0\r\n');
fprintf(fid_raw,'//<Encoding>UTF-8</Encoding>\r\n');
fprintf(fid_raw,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
       '\r\n']);
fprintf(fid_raw,['//WMO ID: ',info.WMO(1,:),'\r\n']);
fprintf(fid_raw,['//MBARI ID: ','ODV',strtrim(info.WMO(1,:)),'\r\n']);
fprintf(fid_raw,['//Data source: ',info.data_file,'\r\n']);
fprintf(fid_raw,['//DAC: ',info.DAC,'\r\n']);
fprintf(fid_raw,['//Float type: ',info.type(1,:),'\r\n']);
fprintf(fid_raw,'//\r\n');

fprintf(fid_raw,'//PLEASE READ:\r\n');
fprintf(fid_raw,['//Data for this file has been extracted from the ',...
    'Synthetic profiles (*SProf.nc)\r\n']);
fprintf(fid_raw,'//Only RAW data parameters were used\r\n');
fprintf(fid_raw,['//', missing_profile_str,'\r\n']);
if ~isempty(interp_pos_str)
    fprintf(fid_raw,['//',interp_pos_str,'\r\n']);
end
if ~isempty(missing_pos_str)
    fprintf(fid_raw,['//',missing_pos_str,'\r\n']);
end

fprintf(fid_raw,'//\r\n');


% fprintf(fid_raw,'//Data was binned according to the following table:\r\n');
% fprintf(fid_raw,'//Start depth\tStop depth\tbin size\r\n');
% for i = 1:size(info.bins,1)
%     if i == size(info.bins,1)
%         fprintf(fid_raw,'//%4.2f\tmax depth\t%4.2f\r\n',info.bins(i,1), ...
%             info.bins(i,2));
%         fprintf(fid_raw,'//\r\n');
%     else
%         fprintf(fid_raw,'//%4.2f\t%4.2f\t%2.2f\r\n',info.bins(i,1), ...
%             info.bins(i+1,1),info.bins(i,2));
%     end
% end

fprintf(fid_raw,['//Missing data value = ',MVI_str,'\r\n']);
fprintf(fid_raw,['//Argo to ODV QF conversion: Argo 0|99 => ODV 1, ', ...
    'Argo 1|2|8 => ODV 0, Argo 3|5 => ODV 4, Argo 4 => ODV 8\r\n']);
fprintf(fid_raw,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
    '1=Missing or not inspected \r\n']);

% NOW PRINT THE RAW DATA HEADER
std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
                  'Lon [°E]' 'Lat [°N]' 'QF'}; % SIZE = 8
std_size = size(std_ODV_vars,2);

for i = 1:std_size % PRINT STANDARD HEADER VARS
    fprintf(fid_raw,'%s\t',std_ODV_vars{1,i}); % std vars
end

for i = 1:raw_var_ct % PRINT FLOAT SPECIFIC HEADER VARS
    if i < raw_var_ct
        fprintf(fid_raw,'%s\t%s\t',ODV_raw{i,1},'QF'); % std vars
    else
        fprintf(fid_raw,'%s\t%s\r\n',ODV_raw{i,1},'QF'); % std vars
    end
end

% BUILD DATA FORMAT STRING & CREATE DUMMY MATRIX TO DUMP DATA TO
% THIS CATCHES FLOATS WHICH DON'T HAVE SENSORS FOR GIVEN VARIABLES
% THESE COLUMNS WILL GET FILLED WITH -1e10 and QC FLAG = 1
dummy_out  = ones(raw_r, raw_var_ct) * NaN;
fill_MVI = ones(raw_r, 1) * -1e10;
fill_QC  = ones(raw_r, 1);
ODV_std_f = '%s\t%0.0f\t%s\t%s\t%s\t%0.3f\t%0.3f\t%0.0f\t'; %std_vars
ODV_raw_f =''; ODV_space_f = '';
for i = 0:raw_var_ct-1
    c_ct = i*2+1; % Need to add data and QC col
    ind = find(strcmp(ODV_raw{i+1,3}, raw_hdr) == 1);
    
    if ~isempty(ind)
        if i < raw_var_ct-1
            ODV_raw_f   = [ODV_raw_f, ODV_raw{i+1,2},'\t%0.0f\t'];
            ODV_space_f = [ODV_space_f,'\t\t'];
        else
            ODV_raw_f   = [ODV_raw_f, ODV_raw{i+1,2},'\t%0.0f\r\n'];
            ODV_space_f = [ODV_space_f,'\t\r\n'];
        end
        dummy_out(:,c_ct:c_ct+1)  = raw_data(:,ind:ind+1);
    else
        if i < raw_var_ct-1
            ODV_raw_f   = [ODV_raw_f, '%1.0E\t%0.0f\t'];
            ODV_space_f = [ODV_space_f,'\t\t'];
        else
            ODV_raw_f   = [ODV_raw_f, '%1.0E\t%0.0f\r\n'];
            ODV_space_f = [ODV_space_f,'\t\r\n'];
        end
        dummy_out(:,c_ct:c_ct+1)  = [fill_MVI,fill_QC]; % No data, but need cols
    end
end


% NOW PRINT DATA LINES TO FILE
cast_num  = 0; %initalize
prof_num  = 0;
line_ct   = 0;
for sample_ct = 1 : raw_r
    if raw_data(sample_ct,2) - cast_num > 0 % build standard cast part
        if sample_ct > 1
            fprintf(fid_raw,[std_str,ODV_space_f]); % add profile spacer line
            line_ct = line_ct+1;
        end
        prof_num  = 0; % reset profile counter
%         cast_num = raw_data(sample_ct,2);
%         date_str = datestr(raw_data(sample_ct,4),'mm/dd/yyyy');
%         time_str = datestr(raw_data(sample_ct,4),'HH:MM');
%         std_str  = sprintf(ODV_std_f, info.WMO(1,:), cast_num, 'C', ...
%             date_str, time_str, raw_data(sample_ct,5:7));
%         std_str = regexprep(std_str,'NaN',MVI_str);
    end
    
    %This changes for every profile
    if raw_data(sample_ct,3) - prof_num > 0 % build standard cast part
        cast_num = raw_data(sample_ct,2);
        prof_num = raw_data(sample_ct,3);
        date_str = datestr(raw_data(sample_ct,4),'mm/dd/yyyy');
        time_str = datestr(raw_data(sample_ct,4),'HH:MM');
        std_str  = sprintf(ODV_std_f, info.WMO(1,:), cast_num, 'C', ...
            date_str, time_str, raw_data(sample_ct,5:7));
        std_str = regexprep(std_str,'NaN',MVI_str);
    end    
    
   
    
    data_str = sprintf(ODV_raw_f, dummy_out(sample_ct,:));
    out_str = [std_str,data_str];
    % replace NaN' w/ missing value indicator
    out_str = regexprep(out_str,'NaN',MVI_str);  % replace NaN' w/ missing value indicator
    fprintf(fid_raw, '%s', out_str);
    line_ct = line_ct+1;
    clear out_str data_str date_str time_str
    if sample_ct == raw_r
        fprintf(fid_raw,[std_str,ODV_space_f]); % add profile spacer line
        line_ct = line_ct+1;
    end    
end

fclose(fid_raw);
clear fid_raw cast_num sample_ct

disp(['DONE printing raw data to: ',outputDIR, 'ODV',strtrim(info.WMO(1,:)),'.TXT']);

% MAKE CONFIG FILE FOR FLOATVIZ
% fid_raw  = fopen([dirs.txt, 'ODV',strtrim(info.WMO(1,:)),'.CFG'],'w');
% fprintf(fid_raw,'//%0.0f\r\n',line_ct);
% fclose(fid_raw);
% clear fid_raw line_ct dummy_out
% 
% copy_dest = '\\atlas\tempbox\Plant\Mprof2ODV\';
% disp(['Copying ODV',strtrim(info.WMO(1,:)), ' to ',copy_dest]);
% copyfile([dirs.txt, 'ODV',strtrim(info.WMO(1,:)),'.*'], copy_dest);


