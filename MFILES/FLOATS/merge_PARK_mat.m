function out = merge_PARK_mat(WMO_ID, dirs, print_ODV)
% ************************************************************************
% PURPOSE: 
%    This function merges all PARK data from *.mat profile files for a
%    given float. This merged dataset can then be used to build ODV
%    compatible text files using the print_ODV flag if desired.
%
% USAGE:
%	d = merge_PARK_mat(WMO_ID, dirs,, print_ODV)
%
% INPUTS:
%   WMO_ID     = WMO ID as a string
%
%   dirs       = Either an empty variable or a structure with directory
%                strings where files are located. It must contain these
%                fields at a minimum or be empty:
%                  dirs.mat = path to mat profile files main directory
%                  dirs.cal = path to cal files and matser MBARI float list
%                  dirs.pk  = path to destination for ODV style text files
%                  dirs.msg = path to msg file repository (to get 000.msg's
%                If dirs is empty, default paths will be used.
%
%   print_ODV  = 0 or 1 flag, 1 generate an ODV style text file
%
% OUTPUTS:
%   out             = A stucture of merged float data & information
%     out.hdr       = cell array of data column parameter names
%     out.data      = a matrix of float park data
%     out.park_vars = desired park variables found in TRAJ
%     out.pos_fixes = float profile position fixes;
%     out.msgs      = Any processing msgs that may get triggered

% ************************************************************************
% CHANGE HISTORY
%
%
% 9/12/2023 EC removed block autoflagging QC to 3
% ************************************************************************

%TESTING
% dirs =[];
% print_ODV = 1;
% WMO_ID = '5906306'; % un1114 
% WMO_ID = '5903593'; % ua7552 
% WMO_ID = '5906226'; % ua18829
% WMO_ID = '5905074';
% WMO_ID = '5905383'; % date str problem NaN'a in park dath date??
% WMO_ID = '5905380'; % date str problem NaN'a in park dath date??
% WMO_ID = '5906517';
%WMO_ID = '5906517';
% WMO_ID = '5906043'; % Fill values in CHL / BBP
% WMO_ID = '5906245'; % float cross 0E a few times!
%WMO_ID = '5905135';
% WMO_ID = '5906224'; % cycle 25 issue


% *************************************************************************
%  SET DATA DIRS AND PATHS
% *************************************************************************
out.data = []; % default function return in not successful
out.hdr  = {};
out.info = [];
out.msgs = {};
msgs     = cell(300,1); % predim error msg log
msg_ct   = 0;

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS. ONLY TWO ARE NEEDED FOR THIS FUNCTION TO OPERATE. ****
if isempty(dirs)
    user_dir  = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir  = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];
    dirs.mat  = [user_dir,'FLOATS\'];
    %dirs.mat  = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATS\'; % TESTING
    dirs.cal  = [user_dir,'CAL\'];
    %dirs.pk   = [user_dir,'FLOATVIZ\PARK\'];
    dirs.pk       = [user_dir,'FLOATVIZ_ParkData\'];
    dirs.msg  = '\\seaecho.shore.mbari.org\floats\';
elseif ~isstruct(dirs)
    disp(['Check "dirs" input. dirs must be an empty variable or a ', ...
        'structure of directory paths as strings'])
    return
end

list_path  = [dirs.cal, 'MBARI_float_list.mat']; % MBARI master list file path
% *************************************************************************
%  GET ALL ID'S ASSOCIATED WITH WMO INPUT
% *************************************************************************

if exist(list_path, 'file') % load list variables
    load(list_path)
else
    disp('BUILDING FLOAT ID LIST ...')
    d = MBARI_float_list([]); %leave dirs structure empty; the function will build it, ensuring all fields are filled appropriately.
end
iWMO  = find(strcmp('WMO',d.hdr) == 1);
iMB   = find(strcmp('MBARI ID',d.hdr) == 1);
iINST = find(strcmp('INST ID',d.hdr) == 1);
iMDIR = find(strcmp('msg dir',d.hdr) == 1);
ind   = find(strcmpi(WMO_ID, d.list(:,iWMO)) == 1);

if isempty(ind)
    str    = sprintf('No WMO_ID match found in lookup table for %s',WMO_ID);
    msg_ct = msg_ct+1;
    msgs{msg_ct} = str;
    disp(str);
    disp('NO MATCH FOUND - EXITING')
    return
end

INST_ID_str  = d.list{ind,iINST};
MBARI_ID_str = d.list{ind,iMB};
clear d

% *************************************************************************
%  MERGE CYCLE PARK DATA IF IT EXISTS
% *************************************************************************
fprintf('Merging park data for float %s (%s)\n',MBARI_ID_str, WMO_ID);
float_fp = fullfile(dirs.mat,WMO_ID,'\');

if isfolder(float_fp) % check for WMO folder, get file listing if it exists
    tmp       = dir(fullfile(float_fp,[WMO_ID,'*.mat']));
    flist     = {tmp.name}';
    if size(flist,1) == 0 %directory exists but no .mat files inside yet.
        str = sprintf('WMO directory found for %s (%s) but no .mat files inside',...
            WMO_ID, MBARI_ID_str);
        msg_ct = msg_ct+1;
        msgs{msg_ct} = str;
        disp(str);
        return
    end

    t000      = contains(flist,'000.mat'); % test for 000.mat file. Remove if exists
    flist = flist(~t000,:);
    clear tmp
else
    str = sprintf('No WMO directory found for %s (%s)', WMO_ID, MBARI_ID_str);
    msg_ct = msg_ct+1;
    msgs{msg_ct} = str;
    disp(str);
    return
end

r_list    = size(flist,1);

% ************************************************************************
% LOOK FOR 000.MSG FILE FOR DEPLOMENT FIX
% *.mat cycle data exist! Try and get 000.msg position for starting fix.
% This will be used for fix interpolation.
%
% If working from home & VPN off Matlab bogs down trying to get to sea echo
% no matter how you test!  Use DOS "ping" to check for existance.
% ************************************************************************
msg0_fix    = [];
[status, ~] = system('ping seaecho.shore.mbari.org -n 1'); % should work for  pc & mac?
if status == 0 % 0 means that seaecho was sucessfully pinged
    cal_fn = sprintf('cal%s.mat',MBARI_ID_str);
    load(fullfile(dirs.cal,cal_fn));
    d    = get_msg_list(cal.info, 'msg');
    t000 = ~cellfun(@isempty, regexp(d.list(:,1),'000\.msg','once'));
    if sum(t000) == 1
        fp000 = fullfile(d.list{t000,2},d.list{t000,1});
        opts  = delimitedTextImportOptions('Delimiter','\n');
        d     = readcell(fp000, opts);
        tFIX  = ~cellfun(@isempty, regexp(d,'Fix:','once'));
        if sum(tFIX) == 1
            tmp = regexp(d{tFIX},'\s+','split');
            % lon, lat, date, time
            msg0_fix = [0, datenum([tmp{4},' ',tmp{5}],'mm/dd/yyyy HHMMSS'), ...
                str2double(tmp(2:3))];
            if msg0_fix(3) < 0
                msg0_fix(3) = msg0_fix(3)+360; % Convert 2 0-360 E
            end
        else
            str = sprintf('No 000.msg file found for %s (%s)',...
                MBARI_ID_str, WMO_ID);
            msg_ct = msg_ct+1;
            msgs{msg_ct} = str;
            disp(str);
        end
    else
        str = sprintf('No 000.msg file found for %s (%s)',...
            MBARI_ID_str, WMO_ID);
        msg_ct = msg_ct+1;
        msgs{msg_ct} = str;
        disp(str);
    end
    clear d t000 cal_fn
else
    fprintf('No network connection to seacho - no 000.msg retrieved\n')
end

% *************************************************************************
% LIST ALL POSSIBLE PARK VARIABLES TO LOOK FOR. THIS CAN BE UPDATED AS NEW
% VARIABLES COME ON LINE OR TO CHANGE DEFAULTS
% This list defines what to look for in TRAJ.PARK, how to name the column
% in the ODV file and what formating to use when printing the variable
% *************************************************************************
park_vars(1,:)  = {'Station'              '%0.0f'  'Cycle'}; 
park_vars(2,:)  = {''                     ''       'SDN'};
park_vars(3,:)  = {'Lon [ºE]'             '%0.3f'  'Lon'}; 
park_vars(4,:)  = {'Lat [ºN]'             '%0.3f'  'Lat'};
park_vars(5,:)  = {'QF'                   '%0.0f'  'SDN_QC'};
park_vars(6,:)  = {'Pressure[dbar]'       '%0.2f'  'PRES'}; 
park_vars(7,:)  = {'QF'                   '%0.0f'  'PRES_QC'}; 
park_vars(8,:)  = {'Temperature[°C]'      '%0.4f'  'TEMP'};
park_vars(9,:)  = {'QF'                   '%0.0f'  'TEMP_QC'};
park_vars(10,:)  = {'Salinity[pss]'        '%0.4f'  'PSAL'}; 
park_vars(11,:) = {'QF'                   '%0.0f'  'PSAL_QC'};      
park_vars(12,:) = {'Oxygen_adj[µmol/kg]'  '%0.2f'  'DOXY_ADJUSTED'};    % some NAVIS floats
park_vars(13,:) = {'QF'                   '%0.0f'  'DOXY_ADJUSTED_QC'}; % some NAVIS floats
park_vars(14,:) = {'Chl_a_adj[mg/m^3]'    '%0.4f'  'CHLA_ADJUSTED'};  
park_vars(15,:) = {'QF'                   '%0.0f'  'CHLA_ADJUSTED_QC'};
park_vars(16,:) = {'b_bp700_adj[1/m]'     '%0.6f'  'BBP700_ADJUSTED'};
park_vars(17,:) = {'QF'                   '%0.0f'  'BBP700_ADJUSTED_QC'};
park_vars(18,:) = {'pHinsitu_adj[Total]'  '%0.4f'  'PH_IN_SITU_TOTAL_ADJUSTED'};    % some NAVIS floats
park_vars(19,:) = {'QF'                   '%0.0f'  'PH_IN_SITU_TOTAL_ADJUSTED_QC'}; % some NAVIS floats

% park_vars(12,:) = {'Oxygen[µmol/kg]'      '%0.2f'  'DOXY'};             % some NAVIS floats
% park_vars(13,:) = {'QF'                   '%0.0f'  'DOXY_QC'}; 
% park_vars(16,:) = {'Chl_a[mg/m^3]'        '%0.4f'  'CHLA'};  
% park_vars(17,:) = {'QF'                   '%0.0f'  'CHLA_QC'};


% ************************************************************************
%  PROCESS *.MAT FILES
% ************************************************************************
data          = [];
tf_vars_found = 0;
pos_fixes     = ones(r_list,4)*NaN; % cyle, sdn, lon, lat % used for pos interp
line_ct       = 1;

disp(' ');
disp(['Merging float profiles for ' ,MBARI_ID_str,':'])

for file_ct = 1 : r_list
    clear TRAJ PARK
    float_file = flist{file_ct};
    fp         = fullfile(float_fp,float_file);
    load(fp)
    fprintf('%0.0f ',INFO.cast)
    
    % CHECK FOR PARK DATA % should be true for all profiles but check each
    % one just in case (use continue instead of return)
    if ~exist('TRAJ','var') || ~isfield(TRAJ,'PARK')
        str = sprintf(['No PARK data found for float %s (%s) ', ...
            'cycle %0.0f'],WMO_ID, MBARI_ID_str, INFO.cast);
        msg_ct = msg_ct+1;
        msgs{msg_ct} = str;
        disp(str);
        continue % or return??
    end
    PARK = TRAJ.PARK;
    
    % ******************* TEMPORARY FIX for 9095.125 ******************
    % jp 07/05/2022 - this is related to stale files on chem
    if size(INFO.gps,2) == 2
        disp(['WARNING: Funky GPS data for ',float_file, ...
            '. Only 2 cols generating a temporary fix']);
        INFO.gps = [INFO.gps(:,1)*0 + INFO.sdn,INFO.gps];
    end
    % ******************* TEMPORARY FIX for 9095.125 ******************
        
    % ONLY WANT 1 GPS POINT
    if size(INFO.gps,1) > 1
        gps = median(INFO.gps); % if multiple fixes take median
    else
        gps = INFO.gps(1,:);
    end
    if gps(2) < 0, gps(2) = gps(2) + 360; end %index changed from 1 to 2 after addition of sdn to gps vector (TM 10/2020)
    pos_fixes(file_ct,:) = [INFO.cast,gps];

    % FIND OUT WHICH DESIRED PARK VARIABLES EXIST
    if tf_vars_found == 0
        pk_fields = fieldnames(PARK);
        ia        = ismember(park_vars(:,3), ...
            [pk_fields;'Cycle';'Lon';'Lat';'SDN_QC']); 
        wrk_vars  = park_vars(ia,3); % these are the float params found on master parm list
        park_vars = park_vars(ia,:);
        %wrk_vars = ['Cycle';'Lon';'Lat';wrk_vars];
        data      = ones(1e5,size(wrk_vars,1))*NaN; %predim
        tf_vars_found = 1;
        %clear pk_fields Lia
    end
 
    % BUILD PARK DATA MATRIX OF MERGED CYCLE DATA
    r_samp    = size(PARK.SDN,1); % get # park smples for cycle
    for i = 1:size(wrk_vars,1)
        if strcmp(wrk_vars{i},'Cycle')
            data(line_ct:line_ct+r_samp-1,i) = INFO.cast;
        elseif regexp(wrk_vars{i},'Lon|Lat|SDN_QC','once')
            continue
        else
            tmp   = PARK.(wrk_vars{i});
            tFill = tmp == 99999;

%             % FORCE MIN QF == 3 - TEMPORARY UNTIL FIXED IN PROCESSING CODE
%             if ~isempty(regexp(wrk_vars{i},'QC$','once')) && ...  
%                 isempty(regexp(wrk_vars{i},'^PRES|^TEMP|^PSAL|^SDN','once'))
%                 tmp(tmp<3) = 3; % Argo QF's here stil
%             end

            
            if any(tFill)
                str = sprintf(['Argo fill values detected for %s on ', ...
                    'cylcle %0.0f. Replacing fill values with NaNs'],...
                    wrk_vars{i},INFO.cast);
                msg_ct = msg_ct+1;
                msgs{msg_ct} = str;
                disp(str);
                tmp(tFill) = NaN;

                



            end
            data(line_ct:line_ct+r_samp-1,i) = tmp;

        end
    end
    line_ct = line_ct+r_samp;
end
fprintf('\n');
data = data(1:line_ct-1,:);

if size(data,1) < 2
    str = sprintf(['Mat cycle files exist but minimal data returned', ...
        ' - skiping %s (%s) for now'],WMO_ID, MBARI_ID_str);
    msg_ct = msg_ct+1;
    msgs{msg_ct} = str;
    disp(str);
    return
end

if ~isempty(msg0_fix)
    pos_fixes = [msg0_fix; pos_fixes];
end

% ************************************************************************
% INTERPOLATE PARK LAT AND LON POSITIONS FROM SURFACE POSITION FIXES
% Since park data for cycle 1 includes park data from deployment to cycle 1
% some park time stamps may occur before the cycle 1 surfce time & result
% in an NaN during interpolation if 000.msg can not be found for a
% deployment postion
iLON = find(strcmp(wrk_vars,'Lon') == 1);
iLAT = find(strcmp(wrk_vars,'Lat') == 1);
iCYC = find(strcmp(wrk_vars,'Cycle') == 1);
iSDN = find(strcmp(wrk_vars,'SDN') == 1);

tNaN   = isnan(pos_fixes(:,3)); % NaN's in Lon?
if sum(~tNaN,1) <2
    str = sprintf(['Not enough valid profile position fixes for %s (%s)',...
        'position interpolation - skipping float'],WMO_ID, MBARI_ID_str);
    msg_ct = msg_ct+1;
    msgs{msg_ct} = str;
    disp(str);
    return
end
clear tg sum_tg

% INTERPOLATE & DEAL WITH 0 EAST CROSSING
fixes_no_nan = pos_fixes(~tNaN,:);
lon_diff     = [0;diff(fixes_no_nan(:,3))];
tX           = abs(lon_diff)  > 300;
iX           = find(abs(lon_diff)  > 300);
%iX           = [iX; size(fixes_no_nan,1)]; % add end point

if ~isempty(iX) % 0 MERIDAIN CROSSING EXISTS!
    str = sprintf(['%s (%s) crosses the 0 meridian. Stepwise position ',...
        'interpolation applied'],WMO_ID, MBARI_ID_str);
    msg_ct = msg_ct+1;
    msgs{msg_ct} = str;
    disp(str);
    LonLat = ones(size(data,1),2) * NaN; % pre dim LON, LAT
    ind1 = 1;
    for ict = 1: size(iX,1)
        ind2 = iX(ict);
        tmp_fix = fixes_no_nan(ind1:ind2,:);
        %lon_diff > 0 = Eastward travel, < 0 = westward travel
        tmp_fix(end,3) = tmp_fix(end,3) - 360*sign(lon_diff(ind2)); 
        
        tmpLL = interp1(tmp_fix(:,2), tmp_fix(:,3:4),data(:,iSDN));
        tg    = ~isnan(tmpLL(:,1)); % interpolated all but just want non NaN values
        LonLat(tg,:) = tmpLL(tg,:);
        
        % Do last part seprately - probably a better way - finish last part
        % with normal interp
        if ict == size(iX,1) && ind2 < size(fixes_no_nan,1)
            tmp_fix = fixes_no_nan(ind2:end,:);
            tmpLL = interp1(tmp_fix(:,2), tmp_fix(:,3:4),data(:,iSDN));
            tg    = ~isnan(tmpLL(:,1)); % interpolated all but just want non NaN values
            LonLat(tg,:) = tmpLL(tg,:);
        end
        
        ind1 = ind2; % update lower bound
    end
    % Bring neg lon & lon > 360 back to reality
    LonLat(LonLat(:,1) < 0, 1) = LonLat(LonLat(:, 1) < 0,1) + 360;
    LonLat(LonLat(:,1)> 360, 1) = LonLat(LonLat(:, 1) > 360,1) - 360;
    data(:,[iLON,iLAT]) = LonLat;
else
    data(:,[iLON,iLAT]) = interp1(fixes_no_nan(:,2),fixes_no_nan(:,3:4),data(:,iSDN));
end

SDN_QC = ones(size(data,1),1)+1; % default = 2 = pretty good
ia    = ismember(data(:,iCYC), pos_fixes(tNaN,1));
SDN_QC(ia) = 3;
data(:,iLAT+1) = SDN_QC;

% REORDER /RENAME SOME ITEMS
out.hdr       = wrk_vars';
out.data      = data;
out.park_vars = park_vars;
out.pos_fixes = pos_fixes;

if msg_ct > 0
    out.msgs = msgs(1:msg_ct); % add processing msgs if they have been generated
end

clearvars -except park_vars dirs out INFO cal WMO_ID print_ODV iCYC iSDN iLAT iLON

% *************************************************************************
% *************************************************************************
%                      PRINT THE DATA T0 FILE
% *************************************************************************
% *************************************************************************
if print_ODV == 0
    return
end

pk_fn           = sprintf('%s_PK.TXT',WMO_ID);
pdata           = out.data;
hdr             = out.hdr;
odv_hdr         = park_vars(:,1)';
[rpdata,cpdata] = size(pdata);


% CONVERT ARGO QFs TO ODV QF SCHEME
for i = 1:cpdata
    if regexp(hdr{i},'\_QC$', 'once')
        tmp2 = pdata(:,i);

        tmp2(tmp2 == 5) = 3; % This is for NPQ corrected data
        tmp2(tmp2 == 4) = 8;
        tmp2(tmp2 == 3) = 4;
        tmp2(tmp2 == 0) = 10; % temporary
        tmp2(tmp2 == 2 | tmp2 == 1) = 0;
        tmp2(tmp2 == 10 | tmp2 == 99) = 1; % NO QC or Missing value
        pdata(:,i) = tmp2;
    end
end

% ************************************************************************
% open file & print meta data
fprintf('Printing park data ODV text file for %s to:\n %s...\n', ...
    WMO_ID, fullfile(dirs.pk, pk_fn));
fid     = fopen(fullfile(dirs.pk, pk_fn),'w');
fprintf(fid,'//0\r\n');
fprintf(fid,'//<Encoding>UTF-8</Encoding>\r\n');
fprintf(fid,['//!!!!! THIS IS AN EXPERIMENTAL DATA SET AND NOT IN A ',...
    'FINAL FORM YET - USER BEWARE !!!!!\r\n']);
fprintf(fid,'//ALL BGC PARAMETER VALUES HAVE QUALITY FLAGS SET TO QUESTIONABLE (4) AT BEST!\\r\n');
fprintf(fid,['//If you have any questions about the data set please ', ...
    'contact Emily Clark (eclark@mbari.org)\r\n']);
fprintf(fid,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'),'\r\n']);
fprintf(fid,['//WMO ID: ',cal.info.WMO_ID,'\r\n']);
fprintf(fid,['//Institution ID: ',cal.info.INST_ID,'\r\n']);
fprintf(fid,['//MBARI ID: ',cal.info.name,'\r\n']);
fprintf(fid,['//Project Name: ',cal.info.Program,'\r\n']);
fprintf(fid,['//Region: ',cal.info.Region,'\r\n']);
fprintf(fid,['//Float type: ',cal.info.float_type,'\r\n']);
fprintf(fid,'//\r\n');

if sum(strcmp(park_vars(:,1), 'Chl_a_adj[mg/m^3]')) == 1 % Park depth chl exist
    chl_str = ['//Park depth Chl_a_adj is calculated in a similar manner to profile data excluding an NPQ correction.\r\n', ...
        '//In situ dark counts (SWDC) are 1st removed from the measured signal if available,\r\n', ...
        '//otherwise factory dark counts (DC) are used before multiplying by the scale factor\r\n', ...
        '//and then dividing by 2 (Roesler global bias correction factor)\r\n'];
    fprintf(fid, chl_str);
    fprintf(fid,'//Chl_a_adj = [counts - DarkCounts] * ScaleFactor / 2\r\n');
    if isfield(cal.CHL,'SWDC')
        fprintf(fid,'//    DarkCounts=%0.0f (factory DC= %0.0f)\r\n', ...
            cal.CHL.SWDC.DC, cal.CHL.ChlDC);
    else
        fprintf(fid,'//    DarkCounts=%0.0f (factory DC= %0.0f)\r\n', ...
            cal.CHL.ChlDC, cal.CHL.ChlDC);
    end
    fprintf(fid,'//    ScaleFactor=%0.4f\r\n',cal.CHL.ChlScale);
    fprintf(fid,['//b_bp700_adj =b_bp700. No adjustements have been made to raw_bp700 and\r\n',...
        '//calculations follow BGC Argo protocols\r\n']);

end

fprintf(fid,'//\r\n');
fprintf(fid,'//Missing data value = -1e10\r\n');
fprintf(fid,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
    '1=Missing or not inspected \r\n']);
fprintf(fid,'//Note: all timestamps are in GMT. \r\n');


% ************************************************************************
% BUILD & PRINT HEADER
phdr    = ['Cruise',odv_hdr(iCYC),'Type','mon/day/yr','hh:mm',odv_hdr(iLON:end)];
hdr_fmt = [repmat('%s\t',1,size(phdr,2)-1),'%s\r\n'];
fprintf(fid, hdr_fmt, phdr{:});

% ************************************************************************
% PRINT DATA

% build fmt str
data_fmt = [];
for i = 1:size(phdr,2)
    tmp = phdr{i};
    t1 = strcmp(park_vars(:,1), tmp);
    
    if i == size(phdr,2) % Will always end with QF
        %data_fmt = [data_fmt,'%0.0f\r\n'];
        data_fmt = [data_fmt,'%0.0f'];
    elseif strcmp(tmp,'QF')
        data_fmt = [data_fmt,'%0.0f\t'];
    elseif regexp(tmp,'Cruise|Type|mon/day/yr|hh:mm','once')
        data_fmt = [data_fmt,'%s\t'];
    else
        data_fmt = [data_fmt, park_vars{t1,2},'\t'];
    end
end

% print data lines - 2 step 
% 1) print data to string & search & replace NaN with -1e10 (ODV fill value)
% 2) then print string to file
for i = 1:rpdata
    if isnan(pdata(i,iSDN)) % chck tp make sure SDN exists
        str = sprintf(data_fmt, WMO_ID, pdata(i, iCYC),'C', ...
            '-1e10', '-1e10', pdata(i,iLON:end));
        fprintf('No park data SDN resolved for:\n%s\n',str);
    else
        str = sprintf(data_fmt, WMO_ID, pdata(i, iCYC),'C', ...
            datestr(pdata(i,iSDN),'mm/dd/yyyy'), datestr(pdata(i,iSDN),'HH:MM'), ...
            pdata(i,iLON:end));
    end
    str = regexprep(str,'NaN','-1e10');
    fprintf(fid,'%s\r\n',str);
end

fclose(fid);

return
    

