function mprof_out = mprof2mat(f_info)

% ************************************************************************
% PURPOSE: 
%    This function should be used prior to "mprofmat2ODV.m".  It is an 
%    intermediary step in generating  ODV compaitble text files from Argo
%    *Mprof NetCDF files. 
%	
%
% INPUTS:
%   f_info: a structure containing the following:
%       f_info.WMO = WMO number (string) of float of interest 
%       f_info.fn = Mprof file name of interest (to convert)
%       f_info.dac_path = path to GDAC where Mprof file lives
%       f_info.local_path = Local path where Mprof file lives
%       f_info.dac = dac that 'owns' the float.
%
%   For example:
% **********************************************************
% TESTING
% f_info.WMO        = '2901562';
% f_info.fn         = '2901562_Mprof.nc';
% f_info.dac_path   = '/ifremer/argo/dac/csio/2901562/';
% f_info.local_path = 'C:/temp/';
% f_info.dac        = regexp(finfo.dac_path,'(?<=dac/)\w+','once','match');
% f_info.dac        = 'csio';
% **********************************************************

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

% FILE NAMES
mprof_file = f_info.fn;
mprof_path = f_info.local_path;

INFO.data_file = mprof_file;
INFO.DAC = f_info.dac;

INFO.bins = [0   0.5; ... % start depth range  and bin size
             5   1; ...
             50  2; ...
             200 5];        

% FILL VALUES
fv.bio = 99999;
fv.QC  = 99;

% RANGE CHECk VALUES
RCR.S     = [26 38]; % from argo parameter list
RCR.T     = [-2.5 40]; % from argo parameter list
RCR.O     = [-5 550]; % from argo parameter list
RCR.OP    = [10 70]; % optode phase, from argo parameter list
RCR.OT    = [-2.5 40]; % optode temperature, from argo parameter list
RCR.CHL   = [-0.1 150]; % argoBGC QC manual 09July2016
RCR.BB700 = [-0.000025 0.1]; % argoBGC QC manual 09July2016
RCR.BB532 = [-0.000005 0.1]; % argoBGC QC manual 09July2016
RCR.CDOM = [-1000 1000]; % Place Holder
RCR.NO3  = [-15 65];

% GET CORE AND BGC VARIABLE NAMES
tmp = get_BGC_param_list(param_list_file);
bgc_param_list = tmp.list(:,2);

% ONLY VALID FOR 
bgc_param_list = [bgc_param_list; 'BISULFIDE'];

clear tmp


% CREATE NetCDF OBJECT
ds = ncdataset([mprof_path, mprof_file]);
flt_param_list = ds.variables; % GET FLOAT PARAMETERS

% CHECK FLOAT FOR BGC VARIABLES
tf_var = ones(size(bgc_param_list,1),1)*0;
for i = 1:size(tf_var,1)
    tf = sum(strcmp(bgc_param_list{i},flt_param_list));
    if sum(tf) == 1
        tf_var(i) = 1;
    end
end
bgc_vars = bgc_param_list(logical(tf_var)); % LIST OF BGC VARS FROM FLOAT
clear bgc_param_list flt_param_list tf tf_var i param_list_file

INFO.params = bgc_vars;

% ************************************************************************
% **********               OK NOW GET FLOAT DATA            **************
% Data col are depth levels & rows are different "data profiles"
% Rows can be an actual profile or it can be a different sampling on the
% same profile (same SDN0
% ************************************************************************
C_NUM = double(ds.data('CYCLE_NUMBER'));
LAT   = double(ds.data('LATITUDE')); 
LON   = double(ds.data('LONGITUDE'));

szWMO = ds.size('PLATFORM_NUMBER');
INFO.WMO   = ds.data('PLATFORM_NUMBER',[1,1],[1,szWMO(2)]); % this is an m x 8 char array to start
INFO.WMO   = strtrim(INFO.WMO); % sometimes a trailing space
WMO = ones(szWMO(1),1) * str2double(INFO.WMO(1,:));
sztype = ds.size('PLATFORM_TYPE');
INFO.type   = strtrim(ds.data('PLATFORM_TYPE',[1,1],[1,sztype(2)])); 
clear szWMO sztype

% CONVERT JULD TO MATLAB SDN - NEED TO GET PIVOT TIME FIRST FROM ATTRIBUTES
JULDa = ds.attributes('JULD'); %  
units = JULDa{strcmp('units',JULDa(:,1)),2}; % string of pivot time
start_sdn = regexp(units,'\d{4}\D\d{2}\D\d{2}\D\d{2}\D\d{2}\D\d{2}',...
    'match', 'once'); % get date string part of string: yyyy-mm-dd HH:MM:SS
ymdhms = str2double(regexp(start_sdn,'\-|\ |\:', 'split')); % get time #'s
SDN = double(ds.data('JULD')) + datenum(ymdhms); % convert to Matlab SDN
clear JULDa UNITS start_sdn ymdhms

flt_info = [WMO, C_NUM, SDN, LON, LAT]; % WMO CYCLE SDN LON LAT
uSDN = unique(SDN); % EACH PROFILE SHOULD HAVE A UNIQUE SDN

% CHECK FOR NaN's IN TIME - DON'T PROCESS FOR NOW
tnan   = isnan(flt_info(:,3));
uCycle = unique(flt_info(tnan,2));
for i = 1:size(uCycle,1)
    if i == 1
        fprintf('No SDN for the folowing cycles: ')
        fprintf('%3.0f ',uCycle(i));
    elseif i == size(uCycle,1)
        fprintf('%3.0f ',uCycle(i));
        fprintf('\nTHESE CYCLES WILL NOT BE PROCESSED!!\n');
    else
        fprintf('%3.0f ',uCycle(i));
    end
end
uSDN(isnan(uSDN)) = [];


% TEST FOR NaN's in SDN - NOTIFY & THEN REMOVE FOR NOW

% PUT ALL DATA MATRICES IN STRUCTURE SO I CAN ACCESS FIELDS INSTEAD OF
% USING EVAL A BUNCH
for i = 1:size(bgc_vars,1)
    str = ['d.',bgc_vars{i},'=double(ds.data(''',bgc_vars{i},'''));'];
    eval(str)
    t_fill = d.(bgc_vars{i}) == fv.bio; % fill values to NaN
    d.(bgc_vars{i})(t_fill) = NaN;
end
clear vdir tAD str t_fill i C_NUM LAT LON WMO

% PUT ALL QC MATRICES IN STRUCTURE SO I CAN ACCESS FIELDS INSTEAD OF
% USING EVAL A BUNCH. QC VALUES ARE STORED AS A CHAR ARRAY ONLY THE LENGTH
% OF PROFILE DATA NOT THE MAX MATRIX COLUMN SIZE (TO SAVE SPACE I THINK)
% EXPAND OUT AND MAKE NUMERIC MATRICES FOR EACH VARIABLE
disp('Converting QC character arrays to numeric matrix ....');
for i = 1:size(bgc_vars,1)
    str = ['stmp=ds.data(''',bgc_vars{i},'_QC',''');']; %get char matrix
    eval(str) % char array of QC values and spaces
    tmp = ones(size(stmp))*NaN; % predim QC matrix
    for j = 1:size(stmp,1)
        qcvals = sscanf(stmp(j,:),'%1f')';
        tmp(j,1:size(qcvals,2)) = qcvals;
    end
    eval(['QC.',bgc_vars{i},'=tmp;'])
    clear str tmp tmps qcvals j
end
clear vdir tAD str t_fill i

% GET CYCLE DIRECTION - ONLY WANT ASCENDING PROFILES
cycle_dir = ds.data('DIRECTION');

% ESTIMATE MASTER MATRIX SIZE FROM PRESSURE MATRIX
r_ct = 0;
for ct = 1:size(uSDN,1)
    t1 = SDN == uSDN(ct);
    P  = d.PRES(t1,:);
    P(isnan(P)) = [];
    P = unique(P(:));
    r_ct   = r_ct + size(P,1);
end
total_rows = r_ct*2; % over ct by *3
clear ct r_ct P

% BUILD HEADERS FOR EVENTUAL OUPUT AND PREDIMENSION OUTPUT MATRIX
info_hdr = {'WMO', 'CYCLE','Profile' 'SDN' ,'Lon', 'Lat'};
var_hdr ={};
for vct = 1:size(bgc_vars,1)
    var_hdr = [var_hdr, bgc_vars{vct}, [bgc_vars{vct},'_QC']];
end

out_hdr = [info_hdr,var_hdr];

% GET SOME INDICES
iWMO = find(strcmp('WMO',info_hdr) ==1);
iCYC = find(strcmp('CYCLE',info_hdr) ==1);
iSDN = find(strcmp('SDN',info_hdr) ==1);
iLON = find(strcmp('Lon',info_hdr) ==1);
iLAT = find(strcmp('Lat',info_hdr) ==1);
iP   = find(strcmp('PRES',var_hdr) ==1);
iO   = find(strcmp('DOXY',var_hdr) ==1);
iS   = find(strcmp('PSAL',var_hdr) ==1);
iT   = find(strcmp('TEMP',var_hdr) ==1);

% ASSUMING EACH UNIQUE PROFILE HAS ITS OWN SDN

% STEP THROUGH PROFILES DEFINED BY UNIQUE SDN (CAN BE MULTIPLE IN CYCLE)
out = ones(total_rows, size(out_hdr,2))* NaN; % overdim & clean up later
profile_ct = 0;
cycle_ct = 0;
line_ct = 0;

% % CHECK FOR NaN's IN TIME - DON'T PROCESS FOR NOW
% tnan   = isnan(flt_info(:,3));
% uCycle = unique(flt_info(tnan,2));
% fprintf('No SDN for the folowing cycles: ')
% for i = 1:size(uCycle,1)
%     if i == size(uCycle,1)
%         fprintf('%3.0f ',uCycle(i));
%         fprintf('\nTHESE CYCLES WILL NOT BE PROCESSED!!\n');
%     else
%         fprintf('%3.0f ',uCycle(i));
%     end
% end

clear i tnan uCycle

for ct = 1:size(uSDN,1) % STEP THROUGH PROFILES
    if strcmp(cycle_dir(ct),'D')
        disp('skipping descending profile')
        continue
    end
    t1 = flt_info(:,3) == uSDN(ct); % flag unique sdn blocks
    
    info = flt_info(t1,:);
    info = info(1,:); % just want first row - all rows the same
    if cycle_ct == info(2)
        profile_ct = profile_ct +1;
    else
        profile_ct = 1;
    end
    cycle_ct = info(2);
    info = [info(1:2),profile_ct,info(3:end)];
    disp(['Cycle ', num2str(info(2)),'   profile ', num2str(profile_ct)])
    
    % ESTIMATE SIZE FOR PREDIM
    P = d.PRES(t1,:);
    % Make one col 
    % acum col by col, order not important if all data treated the same
    P       = P(:);
    Pnan    = isnan(P); % use flags in parameter loop below
    P(Pnan) = [];
    [~,ind] = sort(P,1); % shallow to deep
    %[~,ind] = sort(P,1,'descend'); % deep to shallow
    nP      = size(P,1); 
    mprof   = ones(nP, size(bgc_vars,1)*2)*NaN; % *2 for QC cols
    clear P
    
    % STEP THROUGH BGC PARAMETERS WITHIN A PROFILE & FILL mprof matrix
    for vct = 1:size(bgc_vars,1) 
        %iX = find(strcmp(bgc_vars{vct},var_hdr) ==1);
        
        X  = d.(bgc_vars{vct})(t1,:);
        if isempty(X)
            disp('No parameter data found')
            continue
        end
        X  = X(:);
        X(Pnan) = []; % remove NaN's in P due to "squareness" of netcdf
        
        
        XQC  = QC.(bgc_vars{vct})(t1,:); % QC flags for parameter X
        XQC  = XQC(:);
        XQC(Pnan) = []; % remove NaN's in P due to "squareness" of netcdf
        
%         if strcmp(bgc_vars{vct},'CHLA')
%             Xbad = XQC > 3;
%         else
%             Xbad = XQC > 2; % set data flagged bad & bad flags to NaN
%         end
        
        % UPDATED QC FLAG CRITERIA SO  MBARI RAW DOES GET SET TO NaN
        % 04/16/2017. THIS MAY allow more bad NON MBARI ARGO DATA through
        Xbad = XQC == 4;
        
        X(Xbad)   = NaN;
        XQC(Xbad) = NaN;
        
        mprof(:,vct*2-1) = X(ind); % shallow to deep
        mprof(:,vct*2)   = XQC(ind); % shallow to deep
        
        % Save some QC flag test, remove these line after looping done
        if strcmp(bgc_vars{vct},'PRES') 
            Pbad = Xbad(ind); % sorted shallow to deep
        elseif strcmp(bgc_vars{vct},'PSAL')
            Sbad = Xbad(ind); % sorted shallow to deep
        elseif strcmp(bgc_vars{vct},'TEMP')
            Tbad = Xbad(ind); % sorted shallow to deep
        end
    end
    
    mprof(Pbad|Sbad|Tbad,:) = [];
  
    t_PSAL = ~isnan(mprof(:,iS)); % good PSAL
    if all(isnan(mprof(:))) % CHECK FOR DATA
        fprintf(['No good data retrieved for cycle %0.0f profile ', ...
            '%0.0f. Processing next profile\r\n'],cycle_ct, profile_ct)
        continue
    elseif sum(t_PSAL,1) == 0
        fprintf(['No good salinity data retrieved for cycle %0.0f profile ', ...
            '%0.0f. Processing next profile\r\n'],cycle_ct, profile_ct)
        continue        
    end
   
    % ********************************************************************
    % CHECK IF S & T EXISTS FOR EACH O2 VALUE, IF NOT INTERPOLATE
    % NOT PERFECT BUT PROBABLY GOOD ENOUGH, NEED T, S & O at a P for MLR

    if iO
        t_DOXY = ~isnan(mprof(:,iO)); % good O2
        t_S    = t_DOXY & ~t_PSAL; % O2 but no PSAL
        if sum(t_S,1) > 0 % need to interpolate
            tmp = mprof(t_PSAL,[iP, iT,iT+1 iS, iS+1]); % good S rows only
            % only 1 row of data - can't interp, assign S to O's
            if size(tmp,1) == 1
                mprof(t_S,[iT,iT+1,iS, iS+1]) = ones(sum(t_S),1)* tmp(1,2:5);
            else
                [C,ia,ic] = unique(tmp(:,1));
                tmp = tmp(ia,:); % good S & unique P
                mprof(t_S,[iT,iT+1,iS, iS+1]) = interp1(tmp(:,1),tmp(:,2:5), ...
                    mprof(t_S,iP)); % QC flags are getting interpolated need to deal with this
            end
        end
    else
        disp(['No Oxygen data found in ',mprof_file]);
    end
    
    % BIN DATA ACORDING TO BINS MATRIX
    Pmax = max(mprof(:,1)) + max(INFO.bins(:,2)) *2; % max profile depth ++
    bins = INFO.bins;
    bins(bins(:,1)> Pmax,:) =[]; % remove any bin ranges > max profile depth
    for bin_ct = 1:size(bins,1)
        if bin_ct == size(bins,1)
            tz = mprof(:,1)> bins(bin_ct,1);
            edges = bins(bin_ct,1): bins(bin_ct,2): Pmax;
        else
            tz = mprof(:,1)> bins(bin_ct,1) & mprof(:,1) <= bins(bin_ct+1,1);
            edges = bins(bin_ct,1): bins(bin_ct,2): bins(bin_ct+1,1);
        end
        tmp = mprof(tz,:);
        if isempty(tmp)
            continue
        end
        
        [N,edges,bin_ind] = histcounts(tmp(:,1),edges);
        ubin = unique(bin_ind);
        
        for bct = 1:size(ubin,1)
            t1 = bin_ind == ubin(bct);
            line_ct = line_ct+1;
            out(line_ct,:) = [info, nanmean(tmp(t1,:),1)];
        end
    end
end
t_out = isnan(out);
out(sum(~t_out,2) == 0,:) =[];

clearvars -except out out_hdr RCR INFO mprof_path mprof_file

% ***********************************************************************
% ***********************************************************************
% DO SOME ADDITIONAL QC ON THE DATA BEFORE FINALIZING
% GET SOME INDICES
iWMO = find(strcmp('WMO',out_hdr) ==1);
iCYC = find(strcmp('CYCLE',out_hdr) ==1);
iSDN = find(strcmp('SDN',out_hdr) ==1);
iLON = find(strcmp('Lon',out_hdr) ==1);
iLAT = find(strcmp('Lat',out_hdr) ==1);
iP   = find(strcmp('PRES',out_hdr) ==1);
iT   = find(strcmp('TEMP',out_hdr) ==1);
iS   = find(strcmp('PSAL',out_hdr) ==1);
iO   = find(strcmp('DOXY',out_hdr) ==1);
iN   = find(strcmp('NITRATE',out_hdr) ==1);
iC   = find(strcmp('CHLA',out_hdr) ==1);
iB   = find(strcmp('BBP700',out_hdr) ==1);

% RANGE CHECKS
tS     = out(:,iS) < RCR.S(1) | out(:,iS) > RCR.S(2); % SALINITY
bad_ct = sum(tS);
if bad_ct > 0
    out(tS,iS+1) = 4;
    disp([num2str(bad_ct),' out of range PSAL values found. QC flags set to 4']);
end
clear tS

tT     = out(:,iT) < RCR.T(1) | out(:,iT) > RCR.T(2); % TEMP
bad_ct = sum(tT);
if bad_ct > 0
    out(tT,iT+1) = 4;
    disp([num2str(bad_ct),' out of range TEMP values found. QC flags set to 4']);
end
clear tT

if ~isempty(iO)
    tO     = out(:,iO) < RCR.O(1) | out(:,iO) > RCR.O(2); % DOXY
    bad_ct = sum(tO);
    if bad_ct > 0
        out(tO,iO+1) = 4;
        disp([num2str(bad_ct),' out of range DOXY values found. QC flags set to 4']);
    end
    clear tO
end

if ~isempty(iN)
    tN     = out(:,iN) < RCR.NO3(1) | out(:,iN) > RCR.NO3(2); % NITRATE
    bad_ct = sum(tN);
    if bad_ct > 0
        out(tN,iN+1) = 4;
        disp([num2str(bad_ct),' out of range DOXY values found. QC flags set to 4']);
    end
    clear tN
end

if ~isempty(iC)
    tC     = out(:,iC) < RCR.CHL(1) | out(:,iC) > RCR.CHL(2); % CHLA
    bad_ct = sum(tC);
    if bad_ct > 0
        out(tC,iC+1) = 4;
        disp([num2str(bad_ct),' out of range CHLA values found. QC flags set to 4']);
    end
    clear tC
end

if ~isempty(iB)
    tB     = out(:,iB) < RCR.BB700(1) | out(:,iB) > RCR.BB700(2); % BBP
    bad_ct = sum(tB);
    if bad_ct > 0
        out(tB,iB+1) = 4;
        disp([num2str(bad_ct),' out of range BBP700 values found. QC flags set to 4']);
    end
    clear tB
end

% NOW STEP THROUGH BBP700 PROFILES FOR DETAILED QC
if ~isempty(iB)
    t_bad = out(:,iB+1) == 4;
    cycles = unique(out(:,iCYC));
    for i = 1: size(cycles,1)
        t1 = out(:,iCYC) == cycles(i) & ~t_bad;
        tmp = out(t1,:);
        tz = tmp(:,iP) > 400 ;
        if sum(tz) > 0
            tmpz = tmp(tz,iB);
            bbp_med = median(tmpz,'omitnan');
            bbp_std = std(tmpz,'omitnan');
            t_flyers = tmpz < (bbp_med - 3*bbp_std) | ...
                       tmpz > (bbp_med + 3*bbp_std); % set flyers to NaN
            tmpz(t_flyers) = NaN;  %Set flyers to NaN     
            bbp_med = median(tmpz,'omitnan'); %Redo stats
            bbp_std = std(tmpz,'omitnan');
                         
            tQC = bbp_med > 0.01 | bbp_std > 0.001;
            if sum(tQC) > 0
%                 F1 = figure(1);
%                 plot(tmp(:,iB),tmp(:,iP),'bo-','MarkerSize', 4);
%                 set(gca,'Ydir','Reverse')
                disp(['BBP700 on Cycle ', num2str(cycles(i)),' BBP profile ', ...
                    'failed median (0.01) or std dev test(0.001) - whole ', ...
                    'CHL & BBP profiles will be flagged bad (=4)']);
                out(t1,iB+1) = 4;
                out(t1,iC+1) = 4;
%                 pause
            end
        end
    end
end

% DO A FINAL WIPE BASED ON ODVbad_sensor_list
bad  = parse_bad_ODVsensor_list([dirs.bad_sensor,bad_sensor_file]);
bad_hdr  = bad.hdr;
bad_list = bad.list;

% GET HEADER INDICES
iWMO = find(strcmp('WMO #',bad_hdr) == 1);
iM   = find(strcmp('MBARI ID STR',bad_hdr) == 1);
iSEN = find(strcmp('SENSOR',bad_hdr) == 1);
iBCYC = find(strcmp('CYCLES',bad_hdr) == 1); % idividual bad profiles
iBCB  = find(strcmp('CYCLE BLOCKS',bad_hdr) == 1); % start block of bad

for i = 1:size(bad_list,1)
    tf = strcmpi(bad_list{i,iWMO}, INFO.WMO);
    if sum(tf) ~= 0
        t2 = out(:,iWMO) == str2double(INFO.WMO);
        tmp = out(t2,:);
        
        % BUILD FLAG ARRAY
        cycle_blocks = bad_list{i,iBCB};
        tbad = ones(size(tmp(:,1)))*0; %start with all good
        if ~isempty(cycle_blocks)
            for j = 1:size(cycle_blocks,1)
                tbad = tbad | (tmp(:,iCYC) >= cycle_blocks(j,1) & ...
                    tmp(:,iCYC) <= cycle_blocks(j,2));
            end
            disp('Blocks of bad profiles found on bad sensor list')
        end
        
        cycle_list = bad_list{i,iBCYC};
        if ~isempty(cycle_list)
            for j = 1:size(cycle_list,2)
                tbad = tbad | tmp(:,iCYC) == cycle_list(j);
            end
            disp('Individual bad profiles found on bad sensor list')
        end
        
        % FOR NOW if CHL or BBP is bad, sensor is bad
        if strcmp(bad_list{i,iSEN},'BBP') || strcmp(bad_list{i,iSEN},'CHL')
            tnanB = isnan(tmp(:,iB));
            tmp(tbad & ~tnanB, iB+1) = 4;
            out(t2 ,iB+1) = tmp(:,iB+1);
            
            tnanC = isnan(tmp(:,iC));
            tmp(tbad & ~tnanC, iC+1) = 4;
            out(t2 ,iC+1) = tmp(:,iC+1);
        end
    end
end

mprof_out.hdr = out_hdr;
mprof_out.data = out;
mprof_out.INFO = INFO;

clear ds
% delete([mprof_path, mprof_file]);
clearvars -except mprof_out
                        
                        
                        
                        
                        
                        
                        
                        