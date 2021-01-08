% Float_vs_GLODAPv2
% PURPOSE: 
%   Extract all GLODAPv2 nutrient data within a given distance from each
%   track position if it exists.
fig_switch = [ 0 0 0 0 0];
% ************************************************************************
% CHOOSE FLOAT VARIABLE TO COMPARE {Float Bottle}
% ************************************************************************
%  there's a bug - you have to run pH first to get NO3 or O2 to work
% compare_var = {'Oxygen[µmol/kg]' 'OXYGEN' '[µmol/kg]'};
% compare_var = {'Nitrate[µmol/kg]' 'NITRAT' '[µmol/kg]'};
 compare_var = {'pHinsitu[Total]' 'PH_TOT_INSITU' ' '}; 
%compare_var = {'Chl_a[mg/m^3]' 'CHLORA' '[ug/kg]'}; 

% ************************************************************************
% FLOAT / BOTTLE SETTINGS
QC_flag = 1; % 1 for good data only QC flag must = 1
depth_bnds = [0 2000]; % depth bounds in meters

% ***********************************************************************
%   SET NAMES DIRS AND PATHS, LOAD FLOAT LISt & BOTTLE LOOKUP TABLE
% ***********************************************************************
% returns user path,i.e. 'C:\Users\jplant
 user_dir = getenv('USERPROFILE'); 
 user_dir = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\'];
%user_dir = ['C:\data\isus\argo\Matlab\ARGO_PROCESSING\'];

% LOAD BOTTLE DATA LOOK UP TABLE & STORE IN handles STRUCTURE
% LOOK UP TABLE  HEADER = [UW_ID  WMO   CRUISE   STN   CAST   Data file]
bottle_dir   = [user_dir,'DATA\SHIPBOARD\'];
btable_fn = 'BottleData_lookup_table.txt';
fid       = fopen([bottle_dir,btable_fn]);
btable    = textscan(fid, '%s %s %s %s %f %f %s','Delimiter', '\t','HeaderLines',1);
fclose(fid);

% LOAD FLOAT LIST
load([user_dir,'DATA\CAL\MBARI_float_list.mat']); % MBARI ID, UW ID, WMO, TYPE
FLOAT_LIST = list;

exclude_floats ={'0412HAWAII' '8514HAWAII' '8501CALCURRENT' };
for i = 1:size(exclude_floats,2)
    t1 = strcmp(exclude_floats{i},FLOAT_LIST(:,1));
    FLOAT_LIST(t1,:) = [];
end


% load(['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%       'float_type_lists.mat'])
% jp = regexpi(NON_SOCCOM_list(:,1),'RosSea|SoAtl|SoPac|SoOcn|drake');
% t1 = cellfun(@isempty,jp);
% FLOAT_LIST = [SOCCOM_list; NON_SOCCOM_list(~t1,:)];
[~, ind]   = sort(str2double(FLOAT_LIST(:,2))); % sort by UW ID
FLOAT_LIST = FLOAT_LIST(ind,:);
GOOD_LIST = FLOAT_LIST; % Whittle away at this one


%TESTING
% t1 = strcmp('0570SOOCN',GOOD_LIST(:,1));
% FLOAT_LIST = FLOAT_LIST(t1,:);
% GOOD_LIST = FLOAT_LIST; % Whittle away at this one

clear NON_SOCCOM_list SOCCOM_list jp t1 ind

% ************************************************************************
% NOW LOAD FLOAT DATA     
FV_dir =  [user_dir,'DATA\FLOATVIZ\QC\'];
   


if strcmp(compare_var{1,1},'Oxygen[µmol/kg]')
    hdr_add_str = 'Float O2 [µmole kg^{-1}]';
elseif strcmp(compare_var{1,1},'Nitrate[µmol/kg]')
    hdr_add_str = 'Float NO3 [µmole kg^{-1}]';
elseif strcmp(compare_var{1,1},'pHinsitu[Total]')
    hdr_add_str = 'Float pH';
elseif strcmp(compare_var{1,1},'Chl_a[mg/m^3]')
    hdr_add_str = 'Float CHL[mg m^{-3}]';
    FV_dir =  [user_dir,'DATA\FLOATVIZ\'];
end
    
hdr = {'UW_ID' 'Float Cycle','Cycle Date','Longitude', 'Latitude'}; 
hdr = [hdr, 'Bottle SDN','bottle pres',compare_var{1,2}, ...
       'interp T',['interp ',hdr_add_str], ...
       'lookup T',['lookup ',hdr_add_str], ...
       'Sun elevation','CHL type flag','pH meas ?']; 
data = [];

for flt_ct = 1 : size(FLOAT_LIST,1)
    if strcmp(compare_var{1},'Chl_a[mg/m^3]')
        CHL_str = compare_var{1,2};
        CHL_flag = 0;
        fname = [FLOAT_LIST{flt_ct,1},'.TXT'];
    else
        fname = [FLOAT_LIST{flt_ct,1},'QC.TXT'];
    end
    
    MBARI_ID = FLOAT_LIST{flt_ct,1};
    UW_str = regexp(fname, '^\d{3}\d+', 'once', 'match');
    UW_ID  = str2double(UW_str);
    disp(['Processing float ', fname, ' .....'])
    
    d = get_FloatViz_data([FV_dir,fname]);
    
    
    if isempty(d)
        disp(['Could not find float data for ',fname]);
        tf = strcmp(FLOAT_LIST{flt_ct,1}, GOOD_LIST);
        GOOD_LIST(tf,:) =[];
        continue
    end
    
    iPF  = find(strcmp('Pressure[dbar]',d.hdr)     == 1);
    iTF  = find(strcmp('Temperature[°C]',d.hdr)    == 1);
    iXF  = find(strcmp(compare_var{1,1},d.hdr)     == 1);
    
    if isempty(iXF)
        disp (['No variable data for ', fname, ' moving to next float'])
        tf = strcmp(FLOAT_LIST{flt_ct,1}, GOOD_LIST);
        GOOD_LIST(tf,:) =[];
        continue
    end
    
    % SUBSET - ONLY WANT 1st CAST FOR NOW
    prof1 = min(d.data(:,2));
    ind = find(d.data(:,2) > prof1, 1);
    if isempty(ind)
        prof2 = prof1;
    else
        prof2 = d.data(ind,2);
    end
    
    prof = prof1;
    %prof = prof2;
    
    t1 = d.data(:,2) == prof;
    FD = d.data(t1,:);
    Fhdr = d.hdr;
    clear d
    
    
    [Az,El] = SolarAzElq(FD(1,1),FD(1,4),FD(1,3),0); % sdn lat lon alt 
    
    % GOOD FLOAT DATA ONLY (QC FLAG = 0)
    if QC_flag == 1
        tQC = FD(:,iXF+1) == 0;  % good data only - checking quality flag
        FD(~tQC,iXF) = NaN;
    end
    
    tMIV = FD(:,iXF) == -1e10; % good data only - 2nd check on value
    FD(tMIV,iXF) = NaN;
    
    % NOW DO VARIABLE SPECIFC RANGE CHECKS TOO
    if strcmp(compare_var{1,1},'Oxygen[µmol/kg]')
        tQC = FD(:,iXF) < -5 | FD(:,iXF) > 450;
        FD(tQC,iXF) = NaN;
    elseif strcmp(compare_var{1,1},'Nitrate[µmol/kg]')
        tQC = FD(:,iXF) < -5 | FD(:,iXF) > 55;
        FD(tQC,iXF) = NaN;
    elseif strcmp(compare_var{1,1},'pHinsitu[Total]')
        tQC = FD(:,iXF) < 7.3 | FD(:,iXF) > 8.5;
        FD(tQC,iXF) = NaN;
    end
    
    clear  t1 tQC tMIV
    
    % GET CALIBRATION BOTTLE DATA IF IT EXISTS
    % LOOKUP TABLE  HEADER = [UW_ID  WMO   CRUISE   STN   CAST   Data file]
    ind = strcmp(MBARI_ID, btable{1,1});
    if sum(ind) > 0 ; % float exists in lookup table
        bottle_fname = btable{1,7}{ind};
        stn          = btable{1,5}(ind);
        cst          = btable{1,6}(ind);
        
        % data file & data exist for float
        if ~isempty(bottle_fname) && stn ~= -999 && cst ~= -999
            d = get_shipboard_data([bottle_dir, bottle_fname]);

            iStn  = find(strcmp(d.hdr,'STNNBR') == 1);
            iCast = find(strcmp(d.hdr,'CASTNO') == 1);
            tStn  = d.data(:,iStn)  == stn;
            tCast = d.data(:,iCast) == cst;
            Bcruise  = d.cruise;
            Bhdr     = d.hdr;
            %handles.bdata.units   = d.units;
            d.data(d.data == -999) = NaN;
            BD    = d.data(tStn&tCast,:);
            iPB   = find(strcmp('CTDPRS',Bhdr) == 1);
            tZ = BD(:,iPB) >= depth_bnds(1) & BD(:,iPB) <= depth_bnds(2);
            BD(~tZ,:) =[];
            
        else
            disp(['No bottle data exists for float ', UW_str])
            tf = strcmp(FLOAT_LIST{flt_ct,1}, GOOD_LIST);
            GOOD_LIST(tf,:) =[];
            continue
        end
        clear bottle_fname stn cst d IStn iCast tStn tCast
    else
        disp(['Float ', UW_str, ' was not in bottle data lookup table'])
        tf = strcmp(FLOAT_LIST{flt_ct,1}, GOOD_LIST);
        GOOD_LIST(tf,:) =[];
        continue
    end
    
    
    iXB  = find(strcmp(compare_var{1,2},Bhdr) == 1);
    
    % CHECK CHL FLAVOR
    CHL_flag = NaN;
    if isempty(iXB) && strcmp(compare_var{1,2},'CHLORA')
        CHL_str = 'TOT_CHL_A';
        CHL_flag = 1;
        iXB  = find(strcmp('TOT_CHL_A',Bhdr) == 1); % try alt dfn in file
    end
    
    % CHECK PH FLAVOR if no measured try calculated pH 02/07/2018 - jp
    if strcmp(compare_var{1,2},'PH_TOT_INSITU') && ~isempty(iXB)
        if sum(~isnan(BD(:,iXB))) == 0 % no bottle measured pH
            iXB  = find(strcmp('PH_TOT_INSITU_ALKDIC',Bhdr) == 1);
            disp('No measured pH in bottle data - trying calculated')
            meas_pH = 0;
        else
            meas_pH = 1;
        end
    end
    
    if isempty(iXB)
        disp([compare_var{1,2}, ' does not exist in bottle data ',...
            ' for float ',UW_str]); % dat exits, but not desired var
        tf = strcmp(FLOAT_LIST{flt_ct,1}, GOOD_LIST);
        GOOD_LIST(tf,:) =[];
        continue
    elseif sum(~isnan(BD(:,iXB))) == 0
        disp([compare_var{1,2}, ' exist in bottle data for float ', ...
            UW_str, 'but all data = NaN']); % data exits, but all nan
        tf = strcmp(FLOAT_LIST{flt_ct,1}, GOOD_LIST);
        GOOD_LIST(tf,:) =[];
        continue
    end
    
    % VAR EXISTS SO GET REST OF BOTTLE COL INDICES
    iPB   = find(strcmp('CTDPRS',Bhdr) == 1);
    iSDNB = find(strcmp('DATE',Bhdr)      == 1); % DATASET
    iTB   = find(strcmp('CTDTMP',Bhdr)      == 1); % DATASET
    iSB   = find(strcmp('CTDSAL',Bhdr)      == 1); % DATASET
    
    % RANGE CHECK BOTTLE DATA TOO
    % NOW DO VARIABLE SPECIFC RANGE CHECKS TOO
    if strcmp(compare_var{1,2},'OXYGEN')
        tQC = BD(:,iXB) < -5 | BD(:,iXB) > 450;
        BD(tQC,iXB) = NaN;
    elseif strcmp(compare_var{1,2},'NITRAT')
        tQC = BD(:,iXB) < -5 | BD(:,iXB) > 55;
        FD(tQC,iXB) = NaN;
    elseif strcmp(compare_var{1,2},'PH_TOT_INSITU')
        tQC = BD(:,iXB) < 7.3 | BD(:,iXB) > 8.5;
        BD(tQC,iXB) = NaN;
    end
    
    % PLOT THE DATA AS A CHECK
    if fig_switch(5) == 1
        F5 = figure(20);
        F5.Position = [488.2000 437.8000 709.6000 420.0000];
        subplot(1,2,1)
        plot(FD(:,iXF),FD(:,iPF),'bo')
        set(gca,'Ydir','Reverse')
        title([MBARI_ID,'  Cast ',num2str(prof)])
        xlabel(compare_var{1})
        ylabel('Pressure')
        ylim([0 250])
        hold on
        plot(BD(:,iXB),BD(:,iPB),'ro', 'MarkerfaceColor', 'r')
        hold off
        if exist('CHL_str', 'var')
            L5 = legend(compare_var{1},CHL_str,'Location','SouthEast');
        else
            L5 = legend(compare_var{1},'Location','SouthEast');
        end
        
        L5.Interpreter = 'none';
        
        subplot(1,2,2)
        plot(FD(:,iTF),FD(:,iPF),'bo')
        set(gca,'Ydir','Reverse')
        xlabel('Temperature')
        ylim([0 250])
        hold on
        plot(BD(:,iTB),BD(:,iPB),'ro', 'MarkerfaceColor', 'r')
        hold off
        
        pause
    end
    
    % INTERPOLATE FLOAT DATA ON TO BOTTLE DATA PRESSURE GRID

    if ~isempty(iXF) && ~isempty(iXB)
        tnan = isnan(BD(:,iXB));
        BD(tnan,:) = [];
        
        ncol = size(BD,1);
        % interp float on bottle data P
        
        Xinterp  = interp1(FD(:,iPF),FD(:,[iTF,iXF]), BD(:,iPB));
        
        % OK NOW DO LOOK UP TOO
        tol = 10; % meters, if min diff greater than this don't use
        Xlookup = ones(ncol,2) * NaN;
        for i = 1: ncol
            df = abs(FD(:,iPF) - BD(i,iPB)); % abs(float P - bottle pt P)
            t1 = df == min(df); % find min in absolute diff
            t2 = df(t1) <= tol;
            if sum(t2) == 2
                Xlookup(i,:) = mean(FD(t1,[iTF,iXF]),1);
            elseif sum(t2) == 1
                Xlookup(i,:) = FD(t1,[iTF,iXF]);
            end
        end

        %             if  sum(Xi < -5)  > 0 % sum(isnan(Ni)) > 0 ||
        %                 plot(FD(:,iXF), FD(:,iPF),'bo-', 'MarkerSize', 6)
        %                 hold on
        %                 plot(Xi,GD(t1,iP),'ko', 'MarkerSize', 4, ...
        %                     'MarkerFaceColor', 'r')
        %                 plot(GD(t1,iX),GD(t1,iP),'ko', 'MarkerSize', 6, ...
        %                     'MarkerFaceColor', 'g');
        %
        %                 set(gca, 'Ydir', 'reverse')
        %                 xlabel('campare_var{1,1}')
        %                 ylabel('Pressure')
        %                 str = [num2str(UW_ID),'  Cast: ',num2str(track(cycle,2))];
        %                 title(str)
        %                 legend('Float','interp float', 'GLODAP')
        %                 hold off
        %
        %                 pause
        %             end
        
        
        data = [data; ...
            ones(ncol,1) * UW_ID,   ...         % UW ID
            ones(ncol,1) * FD(1,2), ...         % Cycle #     
            ones(ncol,1) * FD(1,1), ...         % SDN
            ones(ncol,1) * FD(1,3), ...         % LON
            ones(ncol,1) * FD(1,4), ...         % LAT
            BD(:,[iSDNB,iPB,iXB]), ...          % bottle data
            Xinterp, Xlookup, ...               % interp & lookup float data,
            ones(ncol,1) * El, ...              % sun elevation
            ones(ncol,1) * CHL_flag, ...        % CHL type flag CHLORA or TOT_CHL_A
            ones(ncol,1) * meas_pH];
    end
end

% FINAL LOOK for MIV's
tMIV = data == -1e10;
data(tMIV) = NaN;

% GET SOME INDICES
iSDNb  = find(strcmp('Bottle SDN',hdr)    == 1); % bottle SDN
iPb  = find(strcmp('bottle pres',hdr)    == 1); % bottle pressure
iXb  = find(strcmp(compare_var{1,2},hdr) == 1); % bottle ind
iXf  = find(strcmp(['interp ',hdr_add_str],hdr)     == 1); % FLOAT interp
iXf2 = find(strcmp(['lookup ',hdr_add_str],hdr)     == 1); % FLOAT lookup
iMEAS_PH  = find(strcmp('pH meas ?',hdr)    == 1); % bottle SDN

% SET meas pH? flag to -1 for NaN
if strcmp(compare_var{1,2},'PH_TOT_INSITU')
    tnan = isnan(data(:,iXb));
    data(tnan,iMEAS_PH) = -1;
else
    data(:,iMEAS_PH) = -1;
end




% MAKE SHALLOW AND DEEP SUBSETS
surf_z = [0 30];
deep_z = [950 1550];
UW_ID = str2double(GOOD_LIST(:,2));
SURF = ones(size(UW_ID,1), size(data,2))* NaN;
DEEP = ones(size(UW_ID,1), size(data,2))* NaN;
SDhdr = hdr;
SDhdr{10} = 'Bottle - Interp float O2'; 
SDhdr{12} = 'Bottle - Lookup float O2'; 
for i = 1:size(UW_ID,1)
    td = data(:,1) == UW_ID(i);
    tmp = data(td,:);
    
    tzs = tmp(:,iPb) >= surf_z(1) & tmp(:,iPb) <= surf_z(2);
    SURF(i,:)    = nanmean(tmp(tzs,:));
    SURF(i,iXf)  = nanmean(tmp(tzs,iXb) - tmp(tzs,iXf));
    SURF(i,iXf2) = nanmean(tmp(tzs,iXb) - tmp(tzs,iXf2));
    
    tzd = tmp(:,iPb) >= deep_z(1) & tmp(:,iPb) <= deep_z(2);
    DEEP(i,:)    = nanmean(tmp(tzd,:));
    DEEP(i,iXf)  = nanmean(tmp(tzd,iXb) - tmp(tzd,iXf));
    DEEP(i,iXf2) = nanmean(tmp(tzd,iXb) - tmp(tzd,iXf2));    
end
SURF_AVG = nanmean(SURF);
SURF_STD = nanstd(SURF);
DEEP_AVG = nanmean(DEEP);
DEEP_STD = nanstd(DEEP);



% PLOT THE DATA - get indces from output hdr
% iPb  = find(strcmp('bottle pres',hdr)    == 1); % bottle pressure
% iXb  = find(strcmp(compare_var{1,2},hdr) == 1); % bottle ind
% iXf  = find(strcmp(['interp ',hdr_add_str],hdr)     == 1); % FLOAT interp
% iXf2 = find(strcmp(['lookup ',hdr_add_str],hdr)     == 1); % FLOAT lookup

range = [(min([data(:,iXb);data(:,iXf);data(:,iXf2)])), ...
    (max([data(:,iXb);data(:,iXf);data(:,iXf2)]))];

% GET REGRESSIONS
tnan = isnan(data(:,iXb)) | isnan(data(:,iXf));
[m,b,r,sm,sb] = lsqfitma(data(~tnan,iXb), data(~tnan,iXf));
reg_info1   = [m,b,r,sm,sb]; % DIFF REGRESSION OUTPUT
title_str1  = sprintf('Interp Float = %1.3f * bottle + %3.4f  R*R= %1.4f', ...
    reg_info1(1:2), reg_info1(3).^2);


tnan = isnan(data(:,iXb)) | isnan(data(:,iXf2));
[m,b,r,sm,sb] = lsqfitma(data(~tnan,iXb), data(~tnan,iXf2));
reg_info2   = [m,b,r,sm,sb]; % DIFF REGRESSION OUTPUT
title_str2  = sprintf('Lookup Float = %1.3f * bottle + %3.4f  R*R= %1.4f', ...
    reg_info2(1:2), reg_info2(3).^2);

% ************************************************************************
if fig_switch(1) == 1
figure(1)
plot(data(:,iXb),data(:,iXf), 'bo', 'MarkerSize', 2)
xlim(range)
ylim(range)
title({title_str1,title_str2})
ylabel(hdr_add_str)
x_str = ['bottle ',hdr{iXb}];
xlabel(x_str, 'interpreter', 'none')

hold on
plot(data(:,iXb),data(:,iXf2), 'go', 'MarkerSize', 2)
plot(xlim, ylim, 'k--', 'LineWidth', 2)
plot(xlim, xlim .* reg_info1(1) + reg_info1(2),'r-');
legend('interp data','lookup data', '1:1', 'interp model II', 'Location','NorthWest')
hold off
end

% ************************************************************************
if fig_switch(2) == 1
figure(2)
FBdiff  = data(:,iXb)-data(:,iXf); % interp
mFBdiff = nanmean(FBdiff);
sFBdiff = nanstd(FBdiff);

FBdiff2  = data(:,iXb)-data(:,iXf2); % lookup
mFBdiff2 = nanmean(FBdiff2);
sFBdiff2 = nanstd(FBdiff2);


title_str = sprintf(['Interp Mean diff = %2.4f ± %2.4f   ', ...
    'Lookup Mean diff = %2.4f ± %2.4f   '], ...
    mFBdiff,sFBdiff,mFBdiff2,sFBdiff2);

plot(data(:,iXb),FBdiff, 'bo', 'MarkerSize', 2)
xlim(range)
diff_range = max(abs(data(:,iXb)-data(:,iXf)));
ylim([-diff_range diff_range])
title(title_str)
ylabel('Bottle - Float')
x_str = ['bottle ',hdr{iXb}];
xlabel(x_str, 'interpreter', 'none')

hold on
plot(data(:,iXb),FBdiff2, 'go', 'MarkerSize', 2)
plot(xlim, [0 0], 'k--', 'LineWidth', 2)
plot(xlim, [mFBdiff mFBdiff], 'r-');
legend('intep data','lookup data', '0', 'mean diff', 'Location','NorthWest')
hold off
end

% ************************************************************************
if fig_switch(3) == 1
figure(3)
flts    = {'7567' '6091' '9091' '7614' '7557'};
colors = {'k' 'g' 'y' 'm' 'r'};

plot(data(:,iXb),FBdiff, 'bo', 'MarkerSize', 2)
xlim(range)
diff_range = max(abs(data(:,iXb)-data(:,iXf)));
ylim([-diff_range diff_range])
title(title_str)
ylabel('Bottle - Float')
x_str = ['bottle ',hdr{iXb}];
xlabel(x_str, 'interpreter', 'none')

hold on
for i = 1 : size(flts,2)
    t1 = data(:,1) == str2double(flts{i});
    sum(t1)
    plot(data(t1,iXb),FBdiff(t1), 'o', 'MarkerSize', 4, ...
        'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k')
end
hold off
legend(['all', flts])
end    


% ************************************************************************
if fig_switch(4) == 1

surf_str = sprintf('%0.0f to %0.0f m',surf_z(1), surf_z(2));
deep_str = sprintf('%0.0f to %0.0f m',deep_z(1), deep_z(2));
figure(4)
set(gcf, 'Position', [48.2000 164.2000 1.3544e+03 693.6000]);

% ********* SURFACE *************
ax1 = subplot(2,1,1);
x_ct = 1:size(SURF,1);
yrange = max(abs(SURF(:,iXf)));
p41 = plot(x_ct, SURF(:,iXf), 'ko', 'MarkerSize', 4, ...
    'MarkerFaceColor',[255 128 0]./256);
% ax4 = gca;
% ax4.Position = [0.050 0.1100 0.9 0.8150];

xlim([1 size(SURF,1)])
ylim([-yrange yrange])
ylabel('bottle - float, µmol/kg')
set(gca,'XTick', x_ct)
set(gca, 'XtickLabels',GOOD_LIST(:,2),'XTickLabelRotation',90)
title(['Oxygen ',surf_str, '    avg interp diff = ', ...
    num2str(SURF_AVG(iXf),'%0.2f'),' ± ', num2str(SURF_STD(iXf),'%0.2f'), ...
    '    avg lookup diff = ',num2str(SURF_AVG(iXf2),'%0.2f'),' ± ', ...
    num2str(SURF_STD(iXf2),'%0.2f'),])

hold on
p42 = plot(x_ct, SURF(:,iXf2), 'ko', 'MarkerSize', 4, ...
    'MarkerFaceColor',[51 53 255]./256);
% tnan = isnan(SURF(:,iXf));
% p43 = plot(x_ct(tnan), SURF(tnan,iXf2), 'ro', 'MarkerSize', 4, ...
%     'MarkerFaceColor','y');
p43 = plot(xlim, [SURF_AVG(iXf) SURF_AVG(iXf)], 'r--', xlim,  ...
    [SURF_AVG(iXf2) SURF_AVG(iXf2)], 'r--', 'LineWidth', 1);

p44 = plot(xlim, [0 0], 'k--', 'LineWidth', 2);
hold

%legend('interp float', 'lookup float', 'no data')
legend('interp float', 'lookup float', 'avg diff')

% ********* DEEP *************
ax2 = subplot(2,1,2);
x_ct = 1:size(DEEP,1);
yrange = max(abs(DEEP(:,iXf)));
p41 = plot(x_ct, DEEP(:,iXf), 'ko', 'MarkerSize', 4, ...
    'MarkerFaceColor',[255 128 0]./256);
% ax4 = gca;
% ax4.Position = [0.050 0.1100 0.9 0.8150];

xlim([1 size(DEEP,1)])
%ylim([-yrange yrange])
ax2.YLim = ax1.YLim;
ylabel('bottle - float, µmol/kg')
set(gca,'XTick', x_ct)
set(gca, 'XtickLabels',GOOD_LIST(:,2),'XTickLabelRotation',90)
title(['Oxygen ',deep_str, '    avg interp diff = ', ...
    num2str(DEEP_AVG(iXf),'%0.2f'),' ± ', num2str(DEEP_STD(iXf),'%0.2f'), ...
    '    avg lookup diff = ',num2str(DEEP_AVG(iXf2),'%0.2f'),' ± ', ...
    num2str(DEEP_STD(iXf2),'%0.2f'),])

hold on
p42 = plot(x_ct, DEEP(:,iXf2), 'ko', 'MarkerSize', 4, ...
    'MarkerFaceColor',[51 53 255]./256);
% tnan = isnan(DEEP(:,iXf));
% p43 = plot(x_ct(tnan), DEEP(tnan,iXf2), 'ro', 'MarkerSize', 4, ...
%     'MarkerFaceColor','y');

p43 = plot(xlim, [DEEP_AVG(iXf) DEEP_AVG(iXf)],'r--', xlim, ...
    [DEEP_AVG(iXf2) DEEP_AVG(iXf2)], 'r--', 'LineWidth', 1);

p44 = plot(xlim, [0 0], 'k--', 'LineWidth', 2);
hold

%legend('interp float', 'lookup float', 'no data')
legend('interp float', 'lookup float', 'avg diff')
end

% ************************************************************************
% ************************************************************************
%PRINT DATA TO FILE
out = data;
out(:,[3,6]) = excelsdn(out(:,[3,6]));
tmp = regexp(hdr_add_str,' ','split'); % use this to build file name


outname = ['FLOATvsBottle_',tmp{1,2}];
if QC_flag == 1
    outname =[outname,'QC'];
end
    
fid = fopen([bottle_dir, outname,'.txt'], 'w');

fprintf(fid,'File created on %s by %s\r\n', datestr(now),getenv('USERNAME'));
fprintf(fid,'Created with %s\r\n', mfilename);
fprintf(fid,'USING GOOD QC DATA ONLY !!!\r\n');

if QC_flag == 1
    fprintf(fid,'ONLY USING QC DATA WITH FLAG = 0.\r\n', tol);
else
    fprintf(fid,'USING ALL QC (GOOD AND BAD).\r\n', tol);
end

fprintf(fid,'Bottle depth range %4.0f to %4.0f meters\r\n', depth_bnds);
fprintf(fid,'Look up tolerance is %2.0f meters\r\n', tol);
fprintf(fid,['meas_pH ? flag (bottle data): 1 = measured, ', ...
             '0 = calculated, -1 = no data']);
    
% BUILD DATA FORMAT STR
df ='';
for i = 1:size(out,2)
    if i == size(out,2);
        df = [df, '%f\r\n'];
    else
        df = [df,'%f\t'];
    end
end

for i = 1:size(hdr,2)
    if i == size(hdr,2)
        fprintf(fid,'%s\r\n', hdr{i});
    else
        fprintf(fid,'%s\t', hdr{i});
    end
end

disp('printing data to file')
for i = 1: size(out,1)
    str = sprintf(df,out(i,:));    
    fprintf(fid, '%s', strrep(str,'NaN',''));
end

fclose(fid);

