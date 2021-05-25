function Plot_sensor_coverage_phver(PHverNum,save_dir)
% Plot_sensor_coverage.m
%
%--------------------------------------------------------------------------
%
% Program to plot coverage maps showing lat/lon of all SOCCOM profiles,
% color-coded for working vs failing sensors.
%
% INPUTS:  PHverNUM = 0, 1, 1.1, 1.2, 2, 3, 4
%          save_dir ='path\to\image\dir\'
% Tanya Maurer
% MBARI
% 4/25/18
% Modified 5/12/20 to plot pH data by sensor version type
%--------------------------------------------------------------------------
if PHverNum == 1.1
    phplottag = '1pt1';
elseif PHverNum == 1.2
    phplottag = '1pt2';
else
    phplottag = num2str(PHverNum);
end
sensor_type = 'PH'; %This adapted code is pH specific
% save_dir = 'C:\Users\tmaurer\Documents\FLOAT_QC\sensor_coverage_maps\'
% sensor_type = 'CHL'
%
% LOAD PH VERSION NUMBERS
load X:\DuraFET\Diagnostic_data\DATA\pump_offset.mat
UW = out(:,1);
PHver = out(:,9);
[uw,IA,IC] = unique(UW);
ph_ver = PHver(IA);
% 7 cateogries of pH sensor, essentially
% pH version number is ph_ver
    % 1 Ag wire, thinner isfet cover
    % 1.1 Ag wire, thicker isfet cover
    % 1.2 Pt wire, thicker isfet cover
    % 2 Double O-ring, Pt wire, thicker isfet cover
    % 3 Smaller Pt feed holes, double oring, Pt wire, thicker isfet cover
    % 4 Roll pin, smaller Pt feed holes, double oring, Pt wire, thicker isfet cover
    % 0 sbe stem
Xver = find(ph_ver == PHverNum);
uwSub = uw(Xver);


load C:\Users\bgcargo\Documents\MATLAB\ARGO_MBARI2AOML\float_type_lists.mat
SOCfloats = SOCCOM_list(:,1);
NonSOCfloats = NON_SOCCOM_list(:,1);
allflt = [SOCCOM_list;NON_SOCCOM_list];
ALLfloats = allflt(:,1);
mylist = SOCCOM_list(:,1);
myfloats = str2num(char(SOCCOM_list(:,2)));
DEADlist = {'0508SOOCNQC','0509SOOCNQC','0511SOOCNQC','0564SOOCNQC','7557SOOCNQC','7567SOOCNQC','8514SOOCNQC','9031SOOCNQC','9095SOOCNQC','9761SOOCNQC','0693SOOCNQC','0506SOOCNQC','12390SOOCNQC','9260SOOCNQC','9254SOOCNQC'};
DEADlist = {};
K = 1;
NS = 1;
TOTBAD = 0;
TOTGOOD =0;
% for kk = 1:length(mylist)
for kk = 1:length(uwSub)
    thflt = uwSub(kk);
    xf = find(myfloats == thflt);
    if isempty(xf)
        continue
    end
    FLTname = mylist(xf);
% for kk = 1:10
    fname = [char(FLTname),'QC'];
    if strcmp(fname(1),'0') == 1
        flttype = 'NAVIS';
    else
        flttype = 'APEX';
    end
    disp(['EXTRACTING LOCATION DATA AND SENSOR HEALTH INFO FROM ',fname,'...'])
    % Load data------------------------------------------------------------
    d = get_FloatViz_data(fname);
    if isempty(d)
        disp(['get_FloatViz_data returned empty for ',fname])
        continue
    end
    hdr = d.hdr;
    FV  = d.data;
    MVI = -1e10;
    FV(FV == MVI) = NaN; 
    clear d
    
    % Get indices----------------------------------------------------------
    iSTA = find(strcmp('Station',hdr)          == 1);
    iSDN = find(strcmp('SDN',hdr)              == 1);
    iLAT = find(strcmp('Lat [°N]',hdr)              == 1);
    iLON = find(strcmp('Lon [°E]',hdr)              == 1);
    iP   = find(strcmp('Pressure[dbar]',hdr)   == 1);
    iT   = find(strcmp('Temperature[°C]',hdr)  == 1);
    iS   = find(strcmp('Salinity[pss]',hdr)    == 1);
    iC   = find(strcmp('Chl_a[mg/m^3]',hdr)    == 1);
    iB   = find(strcmp('b_bp700[1/m]',hdr)     == 1);
    iO   = find(strcmp('Oxygen[µmol/kg]',hdr)  == 1);
    iN   = find(strcmp('Nitrate[µmol/kg]',hdr) == 1);
    ipH  = find(strcmp('pHinsitu[Total]',hdr)  == 1);
    
    % Get sensor type------------------------------------------------------
    if strcmp(sensor_type,'N') == 1
        ind = iN;
        plottag = 'NO3';
    elseif strcmp(sensor_type,'PH') == 1
        ind = ipH;
        plottag = 'PH';
    elseif strcmp(sensor_type,'O') == 1
        ind = iO;
        plottag = 'O2';
    elseif strcmp(sensor_type,'CHL') == 1
        ind = iC;
        plottag = 'CHL';
    elseif strcmp(sensor_type,'B') == 1
        ind = iB;
        plottag = 'BBP';
    end
    if isempty(ind) %if no "sensor_type" sensor, skip to next float
        NOsensor{NS} = fname;
        NS = NS+1;
        continue
    end
    
    % Extract Data---------------------------------------------------------
    cycles = unique(FV(:,iSTA),'stable');
    k = 1;
    LAT=[];
    STN=[];
    LON=[];
    ISBAD=[];
    for j = 1:length(cycles)
        tmp = FV(FV(:,iSTA)==cycles(j),:);
        n = find(tmp(:,ind+1)>4);
        if length(n)./size(tmp,1)*100 > 50
            isbad = 1;
            totbad = 1;
            totgood = 0;
        else
            isbad = 0;
            totgood = 1;
            totbad = 0;
        end
        LAT(k) = nanmean(tmp(:,iLAT));
        LON(k) = nanmean(tmp(:,iLON));
        STN(k) = j;
        ISBAD(k) = isbad;
        TOTBAD = TOTBAD+totbad;
        TOTGOOD = TOTGOOD+totgood;
        k=k+1;
    end
    FLT = [STN' LAT' LON' ISBAD'];
    FLOATS{K} = FLT;
    FLOATNAME{K} = fname;
    K=K+1;
end
    
% Make plot of full float track--------------------------------------------
H = figure('color','white','Position',[200 -200 1000 700]);
cmp = load('SOCCOM_NCP_cmap.mat');
CMAP = cmp.cmap;
colormap(CMAP);
legend_cell = {};
hold on
m_proj('stereographic','latitude',-90,'radius',60);
m_tbase('contourf','edgecolor','none');
m_gshhs_i('color','k','linewidth',1.5);              % Coastline...
set(gca,'ydir','normal');
m_grid('linewidth',2,'tickdir','out','xaxislocation','top');
for ii = 1:length(FLOATS)
    theflt = FLOATS{ii};
    hold on
    h2=m_line(theflt(:,3),theflt(:,2),'marker','o','color','k',...
              'markersize',3,'markerfacecolor','k');
    if ii == 1;
    hold on
    H2=m_plot(theflt(:,3),theflt(:,2),'o','color','r','markersize',3,'markerfacecolor','r','markeredgecolor','r');
    hold on
    H3=m_plot(theflt(:,3),theflt(:,2),'o','color','g','markersize',3,'markerfacecolor','g','markeredgecolor','g');
    end
    thebad = find(theflt(:,4)==1);
    thegood = find(theflt(:,4)==0);
    hold on
%     if ~isempty(thebad) && ~isempty(thegood)
%         h2a=m_plot(theflt(thebad,3),theflt(thebad,2),'o','color','r','markersize',3,'markerfacecolor','r','markeredgecolor','r');
%         hold on
%         h2b=m_plot(theflt(thegood,3),theflt(thegood,2),'o','color','g','markersize',3,'markerfacecolor','g','markeredgecolor','g');
%     else %for legend, need both h2a and h2b to be un-empty
        HA=m_plot(theflt(thebad,3),theflt(thebad,2),'o','color','r','markersize',3,'markerfacecolor','r','markeredgecolor','r');
        hold on
        HB=m_plot(theflt(thegood,3),theflt(thegood,2),'o','color','g','markersize',3,'markerfacecolor','g','markeredgecolor','g');    
%     end
end
BAD = ['BAD ',plottag,', ',num2str(TOTBAD),' profiles'];
GOOD = ['GOOD ',plottag,', ',num2str(TOTGOOD),' profiles'];
% BAD = [plottag,' Profile, QF=BAD'];
% GOOD = [plottag,' Profile, QF=GOOD'];
legend_cell = [legend_cell, BAD, GOOD];
track_legend = legend([H2 H3],legend_cell,'location','northeast');
title({['SOCCOM Float Tracks, pHv',num2str(PHverNum)],' '})
set(track_legend,'fontsize',14)
set(gca,'fontsize',18)
POS = get(track_legend,'position');
P = POS; %keep original
P(1,1) = P(1,1)*1.35;
set(track_legend,'position',P)
set(gcf,'PaperPositionMode','auto')
print([save_dir,'SOCCOM_floattracks_',plottag,'sensorhealth_pHv',phplottag,'.png'],'-dpng','-r0')

% Make plot of just active floats (last location and if good or bad)-------
H = figure('color','white','Position',[200 -200 1000 700]);
cmp = load('SOCCOM_NCP_cmap.mat');
CMAP = cmp.cmap;
colormap(CMAP);
legend_cell = {};
hold on
m_proj('stereographic','latitude',-90,'radius',60);
m_tbase('contourf','edgecolor','none');
m_gshhs_i('color','k','linewidth',1.5);              % Coastline...
set(gca,'ydir','normal');
m_grid('linewidth',2,'tickdir','out','xaxislocation','top');
whichfirst = [];
k2 = 1;
k3 = 1;
for ii = 1:length(FLOATS)
    tmpflt = FLOATS{ii};
    fltname = FLOATNAME{ii};
    IndexC = strfind(DEADlist, fltname);
    Index = find(not(cellfun('isempty', IndexC)));
    if ~isempty(Index)
        disp([fltname,' is a dead float.  Not plotting on "Active" float map.'])
        continue
    end
    xx = find(~isnan(tmpflt(:,2)));
    theflt = tmpflt(xx,:); % avoid having last position NaN (rare)
    lastloc = theflt(end,:);
    hold on
    HH2 = m_plot(0,0,'o','color','r','markersize',4,'markerfacecolor','r','markeredgecolor','r');
    hold on
    HH3 = m_plot(0,0,'o','color','g','markersize',4,'markerfacecolor','g','markeredgecolor','g');
    hold on
    if lastloc(1,end) == 0 %good
        goodflts{k3} = fltname;
        k3=k3+1;
        if strcmp(FLOATNAME{ii}(1:3),'127')==1 && strcmp(sensor_type,'PH') == 1 %127?? series are newer pH model
%             H2 =
%             m_plot(lastloc(:,3),lastloc(:,2),'^','color','g','markersize',4,'markerfacecolor','g','markeredgecolor','g');%
%             just keep the same symbol for now.  keep code here in case I
%             want to change.
%             
            H2 = m_plot(lastloc(:,3),lastloc(:,2),'o','color','g','markersize',4,'markerfacecolor','g','markeredgecolor','g'); 
        else
           H2 = m_plot(lastloc(:,3),lastloc(:,2),'o','color','g','markersize',4,'markerfacecolor','g','markeredgecolor','g');
        end
        hold on
        m_text(lastloc(:,3),lastloc(:,2),fltname(1:end-7),'horizontalalignment','right','fontsize',7)            
    end
    if lastloc(1,end) == 1 %bad
        badflts{k2} = fltname;
        k2=k2+1;
        if strcmp(FLOATNAME{ii}(1:3),'127')==1 && strcmp(sensor_type,'PH') == 1 %127?? series are newer pH model
%             H2 = m_plot(lastloc(:,3),lastloc(:,2),'^','color','r','markersize',4,'markerfacecolor','r','markeredgecolor','r');
            H2 = m_plot(lastloc(:,3),lastloc(:,2),'o','color','r','markersize',4,'markerfacecolor','r','markeredgecolor','r');
        else
            H2 = m_plot(lastloc(:,3),lastloc(:,2),'o','color','r','markersize',4,'markerfacecolor','r','markeredgecolor','r');
        end
        hold on
        m_text(lastloc(:,3),lastloc(:,2),fltname(1:end-7),'horizontalalignment','right','fontsize',7)            
    end
end
BAD = ['BAD ',plottag];
GOOD = ['GOOD ',plottag];
legend_cell = [legend_cell, BAD, GOOD];
track_legend = legend([HH2 HH3],legend_cell);
% title({'Active SOCCOM Floats, most recent location','(triangles depict floats with new pH sensor design)',' '})
title({['SOCCOM Floats pHv',num2str(PHverNum)],'(most recent location)',' '})
set(track_legend,'fontsize',14)
set(gca,'fontsize',18)
POS = get(track_legend,'position');
P = POS; %keep original
P(1,1) = P(1,1)*1.35;
set(track_legend,'position',P)
set(gcf,'PaperPositionMode','auto')
print([save_dir,'SOCCOM_current_',plottag,'sensorhealth_pHv',phplottag,'.png'],'-dpng','-r0')
    
TOTBAD
TOTGOOD