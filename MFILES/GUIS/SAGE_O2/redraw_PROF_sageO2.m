function redraw_PROF_sageO2(dirs,gui,DATA,inputs)

%UPDATE:  incorporated Josh's changes to bottle lookup-table parsing
%(8/21/17)
%
% LG 9/10/2024: Made a bunch of changes to account for some specific float issues
% - Calculating the std. dev. for plot 1 now omits NaNs to account for floats with
%   mismatched resolutions between variables (1902494 wn1557)
%
% - Created a specific condition for floats with more than 10 years of GLODAP data.
%   Float 5906293 has 11+ years of GLODAP data, which goes over the number of symbols
%   we use to plot the data on plot 1. The new condition checks if there are more than
%   10 years of data, then subsets the YRS variable to only include the most recent decade.
%
% - The bottle-to-float profile plots (axes 3 and 4) would break down if a float 
%   was missing cycles 1-3. A rare case that was broken by wn1557. Plots now pull
%   data from the first two cycles in the dataset regardless of cycle number.
%
% - As an addition to the previous change, plot 3's legend now reflects the true
%   cycle number being plotted (i.e. wn1557 is missing cycles 1-3, so the first
%   available cycles are 4 and 5, which are plotted and labelled accordingly)
%
% LG 9/13/2024: Float 5906217 (ua18771) is a strange case where there IS NO BOTTLE DATA 
%    but there IS CTD DATA. The plot code would crash and not continue with plot #4, but
%    now it displays the CTD residuals and the bottle residuals becomes NaN.
%
%   -Added an if statement for calculating bottle residuals for plot 4. 
%
%   -Also added a logical indexer for the plot 3 legend. In the specific case of float 5906217,
%    the CTD O2 data was being called the bottle O2 data because the bottle O2 does not exist.
%
% LG 11/20/2024: The interp1 function for the first profile plot (axis 1) does not work if the 
%      pressure and O2 arrays being interpolated don't have 2 data points. This condition was 
%      not met when plotting O2 for float 4903026 (ss0001), causing an error. 
%    
%    - Added a test to check if the 2 variables, tmpOx(tu,1) and tmpOx(tu,2), have at least 
%      2 data points each before attempting the interp, otherwise that cycle is skipped and 
%      left as a column of NaNs.
% 
%**************************************************************************
% % Define my colors
%**************************************************************************
mycolors(1,:) = [239 90 17] ./ 255; %orange
mycolors(2,:) = [197 109 39] ./ 255; %burnt orange
mycolors(3,:) = [213 185 31] ./ 255; %mustard
mycolors(4,:) = [21 114 124] ./ 255; %teal
mycolors(5,:) = [46 161 43] ./ 255; %green
mycolors(6,:) = [53 85 40] ./ 255; %forest green
mycolors(7,:) = [24 144 243] ./ 255; %light blue
mycolors(8,:) = [6 62 137] ./ 255; %dark purple/blue

cla(gui.whichAX(1),'reset')
cla(gui.whichAX(2),'reset')
cla(gui.whichAX(3),'reset')
cla(gui.whichAX(4),'reset')
cla(gui.whichAX(5),'reset')

if isfield(DATA,'floatTYPE')
    Type = DATA.floatTYPE;
else
    if strcmp(inputs.floatTYPE,'APEX')==1 %APEX float
        Type = 'APEX';
    else %NAVIS float
        Type = 'NAVIS';
    end
end

%**************************************************************************
% AXIS 1: plot O2 profile data + GLODAP data
%**************************************************************************
D = DATA.RAWorQCprof;
z = [8,11:5:100,111:10:400,450:50:1000,1100:100:1600]';
max_z   = max(D(:,4));
z(z > max_z) =[];
r_z     = size(z,1);
ncasts = unique(D(:,2));
Xi = ones(r_z,length(ncasts)) * NaN; % PREDIMMENSION MATRICES
for icast = 1:length(ncasts)
    TMPO = [D(D(:,2)==ncasts(icast),4) D(D(:,2)==ncasts(icast),7)]; %[P O2] prof data, single cast
    tnan = isnan(TMPO(:,1));
    TMPO=TMPO(~tnan,:);
    tnan = isnan(TMPO(:,2));
    tmpOx=TMPO(~tnan,:);
    [~,tu,~] = unique(tmpOx(:,1));
    
    %LG 11/20/2024 Added a test to check if the 2 variables being interpolated below:
    % tmpOx(tu,1) and tmpOx(tu,2) have at least 2 data points each.
    %The interp1 function does not work if the arrays being interpolated don't have at
    %least 2 data points.

    %Test is the interpolation is possible based on array lengths
    interpPossible = length(tmpOx(tu,1))>=2 & length(tmpOx(tu,2))>=2;

    %If the test passes, interpolate the variables over depths z
    if interpPossible
        Xi(:,icast)  = interp1(tmpOx(tu,1), tmpOx(tu,2), z);
    
    %If interpolation IS NOT possible, then leave the final array column Xi(:,icast) as NaNs
    else
        continue;
    end
end

%LG 9/10/24 find the mean of all profiles omitting NaNs
meanprof = mean(Xi,2,"omitnan");
%LG 9/10/24 do the same for std, this line was running into issues with NAVIS floats
%which may have more NaNs without O2 due to differences in pres axes
stdprof = std(Xi,0,2, "omitnan");
p1(1)=plot(z,meanprof,'b','linewidth',2,'Parent',gui.whichAX(1));
hold(gui.whichAX(1),'on')
p1(2)=plot(z,meanprof-stdprof,'b--','Parent',gui.whichAX(1));
ylabel(gui.whichAX(1),'[O_2] (umol/kg)','fontunits','normalized','fontsize',0.125);
if inputs.rorq == 1 % RAW TAB SELECTED?
    title(gui.whichAX(1),[Type,' Float ',inputs.MBARI_ID,', (',inputs.WMO_ID,')'],'fontunits','normalized','fontsize',0.15)
else
    title(gui.whichAX(1),[Type,' Float ',inputs.MBARI_ID,', (',inputs.WMO_ID,'). QC Adjustments applied.'],'fontunits','normalized','fontsize',0.15)
end
%     XLIMS=get(gui.whichAX(1),'xlim')
XLIMS=[inputs.depthedit(1,1) inputs.depthedit(1,2)];
% Get GLODAP data
track = [DATA.track(:,1:2) DATA.track(:,4) DATA.track(:,3)]; %switch lat/lon columns
G = get_GLODAPv2_local_sO2(track,inputs.GLDPkm,[0 2000],dirs.user_dir);
GLODAP_color = [0   255 255;... % keep color and symbol count equal
                153  51 255;
                255 51 255;
                0 255 0;
                1 1 0;
                1 0 0;
                255 128 0]./255;
GLODAP_symbols = {'^' 'd' 's' '*' '<' 'x' 'p'};
iGcyc = find(strcmp('float cycle',  G.hdr) == 1);
iGSDN = find(strcmp('Date',  G.hdr)        == 1);
iGP   = find(strcmp('G2pressure',G.hdr)    == 1);
iGT   = find(strcmp('G2temperature',G.hdr) == 1);
iGS   = find(strcmp('G2salinity',G.hdr)    == 1);
iGO   = find(strcmp('G2oxygen',G.hdr)      == 1);
iGN   = find(strcmp('G2nitrate',G.hdr)     == 1);
iGPH  = find(strcmp('ph_insitu',G.hdr)     == 1);
GIND = iGO;  %For Oxygen GUI (keep other parameter indices for now)
% mean glodap
GG = [G.data(:,iGP) G.data(:,GIND)];
[c,ia,ic] = unique(G.data(:,iGP));
GG2 = GG(ia,:);
hlegend_cell = {'Mean profile','1 std'};
if ~isempty(GG2) % GLODAP data exists.  plot it
    [~,I] = sort(GG2(:,1));
    G2 = GG2(I,:);
    xg = isnan(G2(:,2));
    G3 = G2(~xg,:);
    Gd = interp1(G3(:,1), G3(:,2),z);
    hlegend_cell = [hlegend_cell, {'Mean GLODAP'}];
    hold(gui.whichAX(1),'on')

    %Plot mean GLODAP
    plot(z,Gd,'r','Parent',gui.whichAX(1),'linewidth',2);
    if ~isempty(GIND)
        [G_YRS,~,~,~,~,~] = datevec(G.data(:,iGSDN));% Years only
        YRS = unique(G_YRS(~isnan(G_YRS)));
        ct = 0; ct2 = 1;

        %LG 9/10/24 rare issue with floats that have more than 10 years of float data (5906293),
        %loop will "run out" of GLODAP_symbols and crash. Adding this condition
        %where if there are more than 10 years, only use the most recent decade.
        if length(YRS)>10
            %If the number of YRS is more than 10, only use the final ten years in the series
            YRS = YRS(end-10:end);
        end
        
        %Iterate over the numbers of years
        for i = 1: size(YRS,1)
            ct = ct+1;
            %Create a year index "t1" to isolate GLODAP data from struct "G"
            t1 = G_YRS == YRS(i);
            hold(gui.whichAX(1),'on')
            %Add GLODAP spot samples
            vp = plot(G.data(t1,iGP),G.data(t1,GIND),'parent',gui.whichAX(1), ...
                'Linestyle', 'none', ...
                'Marker', GLODAP_symbols{ct2}, ...
                'MarkerFaceColor', GLODAP_color(ct,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 4);
            hlegend_cell = [hlegend_cell, ...
                {['GLODAP ',num2str(YRS(i))]}];
            %Once the color counter "ct" gets to 4, start iterating the sybmol
            %counter "ct2" by one.
            if ct == 4
                ct = 0;
                ct2 = ct2+1;
            end
            
        end
    end
else
    disp('NO GLODAP DATA EXISTS IN RANGE DISTANCE.')
end
hold(gui.whichAX(1),'on')
plot(z,meanprof+stdprof,'b--','Parent',gui.whichAX(1)); %plot second errorbar last, to exclude from legend
set(gui.whichAX(1),'xticklabels',[],'xlim',XLIMS)
hlegend = legend(gui.whichAX(1),hlegend_cell);
posL=get(hlegend,'Position');
posL(1) = 1.35*posL(1);
posL(2) = 0.75*posL(2);
set(hlegend,'Position',posL,'fontunits','normalized','fontsize',10);
grid(gui.whichAX(1),'on')

%**************************************************************************
% AXIS 2: plot resid
%**************************************************************************
if ~isempty(GG2)
    FtoG = interp1(z,meanprof,GG2(:,1));
    GLO_resid = GG2(:,2)-FtoG;
    Ggain = GG2(:,2)./FtoG;
    %         meandiffG = nanmean(GLO_resid);
    %         meanGgain = nanmean(Ggain);
    plot(GG2(:,1),GLO_resid,'o','Parent',gui.whichAX(2),...
        'MarkerSize',4, 'MarkerFaceColor','k',...
        'MarkerEdgeColor', 'k');
    set(gui.whichAX(2),'xlim',XLIMS,'xticklabels',[])
    meandiffG = mean(GLO_resid(GG2(:,1)>=XLIMS(1,1) & GG2(:,1)<=XLIMS(1,2)),"omitnan");
    meanGgain = mean(Ggain(GG2(:,1)>=XLIMS(1,1) & GG2(:,1)<=XLIMS(1,2)),"omitnan");
    ylabel(gui.whichAX(2),'GLODAP-Flt (umol/kg)','fontunits','normalized','fontsize',0.125);
    Xls = get(gui.whichAX(2),'xlim');
    Yls = get(gui.whichAX(2),'ylim');
    Xloc = Xls(2)+(Xls(2)-Xls(1))/15;
    Yloc = (Yls(2)-Yls(1))/8;
    text(Xloc,Yls(1)+6*Yloc,['Mean Resid = ',num2str(meandiffG)],'parent',gui.whichAX(2),'fontunits','normalized','fontsize',0.125)
    text(Xloc,Yls(1)+5*Yloc,['Mean Gain = ',num2str(meanGgain)],'parent',gui.whichAX(2),'fontunits','normalized','fontsize',0.125)
else
    text(0.075,0.5, 'NO GLODAP DATA EXISTS WITHIN DEFINED DISTANCE.  EXPAND RANGE!','parent',gui.whichAX(2),'fontunits','normalized','fontsize',0.1,'fontweight','bold','color','r')
end
grid(gui.whichAX(2),'on')
%**************************************************************************
% AXIS 3: plot O2 profile data + bottle data
%**************************************************************************
%LG 9/10/2024 commented this section out
%plot cast 1 and 2 with bottle data
% CAST1 = [D(D(:,2)==1,4) D(D(:,2)==1,7)]; %[P O2] prof data, single cast
% CAST2 = [D(D(:,2)==2,4) D(D(:,2)==2,7)]; %[P O2] prof data, single cast
% CAST3 = [D(D(:,2)==3,4) D(D(:,2)==3,7)]; %[P O2] prof data, single cast

% Float 7567 missing first cast; quick fix, if cast 1 missing --> use
% cast 2 and 3.  Does not fix potential case where cast 2 and/or 3 also missing!
% if isempty(CAST1)
%     CAST1 = CAST2;
%     CAST2 = CAST3;
% end

%LG 9/10/24 Implementing a system that does not rely on having a 1st, 2nd, or 3rd, but instead
% takes the first three cycles based on the unique cast # (ncasts). Float wn1557 is 
% missing first 3, not good but trying out a more dynamic fix.
CAST1 = [D(D(:,2)==ncasts(1),4) D(D(:,2)==ncasts(1),7)]; %[P O2] prof data, single cast
CAST2 = [D(D(:,2)==ncasts(2),4) D(D(:,2)==ncasts(2),7)]; %[P O2] prof data, single cast

plot(CAST1(:,1),CAST1(:,2),'b','linewidth',2,'Parent',gui.whichAX(3));
hold(gui.whichAX(3),'on')
plot(CAST2(:,1),CAST2(:,2),'b--','Parent',gui.whichAX(3));
set(gui.whichAX(3),'xlim',XLIMS,'xticklabels',[])
ylabel(gui.whichAX(3),'[O_2] (umol/kg)','fontunits','normalized','fontsize',0.125);

b=DATA.bdata;
ibSDN = find(strcmp('DATE',  b.hdr) == 1);
ibZ   = find(strcmp('DEPTH', b.hdr) == 1);
ibP   = find(strcmp('CTDPRS',b.hdr) == 1);
ibT   = find(strcmp('CTDTMP',b.hdr) == 1);
ibS   = find(strcmp('CTDSAL',b.hdr) == 1);
ibO   = find(strcmp('OXYGEN',b.hdr) == 1);
ibOCTD = find(strcmp('CTDOXY',b.hdr) == 1);
ibN   = find(strcmp('NITRAT',b.hdr) == 1);
ibPH1  = find(strcmp('PH_TOT_INSITU',b.hdr) == 1);
ibPH2  = find(strcmp('PH_TOT_INSITU_ALKDIC',b.hdr) == 1);

bIND = ibO; %For Oxygen GUI (keep other parameter indices for now)
bIND2 = ibOCTD; % For Oxygen GUI; plot CTD oxygen in addition to Winklers
if ~isempty(bIND)
    hold(gui.whichAX(3),'on')
    plot(b.data(:,ibP),b.data(:,bIND),'o','parent',gui.whichAX(3),...
        'MarkerSize',6, 'MarkerFaceColor',[255 255 0]/255,...
        'MarkerEdgeColor', 'k')
    %LG 9/13/2024 Adding a logical test for whether the float has bottle O2 or not, which is
    %used to create the legend later on.
    hasBottleO2 = true;
else
    %If float doesn't have bottle O2, the test is set to false
    hasBottleO2 = false;
end
if ~isempty(bIND2)
    hold(gui.whichAX(3),'on')
    plot(b.data(:,ibP),b.data(:,bIND2),'o','parent',gui.whichAX(3),...
        'MarkerSize',6, 'MarkerFaceColor',[255 255 0]/255,...
        'MarkerEdgeColor', 'r')
    %Same test as above except this is for CTD O2 data
    hasCTDO2 = true;
else
    %If float doesn't have CTD O2 data, set this test as false
    hasCTDO2 = false;
end
set(gui.whichAX(3),'xticklabels',[],'xlim',XLIMS)
%LG changing the legend labels since first cast may not be cycle 1 or 2.
%Create the cell array of legend labels for plot 3
bot_leg_cell = {"CAST "+ ncasts(1),"CAST "+ ncasts(2),'Bottle O2', 'CTD O2'};

%Create a logical array that indicates whether certain labels are used (i.e. if the float
%doesn't have bottle data but has CTD data, the label "Bottle O2" will not be used, but 
%"CTD O2" will be used.
bot_leg_test = [true true hasBottleO2 hasCTDO2];

%Apply the legend to the plot.
bot_leg = legend(gui.whichAX(3), bot_leg_cell{bot_leg_test});
posL=get(bot_leg,'Position');
posL(1) = 1.25*posL(1);
posL(2) = 0.75*posL(2);
set(bot_leg,'Position',posL,'fontunits','normalized','fontsize',10);
grid(gui.whichAX(3),'on')

%**************************************************************************
% AXIS 4: plot resid
%**************************************************************************
xlabel(gui.whichAX(4),'Depth (m)','fontunits','normalized','fontsize',0.125);
set(gui.whichAX(5),'Visible','Off');
[~,XX,~] = unique(CAST1(:,1));
CAST1=CAST1(XX,:);
nonans = isnan(CAST1(:,2));
C1 = CAST1;
C1(nonans,:) = [];
FtoB = interp1(C1(:,1),C1(:,2),b.data(:,ibP));
bot_resid = b.data(:,bIND)-FtoB;
bot_gain = b.data(:,bIND)./FtoB;
ctd_resid = b.data(:,bIND2)-FtoB;
ctd_gain = b.data(:,bIND2)./FtoB;

%LG 9/13/2024 Adding this as an if statement in case a float DOES NOT have bottle data but
%DOES have CTD data (5906217 ua18771). Before, the plot code would crash and not continue 
%with plot #4, now it displays the CTD residuals and the bottle residuals becomes NaN.
if ~isempty(bot_resid)
    plot(b.data(:,ibP),bot_resid,'o','parent',gui.whichAX(4),...
        'MarkerSize',4, 'MarkerFaceColor','k',...
        'MarkerEdgeColor', 'k');
    meandiffB = mean(bot_resid(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)),'omitnan');
    meanBgain = mean(bot_gain(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)),'omitnan');
    hold(gui.whichAX(4),'on')
else
    meandiffB = NaN;
    meanBgain = NaN;
end
if ~isempty(ctd_resid)
    plot(b.data(:,ibP),ctd_resid,'o','parent',gui.whichAX(4),...
        'MarkerSize',4, 'MarkerFaceColor','r',...
        'MarkerEdgeColor', 'r');
    meandiffCTD = mean(ctd_resid(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)),'omitnan');
    meanCTDgain = mean(ctd_gain(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)),'omitnan');
else
    meandiffCTD = nan;
    meanCTDgain = nan;
end
set(gui.whichAX(4),'xlim',XLIMS)

%LG 9/13/2024 Commented out so these statements could be added to the "if ~isempty(bot_resid)" statement
% meandiffB = nanmean(bot_resid(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)));
% meanBgain = nanmean(bot_gain(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)));

ylabel(gui.whichAX(4),'Bottle-Flt (umol/kg)','fontunits','normalized','fontsize',0.125);
xlabel(gui.whichAX(4),'Pressure (db)','fontunits','normalized','fontsize',0.125);
Xls = get(gui.whichAX(4),'xlim');
Yls = get(gui.whichAX(4),'ylim');
Xloc = Xls(2)+(Xls(2)-Xls(1))/15;
Yloc = (Yls(2)-Yls(1))/8;
text(Xloc,Yls(1)+6*Yloc,[char(956),' Bottle Resid = ',num2str(meandiffB)],'parent',gui.whichAX(4),'fontunits','normalized','fontsize',0.125)
text(Xloc,Yls(1)+5*Yloc,[char(956),' Bottle Gain = ',num2str(meanBgain)],'parent',gui.whichAX(4),'fontunits','normalized','fontsize',0.125)
text(Xloc,Yls(1)+4*Yloc,[char(956),' CTD Resid = ',num2str(meandiffCTD)],'parent',gui.whichAX(4),'fontunits','normalized','fontsize',0.125,'color','r')
text(Xloc,Yls(1)+3*Yloc,[char(956),' CTD Gain = ',num2str(meanCTDgain)],'parent',gui.whichAX(4),'fontunits','normalized','fontsize',0.125,'color','r')
grid(gui.whichAX(4),'on')


end % redraw_PROF_sageO2
