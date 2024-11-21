function redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)


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
mycolors(9,:) = [247 30 237] ./ 255; %hot pink

cla(gui.whichAX(1),'reset')
cla(gui.whichAX(2),'reset')
cla(gui.whichAX(3),'reset')
cla(gui.whichAX(4),'reset')
cla(gui.whichAX(5),'reset')


if isfield(DATA,'floatTYPE')
    Type = DATA.floatTYPE;
else
    if strcmp(inputs.floatTYPE,'APEX')==1; %APEX float
        Type = 'APEX';
    else %NAVIS float
        Type = 'NAVIS';
    end
end

%**************************************************************************
% AXIS 1: plot pO2
%**************************************************************************
cla(gui.whichAX(1),'reset')
plot(DATA.refsub(:,2),DATA.refsub(:,4),'r-','Parent',gui.whichAX(1),...
    'markersize',20,'linewidth',2);
hold(gui.whichAX(1),'on')
if ~isempty(DATA.airsub{1}{1})
    errorbar(DATA.airsub{1}{1}(:,2),DATA.airsub{1}{1}(:,9),DATA.airsub{2}{1}(:,9),...
        'color',mycolors(7,:),'Parent',gui.whichAX(1),'linewidth',2); %AIRold
end
if DATA.howmanyairs >1 && ~strcmp(DATA.floatTYPE,'SOLO')
    hold(gui.whichAX(1),'on')
    errorbar(DATA.airsub{1}{3}(:,2),DATA.airsub{1}{3}(:,9),DATA.airsub{2}{3}(:,9),...
        'color',mycolors(8,:),'Parent',gui.whichAX(1),'linewidth',2);%AIRnew(air)
    lhandle = legend(gui.whichAX(1),DATA.reftag,'AIRold','AIRnew');
elseif DATA.howmanyairs > 1 && strcmp(DATA.floatTYPE,'SOLO')
    errorbar(DATA.airsub{1}{3}(:,2),DATA.airsub{1}{3}(:,9),DATA.airsub{2}{3}(:,9),...
        'color',mycolors(7,:),'Parent',gui.whichAX(1),'linewidth',2); %AIRold
    lhandle = legend(gui.whichAX(1),DATA.reftag,'AIRobs');
else
    lhandle = legend(gui.whichAX(1),DATA.reftag,'AIRobs');
end
posL=get(lhandle,'Position');
posL(1) = 1.25*posL(1);
posL(2) = 0.75*posL(2);
set(lhandle,'Position',posL,'fontunits','normalized','fontsize',12);
axis(gui.whichAX(1),'tight')
xlim(gui.whichAX(1),DATA.xlims{2});
ylabel(gui.whichAX(1),'pO_2 (hPa)','fontunits','normalized','fontsize',0.125);
set(gui.whichAX(1),'xtick',DATA.xticks{2});
if inputs.rorq == 1 % RAW TAB SELECTED?
    title(gui.whichAX(1),[Type,' Float ',inputs.MBARI_ID,', (',inputs.WMO_ID,')'],'fontunits','normalized','fontsize',0.15)
else
    title(gui.whichAX(1),[Type,' Float ',inputs.MBARI_ID,', (',inputs.WMO_ID,'). QC Adjustments applied.'],'fontunits','normalized','fontsize',0.15)
end
grid(gui.whichAX(1),'on')

%**************************************************************************
% AXIS 2: plot gain
%**************************************************************************
cla(gui.whichAX(2))
cla(gui.whichAX(2),'reset')
%     x = DATA.refsub(:,1);
if ~isempty(DATA.GAINS{1})
    x = DATA.GAINStime{1};
    y1 = DATA.GAINS{1};
    meanGAIN1 = nanmean(y1);
    X = x(~isnan(y1));
    Y = y1(~isnan(y1));
    [my1,by1,~,~,~]=lsqfity(X,Y);
    plot(x,y1,'o','parent',gui.whichAX(2),'markersize',5,'markerfacecolor',mycolors(7,:),'markeredgecolor',mycolors(7,:)); %AIRold
    xlim(gui.whichAX(2),DATA.xlims{1});
    Xls = get(gui.whichAX(2),'xlim');
    Yls = get(gui.whichAX(2),'ylim');
    hold(gui.whichAX(2),'on')
    newY = Xls.*my1+by1;
    %     plot(Xls,newY,'--','parent',gui.whichAX(2),'color',mycolors(7,:)) single regression line.  replaced with line showing mean
    plot(Xls,repmat(meanGAIN1,length(Xls)),'--','parent',gui.whichAX(2),'color',mycolors(7,:))
    hold(gui.whichAX(2),'on')
    plot(Xls,newY,'--','parent',gui.whichAX(2),'color',mycolors(7,:),'linewidth',2)
    ysies = get(gui.whichAX(2),'ylim');
        if inputs.rorq==1
        for jt = 1:size(DATA.tableDATA,1)
            %                 Ex = DATA.refsub(DATA.refsub(:,2)>=DATA.tableDATA(jt,1)&DATA.refsub(:,2)<=inputs.cyEND(jt),1);
            Ex = DATA.refdata(DATA.refdata(:,2)>=DATA.tableDATA(jt,1)&DATA.refdata(:,2)<=inputs.cyEND(jt),1);
            Exes = Ex-nanmin(Ex);
            Whys = Exes.*(DATA.tableDATA(jt,3))./365+DATA.tableDATA(jt,2);
            hold(gui.whichAX(2),'on')
            plot(Ex,Whys,'-','parent',gui.whichAX(2),'color',mycolors(9,:),'linewidth',2)
            hold(gui.whichAX(2),'on')
            myEx = DATA.refsub(DATA.refsub(:,2)==DATA.tableDATA(jt,1),1);
            %                 ysies = get(gui.whichAX(2),'ylim');
            if ~isempty(myEx)
                plot([myEx myEx],ysies,'--','parent',gui.whichAX(2),'color',mycolors(3,:))
            end
        end
    
end
if DATA.howmanyairs >1 && ~strcmp(DATA.floatTYPE,'SOLO')
    x2 = DATA.GAINStime{3};
    y2 = DATA.GAINS{3};
    meanGAIN2 = nanmean(y2);
    hold(gui.whichAX(2),'on')
    plot(x2,y2,'o','parent',gui.whichAX(2),'markersize',5,'markerfacecolor',mycolors(8,:),'markeredgecolor',mycolors(8,:));%AIRnew(air)
    X = x2(~isnan(y2));
    Y = y2(~isnan(y2));
    [my2,by2,~,~,~]=lsqfity(X,Y);
    hold(gui.whichAX(2),'on')
    newY = Xls.*my2+by2;
    %         plot(Xls,newY,'--','parent',gui.whichAX(2),'color',mycolors(8,:)) single regression line.  replaced with line showing mean
    plot(Xls,repmat(meanGAIN2,length(Xls)),'--','parent',gui.whichAX(2),'color',mycolors(8,:))
    % Plot Breakpoint analysis:
% % %     if inputs.rorq==1
% % %         for jt = 1:size(DATA.tableDATA,1)
% % %             %                 Ex = DATA.refsub(DATA.refsub(:,2)>=DATA.tableDATA(jt,1)&DATA.refsub(:,2)<=inputs.cyEND(jt),1);
% % %             Ex = DATA.refdata(DATA.refdata(:,2)>=DATA.tableDATA(jt,1)&DATA.refdata(:,2)<=inputs.cyEND(jt),1);
% % %             Exes = Ex-nanmin(Ex);
% % %             Whys = Exes.*(DATA.tableDATA(jt,3))./365+DATA.tableDATA(jt,2);
% % %             hold(gui.whichAX(2),'on')
% % %             plot(Ex,Whys,'--','parent',gui.whichAX(2),'color',mycolors(9,:),'linewidth',2)
% % %             myEx = DATA.refsub(DATA.refsub(:,2)==DATA.tableDATA(jt,1),1);
% % %             ysies = get(gui.whichAX(2),'ylim');
% % %             if ~isempty(myEx)
% % %                 hold(gui.whichAX(2),'on')
% % %                 plot([myEx myEx],ysies,'--','parent',gui.whichAX(2),'color',mycolors(3,:))
% % %             end
% % %         end
% % %     end
    set(gui.whichAX(2),'ylim',ysies)
    Xls = get(gui.whichAX(2),'xlim');
    Yls = get(gui.whichAX(2),'ylim');
    Xloc = Xls(2)+(Xls(2)-Xls(1))/10;
    Yloc = (Yls(2)-Yls(1))/8;
    text(Xloc,Yls(1)+7*Yloc,['mean  = ',num2str(meanGAIN1)],'parent',gui.whichAX(2),'color',mycolors(7,:),'fontunits','normalized','fontsize',0.125)
    text(Xloc,Yls(1)+5.5*Yloc,['mean = ',num2str(meanGAIN2)],'parent',gui.whichAX(2),'color',mycolors(8,:),'fontunits','normalized','fontsize',0.125)
    text(Xloc,Yls(1)+4*Yloc,['BIC = ',num2str(DATA.BIC)],'parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
    if inputs.rorq==1
        text(Xloc-(Xls(2)-Xls(1))/20,Yls(1)+2.5*Yloc,'--- Gain used in QC','parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
    end
elseif   DATA.howmanyairs > 1 && strcmp(DATA.floatTYPE,'SOLO')
    x = DATA.GAINStime{3};
    y1 = DATA.GAINS{3};
    meanGAIN1 = nanmean(y1);
    X = x(~isnan(y1));
    Y = y1(~isnan(y1));
    [my1,by1,~,~,~]=lsqfity(X,Y);
    plot(x,y1,'o','parent',gui.whichAX(2),'markersize',5,'markerfacecolor',mycolors(7,:),'markeredgecolor',mycolors(7,:)); %AIRold
    xlim(gui.whichAX(2),DATA.xlims{1});
    Xls = get(gui.whichAX(2),'xlim');
    Yls = get(gui.whichAX(2),'ylim');
    hold(gui.whichAX(2),'on')
    newY = Xls.*my1+by1;
    %     plot(Xls,newY,'--','parent',gui.whichAX(2),'color',mycolors(7,:)) single regression line.  replaced with line showing mean
    plot(Xls,repmat(meanGAIN1,length(Xls)),'--','parent',gui.whichAX(2),'color',mycolors(7,:))
    hold(gui.whichAX(2),'on')
    plot(Xls,newY,'--','parent',gui.whichAX(2),'color',mycolors(7,:),'linewidth',2)
    ysies = get(gui.whichAX(2),'ylim');
    % Plot Breakpoint analysis:
% % %     if inputs.rorq==1
% % %         for jt = 1:size(DATA.tableDATA,1)
% % %             %                 Ex = DATA.refsub(DATA.refsub(:,2)>=DATA.tableDATA(jt,1)&DATA.refsub(:,2)<=inputs.cyEND(jt),1);
% % %             Ex = DATA.refdata(DATA.refdata(:,2)>=DATA.tableDATA(jt,1)&DATA.refdata(:,2)<=inputs.cyEND(jt),1);
% % %             Exes = Ex-nanmin(Ex);
% % %             Whys = Exes.*(DATA.tableDATA(jt,3))./365+DATA.tableDATA(jt,2);
% % %             hold(gui.whichAX(2),'on')
% % %             plot(Ex,Whys,'-','parent',gui.whichAX(2),'color',mycolors(9,:),'linewidth',2)
% % %             hold(gui.whichAX(2),'on')
% % %             myEx = DATA.refsub(DATA.refsub(:,2)==DATA.tableDATA(jt,1),1);
% % %             %                 ysies = get(gui.whichAX(2),'ylim');
% % %             if ~isempty(myEx)
% % %                 plot([myEx myEx],ysies,'--','parent',gui.whichAX(2),'color',mycolors(3,:))
% % %             end
% % %         end
% % %     end
end
if DATA.howmanyairs == 1 || strcmp(DATA.floatTYPE,'SOLO')
    set(gui.whichAX(2),'ylim',ysies)
    Xls = get(gui.whichAX(2),'xlim');
    Yls = get(gui.whichAX(2),'ylim');
    Xloc = Xls(2)+(Xls(2)-Xls(1))/10;
    Yloc = (Yls(2)-Yls(1))/8;
    text(Xloc,Yls(1)+5.5*Yloc,['mean = ',num2str(meanGAIN1)],'parent',gui.whichAX(2),'color',mycolors(7,:),'fontunits','normalized','fontsize',0.125)
    text(Xloc,Yls(1)+4*Yloc,['BIC = ',num2str(DATA.BIC)],'parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
    if inputs.rorq==1
        text(Xloc-(Xls(2)-Xls(1))/20,Yls(1)+2.5*Yloc,'- Gain used in QC','parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
    end
end
ylabel(gui.whichAX(2),'(Ref pO2)/(Flt pO2)','fontunits','normalized','fontsize',0.125);
set(gui.whichAX(2),'xtick',DATA.xticks{1})
datetick(gui.whichAX(2),'x',2,'keepticks','keeplimits')
grid(gui.whichAX(2),'on')

%**************************************************************************
% AXIS 3: plot O2 data
%**************************************************************************
cla(gui.whichAX(3),'reset')
plot(DATA.O2subset(:,2),DATA.O2subset(:,7),'o','Parent',gui.whichAX(3),...
    'markersize',5,'markerfacecolor',mycolors(5,:),'markeredgecolor','k');
ylabel(gui.whichAX(3),'[O_2] (umol/kg)','fontunits','normalized','fontsize',0.125);
set(gui.whichAX(3),'xtick',DATA.xticks{2},'box','off','ycolor',mycolors(5,:));
xlim(gui.whichAX(3),DATA.xlims{2});
grid(gui.whichAX(3),'on')


%**************************************************************************
% AXIS 4: plot T,S
%**************************************************************************
cla(gui.whichAX(4))
cla(gui.whichAX(5))
%     hold(gui.whichAX(4),'on')
%
xlim(gui.whichAX(5),DATA.xlims{1});
hold(gui.whichAX(5),'on')
plot(DATA.PTSsubset(:,1),DATA.PTSsubset(:,6),'o','Parent',gui.whichAX(5),...
    'markersize',5,'markerfacecolor',mycolors(2,:),'markeredgecolor','k');
ylabel(gui.whichAX(5),'Temperature [C]','fontunits','normalized','fontsize',0.125);
set(gui.whichAX(5),'ycolor',mycolors(2,:));
set(gui.whichAX(5),'yaxislocation','right','color','none','xtick',DATA.xticks{1},'xticklabels',[],...
    'box','off','xcolor','none');
plot(DATA.PTSsubset(:,1),DATA.PTSsubset(:,7),'o','Parent',gui.whichAX(4),...
    'markersize',5,'markerfacecolor',mycolors(5,:),'markeredgecolor','k');
ylabel(gui.whichAX(4),'Salinity [pss]','fontunits','normalized','fontsize',0.125);
xlim(gui.whichAX(4),DATA.xlims{1});
set(gui.whichAX(4),'xtick',DATA.xticks{1},'box','off','ycolor',mycolors(5,:));
datetick(gui.whichAX(4),'x',2,'keepticks','keeplimits')
grid(gui.whichAX(4),'on')

end % redraw_CLIMO_sageO2
