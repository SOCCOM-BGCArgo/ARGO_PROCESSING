 function DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs)

 % ************************************************************************
% redraw_WOA_sageO2Argo.m
% ************************************************************************
%
% Function to plot WOA2013 comparisons within the GUI.
%
%
% INPUTS:
%    gui, inputs, and DATA and dirs area all structures of handles and inputs to the
%    SAGE_O2 gui.
%
% OUTPUTS:
%
%
% AUTHOR: Tanya Maurer
%         MBARI
%         tmaurer@mbari.org
%
% DATE: 03/14/2017
% UPDATES:
% NOTES: 
% ************************************************************************
%
% ************************************************************************

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
        if inputs.floatTYPE == 'A'; %APEX float
            Type = 'APEX';
        else %NAVIS float
            Type = 'NAVIS';
        end    
    end
    
    %**************************************************************************
    % AXIS 1: plot O2sat
    %**************************************************************************
    cla(gui.whichAX(1),'reset')
    plot(DATA.refsub(:,2),DATA.refsub(:,4),'r-','Parent',gui.whichAX(1),...
        'markersize',20,'linewidth',2);
    hold(gui.whichAX(1),'on')
    plot(DATA.refsub(:,2),DATA.SATsubset(:,2),...
        'color',mycolors(7,:),'Parent',gui.whichAX(1),'linewidth',2); %NAVIS mean surface O2 saturation%     errorbar(DATA.refsub(:,2),DATA.SATsubset(:,2),DATA.SATsubset(:,3),...
%         'color',mycolors(7,:),'Parent',gui.whichAX(1),'linewidth',2); %NAVIS mean surface O2 saturation
    lhandle = legend(gui.whichAX(1),'WOA','Float');
    posL=get(lhandle,'Position');
    posL(1) = 1.25*posL(1);
    posL(2) = 0.75*posL(2);
    set(lhandle,'Position',posL,'fontunits','normalized','fontsize',12);
    axis(gui.whichAX(1),'tight');
    xlim(gui.whichAX(1),DATA.xlims{2});
    ylabel(gui.whichAX(1),'Surface O_2 %Sat','fontunits','normalized','fontsize',0.125);
    set(gui.whichAX(1),'xtick',DATA.xticks{2});
    if inputs.rorq == 1; % RAW TAB SELECTED?
        title(gui.whichAX(1),[Type,' Float ',inputs.floatID],'fontunits','normalized','fontsize',0.15)
    else
        title(gui.whichAX(1),[Type,' Float ',inputs.floatID,'. QC Adjustments applied.'],'fontunits','normalized','fontsize',0.15)
    end
    grid(gui.whichAX(1),'on')
    
    %**************************************************************************
    % AXIS 2: plot gain
    %**************************************************************************
    cla(gui.whichAX(2))
    x = DATA.refsub(:,1);
    y1 = DATA.GAINS;
    meanGAIN1 = nanmean(y1);
    DATA.bigG = meanGAIN1; %Save to shared variable for application in Make_Mprof_ODVQC.m
    X = x(~isnan(y1));
    Y = y1(~isnan(y1));
    [my1,by1,~,~,~]=lsqfity(X,Y);
    plot(x,y1,'o','parent',gui.whichAX(2),'markersize',5,'markerfacecolor',mycolors(7,:),'markeredgecolor',mycolors(7,:)); %AIRold
    ysies = get(gui.whichAX(2),'ylim');
    % Plot Breakpoint analysis:
    if inputs.rorq==1;
        for jt = 1:size(DATA.tableDATA,1)
            Ex = DATA.refsub(DATA.refsub(:,2)>=DATA.tableDATA(jt,1)&DATA.refsub(:,2)<=inputs.cyEND(jt),1);
            Exes = Ex-nanmin(Ex);
            Whys = Exes.*(DATA.tableDATA(jt,3))./365+DATA.tableDATA(jt,2);
            hold(gui.whichAX(2),'on')
            plot(Ex,Whys,'-','parent',gui.whichAX(2),'color',mycolors(9,:),'linewidth',2)
            hold(gui.whichAX(2),'on')
            myEx = DATA.refsub(DATA.refsub(:,2)==DATA.tableDATA(jt,1),1);
            if ~isempty(myEx)
                plot([myEx myEx],ysies,'--','parent',gui.whichAX(2),'color',mycolors(3,:))
            end
        end
    end
    set(gui.whichAX(2),'ylim',ysies)
    xlim(gui.whichAX(2),DATA.xlims{1});
    Xls = get(gui.whichAX(2),'xlim');
    Yls = get(gui.whichAX(2),'ylim');
    hold(gui.whichAX(2),'on')
    newY = Xls.*my1+by1;
%     plot(Xls,newY,'--','parent',gui.whichAX(2),'color',mycolors(7,:)) single regression line.  replaced with line showing mean
    plot(Xls,repmat(meanGAIN1,length(Xls)),'--','parent',gui.whichAX(2),'color',mycolors(7,:))
    Xloc = Xls(2)+(Xls(2)-Xls(1))/10;
    Yloc = (Yls(2)-Yls(1))/8;
    text(Xloc,Yls(1)+5.5*Yloc,['mean = ',num2str(meanGAIN1)],'parent',gui.whichAX(2),'color',mycolors(7,:),'fontunits','normalized','fontsize',0.125)
%     text(Xloc,Yls(1)+6*Yloc,['m = ',num2str(my1)],'parent',gui.whichAX(2),'color',mycolors(7,:),'fontunits','normalized','fontsize',0.125)
%     text(Xloc,Yls(1)+5*Yloc,['b = ',num2str(by1)],'parent',gui.whichAX(2),'color',mycolors(7,:),'fontunits','normalized','fontsize',0.125)
        text(Xloc,Yls(1)+4*Yloc,['BIC = ',num2str(DATA.BIC)],'parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
        if inputs.rorq==1
%         text(Xloc-(Xls(2)-Xls(1))/20,Yls(1)+3*Yloc,'--- Gain used in QC','parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
        text(Xloc-(Xls(2)-Xls(1))/20,Yls(1)+2.5*Yloc,'- Gain used in QC','parent',gui.whichAX(2),'color',mycolors(9,:),'fontunits','normalized','fontsize',0.125)
    end
    ylabel(gui.whichAX(2),'(Ref %Sat)/(Flt %Sat)','fontunits','normalized','fontsize',0.125);
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
    xlim(gui.whichAX(3),DATA.xlims{2});
    set(gui.whichAX(3),'box','off','ycolor',mycolors(5,:));
    set(gui.whichAX(3),'xtick',DATA.xticks{2})
    grid(gui.whichAX(3),'on')

    %**************************************************************************
    % AXIS 4: plot T,S
    %**************************************************************************
    cla(gui.whichAX(4))
    cla(gui.whichAX(5))
    plot(DATA.PTSsubset(:,1),DATA.PTSsubset(:,7),'o','Parent',gui.whichAX(4),...
        'markersize',5,'markerfacecolor',mycolors(5,:),'markeredgecolor','k');
    ylabel(gui.whichAX(4),'Salinity [pss]','fontunits','normalized','fontsize',0.125);
    hold(gui.whichAX(4),'on')
    xlim(gui.whichAX(4),DATA.xlims{1});
    set(gui.whichAX(4),'xtick',DATA.xticks{1},'box','off','ycolor',mycolors(5,:));
    plot(DATA.PTSsubset(:,1),DATA.PTSsubset(:,6),'o','Parent',gui.whichAX(5),...
        'markersize',5,'markerfacecolor',mycolors(2,:),'markeredgecolor','k');
    ylabel(gui.whichAX(5),'Temperature [C]','fontunits','normalized','fontsize',0.125);
    set(gui.whichAX(5),'ycolor',mycolors(2,:));
    xlim(gui.whichAX(5),DATA.xlims{1});
    set(gui.whichAX(5),'yaxislocation','right','color','none','xtick',DATA.xticks{1},'xticklabels',[],...
                        'box','off','xcolor','none');
    datetick(gui.whichAX(4),'x',2,'keepticks','keeplimits')
    grid(gui.whichAX(4),'on')

end % redraw_WOA_sageO2Argo