 function redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs)
    
 % ************************************************************************
% redraw_PROF_sageO2Argo.m
% ************************************************************************
%
% Function to plot profile data within the GUI.
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
% UPDATES: incorporated Josh's changes to bottle lookup-table parsing
    %(8/21/17)
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
        Xi(:,icast)  = interp1(tmpOx(tu,1), tmpOx(tu,2), z);
    end
    meanprof = nanmean(Xi,2);
    stdprof = std(Xi,0,2);
    p1(1)=plot(z,meanprof,'b','linewidth',2,'Parent',gui.whichAX(1));
    hold(gui.whichAX(1),'on')
    p1(2)=plot(z,meanprof-stdprof,'b--','Parent',gui.whichAX(1));
    ylabel(gui.whichAX(1),'[O_2] (umol/kg)','fontunits','normalized','fontsize',0.125);
    if inputs.rorq == 1; % RAW TAB SELECTED?
        title(gui.whichAX(1),[Type,' Float ',inputs.floatID],'fontunits','normalized','fontsize',0.15)
    else
        title(gui.whichAX(1),[Type,' Float ',inputs.floatID,'. QC Adjustments applied.'],'fontunits','normalized','fontsize',0.15)
    end
%     XLIMS=get(gui.whichAX(1),'xlim')
    XLIMS=[inputs.depthedit(1,1) inputs.depthedit(1,2)];
    % Get GLODAP data
    track = [DATA.track(:,1:2) DATA.track(:,4) DATA.track(:,3)]; %switch lat/lon columns
    G = get_GLODAPv2_local_sO2Argo(track,inputs.GLDPkm,[0 2000],dirs.user_dir);
    GLODAP_color = [0   255 255;... % keep color and symbol count equal
            153  51 255;
            255 51 255;
            255 128 0]./255;
    GLODAP_symbols = {'^' 'd' 's' '*' '<' 'x'};
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
        plot(z,Gd,'r','Parent',gui.whichAX(1),'linewidth',2);
        if ~isempty(GIND)
                [G_YRS,~,~,~,~,~] = datevec(G.data(:,iGSDN));% Years only
                YRS = unique(G_YRS);
                ct = 0; ct2 = 1;

                for i = 1: size(YRS,1)
                    ct = ct+1;
                    t1 = G_YRS == YRS(i);
                    hold(gui.whichAX(1),'on')
                    vp = plot(G.data(t1,iGP),G.data(t1,GIND),'parent',gui.whichAX(1), ...
                            'Linestyle', 'none', ...
                            'Marker', GLODAP_symbols{ct2}, ...                        
                        'MarkerFaceColor', GLODAP_color(ct,:), ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 4);
                    hlegend_cell = [hlegend_cell, ...
                        {['GLODAP ',num2str(YRS(i))]}];
                    if ct == 4,
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
        meandiffG = nanmean(GLO_resid(GG2(:,1)>=XLIMS(1,1) & GG2(:,1)<=XLIMS(1,2)));
        meanGgain = nanmean(Ggain(GG2(:,1)>=XLIMS(1,1) & GG2(:,1)<=XLIMS(1,2)));
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
    %plot cast 1 and 2 with bottle data
    CAST1 = [D(D(:,2)==1,4) D(D(:,2)==1,7)]; %[P O2] prof data, single cast
    CAST2 = [D(D(:,2)==2,4) D(D(:,2)==2,7)]; %[P O2] prof data, single cast
    CAST3 = [D(D(:,2)==3,4) D(D(:,2)==3,7)]; %[P O2] prof data, single cast
    % Float 7567 missing first cast; quick fix, if cast 1 missing --> use
    % cast 2 and 3.  Does not fix potential case where cast 2 and/or 3 also missing!
    if isempty(CAST1)
        CAST1 = CAST2;
        CAST2 = CAST3;
    end
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
    ibN   = find(strcmp('NITRAT',b.hdr) == 1);
    ibPH  = find(strcmp('PH_TOT_INSITU',b.hdr) == 1);
    bIND = ibO; %For Oxygen GUI (keep other parameter indices for now)
    if ~isempty(bIND);
        hold(gui.whichAX(3),'on')
        plot(b.data(:,ibP),b.data(:,bIND),'o','parent',gui.whichAX(3),...
            'MarkerSize',6, 'MarkerFaceColor',[255 255 0]/255,...
                        'MarkerEdgeColor', 'k')
    end
    set(gui.whichAX(3),'xticklabels',[],'xlim',XLIMS)
    bot_leg_cell = {'CAST 1','CAST 2','Bottle'};
    bot_leg = legend(gui.whichAX(3),bot_leg_cell);       
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
    nonans = find(isnan(CAST1(:,2)))
    C1 = CAST1;
    C1(nonans,:) = [];
    FtoB = interp1(C1(:,1),C1(:,2),b.data(:,ibP));
    bot_resid = b.data(:,bIND)-FtoB;
    bot_gain = b.data(:,bIND)./FtoB;
%     meandiffB = nanmean(bot_resid);
%     meanBgain = nanmean(bot_gain);
    plot(b.data(:,ibP),bot_resid,'o','parent',gui.whichAX(4),...
        'MarkerSize',4, 'MarkerFaceColor','k',...
                        'MarkerEdgeColor', 'k');
    set(gui.whichAX(4),'xlim',XLIMS)
    meandiffB = nanmean(bot_resid(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)));
    meanBgain = nanmean(bot_gain(b.data(:,ibP)>=XLIMS(1,1) & b.data(:,ibP)<=XLIMS(1,2)));
    ylabel(gui.whichAX(4),'Bottle-Flt (umol/kg)','fontunits','normalized','fontsize',0.125);
    xlabel(gui.whichAX(4),'Pressure (db)','fontunits','normalized','fontsize',0.125);
    Xls = get(gui.whichAX(4),'xlim');
    Yls = get(gui.whichAX(4),'ylim');
    Xloc = Xls(2)+(Xls(2)-Xls(1))/15;
    Yloc = (Yls(2)-Yls(1))/8;
    text(Xloc,Yls(1)+6*Yloc,['Mean Resid = ',num2str(meandiffB)],'parent',gui.whichAX(4),'fontunits','normalized','fontsize',0.125)
    text(Xloc,Yls(1)+5*Yloc,['Mean Gain = ',num2str(meanBgain)],'parent',gui.whichAX(4),'fontunits','normalized','fontsize',0.125)
    grid(gui.whichAX(4),'on')


end % redraw_PROF_sageO2Argo
