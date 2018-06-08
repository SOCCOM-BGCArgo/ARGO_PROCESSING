function [] = PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
% PLOT DATA IN float_qc.m GUI
% MAKE PLOTS BASED ON INFO IN GUI handles structure

%   12/09/2016 Salinity and temperature choices have been added
%   10/30/2017 fixed plotting error if no  nitrate sensor exists
%   01/16/2018 modified for use with GUI Layout Toolbox

% ************************************************************************
% FIRST AND FOREMOST: CHECK THAT PARAMETER EXISTS FOR PLOTTING
if isempty(DATA.IND) %parameter doesn't exist
    cla(gui.whichAX(1),'reset')
    cla(gui.whichAX(2),'reset')
    cla(gui.whichAX(3),'reset')
    cla(gui.whichAX(4),'reset')
    set(gui.t1ax(1),'visible','off')
    set(gui.t1ax(2),'visible','off')
    set(gui.t1ax(3),'visible','off')
    set(gui.t1ax(4),'visible','off')
    set(gui.t2ax(1),'visible','off')
    set(gui.t2ax(2),'visible','off')
    set(gui.t2ax(3),'visible','off')
    set(gui.t2ax(4),'visible','off')
    msgbox({'ERROR: PARAMETER DOES NOT EXIST FOR SELECTED FLOAT.','PLEASE CHOOSE ANOTHER PARAMETER TO PLOT.'});
    return
end
 
% ************************************************************************
% GET QC_TABLE CYCLE BREAKS FOR PLOTTING
if ~isempty(DATA.IND) && ~isempty(DATA.tableDATA)
    tdata = DATA.tableDATA(:,1);
else
    tdata = NaN;
end

% ************************************************************************
%                       PREP FOR PLOT MAKING
% ************************************************************************
float_color  = [0 102 204]/255; %light blue
MLR_color    = [1 0 0]; %red
DIFF_color   = [0 153 76]/255; %green
GLODAP_color = [0   255 255;... % keep color and symbol count equal
                153  51 255;
                255 51 255;
                255 128 0]./255;
GLODAP_symbols = {'^' 'd' 's' '*' '<' 'x'};

if strcmp(DATA.paramtag, 'PH')
    txt_units = 'pH units yr^{-1}';
    txt_units1 = 'pH units';
elseif strcmp(DATA.paramtag, 'S')
    txt_units = 'pss yr^{-1}';
    txt_units1 = 'pss';
elseif strcmp(DATA.paramtag, 'T')
    txt_units = 'ºC yr^{-1}';
    txt_units1 = 'ºC';
else 
    txt_units = 'µmol kg^{-1} yr^{-1}';
    txt_units1 = 'µmol kg^{-1}';    
end

% legends = findobj('type', 'legend');
% if ~isempty(legends)
%     delete(legends)
% end
% 
% colorbars = findobj('type', 'colorbar');
% if ~isempty(colorbars)
%     delete(colorbars)
% end

sdn_lim   = [nanmin(DATA.datasub(:,1)) nanmax(DATA.datasub(:,1))];
if ~isempty(sdn_lim)
    sdn_ticks = linspace(sdn_lim(1), sdn_lim(2),4);
else
    hmsg = msgbox({'No Data within auto depth limits.','Expand lower end of range.'},'Warning:');
    sdn_lim = [nanmin(DATA.datatype.data(:,1)) nanmax(DATA.datatype.data(:,1))];
    sdn_ticks = linspace(sdn_lim(1), sdn_lim(2),4);
end


if diff(sdn_lim) == 0
    sdn_lim(2) = sdn_lim(1)+1;
    sdn_ticks = sdn_lim;
end


% *************************************************************************
%                    DEEP SURFACE OR VERTICAL PLOTS
% *************************************************************************
% GET TRENDS

% DATA.datasub_time = DATA.datasub;
% DATA.refsub_time = DATA.refsub;
% DATA.diffsub_time = DATA.diffsub;

DATA.datasub_time = DATA.datasub(DATA.datasub(:,2) >= DATA.xlims{2}(1) & DATA.datasub(:,2) <= DATA.xlims{2}(2),:);
DATA.refsub_time = DATA.refsub(DATA.datasub(:,2) >= DATA.xlims{2}(1) & DATA.datasub(:,2) <= DATA.xlims{2}(2),:);
DATA.diffsub_time = DATA.diffsub(DATA.datasub(:,2) >= DATA.xlims{2}(1) & DATA.datasub(:,2) <= DATA.xlims{2}(2),:);

t_nanF   = ~isnan(DATA.datasub_time(:,DATA.IND)); % float
t_nanMLR = ~isnan(DATA.refsub_time); % MLR
t_nanDIFF = ~isnan(DATA.diffsub_time); % DIFF

if sum(t_nanF) > 1 % more than 2 pts - can do a reg
    %[m,b,r,sm,sb] = lsqfitma(ddata(t_nanF,1), ddata(t_nanF,IND));
    [m,b,r,sm,sb] = lsqfity(DATA.datasub_time(t_nanF,1), DATA.datasub_time(t_nanF,DATA.IND));
    dfloat_reg   = [m,b,r,sm,sb]; % FLOAT REGRESSION OUTPUT
else
    dfloat_reg   = [NaN, NaN, NaN, NaN, NaN];
end

if sum(t_nanMLR) > 1 % more than 2 pts - can do a reg
    t_good = t_nanF & t_nanMLR;
    %[m,b,r,sm,sb] = lsqfitma(ddata(t_good,1), dMLR_X(t_good));
    [m,b,r,sm,sb] = lsqfity(DATA.datasub_time(t_good,1), DATA.refsub_time(t_good));
    dMLR_reg   = [m,b,r,sm,sb]; % MLR REGRESSION OUTPUT
else
    dMLR_reg   = [NaN, NaN, NaN, NaN, NaN];
end

if sum(t_nanDIFF) > 1 % more than 2 pts - can do a reg
    t_good = t_nanF & t_nanDIFF;
    %[m,b,r,sm,sb] = lsqfitma(ddata(t_good,1), dDIFF_X(t_good));
    [m,b,r,sm,sb] = lsqfity(DATA.datasub_time(t_good,1), DATA.diffsub_time(t_good));
    dDIFF_reg   = [m,b,r,sm,sb]; % DIFF REGRESSION OUTPUT

    if ~isempty(tdata) && ~isempty(DATA.diffsub_time(t_good)) % Calc AIC
        aic_K = size(tdata,1)*2+2;
        aic_N = size(DATA.diffsub_time(t_good),1);
        if aic_K *4 > aic_N-1 % Valid data parameters?
            aic = NaN;
            aic_str = sprintf(['AIC = %1.1f (N~>>k)  AIC N ',...
                '= %1.0f AIC K = %1.0f'], aic, aic_N, aic_K);
        else
            aic = mAIC(DATA.diffsub_time(t_good), aic_K, 'owens');
            aic_str = sprintf(['AIC = %1.1f   AIC N = %1.0f',...
                '   AIC K = %1.0f'], aic, aic_N, aic_K);
        end   
    end
else
    dDIFF_reg   = [NaN, NaN, NaN, NaN, NaN];
    aic = NaN;
    aic_N = NaN;
    aic_K = NaN;
    aic_str = sprintf(['AIC = %1.1f  AIC N = %1.0f',...
        'AIC K = %1.0f'], aic, aic_N, aic_K);
end
handles.info.AIC = [aic aic_K aic_N];  

t1 = ~isnan(DATA.diffsub_time);
DATA.RSSN_CT = sum(t1);
DATA.RSSN = sum(DATA.diffsub_time(t1) .* DATA.diffsub_time(t1))./DATA.RSSN_CT;

%%%% MAKE PLOTS %%%%
%**************************************************************************
% AXIS 1: plot by time
%**************************************************************************
cla(gui.whichAX(1),'reset')            
plot(DATA.datasub(:,1),DATA.datasub(:,DATA.IND),'o','Parent',gui.whichAX(1),'MarkerSize',4, ...
    'MarkerEdgeColor','k','MarkerFaceColor', float_color);
xlim(gui.whichAX(1),DATA.xlims{1});
myxlim = get(gui.whichAX(1),'xlim');
set(gui.whichAX(1),'xtick',sdn_ticks);
% ylim auto
datetick(gui.whichAX(1),'x','mm/dd/yy', 'keeplimits', 'keepticks')

ylabel(gui.whichAX(1),inputs.y_label)
hold(gui.whichAX(1),'on')
plot(DATA.datasub(:,1),DATA.refsub,'o','Parent',gui.whichAX(1),'MarkerSize',4, ...
    'MarkerEdgeColor','k','MarkerFaceColor', MLR_color);
hold(gui.whichAX(1),'on')
plot(myxlim, myxlim*dfloat_reg(1) + dfloat_reg(2),'-','Parent',gui.whichAX(1),...
    'LineWidth',2,'Color',float_color);
hold(gui.whichAX(1),'on')
plot(myxlim, myxlim*dMLR_reg(1) + dMLR_reg(2),'-','Parent',gui.whichAX(1),...
    'LineWidth',2,'Color',MLR_color);


s1 = sprintf(['Float trend =   %0.4f  ', txt_units], ...
    dfloat_reg(1)*365);
s2 = sprintf([DATA.reftag,'  trend =   %0.4f  ', txt_units], ...
    dMLR_reg(1)*365);
tstr1 = sprintf('\\fontsize{12}\\color[rgb]{%f, %f, %f}%s', float_color, s1);
tstr2 = sprintf('\\fontsize{12}\\color[rgb]{%f, %f, %f}%s', MLR_color, s2);
title({tstr1,tstr2},'parent',gui.whichAX(1),'FontWeight', 'normal')
% T1 = text(0,0.98,s1,'parent',gui.whichAX(1),'Units', 'Normalized', 'HorizontalAlignment', ...
%     'Left','Color',float_color,'FontSize', 12);
% T2 = text(0,0.96,s2,'parent',gui.whichAX(1),'Units', 'Normalized', 'HorizontalAlignment', ...
%     'Left','Color',MLR_color,'FontSize', 12);

                     
%**************************************************************************
% AXIS 2: plot by cycle
%**************************************************************************
cla(gui.whichAX(2),'reset')
plot(DATA.datasub(:,2),DATA.datasub(:,DATA.IND),'o','Parent',gui.whichAX(2),'MarkerSize',4, ...
    'MarkerEdgeColor','k','MarkerFaceColor', float_color);
xlim(gui.whichAX(2),DATA.xlims{2});
% ylim auto
xlabel(gui.whichAX(2),'Cycle #')
ylabel(gui.whichAX(2),inputs.y_label)
title(gui.whichAX(2),inputs.data_str)

hold(gui.whichAX(2),'on')
plot(DATA.datasub(:,2),DATA.refsub,'o','Parent',gui.whichAX(2),'MarkerSize',4, ...
    'MarkerEdgeColor','k','MarkerFaceColor', MLR_color);
hold(gui.whichAX(2),'on')
            
% PLOT ANY GLODAP CROSSOVER DATA IF IT EXISTS
if ~isempty(DATA.GIND)
    %[handles.deep_figs.UserData(1) handles.deep_figs.UserData(2)]
    %[min(G.data(:,GIND)) max(G.data(:,GIND))]
    if ~isempty(DATA.glosub)
        [G_YRS,~,~,~,~,~] = datevec(DATA.glosub(:,DATA.iGSDN));% Years only
        YRS = unique(G_YRS);
        ct = 0; ct2 = 1;
        legend_cell ={};
        for i = 1: size(YRS,1)
            ct = ct+1;
            t1 = G_YRS == YRS(i);
            %[min(G.data(t1,GIND)) max(G.data(t1,GIND))
            hp(i) = plot(DATA.glosub(t1,DATA.iGcyc),DATA.glosub(t1,DATA.GIND),'Parent',gui.whichAX(2), ...
                'Linestyle', 'none', ...
                'Marker', GLODAP_symbols{ct2}, ...
                'MarkerFaceColor', GLODAP_color(ct,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 6);
            legend_cell{1,i} = ['GLODAP ',num2str(YRS(i))];
            if ct == 4,
                ct  = 0;
                ct2 = ct2 + 1;
                if ct2 > size(GLODAP_symbols,2)
                    ct2 = 1;
                end
            end
        end

        leg_cell = [{[handles.info.data_quality, ' data'], DATA.reftag} legend_cell];
        Glegend = legend(gui.whichAX(2),leg_cell);
        Glegend.Position = [0.775,0.94 0.1 0.05];
        Glegend.Interpreter = 'none';
    end
end
if exist('Glegend')
    posL=get(Glegend,'Position');
    posL(1) = 1.05*posL(1);
    posL(2) = 0.65*posL(2);
    set(Glegend,'Position',posL,'fontunits','normalized','fontsize',10);
end


%**************************************************************************
% AXIS 3: plot differences by time
%**************************************************************************
cla(gui.whichAX(3),'reset')
if strcmp(DATA.paramtag,'NO3')==1 || strcmp(DATA.paramtag,'PH')==1
    plot(DATA.datasub(:,1),DATA.diffsub,'o','MarkerSize',4,'Parent',gui.whichAX(3), ...
    'MarkerEdgeColor','k','MarkerFaceColor', DIFF_color);
    hold(gui.whichAX(3),'on')
    plot(myxlim, myxlim*dDIFF_reg(1) + dDIFF_reg(2),'-','Parent',gui.whichAX(3),...
        'LineWidth',2,'Color',DIFF_color);
%     hold off
    xlim(gui.whichAX(3),DATA.xlims{1});
    myxlim = get(gui.whichAX(3),'xlim');
    set(gui.whichAX(3),'xtick',sdn_ticks);
%     ylim auto
    datetick(gui.whichAX(3),'x','mm/dd/yy', 'keeplimits', 'keepticks')
    ylabel(gui.whichAX(3),inputs.y_label)
    s3 = sprintf(['Anomaly trend =   %0.4f  ', txt_units], ...
        dDIFF_reg(1)*365);
    s4 = sprintf(['Mean Diff (Flt-MLR)=   %0.3f  ',txt_units1], ...
        nanmean(DATA.diffsub_time));

    title({s3,s4},'Color', DIFF_color,'parent',gui.whichAX(3),'FontSize', 12, ...
        'FontWeight', 'normal')
    % T3 = text(0,1.25,s3,'Units', 'Normalized', 'parent',gui.whichAX(3),'HorizontalAlignment', ...
    %     'Left','Color',DIFF_color,'FontSize', 12);
    % T4 = text(0,1.1,s4,'Units', 'Normalized', 'parent',gui.whichAX(3),'HorizontalAlignment', ...
    %     'Left','Color',DIFF_color,'FontSize', 12);
else
    set(gui.t1ax(3),'visible','off')
    set(gui.t2ax(3),'visible','off')
end

%**************************************************************************
% AXIS 4: plot differences by cycle
%**************************************************************************
cla(gui.whichAX(4),'reset')
if strcmp(DATA.paramtag,'NO3')==1 || strcmp(DATA.paramtag,'PH')==1
    plot(DATA.datasub(:,2),DATA.diffsub,'o','MarkerSize',4,'Parent',gui.whichAX(4), ...
    'MarkerEdgeColor','k','MarkerFaceColor', DIFF_color);
    xlim(gui.whichAX(4),DATA.xlims{2});
    myxlim = get(gui.whichAX(4),'xlim');
%     ylim auto
    xlabel(gui.whichAX(4),'Cycle #')
    ylabel(gui.whichAX(4),inputs.y_label)

    s2 = sprintf('RSS = %1.2E  Total N = %1.2E',DATA.RSSN, DATA.RSSN_CT);
    title({aic_str,s2},'Color', DIFF_color,'parent',gui.whichAX(4),'FontSize', 12, ...
        'FontWeight', 'normal')

    hold(gui.whichAX(4),'on')
    yrange = get(gui.whichAX(3),'ylim');
    qc_breaks =[tdata,tdata*0+yrange(1);tdata,tdata*0+yrange(2); ...
        tdata,tdata*NaN];
    [~,I] = sort(qc_breaks(:,1));
    qc_breaks = qc_breaks(I,:);
    plot(qc_breaks(:,1),qc_breaks(:,2),'--','Color', ...
        [255 128 0]/255,'parent',gui.whichAX(4),'LineWidth',2)
    set(gui.whichAX(4),'ylim',yrange)
%     hold off
else
    set(gui.t1ax(4),'visible','off')
    set(gui.t2ax(4),'visible','off')
end
end









