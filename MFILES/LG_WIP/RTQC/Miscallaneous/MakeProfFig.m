function [] = MakeProfFig(s, plotpath)
% Creates figures for the 10 most recent profiles with the last one in
% thick red (crimson). inputs s stuct and plotpath (location to save).
% Ben Werb 5/30/24, based on YT 2020 plot_spraysat.m

% Default plot settings
set(0, 'DefaultAxesFontSize', 20, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontName', 'Times', ...
    'DefaultAxesFontWeight', 'bold', ...
    'DefaultLineMarkerSize', 10, ...
    'DefaultFigureUnit', 'centimeter', ...
    'DefaultFigurePosition', [6 1 25 20]);
% Select graphing display options
fname = 'avenir';
fsize = 12;
newcol = rgb('crimson'); % New profile color settings
newlw = 3; % New line size

% Define number of profiles to plot
if(size(s.lat_,2)<10) % if tot profiles less than 10
    nplot = size(s.lat_,2)-1; % plot n profs - 1
else
    nplot = 9; % how many dives to plot
end
profs2plot = size(s.pres,2)-nplot:size(s.pres,2); % prof idxs
preslimits = [0 round(max(max(s.pres(:,profs2plot))),1)+10]; % pres axis lims

vars = {'tc','psal','doxy','opt','pHin','deltapH','PAR'};
vname = {'Temperature','Salinity','O_2 [umol/kg]','Chl_a', 'pH',...
    'pH_in - pH_CANYONB','PAR'};
xlims = [-inf,inf;-inf,inf;0 350;-inf inf;-inf,inf;-inf,inf;-inf,inf]; % changed PAR limit to 0 to 1500 //YT 6/11/2022
% Need to validate xlims works for all data

s.deltapH = s.pHin - s.pHin_canb; % Define deltapH
s.deltapH_QC = s.pHin_QC; % QC field


% Create a tiled layout with 10 rows and 4 columns
tl = tiledlayout(10, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% title(tl, 'Last 10 Profiles')
% subtitle(tl,'Most recent in red')
cm = colormap(cmocean('dense',11));
cm = cm(2:11,:);
% Define the size of the figure to be full page (8.5 x 11 inches for a standard letter size)
fig = gcf;
fig.Visible = "off";
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.5 11];

% First 4 cols: four unique plots that span 3 rows and 1 col
for i = 1:4
    nexttile([3 1])
    var = char(vars(i));
    var2plot = s.(var);
    var2plot(s.([var '_QC'])~=0) = NaN;
    plot(var2plot(:,profs2plot), s.pres(:,profs2plot));
    hold on
    plot(s.(var)(:,end), s.pres(:,end),'Color', newcol,'Linewidth',newlw);
    set(gca, 'ydir', 'reverse', 'fontsize', fsize, 'FontName', fname);
    xlabel(char(vname(i)));
    xlim(xlims(i,:))
    ylim(preslimits);
    hold off
end

% 5 and 6 span 4 rows and 2 cols
for i = 5:6
    nexttile([4 2])
    var = char(vars(i));
    var2plot = s.(var);
    var2plot(s.([var '_QC'])~=0) = NaN;
    plot(var2plot(:,profs2plot), s.pres(:,profs2plot));
    hold on
    plot(s.(var)(:,end), s.pres(:,end),'Color', newcol,'Linewidth',newlw);
    set(gca, 'ydir', 'reverse', 'fontsize', fsize, 'FontName', fname);
    xlabel(char(vname(i)));
    xlim(xlims(i,:))
    ylim(preslimits);
    hold off
end

% PAR spans 2 rows and 4 cols
for i = 7
    nexttile([2 4])
    var = char(vars(i));
    var2plot = s.(var);
    var2plot(s.([var '_QC'])~=0) = NaN;
    plot(var2plot(:,profs2plot), s.pres(:,profs2plot));
    hold on
    plot(s.(var)(:,end), s.pres(:,end),'Color', newcol,'Linewidth',newlw);
    set(gca, 'ydir', 'reverse', 'fontsize', fsize, 'FontName', fname);
    xlabel(char(vname(i)));
    xlim(xlims(i,:))
    ylim(preslimits);
    hold off
end
% plotname = 'prof';
% SaveFigFolderServer(fig,plotpath,plotname);

end