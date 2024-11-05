% find_new_sensor_floats.m
% ------------------------------------------------------------------------
% % SCRIPT TO USE INFO FROM FLOAT CONFIG FILES & PH LOG FILE TO BUILD LISTING OF NEW
% SENSOR VARIANTS (ie SBE83, OCR, FLBBFL, MBARI GDF, AE GDF)
% PLOT IS ALSO MADE SHOWING COVERAGE OF NEW SENSORS.
%
% ORIGINAL CODE WRITTED BY JOSH PLANT
% UPDATED 4/17/24 BY TANYA MAURER (PLOT OPTIONS MODIFIED, AND MOVED TO VM
% OPERATIONS).
%--------------------------------------------------------------------------

fd             = userpath;
save_fd        = fullfile(fd,'ARGO_PROCESSING\DATA\CAL\NEW_SENSOR_LISTINGS');
new_sensors_fn = 'new_sensors_listing';
if ~isfolder(save_fd)
    mkdir(save_fd)
end

MYcolors = load_mycolors(0);
% ************************************************************************
% load master processing list
fp = fullfile(fd,'ARGO_PROCESSING\DATA\CAL','MBARI_float_list.mat');
load(fp)
mhdr  = d.hdr;
mlist = d.list;

% Get listing of float config filkes
config_fd = fullfile(fd,'ARGO_PROCESSING\DATA\CAL\FLOAT_CONFIG');
%fp = fullfile(fd,'ARGO_PROCESSING\DATA\CAL\FLOAT_CONFIG','*FloatConfig.txt');
tmp     = dir(fullfile(config_fd,'*FloatConfig.txt'));
fnames  = {tmp.name}';

% flter to remove any oddbal config files
tg = ~cellfun(@isempty, regexp(fnames,'^[usw][ans]\d+_FloatConfig.txt','once'));
fnames  = fnames(tg);
rfnames = size(fnames,1);

% ************************************************************************
% SBE83, OCR & FLBBFL floats can be identified from their config files
% PREDIM LOGICAL ARRAYS
tg_83fn  = false(size(fnames));     %  a config file exists
tg_83    = false(size(mlist(:,1))); % on master processing list

tg_OCRfn = false(size(fnames));
tg_OCR   = false(size(mlist(:,1)));

tg_FLBBFLfn = false(size(fnames));
tg_FLBBFL   = false(size(mlist(:,1)));

mbari_id = regexp(fnames,'\w+\d+\w*(?=\_)','match','once');

% STEP THROUGH FLOAT CONFIG FILLES
for fct = 1: rfnames
    fn = fnames{fct};
    id = mbari_id{fct};
    disp(fn)

    fid = fopen(fullfile(config_fd,fn));
    tmp = textscan(fid,'%s','Delimiter','\r\n','CollectOutput',1);
    fclose(fid);

    % SBE83
    tf83 = ~cellfun(@isempty, regexp(tmp{1,1},'SBE83','once'));
    if sum(tf83) >0
        tg_83fn(fct) = 1;
        tf = strcmp(mlist(:,1), id);
        if sum(tf) >0
            tg_83(tf) = 1;
        end
    end

    % OCR
    tfOCR = ~cellfun(@isempty, regexp(tmp{1,1},'^OCR','once'));
    if sum(tfOCR) >0
        tg_OCRfn(fct) = 1;
        tf = strcmp(mlist(:,1), id);
        if sum(tf) >0
            tg_OCR(tf) = 1;
        end
    end

    % FLBBFL
    %tfFLBBFL = ~cellfun(@isempty, regexp(tmp{1,1},'Chl435','once'));
    tfFLBBFL = ~cellfun(@isempty, regexp(tmp{1,1},'FLBBFL','once'));
    if sum(tfFLBBFL) >0
        tg_FLBBFLfn(fct) = 1;
        tf = strcmp(mlist(:,1), id);
        if sum(tf) >0
            tg_FLBBFL(tf) = 1;
        end
    end


end
fn83    = fnames(tg_83fn);
mlist83 = mlist(tg_83,:);

fnOCR    = fnames(tg_OCRfn);
mlistOCR = mlist(tg_OCR,:);

fnFLBBFL    = fnames(tg_FLBBFLfn);
mlistFLBBFL = mlist(tg_FLBBFL,:);


% ************************************************************************
% pH sensors need to be destinguished from the pH-log excel file
% get GDF FLOATS
d = get_ph_version;
tGDFM     = cell2mat(d.data(:,23)) == 7; % MBARI
tGDFA     = cell2mat(d.data(:,23)) == 8; % AE
tGDF      = tGDFM | tGDFA;
GDF_list = d.data(tGDF,:);
GDFAE_list = d.data(tGDFA,:);
GDFM_list = d.data(tGDFM,:);
Lia      = ismember( str2double(mlist(:,2)), cell2mat(GDF_list(:,2)) );
mlistGDF = mlist(Lia,:);

Lia      = ismember( str2double(mlist(:,2)), cell2mat(GDFAE_list(:,2)) );
mlistGDFAE = mlist(Lia,:);

Lia      = ismember( str2double(mlist(:,2)), cell2mat(GDFM_list(:,2)) );
mlistGDFM = mlist(Lia,:);

% Set output
sensors.list_hdr  = mhdr;
sensors.SBE83     = mlist83;
sensors.OCR       = mlistOCR;
sensors.FLBBFL    = mlistFLBBFL;
sensors.AE_GDF    = mlistGDFAE;
sensors.MBARI_GDF = mlistGDFM;

save(fullfile(save_fd,[new_sensors_fn,'.mat']),'sensors');

% ***********************************************************************
% ***********************************************************************
%                             PRINT LISTINGS
% ***********************************************************************
% ***********************************************************************
flds  = fieldnames(sensors);
tDATA = ~strcmp(flds,'list_hdr');
dflds = flds(tDATA);
phdr  = sensors.list_hdr;


I.MBA = strcmp('MBARI ID', phdr);  % DEFINE INDICES
I.INS = strcmp('INST ID', phdr);
I.WMO = strcmp('WMO', phdr);
I.TYP = strcmp('float type', phdr);
I.LON = strcmp('1st lon', phdr);
I.LAT = strcmp('1st lat', phdr);
I.SDN = strcmp('1st date', phdr);
% I.MSG = strcmp('msg dir', phdr);
% I.NC  = strcmp('NC template', phdr);
% I.PRG = strcmp('Program', phdr);
% I.REG = strcmp('Region', phdr);
% I.BFL = strcmp('tf Bfile', phdr);
% I.D   = strcmp('tf Dead', phdr);
% I.NO3 = strcmp('tf NO3', phdr);
% I.PH  = strcmp('tf pH', phdr);
% I.CHL = strcmp('tf Chl', phdr);
% I.OCR = strcmp('tf OCR', phdr);

I.MMF  = strcmp('max cycle msg file', phdr);
I.MMFD = strcmp('max cycle file date', phdr);
I.LMF  = strcmp('latest cycle msg file', phdr);
I.LMFD = strcmp('latest cycle file date', phdr);
I.MPC  = strcmp('max cycle proc', phdr);
I.MPCD = strcmp('max cycle proc date', phdr);

fid = fopen(fullfile(save_fd,[new_sensors_fn,'.TXT']),'w');
for ct = 1:size(dflds,1)
    fld  = dflds{ct};
    ptmp = sensors.(fld);

    % FOR TEXT FILE CONVERT SDN TO DATE STR
    % 1st profile SDN
    plist = ptmp;
    rplist = size(plist,1);
    tmp   = cell2mat(plist(:, I.SDN));
    tnan  = isnan(tmp);
    dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
    plist(~tnan, I.SDN) = cellstr(dstr);
    plist(tnan, I.SDN)  = cellstr('');

    % Max msg cycle file mod SDN
    tmp   = cell2mat(plist(:, I.MMFD));
    tnan  = isnan(tmp);
    dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
    plist(~tnan, I.MMFD) = cellstr(dstr);
    plist(tnan, I.MMFD)  = cellstr('');

    % Last msg cycle file mod SDN
    tmp   = cell2mat(plist(:, I.LMFD));
    tnan  = isnan(tmp);
    dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
    plist(~tnan, I.LMFD) = cellstr(dstr);
    plist(tnan, I.LMFD)  = cellstr('');

    % Last processed cycle date
    tmp   = cell2mat(plist(:, I.MPCD));
    tnan  = isnan(tmp);
    dstr  = datestr(tmp(~tnan),'mm/dd/yyyy HH:MM');
    plist(~tnan, I.MPCD) = cellstr(dstr);
    plist(tnan, I.MPCD)  = cellstr('');

    % FOR TEXT FILE PRINTING, CONVERT NaN's in LAT & LON & TO EMPTY
    mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.LON));
    plist(mask,I.LON) = {[]};
    plist(mask,I.LAT) = {[]};

    mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.MMF));
    plist(mask,I.MMF) = {[]};

    mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.LMF));
    plist(mask,I.LMF) = {[]};

    mask = cellfun(@(C) isnumeric(C) && isscalar(C) && isnan(C), plist(:,I.MPC));
    plist(mask,I.MPC) = {[]};


    % PRINT HEADER  & BUILD FORMAT STR
    if ct == 1
        fmt = '';
        for i = 1:size(phdr,2)-1
            fprintf(fid,'%s\t', phdr{i});
            if regexp(phdr{i},'^NC|^tf|cycle msg|cycle proc(?!\s)','once')
                fmt = [fmt,'%0.0f\t'];
            elseif regexp(phdr{i},'^1st l','once')
                fmt = [fmt,'%0.3f\t'];
            else
                fmt = [fmt,'%s\t'];
            end
        end
        fmt = [fmt,'%s\r\n'];
        fprintf(fid,'%s\r\n',phdr{i+1});
    end

    fprintf(fid,'\r\n ********   %s  ********\r\n',fld);
    D = plist';
    fprintf(fid, fmt, D{:});
end
fclose(fid);


% MAKE MAPS

    F1 = figure(1);
    F1.Units =  'inches';
    F1.Position = [8.5938 3.3667 11.5395 7.7750]; % work wide screen, inches

    % SBS83 *********************************************************
    ax(1) = subplot(2,2,1);
    m_proj('Miller Cylindrical','lat',[-80 80]);
    m_elev('contourf',30, 'edgecolor','none');
    m_grid('xaxislocation','bottom','box','on',  ...
        'ytick',[-45 0 45],'xtick',[-120 0 120], 'fontsize',16);
    hp = m_coast('patch',[.7 .7 .7],'edgecolor','k');

    % LOAD COLOR MAP
%     load('C:\Users\jplant\Documents\MATLAB\COLORMAPS\SOCCOM_NCP_cmap.mat')
%     colormap(cmap)
      colorlist = slanCM('bone');
      colormap(colorlist)


    X = cell2mat(mlist83(:,5));
    Y = cell2mat(mlist83(:,6));
    tnan = isnan(X);

    p1 = m_plot(X(~tnan), Y(~tnan),'ko', 'MarkerfaceColor',MYcolors(11,:));
    tstr = sprintf('SBE 83 oxygen floats (%0.0f)', sum(~tnan));
    title(tstr, 'FontSize', 16)


    % OCR *********************************************************
    ax(2) = subplot(2,2,2);
    m_proj('Miller Cylindrical','lat',[-80 80]);
    m_elev('contourf',30, 'edgecolor','none');
    m_grid('xaxislocation','bottom','box','on',  ...
        'ytick',[-45 0 45],'xtick',[-120 0 120], 'fontsize',16);
    hp = m_coast('patch',[.7 .7 .7],'edgecolor','k');

    X = cell2mat(mlistOCR(:,5));
    Y = cell2mat(mlistOCR(:,6));
    tnan = isnan(X);

    p2 = m_plot(X(~tnan), Y(~tnan),'ko', 'MarkerfaceColor',MYcolors(11,:));
    tstr = sprintf('OCR radiometer floats (%0.0f)', sum(~tnan));
    title(tstr, 'FontSize', 16)


    % GDF *********************************************************
    ax(3) = subplot(2,2,3);
    m_proj('Miller Cylindrical','lat',[-80 80]);
    m_elev('contourf',30, 'edgecolor','none');
    m_grid('xaxislocation','bottom','box','on',  ...
        'ytick',[-45 0 45],'xtick',[-120 0 120], 'fontsize',16);
    hp = m_coast('patch',[.7 .7 .7],'edgecolor','k');

    hold(ax(3), 'on')

    X1 = cell2mat(mlistGDFM(:,5));
    Y1 = cell2mat(mlistGDFM(:,6));
    tnan1 = isnan(X1);

    X2 = cell2mat(mlistGDFAE(:,5));
    Y2 = cell2mat(mlistGDFAE(:,6));
    tnan2 = isnan(X2);

    % p2 = m_plot(X(~tnan), Y(~tnan),'ko', 'MarkerfaceColor','r');
    % tstr = sprintf('GDF  pH on floats (%0.0f)', sum(~tnan));
    % title(tstr, 'FontSize', 16)

    p21 = m_plot(X1(~tnan1), Y1(~tnan1),'ko', 'MarkerfaceColor',MYcolors(11,:));
    tstr = sprintf('GDF pH floats (%d)', sum(~tnan1) + sum(~tnan2));
%     tstr2 = sprintf('MBARI (%d) & AE (%d)', sum(~tnan1),sum(~tnan2));
%     title({tstr, tstr2}, 'FontSize', 16)
    title(tstr,'Fontsize',16)

%     t = title([ ...
%     '\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', Mycolors(11,:)) '} GDF pH floats ' sum(~tnan1) + sum(~tnan2), ...
%     ',\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', Mycolors(6,:)) '} MBARI (%d) & AE (%d) ' sum(~tnan1),sum(~tnan2)], ...
%     'fontsize',16, 'Interpreter', 'tex');


    p22 = m_plot(X2(~tnan2), Y2(~tnan2),'ko', 'MarkerfaceColor',MYcolors(6,:));
    legend([p21 p22],{['MBARI (',num2str(sum(~tnan1)),')'],['AE (',num2str(sum(~tnan2)),')']},'fontsize',12)

    hold(ax(3), 'off')
    %L3 = m_legend([p21,p22],'MBARI', 'AE');
    

    % FLBBFL *********************************************************
    ax(4) = subplot(2,2,4);
    m_proj('Miller Cylindrical','lat',[-80 80]);
    m_elev('contourf',30, 'edgecolor','none');
    m_grid('xaxislocation','bottom','box','on',  ...
        'ytick',[-45 0 45],'xtick',[-120 0 120], 'fontsize',16);
    hp = m_coast('patch',[.7 .7 .7],'edgecolor','k');

    X = cell2mat(mlistFLBBFL(:,5));
    Y = cell2mat(mlistFLBBFL(:,6));
    tnan = isnan(X);

    p2 = m_plot(X(~tnan), Y(~tnan),'ko', 'MarkerfaceColor',MYcolors(11,:));
    tstr = sprintf('FLBBFL floats (%0.0f)', sum(~tnan));
    title(tstr, 'FontSize', 16)

    print(F1, fullfile(save_fd,[new_sensors_fn,'.png']), '-dpng');

    close all
	
	% TM add copy block to sirocco.
	siroccodir = '\\sirocco\wwwroot\gobgc\images\NEW_SENSOR_LISTING\';
	try
	status = copyfile(save_fd,siroccodir);
	catch
	disp('Error copying new-sensor-listing files to sirocco.')
	end
    




