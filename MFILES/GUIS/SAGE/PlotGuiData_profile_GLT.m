function [] = PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)

% ***************************************************************G*********
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
%y_lim   = [0 1600];
y_lim       = [inputs.depthedit(1) inputs.depthedit(2)];

% gui
% inputs
% DATA.datatype.hdr

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

cmp_str = DATA.reftag;

profile_lim = [inputs.profedit(1) inputs.profedit(2)];
iSTA        = strcmp(DATA.datatype.hdr,'Station');
tSTA        = DATA.datatype.data(:,iSTA) >= profile_lim(1) & ...
              DATA.datatype.data(:,iSTA) <= profile_lim(2);
profdata    = DATA.datatype.data(tSTA,:);
vMLRtmp = DATA.refdata(tSTA);

tvert   = profdata(:,DATA.iP) >= y_lim(1) & profdata(:,DATA.iP) <= y_lim(2);
vdata   = profdata(tvert,:); % subset
vMLR_X  = vMLRtmp(tvert);
vDIFF_X = DATA.DIFF_X(tvert);
sdn_lim   = [nanmin(profdata(:,1)) nanmax(profdata(:,1))];
sdn_ticks = linspace(sdn_lim(1), sdn_lim(2),4);

vNaN    = isnan(vdata(:,DATA.IND));
vdata   = vdata(~vNaN,:);
if ~isempty(vdata)
    vMLR_X  = vMLR_X(~vNaN);
    vDIFF_X = vDIFF_X(~vNaN);

    [~,ia,~] = unique(vdata(:,2));
    casts   = vdata(ia,1:4);
    r_casts = size(casts,1);
    clear ia

    % STEP THROUGH AND SORT ALL THE PROFILES. MPROFs CAN HAVE
    % MULTIPLE PROFILES IN A CYCLE
     for i = 1:r_casts
        t1  = vdata(:,2) == casts(i,2);
        tmp = vdata(t1,:); % cast subset
        [~,SX]     = sort(tmp(:,DATA.iP));
        MLR_X_tmp  = vMLR_X(t1);
        DIFF_X_tmp = vDIFF_X(t1);

        vdata(t1,:) = tmp(SX,:);
        vMLR_X(t1)  = MLR_X_tmp(SX);
        vDIFF_X(t1) =  DIFF_X_tmp(SX);
     end
     clear i tmp t1 SX MLR_X_tmp DIFF_X_tmp

    %z      = [8,11:5:100,111:10:400,450:50:1000,1100:100:1600]';
    z      = [8,11:5:100,111:10:400,450:50:1000,1100:100:2000]';
    max_z  = max(vdata(:,DATA.iP));
    z(z > max_z) =[];
    r_z     = size(z,1);
    Xi = ones(r_z,r_casts) * NaN; % PREDIMMENSION MATRICES
    MLR_Xi = Xi;
    DIFF_Xi = Xi;

    for i = 1:r_casts
        t1  = vdata(:,2) == casts(i,2);
        tmp = vdata(t1,:); % cast subset
        MLR_X_tmp  = vMLR_X(t1);
        DIFF_X_tmp = vDIFF_X(t1); 

        % SORT & ORDER BEFORE INTERP
        [~,SX,~] = unique(tmp(:,DATA.iP));
        %[~,SX] = sort(tmp(:,iP),1,'ascend');
        tmp = tmp(SX,:);
        MLR_X_tmp = MLR_X_tmp(SX);
        DIFF_X_tmp = DIFF_X_tmp(SX);



%                     tmp = flipud((vdata(t1,:)));
%                     MLR_X_tmp  = flipud(vMLR_X(t1));
%                     DIFF_X_tmp = flipud(vDIFF_X(t1));

        if size(tmp,1) > 3 % there
            % Interpolate for contour plotting & density grabs
            % pressure needs to be monotonic - quick clean

%                         tpress = 1; % prime loop
%                         while sum(tpress) > 0
%                             tpress = tmp(2:end,iP) - tmp(1:end-1,iP) <= 0;
%                             tmp(tpress,:) = [];
%                             MLR_X_tmp(tpress) = [];
%                             DIFF_X_tmp(tpress) = [];
%                         end

            Xi(:,i)     = interp1(tmp(:,DATA.iP), tmp(:,DATA.IND), z);

            gMLR_X_tmp = ~isnan(MLR_X_tmp);
            if sum(gMLR_X_tmp) > 3
                MLR_Xi(:,i) = interp1(tmp(gMLR_X_tmp,DATA.iP), ...
                              MLR_X_tmp(gMLR_X_tmp), z);
            else
                MLR_Xi(:,i)  = z* NaN;
            end

            gDIFF_X_tmp = ~isnan(DIFF_X_tmp);
            if sum(gDIFF_X_tmp) > 3
                DIFF_Xi(:,i) = interp1(tmp(gDIFF_X_tmp,DATA.iP), ...
                               DIFF_X_tmp(gDIFF_X_tmp), z);
            else
                DIFF_Xi(:,i) = z* NaN;
            end

        else
            Xi(:,i)      = z* NaN;
            MLR_Xi(:,i)  = z* NaN;
            DIFF_Xi(:,i) = z* NaN;
        end
    end

    mXi      = nanmean(Xi,2); % Get "average" profile
    mMLR_Xi  = nanmean(MLR_Xi,2);
    mDIFF_Xi = nanmean(DIFF_Xi,2);

%%% MAKE PLOTS %%%
%**************************************************************************
% AXIS 1: parameter vs cycle
%**************************************************************************    
    cla(gui.whichAX(1),'reset')     
    if strcmp(DATA.paramtag, 'O2')
        set(gui.whichAX(1),'visible','off') %don't plot contour plot for O2
    else
        if ~all(isnan(DIFF_Xi(:)))
            c_avg    = nanmean(DIFF_Xi(:));
            c_std    = nanstd(DIFF_Xi(:));
            c_limits = [c_avg-(3*c_std) c_avg+(3*c_std)];

            contourf(casts(:,1),z,DIFF_Xi,20,'parent',gui.whichAX(1),'linecolor','none')
            xlim(gui.whichAX(1),sdn_lim);
            set(gui.whichAX(1),'xtick',sdn_ticks);
            ylim(gui.whichAX(1),y_lim);
            caxis(gui.whichAX(1),c_limits);
            set(gui.whichAX(1),'Ydir','reverse');
            datetick(gui.whichAX(1),'x','mm/dd/yy', 'keeplimits', 'keepticks')
            ylabel(gui.whichAX(1),'Depth')
            title(inputs.data_str,'parent',gui.whichAX(1))
            %title({'Depth interpolated Float - MLR',['(',MLR_flavor,')']})
            cb = colorbar(gui.whichAX(1),'NorthOutside');
            ax_pos = gui.whichAX(1).Position;
            set(cb, 'Position', [ax_pos(1) .85 ax_pos(3) .05]);
        end
    end

%**************************************************************************
% AXIS 2: vertical profile plots
%**************************************************************************
    cla(gui.whichAX(2),'reset')            
    plot(mXi,z,'-','Color', float_color, 'parent',gui.whichAX(2),'LineWidth', 2)
    ylim(gui.whichAX(2),y_lim);
    set(gui.whichAX(2),'Ydir','reverse');
    xlabel(gui.whichAX(2),['Mean ',inputs.y_label])
    ylabel(gui.whichAX(2),'Pressure, dbar')
    hold(gui.whichAX(2),'on')
    tnan = ~isnan(mMLR_Xi);
    if sum(tnan) > 0
        plot(mMLR_Xi(tnan), z(tnan),'-','Color', MLR_color, ...
            'parent',gui.whichAX(2),'LineWidth', 2)
    else
        plot(mMLR_Xi,z,'-','Color', MLR_color, 'parent',gui.whichAX(2),'LineWidth', 2)
    end

    %title({fname, 'average vertical profiles'})
    hlegend_cell = {[handles.info.data_quality, ' data'], cmp_str};


    % PLOT ANY GLODAP CROSSOVER DATA IF IT EXISTS
    if ~isempty(DATA.GIND)
        [G_YRS,~,~,~,~,~] = datevec(DATA.G.data(:,DATA.iGSDN));% Years only
        YRS = unique(G_YRS);
        ct = 0; ct2 = 1;

        for i = 1: size(YRS,1)
            ct = ct+1;
            t1 = G_YRS == YRS(i);
            hold(gui.whichAX(1),'on')
            vp(i) = plot(DATA.G.data(t1,DATA.GIND),DATA.G.data(t1,DATA.iGP), ...
                'Linestyle', 'none', ...
                'Marker', GLODAP_symbols{ct2}, ...
                'MarkerFaceColor', GLODAP_color(ct,:), ...
                'MarkerEdgeColor', 'k', 'parent',gui.whichAX(2),'MarkerSize', 4);
            hlegend_cell = [hlegend_cell, ...
                ['GLODAP ',num2str(YRS(i))]];
            if ct == 4
                ct = 0;
                ct2 = ct2+1;
            end
        end
    end

    %hlegend = legend([handles.info.data_quality, ' data'], cmp_str);
    hlegend = legend(gui.whichAX(2),hlegend_cell);
    hlegend.Position = [0.8,0.84 0.1 0.05];
    hlegend.Orientation = 'Vertical';
    hlegend.Interpreter = 'none';
    posL=get(hlegend,'Position');
    posL(1) = 1.05*posL(1);
    posL(2) = 0.65*posL(2);
    set(hlegend,'Position',posL,'fontunits','normalized','fontsize',10);
    

%**************************************************************************
% AXIS 3: diff vs time
%**************************************************************************
    cla(gui.whichAX(3),'reset')            
    if ~strcmp(DATA.paramtag, 'O2')
        tnan = ~isnan(mDIFF_Xi);
        if sum(tnan) > 0
            plot(mDIFF_Xi(tnan), z(tnan),'k-', 'parent',gui.whichAX(3),'LineWidth', 2)
        else
            plot(mDIFF_Xi,z,'k-', 'parent',gui.whichAX(3),'LineWidth', 2)
        end

        ylim(gui.whichAX(3),y_lim);
        set(gui.whichAX(3),'Ydir','reverse');
        xlabel(gui.whichAX(3),['Mean diff,',inputs.y_label])
        ylabel(gui.whichAX(3),'Pressure, dbar')
    else
        % GET 1st CAST
        casts       = unique(DATA.datatype.data(:,2));
        start_ind  = find(casts >0, 1, 'first');
        if size(start_ind,2) == 1
            start_cast = casts(start_ind(1));
        end
        c1 = vdata(:,2) == start_cast(1);
        if sum(~isnan(vdata(c1,DATA.IND))) > 0 && ~isempty(DATA.bIND)% DATA?
            IND_interp = interp1(vdata(c1,DATA.iP), vdata(c1,DATA.IND), ...
                DATA.b.data(:,DATA.ibP));

            plot(DATA.b.data(:,DATA.bIND) -IND_interp ,DATA.b.data(:,DATA.ibP), ...
                '-','Color', float_color, 'parent',gui.whichAX(3),'LineWidth', 2)
            set(gui.whichAX(3),'Ydir','reverse');
            xlabel(gui.whichAX(3),'Bottle - float')
        end
    end

%**************************************************************************
% AXIS 4: diff vs cycle
%**************************************************************************    
    cla(gui.whichAX(4),'reset')       
    % GET 1st TWO CASTS
    casts       = unique(profdata(:,2));
    start_inds  = find(casts >0, 2, 'first');
    if size(start_inds,1) == 1
        start_casts = [casts(start_inds(1)) casts(start_inds(1))];
    else
        start_casts = [casts(start_inds(1)) casts(start_inds(2))];
    end


    c1 = vdata(:,2) == start_casts(1);
    c2 = vdata(:,2) == start_casts(2);

    if sum(~isnan(vdata(c1|c2,DATA.IND))) > 0 % DATA TO PLOT?

        plot(vdata(c1,DATA.IND), vdata(c1,DATA.iP),'-','Color', float_color, ...
            'parent',gui.whichAX(4),'LineWidth', 2)
        set(gui.whichAX(4),'Ydir','reverse');
        ylim(gui.whichAX(4),y_lim);

        hold(gui.whichAX(4),'on')
        plot(vdata(c2,DATA.IND), vdata(c2,DATA.iP),'--','Color', float_color, ...
            'parent',gui.whichAX(4),'LineWidth', 2)
        legend_str ={['cycle ',num2str(start_casts(1))], ...
            ['cycle ',num2str(start_casts(2))]};

        % PLOT MLR OR WOA PROFILES
        c1_data = [vMLR_X(c1), vdata(c1,DATA.iP)]; % REMOVE NaN's for Plotting
        tNaN = isnan(c1_data(:,1));
        c1_data(tNaN,:) = [];

        c2_data = [vMLR_X(c2), vdata(c2,DATA.iP)];
        tNaN = isnan(c2_data(:,1));
        c2_data(tNaN,:) = [];

        plot(c1_data(:,1), c1_data(:,2),'-','Color', MLR_color, ...
            'parent',gui.whichAX(4),'LineWidth', 2);
        tNaN = isnan(vMLR_X(c2));
        hold(gui.whichAX(4),'on')
        plot(c2_data(:,1), c2_data(:,2),'--','Color', MLR_color, ...
            'parent',gui.whichAX(4),'LineWidth', 2);
        if sum(~isnan(c1_data(:,1))) > 0
            legend_str = [legend_str, [cmp_str, ' ', ...
                num2str(start_casts(1))],[cmp_str, ' ', ...
                num2str(start_casts(2))]];
        end

        % PLOT BOTTLE DATA IF IT EXISTS
        if ~isempty(DATA.bIND)
            plot(DATA.b.data(:,DATA.bIND), DATA.b.data(:,DATA.ibP), 'ko', ...
                'MarkerSize',6, 'MarkerFaceColor',[255 255 0]/255,...
                'parent',gui.whichAX(4),'MarkerEdgeColor', 'k')
            
            if strcmp(DATA.paramtag, 'O2') && ~isempty(DATA.bIND2)
                plot(DATA.b.data(:,DATA.bIND2), DATA.b.data(:,DATA.ibP), 'ko', ...
                    'MarkerSize',6, 'MarkerFaceColor',[255 255 0]/255,...
                    'parent',gui.whichAX(4),'MarkerEdgeColor', 'r')
                legend_str = [legend_str, 'bottle O2','CTD O2'];
            elseif strcmp(DATA.paramtag, 'PH') && ~isempty(DATA.bIND2)
                plot(DATA.b.data(:,DATA.bIND2), DATA.b.data(:,DATA.ibP), 'ko', ...
                    'MarkerSize',6, 'MarkerFaceColor',[255 255 0]/255,...
                    'parent',gui.whichAX(4),'MarkerEdgeColor', 'r')
                legend_str = [legend_str, 'bottle pH1','bottle pH2'];
            else
                legend_str = [legend_str, 'bottle'];
            end
        end

        L1 = legend(gui.whichAX(4),legend_str, 'Orientation', 'Vertical','interpreter','none');
        L1.Position = [0.8 0.84 0.1 0.05];
%         L1.Units = 'characters';
        L1.Interpreter = 'none';    
        posL=get(L1,'Position');
        posL(1) = 1.05*posL(1);
        posL(2) = 0.825*posL(2);
        set(L1,'Position',posL,'fontunits','normalized','fontsize',10);
        

    end
end

end
