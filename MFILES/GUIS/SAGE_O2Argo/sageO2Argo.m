function varargout = sageO2Argo()
% *************************************************************************
% *************************************************************************
% ** THIS MFILE AND ASSOCIATED SOFTWARE IS COPYRIGHT UNDER MIT LICENSE.  **
% **      PLEASE SEE SAGEO2Argo_MITLisence.txt FOR MORE INFORMATION      **
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% *************************************************************************
%
% ************************************************************************
% sageO2Argo.m
% SOCCOM Assessment and Graphical Evaluation for Oxygen GUI - for use with
% BGC-Argo floats (utilizes BRtraj files for AirO2 data).
% ************************************************************************
%
%
% USE AS:  sageO2Argo
%
% INPUTS:  
%
% OUTPUTS: 
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 9/20/17
% UPDATES: 11/14/17 MIT license added.
%          04/16/18 Added pushbutton to write ODV*QC.TXT file, which
%          includes adjusted oxygen data (average gain viewed in software
%          gets applied to the data and written to file).
% NOTES: Adapted from sageO2 (original version made for MBARI)
% ************************************************************************
%
% ************************************************************************


% Declare shared variables
gui = createInterface();
dirs = SetupDirs_sO2Argo;
DATA = [];
inputs = [];
handlesODVQC = [];


% % Define my colors
% 
mycolors(1,:) = [239 90 17] ./ 255; %orange
mycolors(2,:) = [197 109 39] ./ 255; %burnt orange
mycolors(3,:) = [213 185 31] ./ 255; %mustard
mycolors(4,:) = [21 114 124] ./ 255; %teal
mycolors(5,:) = [46 161 43] ./ 255; %green
mycolors(6,:) = [53 85 40] ./ 255; %forest green
mycolors(7,:) = [24 144 243] ./ 255; %light blue
mycolors(8,:) = [6 62 137] ./ 255; %dark purple/blue


%-------------------------------------------------------------------------%
function gui = createInterface( ~ )
        % Create the user interface for the application and return a
        % structure of handles for global use.
        gui = struct();
        % Open a window and add some menus
        scrsz = get(0,'ScreenSize');
        b=0.1;
        l=0.15;
        w=0.7;
        h=0.8;
        gui.Window = figure( ...
            'Name', 'SAGE-O2 (SOCCOM Assessment and Graphical Evaluation for Oxygen)', ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'Toolbar', 'none', ...
            'HandleVisibility', 'off' ,...
            'Position',[scrsz(3)*l scrsz(4)*b scrsz(3)*w scrsz(4)*h]);
        
        % + Menu items
        helpMenu = uimenu( gui.Window, 'Label', 'Help' );
        uimenu( helpMenu, 'Label', 'Documentation', 'Callback', @onhelp );
        ackMenu = uimenu( gui.Window,'Label','Acknowledgements');
        uimenu(ackMenu, 'Label','Documentation','Callback', @onack);
        
        % ARRANGE THE MAIN INTERFACE
        mainLayout = uiextras.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );
        
        % + Create the panels
        controlPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Define Selections:' );
        gui.ViewPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'View Data: ');
        gui.ViewContainer = uicontainer( ...
            'Parent', gui.ViewPanel );        
        
        % + Adjust the main layout HBoxFlex widths
        set( mainLayout, 'Widths', [-0.5,-2.5]  );
        
        
        % + Create the controls
        BGC = [0.5 0.5 0.5]; %overarching background color

        controlLayout = uix.VBoxFlex( 'Parent', controlPanel, ...
            'Padding', 3, 'Spacing', 3 );
            
            % Control Box 1 (float spec)
            VB1 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);         
            VP1 = uix.Panel('Parent',VB1,'Padding',5,'title','Float Data Specs');
                vb = uix.VBox('Parent',VP1,'Padding',5,'Spacing',3);
                    gui.Fbutton = uicontrol( 'Parent', vb,'Style','pushbutton',...
                    'String','Select Float','fontsize',14,'callback',@on_selectfloat);
                    gui.Mbutton = uicontrol( 'Parent', vb,'Style','pushbutton',...
                        'Enable','off','String','Show Map','fontsize',14,...
                        'callback',@on_showmap);
                    P1 = uix.Panel('Parent',vb,'title','Profile Window');
                        E1 = uix.HBox( 'Parent', P1,'Padding', 5, 'Spacing', 5);
                            gui.profmin = uicontrol('Parent',E1,'Style','edit',...
                                'tag','profmin','callback',@on_profedit);
                            gui.profmax = uicontrol('Parent',E1,'Style','edit',...
                                'tag','profmax','callback',@on_profedit);
                    P2 = uix.Panel('Parent',vb,'title','Pressure Range (dbar)');
                        E2 = uix.HBox( 'Parent', P2, 'Padding', 5, 'Spacing', 5 );
                            gui.Pmin = uicontrol('Parent',E2,'Style','edit',...
                                'tag','Pmin','callback',@on_depthedit);
                            gui.Pmax = uicontrol('Parent',E2,'Style','edit',...
                                'tag','Pmax','callback',@on_depthedit);
                    P3 = uix.Panel('Parent',vb,'title','GLODAP cross-over range (km)');
                        E3 = uix.HBox( 'Parent', P3, 'Padding', 5, 'Spacing', 5 );
                            gui.GLDPkm = uicontrol('Parent',E3,'Style','edit',...
                                'callback',@on_GLODAP);  
           
            % Control Box 2 (Plot Type)
           VB2 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);         
           VP2 = uix.Panel('Parent',VB2,'Padding',5,'title','Plot Type');        
                bbox = uix.VButtonBox('Parent',VP2,'Padding',5,...
                    'ButtonSize',[100,20],'HorizontalAlignment','left');
                    rb2(1) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Surface','tag','1','Value',1,'Callback',@plottype_onClicked ); 
                    rb2(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Deep','tag','2','Value',0,'Callback',@plottype_onClicked ); 
                    rb2(3) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Profile','tag','3','Value',0,'Callback',@plottype_onClicked ); 
                gui.rb2 = rb2;
                            
           % Control Box 3 (Reference Data)
           VB3 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);         
           VP3 = uix.Panel('Parent',VB3,'Padding',5,'title','Reference Data');        
                bbox = uix.VButtonBox('Parent',VP3,'Padding',5,...
                    'ButtonSize',[100,20],'HorizontalAlignment','left');
                    rb3(1) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','NCEP','tag','NCEP','Value',1,'Callback',@ref_onClicked ); 
%                     rb3(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
%                         'String','ERA-INT','tag','ERA_INT','Value',0,'Callback',@ref_onClicked ); 
                    rb3(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','WOA2013','tag','woa','Value',0,'Callback',@ref_onClicked ); 
                gui.rb3 = rb3;
                               
            % Control Box 4 (Oxygen Gain Adjustments)
            VP4 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);        
                P4 = uix.Panel('Parent',VP4,'title','Oxygen Gain Adjustments:');
                    P4v = uix.VBox('Parent',P4,'Padding',3,'Spacing',3);
                        P4mean = uicontrol('Parent',P4v,'Style','pushbutton','string',...
                                'CALC MEAN GAIN','Enable','on','Callback',@on_calcmeanG);
                        P4h = uix.HBox('Parent',P4v,'Padding',5,'Spacing',5);
                            gui.addrow = uicontrol('Parent',P4h,'Style','pushbutton','string',...
                                'ADD ROW','Enable','on','Callback',@on_addrow);
                            gui.removerow = uicontrol('Parent',P4h,'Style','pushbutton','string',...
                                'REMOVE ROW','Enable','on','Callback',@on_removerow);
                        gui.tbl=uitable('Parent',P4v,'Data',[1 1 0],'Enable','on',...
                            'CellEditCallback',@on_celledit);
                        gui.tbl.ColumnName = {'Cycle','Gain','Drift'};
                        gui.tbl.ColumnEditable = [true true true];
                        gui.tbl.ColumnWidth = {50 50 50};
%                         P4b = uicontrol('parent',P4v,'style','pushbutton','string',...
%                             'Reload DAC QC','Callback',@on_reloadQC,'Enable','off');

            % Control Box 5 (APPLY Oxygen Gain Adjustment)
            VP5 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);  
                P5 = uix.Panel('Parent',VP5,'title','Apply Oxygen Gain Adjustment:');
                    gui.doQCbutton = uicontrol('parent',P5,'style','pushbutton','string',...
                            'Write ODV*QC.TXT','fontsize',14,'Callback',@on_applyO2gain,'Enable','on');



                
        set( controlLayout, 'Heights', [-5.5 -2.5 -2.5 -5 -2] ); %Main control vbox heights
        set( vb, 'Heights', [-1 -1 -1.75 -1.75 -1.75] ); %Float spec box heights
        set( P4v, 'Heights', [-1 -1 -3] ); %QC adj box heights
        
        % + Create the view
        p = gui.ViewContainer;
        t = uiextras.TabPanel('Parent',p,'BackgroundColor','w');
                t1 = uix.VBox('Parent',t,'padding',10,'spacing',0.5);
                    t1a = uicontainer('Parent',t1);
                        t1ax(1) = axes('Parent',t1a,'Visible','off');
                    t1b = uicontainer('Parent',t1);
                        t1ax(2) = axes('Parent',t1b,'Visible','off');
                    t1c = uicontainer('Parent',t1);
                        t1ax(3) = axes('Parent',t1c,'Visible','off');
                    t1d = uicontainer('Parent',t1);
                        t1ax(4) = axes('Parent',t1d,'Visible','off');
                        t1ax(5) = axes('Parent',t1d,'Visible','off');
            gui.t1 = t1;
            t2 = uix.VBox('Parent',t,'padding',10,'spacing',0.5);
                t2a = uicontainer('Parent',t2);
                    t2ax(1) = axes('Parent',t2a,'Visible','off');
                t2b = uicontainer('Parent',t2);
                    t2ax(2) = axes('Parent',t2b,'Visible','off');
                t2c = uicontainer('Parent',t2);
                    t2ax(3) = axes('Parent',t2c,'Visible','off');
                t2d = uicontainer('Parent',t2);
                    t2ax(4) = axes('Parent',t2d,'Visible','off');
                    t2ax(5) = axes('Parent',t2d,'Visible','off');
            gui.t2=t2;
            t.TabNames = {'Raw','QC'};
            t.TabSize = 75;
            t.FontSize = 16;
            t.SelectedChild = 1;
            gui.t = t;
            gui.t1ax = t1ax;
            gui.t2ax = t2ax;
        % POSITION AXES, MOVE PLOTS LEFT SO ROOM ON RIGHT FOR STATS
        % RAW DATA TAB POSITIONING:
        pos1 = get(gui.t1ax(1),'Position');
        pos1(3) = 0.8*pos1(3);
        pos1(4) = 0.91*pos1(4);
        pos1(2) = 0.9*pos1(2);
        set(gui.t1ax(1),'Position',pos1);
        pos2 = get(gui.t1ax(2),'Position');
        pos2(3) = 0.8*pos2(3);
        pos2(4) = 0.91*pos2(4);
        set(gui.t1ax(2),'Position',pos2);
        pos3 = get(gui.t1ax(3),'Position');
        pos3(3) = 0.8*pos3(3);
        pos3(2) = 1.25*pos3(2);
        pos3(4) = 0.91*pos3(4);
        set(gui.t1ax(3),'Position',pos3);
        pos4 = get(gui.t1ax(4),'Position');
        pos4(3) = 0.8*pos4(3);
        pos4(2) = 1.8*pos4(2);
        pos4(4) = 0.91*pos4(4);
        set(gui.t1ax(4),'Position',pos4);
        pos5 = get(gui.t1ax(5),'Position');
        pos5(3) = 0.8*pos5(3);
        pos5(2) = 1.8*pos5(2);
        pos5(4) = 0.91*pos5(4);
        set(gui.t1ax(5),'Position',pos5);
        gui.pos2 = pos2;
        % QC DATA TAB POSITIONING
        pos1 = get(gui.t2ax(1),'Position');
        pos1(3) = 0.8*pos1(3);
        pos1(4) = 0.91*pos1(4);
        pos1(2) = 0.9*pos1(2);
        set(gui.t2ax(1),'Position',pos1);
        pos2 = get(gui.t2ax(2),'Position');
        pos2(3) = 0.8*pos2(3);
        pos2(4) = 0.91*pos2(4);
        set(gui.t2ax(2),'Position',pos2);
        pos3 = get(gui.t2ax(3),'Position');
        pos3(3) = 0.8*pos3(3);
        pos3(2) = 1.25*pos3(2);
        pos3(4) = 0.91*pos3(4);
        set(gui.t2ax(3),'Position',pos3);
        pos4 = get(gui.t2ax(4),'Position');
        pos4(3) = 0.8*pos4(3);
        pos4(2) = 1.8*pos4(2);
        pos4(4) = 0.91*pos4(4);
        set(gui.t2ax(4),'Position',pos4);
        pos5 = get(gui.t2ax(5),'Position');
        pos5(3) = 0.8*pos5(3);
        pos5(2) =1.8*pos5(2);
        pos5(4) = 0.91*pos5(4);
        set(gui.t2ax(5),'Position',pos5);
        gui.pos2 = pos2;
        
    end % createInterface
%-------------------------------------------------------------------------%
    function updateInterface()
        % Update various parts of the interface in response to the inputs
        % being changed.
        if isfield(DATA,'GAINS')
           DATA=rmfield(DATA,'GAINS');
        end
        set(gui.profmin,'String',num2str(inputs.profedit(1)));
        if max(DATA.PTSdata(:,2)) < inputs.profedit(2);
            set(gui.profmax,'String',num2str(max(DATA.PTSdata(:,2))));
        else
            set(gui.profmax,'String',num2str(inputs.profedit(2))); %last cast #
        end
        set(gui.Pmin,'String',num2str(inputs.depthedit(1))); 
        set(gui.Pmax,'String',num2str(inputs.depthedit(2))); 
        if max(DATA.PTSdata(:,5)) < inputs.depthedit(2);
            set(gui.Pmax,'String',num2str(floor(max(DATA.PTSdata(:,5)))));
        end
        if max(DATA.PTSdata(:,5)) < inputs.depthedit(1);
            set(gui.Pmin,'String',num2str(floor(max(DATA.PTSdata(:,5)))-100));
        end

        set(gui.GLDPkm,'String',num2str(inputs.GLDPkm));
        depthmin = str2double(get(gui.Pmin,'String'));
        depthmax = str2double(get(gui.Pmax,'String'));
        PROFmin = str2double(get(gui.profmin,'String'));
        PROFmax = str2double(get(gui.profmax,'String'));
        % SET X LIMITS AND TICK LOCATIONS FOR PLOTTING (AS DATE AND CYCLE#)
        DATA.xlims{1} = [DATA.track(DATA.track(:,2)==PROFmin,1) DATA.track(DATA.track(:,2)==PROFmax,1)];
        DATA.xlims{2} = [DATA.track(DATA.track(:,2)==PROFmin,2) DATA.track(DATA.track(:,2)==PROFmax,2)];
        if nanmax(DATA.xlims{2}(2)) >= 6 % at least 6 cycles
            xtckfac{2} = floor((DATA.xlims{2}(2)-DATA.xlims{2}(1))/5);
        else
            xtckfac{2} = 1; %else increment x-axis by 1 cycle
        end
        DATA.xticks{2} = DATA.xlims{2}(1):xtckfac{2}:DATA.xlims{2}(2);
        [~,xti,~] = intersect(DATA.track(:,2),DATA.xticks{2});
        DATA.xticks{1} = DATA.track(xti,1);
        % GET DATA SUBSETS
        DATA.O2subset = DATA.RAWorQCprof(DATA.RAWorQCprof(:,4) >= depthmin & ...
            DATA.RAWorQCprof(:,4) <= depthmax & DATA.RAWorQCprof(:,2) >= PROFmin &...
            DATA.RAWorQCprof(:,2) <= PROFmax,:);
        DATA.PTSsubset = DATA.PTSdata(DATA.PTSdata(:,5) >= depthmin & ...
            DATA.PTSdata(:,5) <= depthmax & DATA.PTSdata(:,2) >= PROFmin &...
            DATA.PTSdata(:,2) <= PROFmax,:);
        DATA.refsub = DATA.refdata(DATA.refdata(:,2) >= PROFmin & ...
            DATA.refdata(:,2) <= PROFmax,:);
        DATA.SATsubset = DATA.RAWorQCwoa(DATA.RAWorQCwoa(:,1) >= PROFmin & ... 
            DATA.RAWorQCwoa(:,1) <= PROFmax,:);
        DATA.iswoaref = get(gui.rb3(2),'Value');
        if DATA.iswoaref == 1
            DATA.GAINS = DATA.refsub(:,4)./DATA.SATsubset(:,2);
        else
            tmp = find(~cellfun(@isempty,DATA.O2air{1}));
            DATA.howmanyairs = length(tmp);
            for  s = 1:DATA.howmanyairs
                DATA.airsub{1}{s} = DATA.RAWorQCair{1}{s}(DATA.RAWorQCair{1}{s}(:,2) >=PROFmin & ...
                    DATA.RAWorQCair{1}{s}(:,2) <= PROFmax,:);
                DATA.airsub{2}{s} = DATA.RAWorQCair{2}{s}(DATA.RAWorQCair{2}{s}(:,2) >=PROFmin & ...
                    DATA.RAWorQCair{2}{s}(:,2) <= PROFmax,:);
%                 DATA.GAINS{s} = DATA.refsub(:,4)./DATA.airsub{1}{s}(:,9);
                [cc,iia,iib] = intersect(DATA.refsub(:,2),DATA.airsub{1}{s}(:,2));
                DATA.GAINS{s} = DATA.refsub(iia,4)./DATA.airsub{1}{s}(iib,9);
                DATA.GAINStime{s} = DATA.refsub(iia,1); %for use in plotting climos
            end
        end
        % Now further subset the data if plotting only night-time cals;
        if strcmp(inputs.AMPM,'all') ~= 1 % if inputs.AMPM = 'PM'
            [~,El] = SolarAzElq(DATA.O2subset(:,1),DATA.O2subset(:,3),DATA.O2subset(:,4),zeros(size(DATA.O2subset,1),1));
            XEl = find(El<5);
            DATA.O2subset = DATA.O2subset(XEl,:);
            Xunq = unique(DATA.O2subset(:,2));
            DATA.PTSsubset = DATA.PTSsubset(XEl,:);
            [~,ia,~] = intersect(DATA.SATsubset(:,1),Xunq);
            DATA.SATsubset = DATA.SATsubset(ia,:);
            [~,ia,~] = intersect(DATA.refsub(:,2),Xunq);
            DATA.refsub = DATA.refsub(ia,:);
            if isfield(DATA,'airsub') == 1
                for ai = 1:length(DATA.airsub)
                    DATA.airsub{ai}{1}= DATA.airsub{ai}{1}(ia,:);
                end
            end
            if iscell(DATA.GAINS) == 0
                DATA.GAINS = DATA.GAINS(ia,:);
            else
                for gi = 1:length(DATA.GAINS)
                    DATA.GAINS{gi} = DATA.GAINS{gi}(ia,:);
                end
            end
        end
    end % updateInterface

%-------------------------------------------------------------------------%
    function onhelp( ~, ~ )
        % User has asked for the documentation
        open([dirs.mfiles,'GUIS',filesep,'SAGE_O2Argo',filesep,'README_SageO2Argo.txt'])
    end % onHelp

%-------------------------------------------------------------------------%
    function onack( ~, ~ )
        % User has asked for the documentation
        open([dirs.mfiles,'GUIS',filesep,'SAGE_O2Argo',filesep,'acknowledgements_sO2Argo.txt'])
    end % onHelp

%-------------------------------------------------------------------------%
    function on_selectfloat( ~, ~ )
        % select float data dir from dialog box
        inputs.rorq = 1;
        x=1;
        while x==1
            [pn] = uigetdir([dirs.Argo],'SELECT FLOAT');
            x=0;
            pattern='\w*';
            [Idx1,debut,fin] = regexp(pn, pattern,'match','start','end');
            inputs.floatID = cell2mat(Idx1(end));
        end
        fp = filesep; % File separator for current platform
        thedatadir = [pn,fp];
        handlesODVQC.info.file_path = thedatadir;
        set( gui.Fbutton,'String','Loading Data ...');
        gui.wrk_color = gui.Fbutton.BackgroundColor;
        set(gui.Fbutton,'BackgroundColor','y');
        set( gui.doQCbutton,'String','Write ODV*QC.TXT','fontsize',16)
        set(gui.doQCbutton,'BackgroundColor',gui.wrk_color);
        gui.rb2(1).Value = 1; % reset to surface plot when float opens
        gui.rb2(2).Value = 0;
        gui.rb2(3).Value = 0;
        drawnow
        set( gui.Fbutton,'String',inputs.floatID);
        DATA = getall_floatdata_sO2Argo(thedatadir,inputs.floatID);
        set(gui.Mbutton,'Enable','on');
        %inputs.depthedit = [1480 1520]; %default pressure range (deep)
        inputs.depthedit = [0 30]; %default pressure range (deep)
        inputs.GLDPkm = 30;
        inputs.profedit = [1 DATA.track(end,2)]; 
        inputs.AMPM = 'all';
        gui.TP.SelectionChangedFcn = @on_PlottypeChanged;
        gui.t.SelectedChild = 1; %default to Raw upon float selection
        gui.t.SelectionChangedFcn = @on_RAWorQC;
%         gui.TP.Selection = 2;   %default to Deep tab 
        %set(gui.rb2(2),'Value',1);
%         set(gui.rb5,'Value',0);
        DATA.RAWorQCprof = DATA.O2data{1}; %profile data  
        DATA.RAWorQCair = DATA.O2air;
        gui.whichAX = gui.t1ax;
        if iscell(DATA.O2air{1})==1
            mytemp = find(~cellfun(@isempty,DATA.O2air{1}));
            DATA.howmanyairs = length(mytemp);
        else
            DATA.howmanyairs = 0;
        end
        %get QC adjustment data
        inputs.qc_path = [dirs.QCadj,DATA.floatNAME,'_FloatQCList.txt'];
        DATA.QCA = get_QCA(inputs.qc_path,DATA.floatNAME);
        %Calculate WOA data and surface o2 sat for all floats
        set( gui.Fbutton,'String','Loading WOA ...');
        set(gui.Fbutton,'BackgroundColor','y');
        drawnow
        [~,inputs.intersect_cycles_WOA,~] = intersect(DATA.track(:,2),DATA.O2data{1}(:,2));
        Wtrack = [DATA.track(inputs.intersect_cycles_WOA,1) DATA.track(inputs.intersect_cycles_WOA,4) DATA.track(inputs.intersect_cycles_WOA,3)];
        try
            Wdata = get_WOA2013_local_sO2(Wtrack, [0 2000], 'O2sat',dirs.user_dir);
            zsurf = Wdata.Z<=25;
            WOA_surf = Wdata.d(zsurf,:);
            DATA.WOAsurf = nanmean(WOA_surf,1);
        catch
            msgbox({'ERROR: getWOA failed.',...
                'Must use NCEP or bottle data for reference.'})
            set(gui.rb3(1),'Value',0,'Enable','off');
            DATA.WOAsurf = [];
        end
        FLT = DATA.O2data{1} ;% sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat
        FLT_surf = FLT(FLT(:,4)<=25,:);
        cst = unique(FLT(:,2));
        DATA.SURF_SAT=[];
        for i = 1:length(cst)
            cstnum = cst(i);
            FLT_tmp = FLT_surf(FLT_surf(:,2)==cstnum,11);
            DATA.SURF_SAT(i,1) = cstnum;
            DATA.SURF_SAT(i,2) = nanmean(FLT_tmp);
            DATA.SURF_SAT(i,3) = nanstd(FLT_tmp);
        end
        DATA.RAWorQCwoa = DATA.SURF_SAT;
        if isempty(DATA.WOAsurf)
            DATA.rawGAINS_WOA{1} = [];
        else
            DATA.rawGAINS_WOA{1} = DATA.WOAsurf'./DATA.RAWorQCwoa(:,2);
        end
        DATA.rawGAINS_WOA{3} = [];
%         DATA.O2air{1}
        if ~isempty(DATA.O2air{1}) && sum(~isnan(DATA.O2air{1}{1}(:,10)))>0 %AIRCAL EXISTS SO USE NCEP OR ERA
            [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2air{1}{1}(:,2));
            %______________________________________________________________
            %TRY NCEP FIRST (REALTIME).  DATA IS GRABBED FROM THE WEB.  IF WEBSITE IS
            %DOWN FOR MAINTENANCE ETC WILL TRY ERA NEXT.
            try
                set( gui.Fbutton,'String','Loading NCEP ...');
                set(gui.Fbutton,'BackgroundColor','y');
                drawnow
                DATA.NCEP = getNCEP(DATA.track(inputs.intersect_cycles,1),DATA.track(inputs.intersect_cycles,3),DATA.track(inputs.intersect_cycles,4),dirs);
            catch
                msgbox({'ERROR: getNCEP failed.',...
                    'Must use ERA or WOA for reference.'})
                set(gui.rb3(1),'Value',0,'Enable','off');
                DATA.NCEP.PRES = [];
            end
            if ~isempty(DATA.NCEP.PRES)
                set(gui.rb3(1),'Value',1,'Enable','on');
%                 set(gui.rb3(2),'Value',0,'Enable','on');
                set(gui.rb3(2),'Value',0);
                DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
                DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
                DATA.refdata(:,3) = DATA.NCEP.PRES;
                DATA.reftag = 'NCEP';
                %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
                if ~isempty(DATA.O2air{1}{3})
                    DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - mean([DATA.RAWorQCair{1}{1}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2)).*0.20946;
                else
                    DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
                end
%                 DATA.refdata(1:5,3)
%                 DATA.RAWorQCair{1}{1}(1:5,10)
                inputs.cyEND = nanmax(DATA.refdata(:,2));
                %Maintain raw gains for adjustment calcs in
                %calc_gain_ADJtable_sO2Argo.m and for toggling between reference
                %datasets
                if DATA.howmanyairs > 1
                    DATA.rawGAINS_NCEP{1} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{1}(:,9);
                    DATA.rawGAINS_NCEP{3} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{3}(:,9);
                    useGAIN = DATA.rawGAINS_NCEP{3};
                else
                    DATA.rawGAINS_NCEP{1} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{1}(:,9);
                    DATA.rawGAINS_NCEP{3} = [];
                    useGAIN = DATA.rawGAINS_NCEP{1};
                end
                DATA.GAINS=DATA.rawGAINS_NCEP;
                DATA.rawGAINS=DATA.rawGAINS_NCEP;
                inputs.isprof = 0;
                %populate gain adjustment table
                cyST = 1;
                DATA.iswoaref = get(gui.rb3(2),'Value');
                if isempty(DATA.QCA.O2)
                    DATA.tableDATA = [cyST nanmean(useGAIN) 0];
                    [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
                    set(gui.tbl,'Data',DATA.tableDATA)
                    msgbox({'FYI: AVG GAIN POPULATES ADJUSTMENT TABLE BY DEFAULT.'});
                else
%                      DATA.tableDATA = [cyST DATA.QCA.O2(2) 0];
                    DATA.tableDATA = [DATA.QCA.O2(:,1) DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
                    cyST = DATA.QCA.O2(:,1);
                    inputs.cyEND= nanmax(DATA.refdata(:,2));
                    if length(cyST)>1
                        inputs.cyEND = [cyST(2:end);inputs.cyEND];
                    end
                    [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
                    set(gui.tbl,'Data',DATA.tableDATA)
                end
                DATA.BIC = calcBIC_sO2Argo(gui,DATA);
                updateInterface()
                DATA = redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
            end
            
            %______________________________________________________________
%***************************************************************************************
% %         ERA FUNCTIONALITY HAS BEEN DEPRECATED; ECMWF HAS MADE IT MORE
%           DIFFICULT TO GRAB LARGE AMOUNTS OF ERA-INTERIM PRODUCTS IN AUTO/REAL TIME.  WE
%           NOW FOCUS ONLY ON NCEP FOR IN-AIR COMPARISONS.  FEEL FREE TO REINSTATE
%           THIS FUNCTIONALITY BUT WILL REQUIRE EFFORT IN DATA RETRIEVAL.
%***************************************************************************************
% %             %TRY ERA DATA NEXT.  DATA IS DOWNLOADED; WILL WORK UNLESS
% %             %FAIRLY NEW FLOAT (ERA DATASET NOT REALTIME).
% %             try
% %                 DATA.ERA = getERA_sO2Argo(dirs.ERA,DATA.track(inputs.intersect_cycles,1),DATA.track(inputs.intersect_cycles,3),DATA.track(inputs.intersect_cycles,4));
% %             catch
% %                 msgbox({'ERROR: getERA failed.',...
% %                     'Possibly no data for time period of interest. Must use NCEP or WOA for reference.'})
% % %                 set(gui.rb3(2),'Value',0,'Enable','off')
% %                 DATA.ERA.PRES=[];
% %             end
% %             if ~isempty(DATA.ERA.PRES)
% %                  %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
% %                 if ~isempty(DATA.O2air{1}{3})
% %                     ERAG = (DATA.ERA.PRES./100 - mean([DATA.RAWorQCair{1}{1}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2)).*0.20946;
% %                 else
% %                     ERAG = (DATA.ERA.PRES./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
% %                 end
% %                 %Maintain raw gains for adjustment calcs in
% %                 %calc_gain_ADJtable_sO2Argo.m and for toggling between reference
% %                 %datasets
% %                 if DATA.howmanyairs > 1
% %                     DATA.rawGAINS_ERA{1} = ERAG./DATA.RAWorQCair{1}{1}(:,9);
% %                     DATA.rawGAINS_ERA{3} = ERAG./DATA.RAWorQCair{1}{3}(:,9);
% %                     useGAIN=DATA.rawGAINS_ERA{3};
% %                 else
% %                     DATA.rawGAINS_ERA{1} = ERAG./DATA.RAWorQCair{1}{1}(:,9);
% %                     DATA.rawGAINS_ERA{3} = [];
% %                     useGAIN=DATA.rawGAINS_ERA{1};
% %                 end
% %                 if isempty(DATA.NCEP.PRES) %SET DEFAULT TO ERA ONLY IF NCEP DIDN'T WORK
% %                     set(gui.rb3(2),'Value',1,'Enable','on');
% %                     set(gui.rb3(3),'Value',0);
% %                     DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
% %                     DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
% %                     DATA.refdata(:,3) = DATA.ERA.PRES;
% %                     DATA.reftag = 'ERA';
% %                     DATA.refdata(:,4)=ERAG;
% %                     DATA.rawGAINS=DATA.rawGAINS_ERA;
% %                     DATA.GAINS = DATA.rawGAINS_ERA;
% %                     inputs.isprof = 0;
% %                     %populate gain adjustment table
% %                     cyST = 1;
% %                     inputs.cyEND = max(DATA.refdata(:,2));
% %                     if isempty(DATA.QCA.O2)
% %                         DATA.tableDATA = [cyST nanmean(useGAIN) 0];
% %                         [DATA.tableDATA,DATA.brkRESIDS] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
% % %                         set(gui.tbl,'Data',DATA.tableDATA)
% %                         msgbox({'FYI: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED','POPULATING ADJUSTMENT TABLE WITH DEFAULT.'});
% %                     else
% %                         DATA.tableDATA = [cyST DATA.QCA.O2(2) 0];
% %                         [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
% % %                         set(gui.tbl,'Data',DATA.tableDATA)
% %                     end
% %                     DATA.AIC = calcAIC_sO2Argo(gui,DATA);
% %                     updateInterface()
% %                     redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs)
% %                 end
% %             end
        end %end if aircal exists
        %__________________________________________________________________
        % NOW SET DEFAULT TO WOA IF NO AIRCAL OR IF BOTH NCEP AND ERA
        % FAILED (WEB-GET NOT WORKING, AND NEWER FLOAT)
%         if isempty(DATA.O2air{1}) || (isempty(DATA.NCEP.PRES) && isempty(DATA.ERA.PRES))
        if isempty(DATA.O2air{1}) || (isfield(DATA,'NCEP') && isempty(DATA.NCEP.PRES)) || sum(~isnan(DATA.O2air{1}{1}(:,10)))==0
            if isempty(DATA.O2air{1}) || sum(~isnan(DATA.O2air{1}{1}(:,10)))==0
                msgbox('NO USABLE IN-AIR DATA.  CALIBRATE TO WOA.');
            end
            set(gui.rb3(1),'Value',0,'Enable','off');
%             set(gui.rb3(2),'Value',0,'Enable','off');
            set(gui.rb3(2),'Value',1);
            DATA.iswoaref = get(gui.rb3(2),'Value');
            if isfield(DATA,'refdata')
                DATA=rmfield(DATA,'refdata');
            end
            DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles_WOA,1);
            DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles_WOA,2);
            DATA.refdata(:,3) = nan(length(DATA.refdata(:,2)),1);
            DATA.refdata(:,4) = DATA.WOAsurf;
            inputs.cyEND = nanmax(DATA.refdata(:,2));
            DATA.rawGAINS=DATA.rawGAINS_WOA;
            DATA.GAINS=DATA.rawGAINS_WOA; %maintain raw gains for adjustment calcs in calc_gain_ADJtable_sO2Argo.m
            inputs.isprof = 0;
            %populate gain adjustment table
            
            %cyST = DATA.QCA.O2(:,1); %This will only work if adjustment table already exists -jp 12/18/18
            inputs.cyEND= nanmax(DATA.refdata(:,2));
            
            if isempty(DATA.QCA.O2) % NO QC yet
                cyST = 1;
                DATA.tableDATA = [cyST nanmean(DATA.GAINS{1}) 0];
                [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
                set(gui.tbl,'Data',DATA.tableDATA)
                msgbox({'FYI: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED','POPULATING ADJUSTMENT TABLE WITH DEFAULT.'});
            else
                cyST = DATA.QCA.O2(:,1);
                if length(cyST)>1 % this block seems like it could put in end value twice?
                    inputs.cyEND = [cyST(2:end);inputs.cyEND];
                end
                DATA.tableDATA = [cyST DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
                disp(DATA.tableDATA)
                [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
                set(gui.tbl,'Data',DATA.tableDATA)
            end
            DATA.BIC = calcBIC_sO2Argo(gui,DATA);
            updateInterface()
            DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);
        end % end if isempty aircal
        
    %__________________________________________________________________
    % PRELOAD BOTTLE DATA FOR PROFILE PLOTS
    % Get bottle data
    Blookup = [dirs.bottle,'BottleData_lookup_table.txt'];
    % LOAD BOTTLE DATA LOOK UP TABLE & STORE IN handles STRUCTURE
    % LOOK UP TABLE  HEADER = [UW_ID  WMO   CRUISE   STN   CAST   Data file]
    if exist(Blookup,'file')
        fid = fopen(Blookup);
        d   = textscan(fid, '%s %s %s %s %f %f %s','Delimiter', '\t','HeaderLines',1);
        fclose(fid);
        bottle_lookup = d;
        inputs.floatID
        ind = strcmp(num2str(str2num(inputs.floatID)), bottle_lookup{1,2});
        if sum(ind) > 0; % float exists in lookup table
            bottle_fname = bottle_lookup{1,7}{ind};
            stn = bottle_lookup{1,5}(ind);

            cst = bottle_lookup{1,6}(ind);

            % data file & data exist for float
            if ~isempty(bottle_fname) && stn ~= -999 && cst ~= -999
                tic
                d = get_shipboard_data([dirs.bottle,bottle_fname]);
                toc
                iStn  = find(strcmp(d.hdr,'STNNBR') == 1);
                iCast = find(strcmp(d.hdr,'CASTNO') == 1);
                tStn  = d.data(:,iStn)  == stn;
                tCast = d.data(:,iCast) == cst;
                DATA.bdata.cruise  = d.cruise;
                DATA.bdata.hdr     = d.hdr;
                %handles.bdata.units   = d.units;
                d.data(d.data == -999) = NaN;
                DATA.bdata.data    = d.data(tStn&tCast,:);
            else
                DATA.bdata.Cruise = '';
                DATA.bdata.hdr    = '';
                DATA.bdata.data   = [];
            end
            clear bottle_fname stn cst d IStn iCast tStn tCast
        else
            DATA.bdata.Cruise = '';
            DATA.bdata.hdr    = '';
            DATA.bdata.data   = [];
        end
    end
    
    %set( gui.Fbutton,'String',inputs.floatID)
    set( gui.Fbutton,'String',DATA.floatNAME)
    set(gui.Fbutton,'BackgroundColor',gui.wrk_color);
    drawnow
    end %end selectfloat

%-------------------------------------------------------------------------%
    function on_showmap( ~,~ )
        track = DATA.track; % float track
        tnan = track(:,3) == 1e-10;
        track(tnan,3:4) = NaN;
        lon_limits = [min(track(:,3))-10 max(track(:,3))+10];
        lat_limits = [min(track(:,4))-6 max(track(:,4))+6];
        handles.float_track = figure;
        set(gcf,'Position', [350 50 860 620]);
        DATA.CMAP = load('SOCCOM_NCP_cmap.mat');
        DATA.CMAP = DATA.CMAP.cmap;
        colormap(handles.float_track,DATA.CMAP);
        if max(track(:,4)) > -30 % Don't use stereographic projection
            m_proj('lambert','lon',lon_limits,'lat',lat_limits);
        else
            m_proj('stereographic','latitude',-90,'radius',60,'rotangle',45);
        end
        m_tbase('contourf','edgecolor','none');
        legend_cell = {};
        hold on
        m_plot(track(:,3),track(:,4),'ko-', 'MarkerSize',3)%     xlim(lon_limits);
        title(['FLOAT ',inputs.floatID],'FontSize', 16)
        hold on
        Hm1 = m_plot(track(1,3),track(1,4),'o', 'MarkerSize',10, 'MarkerFaceColor','g', ...
            'MarkerEdgeColor','k');
        Hm2 = m_plot(track(end,3),track(end,4),'o', 'MarkerSize',10, 'MarkerFaceColor','r', ...
            'MarkerEdgeColor','k');
        whos Hm
        legend_cell = [legend_cell, 'first', 'last'];
        colorbar

        track_legend = legend([Hm1 Hm2],legend_cell);
    %     track_legend.Orientation = 'Horizontal';
        track_legend.Location = 'northeastoutside';
        set(gca,'ydir','normal');
        m_grid('linewidth',2,'tickdir','out','xaxislocation','top');
    end

%-------------------------------------------------------------------------%
    function on_RAWorQC( source,~ )
       inputs.rorq = get(source,'selectedchild'); %'1' is raw, '2' is QC
       [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
       updateInterface()
       if   inputs.isprof == 1; %profile selected?
            redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
       elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1;
           DATA = redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
       else
           DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);  
       end       
    end
 
%-------------------------------------------------------------------------%        
   function on_profedit( source, ~ ) 
       PEtag = get(source,'tag');
       PE = get(source,'String');
       if (strcmp(PEtag,'profmin')) == 1
           inputs.profedit(1,1) = str2double(PE);
       else
           inputs.profedit(1,2) = str2double(PE);
       end
       updateInterface()
       if   inputs.isprof == 1; %profile selected?
            redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
       elseif ~isempty(DATA.O2air{1})&& get(gui.rb3(2),'Value')~=1;
           DATA = redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
       else
           DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);   
       end
   end

%-------------------------------------------------------------------------%        
   function on_depthedit( source, ~ ) 
       DEtag = get(source,'tag');
       DE = get(source,'String');
       if (strcmp(DEtag,'Pmin')) == 1
           inputs.depthedit(1,1) = str2double(DE);
       else
           inputs.depthedit(1,2) = str2double(DE);
       end
       updateInterface()
       if   inputs.isprof == 1; %profile selected?
           redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
       elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
           DATA = redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
       else
           DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);   
       end
   end

%-------------------------------------------------------------------------%        
   function on_GLODAP( source, ~ ) 
       GDkm = get(source,'String');
       inputs.GLDPkm = str2double(GDkm);
       updateInterface()
       if gui.TP.Selection == 3
           redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
       end
   end

%-------------------------------------------------------------------------%
   function plottype_onClicked( source, ~ ) 
        source.Value = 1; % select this 
        RB2 = gui.rb2;
        set( RB2(RB2~=source), 'Value', 0 ) % unselect others 
        % get track data
        reftag = get(source,'tag');
        SC = str2num(reftag);
%         SC = get(source,'Selection'); %tabs: (1)surf, (2)depth, (3)profile
        Dedit = [0 30; 1480 1520; 0 1600];
        inputs.depthedit = Dedit(SC,:);
        if get(gui.rb2(3),'Value')==1
            inputs.isprof = 1; %profile selected?
            updateInterface()
            redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
        else
            inputs.isprof = 0;
            updateInterface()
            if ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1;
               DATA =  redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
            else
                DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);   
            end
        end
   end

%-------------------------------------------------------------------------%
%     function DD_onClicked( source, ~ ) 
%         tabnum = gui.TP.SelectedChild;
%         r2=gui.rb2;
%         source.Value = 1; % select this 
%         set( r2(tabnum,r2(tabnum,:)~=source), 'Value', 0 ) %unselect others 
%         % do something useful 
%     end

%-------------------------------------------------------------------------%
    function ref_onClicked( source, ~ ) 
        r3=gui.rb3;
        source.Value = 1; % select this 
        set( r3(r3~=source), 'Value', 0 ) % unselect others 
        % get track data
        reftag = get(source,'tag');
        % NCEP_____________________________________________________________ 
        [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2air{1}{1}(:,2));
        if (strcmp(reftag,'NCEP')) == 1
            if isfield(DATA,'refdata')
                DATA=rmfield(DATA,'refdata');
            end
            DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
            DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
            DATA.refdata(:,3) = DATA.NCEP.PRES;
            DATA.reftag = 'NCEP';
            %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
            if inputs.rorq == 2 %QC tab
               DATA.RAWorQCair = DATA.O2air;
            end
            if ~isempty(DATA.O2air{1}{3});
                DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - mean([DATA.RAWorQCair{1}{1}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2)).*0.20946;
            else
                DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
            end
            if inputs.rorq == 2 %QC tab
                [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
            end
            updateInterface()
            if   inputs.isprof == 1; %profile selected?
                redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
            else
                DATA = redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
            end
% %         % ERA______________________________________________________________    
% %         elseif (strcmp(reftag,'ERA_INT')) == 1
% %             [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2air{1}{1}(:,2));
% %             if isempty(DATA.ERA.PRES)
% %                 set(gui.rb3(1),'Value',1)
% %                 set(gui.rb3(2),'Value',0,'Enable','off')
% %                 msgbox({'WARNING: No ERA data exists for time period of interest.',...
% %                     'Must use NCEP or WOA as refernce.'});
% %             else
% %                 if isfield(DATA,'refdata')
% %                     DATA=rmfield(DATA,'refdata');
% %                 end
% %                 DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
% %                 DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
% %                 DATA.refdata(:,3) = DATA.ERA.PRES;
% %                 DATA.reftag = 'ERA';
% %                 %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
% %                 if inputs.rorq == 2 %QC tab
% %                     DATA.RAWorQCair = DATA.O2air;
% %                 end
% %                 if ~isempty(DATA.O2air{1}{3});
% %                     DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - mean([DATA.RAWorQCair{1}{1}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2)).*0.20946;
% %                 else
% %                     DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
% %                 end
% %                 if inputs.rorq == 2 %QC tab
% %                     [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
% %                 end
% %                 updateInterface()
% %                 if   inputs.isprof == 1; %profile selected?
% %                     redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs)
% %                 else
% %                     redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs)
% %                 end
% %             end
        % WOA______________________________________________________________    
        else 
            if isfield(DATA,'refdata')
                DATA=rmfield(DATA,'refdata');
            end
            [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2data{1}(:,2));
            DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
            DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
            Wtrack = [DATA.track(inputs.intersect_cycles,1) DATA.track(inputs.intersect_cycles,4) DATA.track(inputs.intersect_cycles,3)];
            Wdata = get_WOA2013_local_sO2Argo(Wtrack, [0 2000], 'O2sat',dirs.user_dir);
            zsurf = Wdata.Z<=25;
            WOA_surf = Wdata.d(zsurf,:);
            DATA.WOAsurf = nanmean(WOA_surf,1);
            DATA.refdata(:,3) = nan(length(DATA.refdata(:,2)),1);
            DATA.refdata(:,4) = DATA.WOAsurf;
            updateInterface()
            if   inputs.isprof == 1; %profile selected?
                redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
            else
                DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);
            end
        end
    end
    
 %-------------------------------------------------------------------------%
    function on_calcmeanG( source, ~ ) 
        if DATA.iswoaref == 0
            if ~isempty(DATA.rawGAINS{3})
                useG = nanmean(DATA.rawGAINS{3});
            else
                useG = nanmean(DATA.rawGAINS{1});
            end
        elseif DATA.iswoaref == 1 % WOA button is toggled
            useG = nanmean(DATA.rawGAINS_WOA{1});
        end
        DATA.tableDATA = [1 useG 0];
        cyST = DATA.tableDATA(1);
        inputs.cyEND= nanmax(DATA.refdata(:,2));
        [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
%             DATA.tableDATA = get(gui.tbl,'Data');
        set(gui.tbl,'Data',DATA.tableDATA);
        [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
        DATA.BIC = calcBIC_sO2Argo(gui,DATA);
        updateInterface()
        if   inputs.isprof == 1 %profile selected?
             redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs)
        elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
            redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs)
        else
            redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs)   
        end
    end 

%-------------------------------------------------------------------------%
    function on_addrow( source, ~ ) 
        celldata = gui.tbl.Data;
        newstrt = celldata(end,1)+1;
        newrow = [newstrt 0 0];
        new_celldata = [celldata;newrow];
        new_ends = [new_celldata(2:end,1);inputs.cyEND];
        DATA.tableDATA=new_celldata;
        set(gui.tbl,'Data',DATA.tableDATA);
    end 

%-------------------------------------------------------------------------%
    function on_removerow( source, ~ ) 
        celldata = gui.tbl.Data;
        if size(celldata,1) == 1
            new_celldata = celldata;
            new_celldata(1,2:end) = 0;
        else
            new_celldata = celldata(1:end-1,:);
        end
        set(gui.tbl,'Data',new_celldata);
    end

%-------------------------------------------------------------------------%
    function on_celledit( source, callbackdata ) 
        i1 = callbackdata.Indices(1); %row
        i2 = callbackdata.Indices(2); %column
        cyST = source.Data(:,1);
        if i2 == 1 %column 1 has been edited (cycle)
            cyST = source.Data(:,1);
            %This if block forces a final node (ends the drift at last
            %available cycle).
            A = DATA.refdata(:,4);
            lastN = ~isnan(A); %find last index to non-nan gain
            Ind = arrayfun(@(x) find(lastN(:, x), 1, 'last'), 1:size(A, 2));
            Val = DATA.refdata(Ind,2);
            if cyST(end) < nanmax(DATA.refdata(:,2))  && cyST(end) ~= Val
                    celldata = gui.tbl.Data;
                    newrow = [Val 0 0];
                    new_celldata = [celldata;newrow];
                    cyST = [cyST;Val];
                    DATA.tableDATA=new_celldata;
                    set(gui.tbl,'Data',DATA.tableDATA);
            elseif cyST(end) > Val
                cyST(end) = Val;     
            end
            if i1 >= 2 && i1> length(cyST)-1 %second row or beyond is edited
                cyST(i1+1:end) = cyST(i1+1:end)+1;
            end
            inputs.cyEND = [cyST(2:end,1);nanmax(DATA.refdata(:,2))];
% %             refN = get(gui.rb3(1),'Value');
% % %             refE = get(gui.rb3(2),'Value');
% %             refW = get(gui.rb3(2),'Value');
% %             if refN == 1
% %                 DATA.rawGAINS = DATA.rawGAINS_NCEP;
% % %             elseif refE ==1
% % %                 DATA.rawGAINS = DATA.rawGAINS_ERA;
% %             else
% %                 DATA.rawGAINS = DATA.rawGAINS_WOA;
% %             end
            [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,0);
            set(gui.tbl,'Data',DATA.tableDATA)
        else % else slope and intercept has been manually edited.  DO NOT AUTO-CALCULATE.
            cyEND= nanmax(DATA.refdata(:,2));
            if length(cyST)>1
                cyEND = [cyST(2:end);cyEND];
            end
            DATA.tableDATA = get(gui.tbl,'Data');
%             [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
        end
        [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
        DATA.BIC = calcBIC_sO2Argo(gui,DATA);
        % Prep for the "reprocess" step, if user goes that route.  Need
        % "QCA"
        DATA.QCA.O2 = []; %clear it out, then replace with new
        DATA.QCA.O2(:,1) = DATA.tableDATA(:,1); %cycle
        DATA.QCA.O2(:,2) = DATA.tableDATA(:,2); %O2 gain
        DATA.QCA.O2(:,3) = zeros(size(DATA.tableDATA,1),1); %insert column of zeros for "offset" (only for compatability)
        DATA.QCA.O2(:,4) = DATA.tableDATA(:,3); %drift
        updateInterface()
        if   inputs.isprof == 1 %profile selected?
             redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs);
        elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
            DATA = redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs);
        else
            DATA = redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs);  
        end
    end

%-------------------------------------------------------------------------%
    function on_applyO2gain( source, ~ ) 
         handlesODVQC.info.Mprof = 1;
         handlesODVQC.info.float_name = ['ODV',inputs.floatID];
         % Calling Make_Mprof_ODVQC which also gets used in Sage. 
         %Maintain QCA structure, as called in this function
         handlesODVQC.QCA.NO3 = [];
         handlesODVQC.QCA.PH_OFFSET = [];
         handlesODVQC.QCA.PH = [];
         
         handlesODVQC.QCA.O2 = [DATA.tableDATA(:,1:2) zeros(size(DATA.tableDATA,1),1) DATA.tableDATA(:,3)]; %add placeholder for offset (not used in O2 QC, only gain and drift for now)
         gui.wrk_color = gui.doQCbutton.BackgroundColor;
         if sum(isnan(handlesODVQC.QCA.O2)) == 0
             set( gui.doQCbutton,'String','Writing QC to file...','fontsize',14);
             set(gui.doQCbutton,'BackgroundColor','y');
             drawnow
             tf = Make_Mprof_ODVQC(handlesODVQC);
             set( gui.doQCbutton,'String','Write ODV*QC.TXT','fontsize',16)
             set(gui.doQCbutton,'BackgroundColor',gui.wrk_color);
         else
             msgbox('ERROR: AVERAGE GAIN ON SCREEN IS NAN.  ODV*QC.TXT WAS NOT WRITTEN.');
         end
    end

%-------------------------------------------------------------------------%
    function on_reloadQC( source, ~ ) 
        DATA.QCA = get_QCA(inputs.qc_path,DATA.floatNAME);
        cyST = DATA.QCA.O2(:,1);
        inputs.cyEND= nanmax(DATA.refdata(:,2));
        if length(cyST)>1
            inputs.cyEND = [cyST(2:end);inputs.cyEND];
        end
        if ~isempty(DATA.QCA.O2)
            DATA.tableDATA = [DATA.QCA.O2(:,1) DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
            [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
            set(gui.tbl,'Data',DATA.tableDATA);
            [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
            DATA.BIC = calcBIC_sO2Argo(gui,DATA);
           updateInterface()
           if   inputs.isprof == 1 %profile selected?
                redraw_PROF_sageO2(dirs,gui,DATA,inputs)
           elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
               redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
           else
               redraw_WOA_sageO2(dirs,gui,DATA,inputs)   
           end     
        else
            msgbox('ERROR: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED');
        end
    end
%-------------------------------------------------------------------------%
% THIS PUSHBUTTON/FUNCTION WAS REMOVED FROM GUI INTERFACE.  LEFT IN
% COMMENTS IN CASE USER WANTS TO DEVELOP GAIN DRIFT FUNCTIONALITY
% (THUS RELOADING ORIGINAL STORED GAINS WOULD THEN BE USEFUL.)
%     function on_reloadQC( source, ~ ) 
%         DATA.QCA = get_QCA(inputs.qc_path,DATA.floatNAME);
%         cyST = 1;
%         inputs.cyEND= max(DATA.refdata(:,2));
%         if ~isempty(DATA.QCA.O2)
%             DATA.tableDATA = [cyST DATA.QCA.O2(2) 0];
%             [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2Argo(DATA,cyST,inputs.cyEND,1);
%             set(gui.tbl,'Data',DATA.tableDATA);
%             [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA);
%             DATA.AIC = calcAIC_sO2Argo(gui,DATA);
%            updateInterface()
%            if   inputs.isprof == 1 %profile selected?
%                 redraw_PROF_sageO2Argo(dirs,gui,DATA,inputs)
%            elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
%                redraw_CLIMO_sageO2Argo(dirs,gui,DATA,inputs)
%            else
%                redraw_WOA_sageO2Argo(dirs,gui,DATA,inputs)   
%            end     
%         else
%             msgbox('ERROR: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED');
%         end
%     end

%-------------------------------------------------------------------------%
end % end sageO2.m
