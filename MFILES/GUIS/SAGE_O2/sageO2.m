function varargout = sageO2()
% *************************************************************************
% *************************************************************************
% ** THIS MFILE AND ASSOCIATED SOFTWARE IS COPYRIGHT UNDER MIT LICENSE.  **
% **       PLEASE SEE SAGEO2_MITLisence.txt FOR MORE INFORMATION         **
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
% IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% *************************************************************************
% *************************************************************************
%
% ************************************************************************
% sageO2.m
% SOCCOM Assessment and Graphical Evaluation for Oxygen GUI.
% ************************************************************************
%
%
% USE AS:  sageO2
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
% DATE: 12/13/16
% UPDATES: 8/17/17 added global map enhancement.
%          11/14/17 MIT license added.
%		   Dec2018 Added modifications to drift in gain.
%          04/15/19 changed Make_Mprof_ODVQC to Make_Sprof_ODVQC.
%          07/22/19 updated to WOA2018 and GLODAP2019
%		   11/10/21 small bug fix to plottype_onClicked
%          02/15/22 TM, integration of ERA5 reference data.
%          10/31/22 TM, modification to ERA5 data ingestion; now working
%          off of a pre-processed matfile repository on \\atlas\chem\ARGO_PROCESSING\DATA\ERA5\FLOAT_REF\ (much faster to load!)
%           4/5/23 TM, bug fix for matching up cycle numbers on ERA vs full
%           track.
%		   02/08/2024 TM, changes to fx calls in support of Navis Nautilus with SBS83.
%          04/29/24 TM, modification to axes limits when missing track
%          info for first cycle(s) (ie error in WHOI 1557!)
%          09/10/2024 LG, added new map code that uses geoaxes.
%          09/16/2024 LG, resolved a bug where loading a SOLO float without ERA data would permanently grey out the "profile" view radiobutton.
%                     This was resolved by deleting the line that would grey out the button, so now if there is no ERA data and the button is
%                     pressed, it will bring up a message saying there is no ERA data, but the profile view can still be looked at.
% NOTES:
% ************************************************************************
%
% ************************************************************************


% Declare shared variables
gui = createInterface();
dirs = SetupDirs_sO2;
DATA = [];
inputs = [];
handles = [];


% Set some paths
tf = set_sage_paths(dirs.user_dir);

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
        rb3(3) = uicontrol('Parent',bbox,'Style','radiobutton',...
            'String','ERA5','tag','ERA5','Value',0,'Callback',@ref_onClicked );
        rb3(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
            'String','WOA2023','tag','woa','Value',0,'Callback',@ref_onClicked );
        gui.rb3 = rb3;
        
        %             % Control Box 5 (O2 Corrections) - This is really control box 4
        %             % as you go down the layout, but added last (10.11.17), so
        %             % calling it 5
        %            VB5 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);
        %            VP5 = uix.Panel('Parent',VB5,'Padding',5,'title','Bittig Corrections');
        %                 bbox = uix.VButtonBox('Parent',VP5,'Padding',5,...
        %                     'ButtonSize',[100,20],'HorizontalAlignment','left');
        %                     rb5(1) = uicontrol('Parent',bbox,'Style','radiobutton',...
        %                         'String','Night Cals Only','tag','PM','Value',0,'Callback',@O2corr_onClicked );
        %                     rb5(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
        %                         'String','Apply Carryover','tag','carryover','Value',0,'Enable','off','Callback',@O2corr_onClicked );
        %                     rb5(3) = uicontrol('Parent',bbox,'Style','radiobutton',...
        %                         'String','Apply T Correction','tag','tcorr','Value',0,'Enable','off','Callback',@O2corr_onClicked );
        %                 gui.rb5 = rb5;
        
        % Control Box 4 (QC Adjustments)
        VP4 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);
        P4 = uix.Panel('Parent',VP4,'title','Oxygen Gain adjustments:');
        P4v = uix.VBox('Parent',P4,'Padding',3,'Spacing',3);
        P4mean = uicontrol('Parent',P4v,'Style','pushbutton','string',...
            'CALC MEAN GAIN','Enable','on','Callback',@on_calcmeanG);
        P4diag = uicontrol('Parent',P4v,'Style','pushbutton','string',...
            'AirCal Diag','Enable','on','Callback',@on_AirCal_diag);
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
        P4b = uicontrol('parent',P4v,'style','pushbutton','string',...
            'Reload FloatQCList.txt','Callback',@on_reloadQC,'Enable','on');
        P4c = uicontrol('parent',P4v,'style','pushbutton','backgroundcolor','red','string',...
            'REPROCESS','fontsize',18,'Callback',@on_reprocess,'Enable','on');
        
        
        
        
        %         set( controlLayout, 'Heights', [-6 -3 -3 -6 -1] ); %Main control vbox heights
        set( controlLayout, 'Heights', [-5.5 -2.5 -2.5 -5] ); %Main control vbox heights
        set( vb, 'Heights', [-1 -1 -1.75 -1.75 -1.75] ); %Float spec box heights
        set( P4v, 'Heights', [-1 -1 -1 -3 -1 -1] ); %QC adj box heights
        
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
        
        sdnprofmin = DATA.track(DATA.track(:,2)==PROFmin,1);
        cycprofmin = DATA.track(DATA.track(:,2)==PROFmin,2);
        if ~isempty(sdnprofmin) %if missing track info at start, use minimum of available float track
            DATA.xlims{1} = [sdnprofmin DATA.track(DATA.track(:,2)==PROFmax,1)];
            DATA.xlims{2} = [cycprofmin DATA.track(DATA.track(:,2)==PROFmax,2)];
        else
            DATA.xlims{1} = [nanmin(DATA.track(:,1)) DATA.track(DATA.track(:,2)==PROFmax,1)];
            DATA.xlims{2} = [nanmin(DATA.track(:,2)) DATA.track(DATA.track(:,2)==PROFmax,2)];
        end

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
                s2 = tmp(s);
                DATA.airsub{1}{s2} = DATA.RAWorQCair{1}{s2}(DATA.RAWorQCair{1}{s2}(:,2) >=PROFmin & ...
                    DATA.RAWorQCair{1}{s2}(:,2) <= PROFmax,:);
                DATA.airsub{2}{s2} = DATA.RAWorQCair{2}{s2}(DATA.RAWorQCair{2}{s2}(:,2) >=PROFmin & ...
                    DATA.RAWorQCair{2}{s2}(:,2) <= PROFmax,:);
                %                 DATA.GAINS{s} = DATA.refsub(:,4)./DATA.airsub{1}{s}(:,9);
                [cc,iia,iib] = intersect(DATA.refsub(:,2),DATA.airsub{1}{s2}(:,2));
                DATA.GAINS{s2} = DATA.refsub(iia,4)./DATA.airsub{1}{s2}(iib,9);
                DATA.GAINStime{s2} = DATA.refsub(iia,1); %for use in plotting climos
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
        open([dirs.mfiles,'GUIS\SAGE_O2\README_sO2.txt'])
    end % onHelp

%-------------------------------------------------------------------------%
    function onack( ~, ~ )
        % User has asked for the documentation
        open([dirs.mfiles,'GUIS\SAGE_O2\acknowledgements_sO2.txt'])
    end % onHelp

%-------------------------------------------------------------------------%
    function on_selectfloat( ~, ~ )
        % first load MBARI float list for cross-referencing
        load([dirs.cal,'MBARI_Float_list.mat']);
        iWMO = find(strcmp('WMO',d.hdr) == 1);
        iMB  = find(strcmp('MBARI ID',d.hdr) == 1);
        iINST = find(strcmp('INST ID',d.hdr) == 1);
        iFLT  = find(strcmp('float type',d.hdr) == 1);
        iPRG  = find(strcmp('Program',d.hdr) == 1);
        iREG  = find(strcmp('Region',d.hdr) == 1);
        % select float by wmo from dialog box
        inputs.rorq = 1;
        %         x=1;
        %         while x==1
        [pn] = uigetdir([dirs.mat],'SELECT FLOAT'); %Select By WMO in the FLOATS dir.
        wmonum = regexp(pn,'\d{7}','match','once');
        if isempty(wmonum)
            msgbox('ERROR: Could not extract any float WMO from path selected.')
            return
        end
        tmi = regexpi(d.list(:,iWMO),wmonum);
        tmi2 = cellfun(@isempty,tmi);
        inputs.MBARI_ID = d.list{~tmi2,iMB};
        inputs.INST_ID = d.list{~tmi2,iINST};
        inputs.WMO_ID = wmonum;
        handles.info.WMO_ID = inputs.WMO_ID; %needed for Process_GUI_float...
        inputs.floatTYPE = d.list{~tmi2,iFLT};
        handles.info.float_type = inputs.floatTYPE; %needed for Process_GUI_float...
        inputs.PROGRAM = d.list{~tmi2,iPRG};
        inputs.REGION = d.list{~tmi2,iREG};
        %         end
        %set( gui.Fbutton,'String',inputs.floatID);
        
        set( gui.Fbutton,'String','Loading msg files ...');
        wrk_color = gui.Fbutton.BackgroundColor;
        set(gui.Fbutton,'BackgroundColor','y');
        gui.rb2(1).Value = 1; % reset to surface plot when float opens
        gui.rb2(2).Value = 0;
        gui.rb2(3).Value = 0;
        
        drawnow
        
        if strcmp(handles.info.float_type,'SOLO')
            DATA = getall_SOLOdata_sO2(dirs,inputs.MBARI_ID);
        else
            DATA = getall_floatdata_sO2(dirs,inputs.MBARI_ID,handles.info.float_type );
        end
%         keyboard
        DATA.MBARI_ID = inputs.MBARI_ID; %why not
        DATA.WMO_ID = inputs.WMO_ID;
        DATA.INST_ID = inputs.INST_ID;
        
        set(gui.Mbutton,'Enable','on');
        %inputs.depthedit = [1480 1520]; %default pressure range (deep)
        inputs.depthedit = [0 30]; %default pressure range (deep)
        inputs.GLDPkm = 30;
        inputs.profedit = [1 DATA.track(end,2)];
        inputs.AMPM = 'all';
        gui.TP.SelectionChangedFcn = @on_PlottypeChanged; % I don't think this is used -jp?
        gui.t.SelectedChild = 1; %default to Raw upon float selection
        gui.t.SelectionChangedFcn = @on_RAWorQC;
        %         gui.TP.Selection = 2;   %default to Deep tab
        %set(gui.rb2(2),'Value',1);
        %         set(gui.rb5,'Value',0);
        DATA.RAWorQCprof = DATA.O2data{1}; %profile data: sdn, cast, S, P, T, Phase,  O2, O2sol, pO2, pH20, O2Sat
        DATA.RAWorQCair = DATA.O2air;
        gui.whichAX = gui.t1ax;
        if iscell(DATA.O2air{1})==1
            mytemp = find(~cellfun(@isempty,DATA.O2air{1}));
            DATA.howmanyairs = length(mytemp);
        else
            DATA.howmanyairs = 0;
        end
        %get QC adjustment data
        inputs.qc_path = [dirs.QCadj,inputs.WMO_ID,'_FloatQCList.txt'];
        DATA.QCA = get_QCA(inputs.qc_path);
        
        %Calculate WOA data and surface o2 sat for all floats
        set( gui.Fbutton,'String','Loading WOA ...');
        set(gui.Fbutton,'BackgroundColor','y');
        drawnow
        [~,inputs.intersect_cycles_WOA,~] = intersect(DATA.track(:,2),DATA.O2data{1}(:,2));
        Wtrack = [DATA.track(inputs.intersect_cycles_WOA,1) DATA.track(inputs.intersect_cycles_WOA,4) DATA.track(inputs.intersect_cycles_WOA,3)];
        try
            Wdata = get_WOA_local(dirs.woa,Wtrack, [0 2000], 'O2sat');
            zsurf = Wdata.Z<=25;
            WOA_surf = Wdata.d(zsurf,:);
            DATA.WOAsurf = nanmean(WOA_surf,1);
            DATA.reftag = 'WOA';
        catch
            msgbox({'ERROR: getWOA failed.',...
                'Must use NCEP or bottle data for reference.'})
            set(gui.rb3(2),'Value',0,'Enable','off');
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
        if strcmp(handles.info.float_type,'SOLO') %SOLO floats don't have air sample flavor associated with telemetry (only the 'sequence' of measurements)
            DATA.inairIND = 3;
        else
            DATA.inairIND = 1;
        end
        if ~isempty(DATA.O2air{1})  %AIRCAL EXISTS SO USE NCEP OR ERA
            [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2air{1}{DATA.inairIND}(:,2));
            %______________________________________________________________
            %TRY NCEP FIRST (REALTIME).  DATA IS GRABBED FROM THE WEB.  IF WEBSITE IS
            %DOWN FOR MAINTENANCE ETC WILL TRY ERA NEXT.
            try
                set( gui.Fbutton,'String','Loading NCEP ...');
                set(gui.Fbutton,'BackgroundColor','y');
                drawnow
%                 save('temp.mat','DATA','inputs','dirs')
%                 pause
                DATA.NCEP = getNCEP(DATA.track(inputs.intersect_cycles,1),DATA.track(inputs.intersect_cycles,3),DATA.track(inputs.intersect_cycles,4),dirs);
                %                 DATA.ERA = getERA_sO2(dirs.ERA,DATA.track(inputs.intersect_cycles,1),DATA.track(inputs.intersect_cycles,3),DATA.track(inputs.intersect_cycles,4));
            catch
                msgbox({'ERROR: getNCEP failed.',...
                    'Must use ERA or WOA for reference.'})
                set(gui.rb3(1),'Value',0,'Enable','off');
                DATA.NCEP.PRES = [];
            end
            try
                set( gui.Fbutton,'String','Loading ERA5 ...');
                set(gui.Fbutton,'BackgroundColor','y');
                drawnow
                %DATA.ERA = getERA_sO2(dirs.ERA,DATA.track(inputs.intersect_cycles,1),DATA.track(inputs.intersect_cycles,3),DATA.track(inputs.intersect_cycles,4));
				myERA = load([dirs.ERA,inputs.WMO_ID,'_ERA5ref.mat']);
                % QUICK & DIRTY FIX if size of intersect_cycles > size PRES
                % need to check that intersect cycle index values are <
                % size of ERA arrays; TM 4.5.23 - fix to the fix (logic not
                % quite comprehensive enough)
                if max(inputs.intersect_cycles,[],'omitnan') > size(myERA.FLT.CYC,1)
                    fprintf(['WARNING: ERA cycle count < track cycle count.', ...
                        'Padding ERA.PRES with NaNs!\n'])
                    [C,IA,IB] = intersect(inputs.intersect_cycles,myERA.FLT.CYC);
%                     tg = inputs.intersect_cycles <= size(myERA.FLT.CYC,1);
%                     ERA_intersect = inputs.intersect_cycles(tg);
                    %ERA_intersect = inputs.intersect_cycles(1:size(myERA.ERA.PRES,1));
                    DATA.ERA.PRES = ones(size(inputs.intersect_cycles))*NaN;
                    DATA.ERA.PRES(IA) = myERA.ERA.PRES(IB);
                    clear C IA IB
                else
                    DATA.ERA.PRES = myERA.ERA.PRES(inputs.intersect_cycles);
                end
%                 DATA.ERA.PRES = myERA.ERA.PRES(inputs.intersect_cycles);
            catch
                msgbox({'ERROR: getERA failed.',...
                    'Must use NCEP or WOA for reference.'})
                DATA.ERA.PRES = [];
            end
            if ~isempty(DATA.ERA.PRES)
                set(gui.rb3(3),'Value',1,'Enable','on');
                %                 set(gui.rb3(2),'Value',0,'Enable','on');
                set(gui.rb3(2),'Value',0);
                set(gui.rb3(1),'Value',0);
                DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
                DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
                DATA.refdata(:,3) = DATA.ERA.PRES;
                DATA.reftag = 'ERA';
                %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
%                 save('tmp2.mat','DATA')
%                 pause
                if ~isempty(DATA.O2air{1}{3})
                    DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - mean([DATA.RAWorQCair{1}{DATA.inairIND}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2,'omitnan')).*0.20946;
                else
                    DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
                end
                %                 DATA.refdata(1:5,3)
                %                 DATA.RAWorQCair{1}{1}(1:5,10)
                inputs.cyEND = nanmax(DATA.refdata(:,2));
                %Maintain raw gains for adjustment calcs in
                %calc_gain_ADJtable_sO2.m and for toggling between reference
                %datasets
                if DATA.howmanyairs > 1
                    DATA.rawGAINS_ERA{1} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{DATA.inairIND}(:,9);
                    DATA.rawGAINS_ERA{3} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{3}(:,9);
                    useGAIN = DATA.rawGAINS_ERA{3};
                else
                    DATA.rawGAINS_ERA{1} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{1}(:,9);
                    DATA.rawGAINS_ERA{3} = [];
                    useGAIN = DATA.rawGAINS_ERA{1};
                end
                DATA.GAINS=DATA.rawGAINS_ERA;
                DATA.rawGAINS=DATA.rawGAINS_ERA;
                inputs.isprof = 0;
                %populate gain adjustment table
                cyST = 1;
                DATA.iswoaref = get(gui.rb3(2),'Value');
                if isempty(DATA.QCA.O2)
                    DATA.tableDATA = [cyST nanmean(useGAIN) 0];
                    [DATA.tableDATA,DATA.brkRESIDS] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
                    set(gui.tbl,'Data',DATA.tableDATA)
                    msgbox({'FYI: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED','POPULATING ADJUSTMENT TABLE WITH DEFAULT.'});
                else
                    %                      DATA.tableDATA = [cyST DATA.QCA.O2(2) 0];
                    DATA.tableDATA = [DATA.QCA.O2(:,1) DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
                    cyST = DATA.QCA.O2(:,1);
                    inputs.cyEND= nanmax(DATA.refdata(:,2));
                    if length(cyST)>1
                        inputs.cyEND = [cyST(2:end);inputs.cyEND];
                    end
                    [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
                    set(gui.tbl,'Data',DATA.tableDATA)
                end
                DATA.BIC = calcBIC_sO2(gui,DATA);
                updateInterface()
                redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
            elseif ~isempty(DATA.NCEP.PRES)
                set(gui.rb3(1),'Value',1,'Enable','on');
                set(gui.rb3(3),'Value',0,'Enable','on');
                set(gui.rb3(2),'Value',0,'Enable','on');
                DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
                DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
                DATA.refdata(:,3) = DATA.NCEP.PRES;
                DATA.reftag = 'NCEP';
                %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
                if ~isempty(DATA.O2air{1}{3})
                    DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - nanmean([DATA.RAWorQCair{1}{DATA.inairIND}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2)).*0.20946;
                else
                    DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
                end
                %                 DATA.refdata(1:5,3)
                %                 DATA.RAWorQCair{1}{1}(1:5,10)
                inputs.cyEND = nanmax(DATA.refdata(:,2));
                %Maintain raw gains for adjustment calcs in
                %calc_gain_ADJtable_sO2.m and for toggling between reference
                %datasets
                if DATA.howmanyairs > 1
                    DATA.rawGAINS_NCEP{1} = DATA.refdata(:,4)./DATA.RAWorQCair{1}{DATA.inairIND}(:,9);
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
                    [DATA.tableDATA,DATA.brkRESIDS] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
                    set(gui.tbl,'Data',DATA.tableDATA)
                    msgbox({'FYI: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED','POPULATING ADJUSTMENT TABLE WITH DEFAULT.'});
                else
                    %                      DATA.tableDATA = [cyST DATA.QCA.O2(2) 0];
                    DATA.tableDATA = [DATA.QCA.O2(:,1) DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
                    cyST = DATA.QCA.O2(:,1);
                    inputs.cyEND= nanmax(DATA.refdata(:,2));
                    if length(cyST)>1
                        inputs.cyEND = [cyST(2:end);inputs.cyEND];
                    end
                    [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
                    set(gui.tbl,'Data',DATA.tableDATA)
                end
                DATA.BIC = calcBIC_sO2(gui,DATA);
                updateInterface()
                redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
            end
            
        end %end if aircal exists
        %__________________________________________________________________
        % NOW SET DEFAULT TO WOA IF NO AIRCAL OR IF BOTH NCEP AND ERA
        % FAILED (WEB-GET NOT WORKING, AND NEWER FLOAT)
        if isempty(DATA.O2air{1}) || (isempty(DATA.NCEP.PRES) && isempty(DATA.ERA.PRES))
            set(gui.rb3(1),'Value',0,'Enable','off');
            set(gui.rb3(3),'Value',0,'Enable','off');
            set(gui.rb3(2),'Value',1);
            if isfield(DATA,'refdata')
                DATA=rmfield(DATA,'refdata');
            end
            DATA.iswoaref = get(gui.rb3(2),'Value');
            DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles_WOA,1);
            DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles_WOA,2);
            DATA.refdata(:,3) = nan(length(DATA.refdata(:,2)),1);
            DATA.refdata(:,4) = DATA.WOAsurf;
            inputs.cyEND = nanmax(DATA.refdata(:,2));
            DATA.rawGAINS=DATA.rawGAINS_WOA;
            DATA.GAINS=DATA.rawGAINS_WOA; %maintain raw gains for adjustment calcs in calc_gain_ADJtable_sO2.m
            inputs.isprof = 0;
            %populate gain adjustment table
            
            %cyST = DATA.QCA.O2(:,1); %This will only work if adjustment table already exists -jp 12/18/18
            inputs.cyEND= nanmax(DATA.refdata(:,2));
            
            if isempty(DATA.QCA.O2) % NO QC yet
                cyST = 1;
                DATA.tableDATA = [cyST nanmean(DATA.GAINS{1}) 0];
                [DATA.tableDATA,DATA.brkRESIDS] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
                set(gui.tbl,'Data',DATA.tableDATA)
                msgbox({'FYI: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED','POPULATING ADJUSTMENT TABLE WITH DEFAULT.'});
            else
                cyST = DATA.QCA.O2(:,1);
                if length(cyST)>1 % this block seems like it could put in end value twice?
                    inputs.cyEND = [cyST(2:end);inputs.cyEND];
                end
                DATA.tableDATA = [cyST DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
                disp(DATA.tableDATA)
                [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
                set(gui.tbl,'Data',DATA.tableDATA)
            end
            DATA.BIC = calcBIC_sO2(gui,DATA);
            updateInterface()
            redraw_WOA_sageO2(dirs,gui,DATA,inputs)
        end % end if isempty aircal
        
        %__________________________________________________________________
        % PRELOAD BOTTLE DATA FOR PROFILE PLOTS
        % Get bottle data
        Blookup = [dirs.bottle,'BottleData_lookup_table.txt'];
        % LOAD BOTTLE DATA LOOK UP TABLE & STORE IN handles STRUCTURE
        % LOOK UP TABLE  HEADER = [UW_ID  WMO   CRUISE   STN   CAST   Data file]
        fid = fopen(Blookup);
        d   = textscan(fid, '%s %s %s %s %f %f %s','Delimiter', '\t','HeaderLines',1);
        fclose(fid);
        bottle_lookup = d;
        IndexB = strfind(bottle_lookup{1,3},inputs.WMO_ID);
        ind = find(not(cellfun('isempty', IndexB)));
        if ~isempty(ind) % float exists in lookup table
            bottle_fname = bottle_lookup{1,7}{ind};
            stn = bottle_lookup{1,5}(ind);
            cst = bottle_lookup{1,6}(ind);
            % data file & data exist for float
            if ~isempty(bottle_fname) && (stn ~= -999) && (cst ~= -999)
                tic
                set( gui.Fbutton,'String','Loading bottle data ...');
                set(gui.Fbutton,'BackgroundColor','y');
                drawnow
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
        
        %set( gui.Fbutton,'String',inputs.floatID)
        set( gui.Fbutton,'String',inputs.WMO_ID)
        set(gui.Fbutton,'BackgroundColor',wrk_color);
        drawnow
        %     saveas(gui.Window,'sageO2.png')
        %     set(gui.Window,'PaperPositionMode','auto')
        %     print(gui.Window,'sageO2interface.png','-dpng','-r300')
    end %end selectfloat
    
%-------------------------------------------------------------------------%
    function on_showmap( ~,~ )


        track = DATA.track; % float track
        tnan = track(:,3) == 1e-10;
        track(tnan,3:4) = NaN;
        lon_limits = [min(track(:,3))-10 max(track(:,3))+10];
        lat_limits = [min(track(:,4))-6 max(track(:,4))+6];
        
        %================ NEW MAP CODE =====================
        figure;
        gax = geoaxes; % initiate the geographic axes
        geobasemap(gax, 'landcover'); % several options to choose from
        geoplot(gax, track(:,4), track(:,3),'k.-','MarkerSize',10); % lat, lon
        hold(gax,'on')
        Hm1 = geoplot(gax, track(1,4), track(1,3),'ko','MarkerSize',7,...
            'MarkerFaceColor', 'g');
        Hm2 = geoplot(gax, track(end,4), track(end,3),'ko','MarkerSize',7,...
            'MarkerFaceColor', 'r');
        hold(gax,'off')
        gax.Title.String = ['FLOAT ',inputs.WMO_ID,' (',inputs.MBARI_ID,')'];%app.finfo.wmo;
        gax.FontSize = 14;
        track_legend = legend([Hm1 Hm2],'first', 'last');
        track_legend.Location = 'northeastoutside';
        %================ END NEW MAP CODE =====================

        % track = DATA.track; % float track
        % tnan = track(:,3) == 1e-10;
        % track(tnan,3:4) = NaN;
        % tmplong = track(:,3);
        % tmplong(tmplong>180) = tmplong(tmplong>180)-360;
        % lon_limits = [min(tmplong)-20 max(tmplong)+20];
        % 
        % 
        % lat_limits = [min(track(:,4))-6 max(track(:,4))+6];
        % handles.float_track = figure;
        % set(gcf,'Position', [350 50 860 620]);
        % DATA.CMAP = load('SOCCOM_NCP_cmap.mat');
        % DATA.CMAP = DATA.CMAP.cmap;
        % colormap(handles.float_track,DATA.CMAP);
        % 
        % if max(track(:,4)) > -30 % Don't use stereographic projection
        %     m_proj('lambert','lon',lon_limits,'lat',lat_limits);
        % else
        %     m_proj('stereographic','latitude',-90,'radius',60,'rotangle',45);
        % end
        % m_tbase('contourf','edgecolor','none');
        % hold on
        % tmplong = track(:,3);
        % tmplong(tmplong>180) = tmplong(tmplong>180)-360;
        % m_plot(tmplong,track(:,4),'ko-', 'MarkerSize',3)%     xlim(lon_limits);
        % title(['FLOAT ',inputs.WMO_ID,' (',inputs.MBARI_ID,')'],'FontSize', 16)
        % hold on
        % Hm1 = m_plot(track(1,3),track(1,4),'o', 'MarkerSize',10, 'MarkerFaceColor','g', ...
        %     'MarkerEdgeColor','k');
        % Hm2 = m_plot(track(end,3),track(end,4),'o', 'MarkerSize',10, 'MarkerFaceColor','r', ...
        %     'MarkerEdgeColor','k');
        % 
        % colorbar
        % 
        % set(gca,'ydir','normal');
        % m_grid('linewidth',2,'tickdir','out','xaxislocation','top');
        % hold off
        % track_legend = legend([Hm1 Hm2],'first', 'last');
        % track_legend.Location = 'northeastoutside';
        %================ END ORIGINAL SAGEO2 CODE =====================
    end

%-------------------------------------------------------------------------%
    function on_RAWorQC( source,~ )
        inputs.rorq = get(source,'selectedchild'); %'1' is raw, '2' is QC
        [gui,DATA] = apply_QCadj_sO2(gui,inputs,DATA);
        updateInterface()
        if   inputs.isprof == 1; %profile selected?
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
        elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1;
            redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
        else
            redraw_WOA_sageO2(dirs,gui,DATA,inputs)
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
        if   inputs.isprof == 1 %profile selected?
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
        elseif ~isempty(DATA.O2air{1})&& get(gui.rb3(2),'Value')~=1
            redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
        else
            redraw_WOA_sageO2(dirs,gui,DATA,inputs)
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
        if   inputs.isprof == 1 %profile selected?
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
        elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
            redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
        else
            redraw_WOA_sageO2(dirs,gui,DATA,inputs)
        end
    end

%-------------------------------------------------------------------------%
    function on_GLODAP( source, ~ )
        GDkm = get(source,'String');
        inputs.GLDPkm = str2double(GDkm);
        updateInterface()
%         if gui.TP.Selection == 3
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
%         end
    end

%-------------------------------------------------------------------------%
    function plottype_onClicked( source, ~ )
        source.Value = 1; % select this
        RB2 = gui.rb2;
        set( RB2(RB2~=source), 'Value', 0 ) % unselect others (NOT SURE if THIS works? -jp)
        
        % get track data
        reftag = get(source,'tag');
        SC = str2num(reftag);
        %         SC = get(source,'Selection'); %tabs: (1)surf, (2)depth, (3)profile
        Dedit = [0 30; 1480 1520; 0 1600];
        inputs.depthedit = Dedit(SC,:);
        if get(gui.rb2(3),'Value')==1
            inputs.isprof = 1; %profile selected?
            updateInterface()
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
        else
            inputs.isprof = 0;
            updateInterface()
            if ~isempty(DATA.O2air{1}) && get(gui.rb2(3),'Value')~=1 && ...
                    get(gui.rb3(2),'Value')~=1
                redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
            else
                redraw_WOA_sageO2(dirs,gui,DATA,inputs)
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
        [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2air{1}{DATA.inairIND}(:,2));
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
            if ~isempty(DATA.O2air{1}{3})
                DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - nanmean([DATA.RAWorQCair{1}{DATA.inairIND}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2)).*0.20946;
            else
                DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
            end
            if inputs.rorq == 2 %QC tab
                [gui,DATA] = apply_QCadj_sO2(gui,inputs,DATA);
            end
            updateInterface()
            if   inputs.isprof == 1 %profile selected?
                redraw_PROF_sageO2(dirs,gui,DATA,inputs)
            else
                redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
            end
            % %         % ERA______________________________________________________________
                    elseif (strcmp(reftag,'ERA5')) == 1
                        [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2air{1}{DATA.inairIND}(:,2));
                        if isempty(DATA.ERA.PRES)
                            set(gui.rb3(1),'Value',1)
                            set(gui.rb3(3),'Value',0,'Enable','off')
                            msgbox({'WARNING: No ERA data exists for time period of interest.',...
                                'Must use NCEP or WOA as refernce.'});
                        else
                            if isfield(DATA,'refdata')
                                DATA=rmfield(DATA,'refdata');
                            end
                            DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
                            DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
                            DATA.refdata(:,3) = DATA.ERA.PRES;
                            DATA.reftag = 'ERA';
                            %fourth column of ref data = pO2 = (ncep.pres - pH20)*airmix
                            if inputs.rorq == 2 %QC tab
                                DATA.RAWorQCair = DATA.O2air;
                            end
                            if ~isempty(DATA.O2air{1}{3})
                                DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - mean([DATA.RAWorQCair{1}{DATA.inairIND}(:,10) DATA.RAWorQCair{1}{3}(:,10)],2,'omitnan')).*0.20946;
                            else
                                DATA.refdata(:,4) = (DATA.refdata(:,3)./100 - DATA.RAWorQCair{1}{1}(:,10)).*0.20946;
                            end
                            if inputs.rorq == 2 %QC tab
                                [gui,DATA] = apply_QCadj_sO2(gui,inputs,DATA);
                            end
                                        updateInterface()

                            if   inputs.isprof == 1 %profile selected?
                                redraw_PROF_sageO2(dirs,gui,DATA,inputs)
                            else
                                redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
                            end
                        end
            % WOA______________________________________________________________
        else
            if isfield(DATA,'refdata')
                DATA=rmfield(DATA,'refdata');
            end
            [~,inputs.intersect_cycles,~] = intersect(DATA.track(:,2),DATA.O2data{1}(:,2));
            DATA.refdata(:,1) = DATA.track(inputs.intersect_cycles,1);
            DATA.refdata(:,2) = DATA.track(inputs.intersect_cycles,2);
            Wtrack = [DATA.track(inputs.intersect_cycles,1) DATA.track(inputs.intersect_cycles,4) DATA.track(inputs.intersect_cycles,3)];
            Wdata = get_WOA_local(dirs.woa,Wtrack, [0 2000], 'O2sat');
            zsurf = Wdata.Z<=25;
            WOA_surf = Wdata.d(zsurf,:);
            DATA.WOAsurf = nanmean(WOA_surf,1);
            DATA.refdata(:,3) = nan(length(DATA.refdata(:,2)),1);
            DATA.refdata(:,4) = DATA.WOAsurf;
            DATA.reftag = 'WOA';
            updateInterface()
            if   inputs.isprof == 1 %profile selected?
                redraw_PROF_sageO2(dirs,gui,DATA,inputs)
            else
                redraw_WOA_sageO2(dirs,gui,DATA,inputs)
            end
        end
    end

% %-------------------------------------------------------------------------%
%     function O2corr_onClicked( source, ~ )
%         s_string = source.String;
%         s_value = source.Value;
%         F = strfind(s_string,'PM');
%         if F>0
%             if s_value == 1
%                 inputs.AMPM = 'PM';
%             else
%                 inputs.AMPM = 'all';
%             end
%         end
%         updateInterface()
%         if   inputs.isprof == 1 %profile selected?
%              redraw_PROF_sageO2(dirs,gui,DATA,inputs)
%         elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
%             redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
%         else
%             redraw_WOA_sageO2(dirs,gui,DATA,inputs)
%         end
%     end

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
        %         useG
        DATA.tableDATA = [1 useG 0];
        cyST = DATA.tableDATA(1);
        inputs.cyEND= nanmax(DATA.refdata(:,2));
        [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
        %             DATA.tableDATA = get(gui.tbl,'Data');
        set(gui.tbl,'Data',DATA.tableDATA)
        [gui,DATA] = apply_QCadj_sO2(gui,inputs,DATA);
        DATA.BIC = calcBIC_sO2(gui,DATA);
        updateInterface()
        if   inputs.isprof == 1 %profile selected?
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
        elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
            redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
        else
            redraw_WOA_sageO2(dirs,gui,DATA,inputs)
        end
    end

%-------------------------------------------------------------------------%
    function on_AirCal_diag(source, ~)
        tempthresh = 0.1;
        F1 = myfig(0.1,0.1,0.7,0.7);
        airG = (1./DATA.rawGAINS_NCEP{3}).*100;
        airCYC = DATA.refdata(:,2);
        O2prof = DATA.O2data{1}; %profile data
        O2Air = DATA.O2data{4}; %air sequence in-air
        allcast = unique(O2prof(:,2)); %unique casts
        k=1;
        for ii = 1:length(allcast)
            %              if k == 65
            %                  surfsat(k) = nan;
            %                  k = k+1;
            %                  continue
            %              end
            tmpO = O2prof(O2prof(:,2)==allcast(ii),:);
            surfpresind = find(tmpO(:,4) == min(tmpO(:,4),[],'omitnan')); %find index of shallowest non-nan pressure
            surfsattmp = tmpO(surfpresind,11);
            surftemptmp = tmpO(surfpresind,5);
            surfsat(k) = surfsattmp;
            surftemp(k) = surftemptmp;
            
            tmpAO = O2Air(O2Air(:,2)==allcast(ii),:);
            tmpTEMP = tmpAO(:,5);
            airtemp(k) = median(tmpTEMP,'omitnan');
            xx = find(airCYC==allcast(ii));
            if isempty(xx)
                continue
            end
            AIRG(k) = airG(xx);
            k = k+1;
        end
        subplot(2,3,1)
        plot(surfsat,AIRG,'bo','markeredgecolor','k','markerfacecolor','b')
        %         xlim([89 98])
        %         ylim([89 98])
        X = surfsat;
        Y = AIRG;
        xx = find(isnan(X));
        X(xx) = [];
        Y(xx) = [];
        yy = find(isnan(Y));
        X(yy) = [];
        Y(yy) = [];
        [m,b,r,sm,sb] = lsqfity(X,Y);
        newX = get(gca,'xlim');
        newY = newX.*m+b;
        myY = get(gca,'ylim');
        hold on
        plot(newX,newY,'r-')
        hold on
        plot(newX,newX,'k--')
        title({['All data; m = ',num2str(m,'%3.2f')],['Mean Gain = ',num2str(nanmean((1./Y).*100),'%5.4f')]}) %Y = ',num2str(m,'%3.2f'),'*X + ',num2str(b,'%3.2f')])
        xlabel('%O2 saturation seawater')
        ylabel('%O2 saturation air')
        
        %Now restrict to data where Temp differs by more than X amount.  X=  0.3
        subplot(2,3,2)
        TEMPdiff = abs(surftemp-airtemp);
        xtempdiff = find(TEMPdiff>=tempthresh);
        plot(surfsat(xtempdiff),AIRG(xtempdiff),'bo','markeredgecolor','k','markerfacecolor','b')
        %         xlim([89 98])
        %         ylim([89 98])
        X = surfsat(xtempdiff);
        Y = AIRG(xtempdiff);
        xx = find(isnan(X));
        X(xx) = [];
        Y(xx) = [];
        yy = find(isnan(Y));
        X(yy) = [];
        Y(yy) = [];
        [m,b,r,sm,sb] = lsqfity(X,Y);
        newX = get(gca,'xlim');
        newY = newX.*m+b;
        myY = get(gca,'ylim');
        hold on
        plot(newX,newY,'r-')
        hold on
        plot(newX,newX,'k--')
        title({['|',char(916),'T| >=',num2str(tempthresh,'%2.1f'),'; m = ',num2str(m,'%3.2f')],['Mean Gain = ',num2str(nanmean((1./Y).*100),'%5.4f')]}) %Y = ',num2str(m,'%3.2f'),'*X + ',num2str(b,'%3.2f')])
        xlabel('%O2 saturation seawater')
        ylabel('%O2 saturation air')
        
        %Now restrict to data where %sat difference is less than 2%
        subplot(2,3,3)
        SATdiff = abs(AIRG - surfsat);
        xsatdiff = find(SATdiff<2);
        plot(surfsat(xsatdiff),AIRG(xsatdiff),'bo','markeredgecolor','k','markerfacecolor','b')
        %         xlim([89 98])
        %         ylim([89 98])
        X = surfsat(xsatdiff);
        Y = AIRG(xsatdiff);
        xx = find(isnan(X));
        X(xx) = [];
        Y(xx) = [];
        yy = find(isnan(Y));
        X(yy) = [];
        Y(yy) = [];
        [m,b,r,sm,sb] = lsqfity(X,Y);
        newX = get(gca,'xlim');
        newY = newX.*m+b;
        myY = get(gca,'ylim');
        hold on
        plot(newX,newY,'r-')
        hold on
        plot(newX,newX,'k--')
        title({['|',char(916),'%sat| < 2% ','; m = ',num2str(m,'%3.2f')],['Mean Gain = ',num2str(nanmean((1./Y).*100),'%5.4f')]}) %Y = ',num2str(m,'%3.2f'),'*X + ',num2str(b,'%3.2f')])
        xlabel('%O2 saturation seawater')
        ylabel('%O2 saturation air')
        
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
            [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,0);
            set(gui.tbl,'Data',DATA.tableDATA);
        else % else slope and intercept has been manually edited.  DO NOT AUTO-CALCULATE.
            cyEND= nanmax(DATA.refdata(:,2));
            if length(cyST)>1
                cyEND = [cyST(2:end);cyEND];
            end
            DATA.tableDATA = get(gui.tbl,'Data');
            %             [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
        end
        [gui,DATA] = apply_QCadj_sO2(gui,inputs,DATA);
        DATA.BIC = calcBIC_sO2(gui,DATA);
        % Prep for the "reprocess" step, if user goes that route.  Need
        % "QCA"
        DATA.QCA.O2 = []; %clear it out, then replace with new
        DATA.QCA.O2(:,1) = DATA.tableDATA(:,1); %cycle
        DATA.QCA.O2(:,2) = DATA.tableDATA(:,2); %O2 gain
        DATA.QCA.O2(:,3) = zeros(size(DATA.tableDATA,1),1); %insert column of zeros for "offset" (only for compatability)
        DATA.QCA.O2(:,4) = DATA.tableDATA(:,3); %drift
        updateInterface()
        if   inputs.isprof == 1 %profile selected?
            redraw_PROF_sageO2(dirs,gui,DATA,inputs)
        elseif ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1
            redraw_CLIMO_sageO2(dirs,gui,DATA,inputs)
        else
            redraw_WOA_sageO2(dirs,gui,DATA,inputs)
        end
    end

%-------------------------------------------------------------------------%
    function on_reloadQC( source, ~ )
        DATA.QCA = get_QCA(inputs.qc_path);
        cyST = DATA.QCA.O2(:,1);
        inputs.cyEND= nanmax(DATA.refdata(:,2));
        if length(cyST)>1
            inputs.cyEND = [cyST(2:end);inputs.cyEND];
        end
        if ~isempty(DATA.QCA.O2)
            DATA.tableDATA = [DATA.QCA.O2(:,1) DATA.QCA.O2(:,2) DATA.QCA.O2(:,4)];
            [DATA.tableDATA,DATA.brkRES] = calc_gainADJtable_sO2(DATA,cyST,inputs.cyEND,1);
            set(gui.tbl,'Data',DATA.tableDATA);
            [gui,DATA] = apply_QCadj_sO2(gui,inputs,DATA);
            DATA.BIC = calcBIC_sO2(gui,DATA);
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
    function on_reprocess( source, ~ )
        % %              set(gui.Window,'PaperPositionMode','auto')
        % %      print(gui.Window,'sageO2interfaceQC','-djpeg','-r300')
        % %      return
        update_flag = 1;
        USER       = getenv('USERNAME');
        sdn        = datestr(now,'mm/dd/yy HH:MM:SS');
        title_txt  = ' ';
        user_input = inputdlg('QC adjustment comments',title_txt, ...
            [1, length(title_txt)+150])  ;
        
        %IF NO COMMENT  OR CANCELED DON'T DO ANYTHING
        if isempty(user_input) % Canceled or no info - either way stop process
            msgbox('NO COMMENT ENTERED - PROCESSING UPDATE CANCELED');
        else
            set( gui.Fbutton,'String','Reprocessing ...');
            wrk_color = gui.Fbutton.BackgroundColor;
            set(gui.Fbutton,'BackgroundColor','y');
            stat_str = sprintf('BIC = %0.3f',DATA.BIC);
            mymsg = figure('Name','UPDATING QC AND REPROCESSING...','NumberTitle','off','units','pixels','position',[500 500 200 50],'windowstyle','modal');
            uicontrol('style','text','string','PLEASE WAIT.','units','pixels','position',[75 10 50 30]);
            drawnow
            %             mymsg = msgbox('UPDATING QC AND REPROCESSING....');
            handles.info.QCadj_log    = 'FloatQCList_log.txt';
            fid = fopen([dirs.cal, handles.info.QCadj_log],'a');
            handles.info.float_name = DATA.MBARI_ID; %for consistency with "sage" mfiles...
            fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\r\n',DATA.WMO_ID,handles.info.float_name,sdn,USER, ...
                stat_str, user_input{1});
            fclose(fid);
            clear USER sdn title_txt user_input fid
            % MAKE NEW QC ADJUSTMENT LIST
            % First, have to define some variables that the function
            % NewFloatQCList_GLT.m uses.  This is the same function used by
            % sage, and should also be used for sageO2 for compatability.
            handles.info.Mprof = 0;
            DATA.QCA.O2 = []; %clear it out, then replace with new
            DATA.QCA.O2(:,1) = DATA.tableDATA(:,1); %cycle
            DATA.QCA.O2(:,2) = DATA.tableDATA(:,2); %O2 gain
            DATA.QCA.O2(:,3) = zeros(size(DATA.tableDATA,1),1); %insert column of zeros for "offset" (only for compatability)
            DATA.QCA.O2(:,4) = DATA.tableDATA(:,3); %drift
            handles.QCA = DATA.QCA;
            handles.info.QCadj_file = [DATA.WMO_ID,'_FloatQCList.txt'];
            handles.float_IDs = MBARI_float_list(dirs);
            Draw = get_FloatViz_data(DATA.WMO_ID);
            if isempty(Draw)
                floatvM = msgbox('get_FloatViz Error.  Check existance of floatViz file.','Error');
            end
            %Now define the handles.info.XX_sensor variables.  Using the
            %floatviz files is the easiest way to do this.
            DATA.iO    = find(strcmp('Oxygen[µmol/kg]', Draw.hdr)  == 1);
            DATA.iN    = find(strcmp('Nitrate[µmol/kg]', Draw.hdr) == 1);
            DATA.iPH   = find(strcmp('pHinsitu[Total]', Draw.hdr)  == 1);
            DATA.iCHL  = find(strcmp('Chl_a[mg/m^3]', Draw.hdr)  == 1);
            if ~isempty(DATA.iO) && all(isnan(Draw.data(:,DATA.iO))); % NO DATA!!
                handles.info.O2_sensor = 0;
            else
                handles.info.O2_sensor = 1;
            end
            if ~isempty(DATA.iN) && all(isnan(Draw.data(:,DATA.iN))); % NO DATA!!
                handles.info.NO3_sensor = 0;
            else
                handles.info.NO3_sensor = 1;
            end
            if ~isempty(DATA.iPH) && all(isnan(Draw.data(:,DATA.iPH))); % NO DATA!!
                handles.info.PH_sensor = 0;
            else
                handles.info.PH_sensor = 1;
            end
            if ~isempty(DATA.iCHL) && all(isnan(Draw.data(:,DATA.iCHL))); % NO DATA!!
                handles.info.CHL_sensor = 0;
                %                handles.info.cal.CHL = []; %TM 12/10/20; not needed anymore.
            else
                handles.info.CHL_sensor = 1;
                %                handles.info.cal.CHL = 1; %TM 12/10/20; not needed anymore.
            end
            
            tf = NewFloatQCList_GLT(handles,dirs); % Make New QC list
            
 
             % SAVE OXYGEN ADJUSTMENT INFO (BY CYCLE) TO CHEM (per Yui request.  10/29/24)
			O2gain_basedir = '\\atlas\Chem\ARGO_PROCESSING\DATA\DOXY_GAIN_INFO\';
            O2gain_savedir = [O2gain_basedir,handles.info.WMO_ID,'\'];
            if ~exist(O2gain_savedir,'dir')
                mkdir(O2gain_savedir);
            end
			DOXADJ.refdata_series.data = DATA.refdata;
            DOXADJ.refdata_series.hdr = {'MATLABdate','Cycle num','Reanalysis spres (pascal)','Reanalysis pO2 (mb) -OR- surfacd % saturation'};
            if strcmp(DATA.reftag,'WOA')
                 DOXADJ.doxygain_series.data = [DATA.refdata(:,1:2) DATA.GAINS];
            else
                 DOXADJ.doxygain_series.data = [DATA.refdata(:,1:2) DATA.GAINS{1}];
            end
            DOXADJ.doxygain_series.hdr = {'MATLABdate','Cycle num','computed gain used in sageO2'};
			DOXADJ.reftag = DATA.reftag;
			DOXADJ.WMO = handles.info.WMO_ID;
            DOXADJ.Date_assessed = now;
            DOXADJ.gainadj_applied.data = DATA.QCA.O2;
            DOXADJ.gainadj_applied.hdr = {'cycle','gain','offset','drift'};
            DOXADJ.inair_series.means = DATA.O2air{1};
            DOXADJ.inair_series.stds = DATA.O2air{2};
            DOXADJ.inair_series.hdr = {'cell1 = telemetry aircal stats','cell2 = in-water air-series stats','cell3 = in-air air-series stats'};
			
            save([O2gain_savedir,handles.info.WMO_ID,'.mat'],'DOXADJ')

            % REPROCESS FLOATVIZ QC DATAFILE
            if handles.info.Mprof == 0
                tf = Process_GUI_float_GLT(handles,dirs);
            elseif handles.info.Mprof == 1
                tf = Make_Sprof_ODVQC(handles);
            end
            %set( gui.Fbutton,'String',inputs.floatID)
            set( gui.Fbutton,'String',DATA.WMO_ID)
            set(gui.Fbutton,'BackgroundColor',wrk_color);
            drawnow
        end
        close(mymsg)
    end
%-------------------------------------------------------------------------%
end % end sageO2.m
