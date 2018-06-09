function varargout = sage()
% *************************************************************************
% *************************************************************************
% ** THIS MFILE AND ASSOCIATED SOFTWARE IS COPYRIGHT UNDER MIT LICENSE.  **
% **       PLEASE SEE SAGE_MITLisence.txt FOR MORE INFORMATION         **
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
% sage.m
% SOCCOM Assessment and Graphical Evaluation GUI.
% ************************************************************************
%
%
% USE AS:  sage
%
% INPUTS:  
%
% OUTPUTS: 
%
%
% AUTHOR: Josh Plant & Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         jplant@mbari.org
%         tmaurer@mbari.org
%
% DATE: 01/10/2018
% UPDATES: %   11/02/2016 qc_table_CellEditCallback now checks for valid entry in
%              table before updating QC data set: cycle # must be monotonic
%              and increasing, only #'s will be accepted, no neg #'s in
%              cycle col. -JP
%   11/02/2016 Fixed REMOVE ROW callback so no crash if all rows removed
%              from table
%   11/02/2016 Changed behavior of Gain col in table. Only first value
%              can be changed - subsequent values will revert to first
%              value if changed.
%   11/07/2016 Fixed comp_data choice crash when going from NO3 WOA TO PH
%              Now GUI "remembers" last choice from NO3 , O2 or PH
%   11/16/2016 updated table behavior so table (and qc adjustment data) is
%              inactive if the sensor does not exist on float
%   11/30/2016 added ph_pump offset text edit box
%   12/01/2016 Fixed WOA data injestion bug - interpolate now vs table look
%   12/02/2016 added section to load float ID & type table used in
%              reprocessing to determine float type (APEX or NAVIS)
%   12/06/2016 added dialog box to accept QC adjustmets comments and add
%              comments to "FloatQCList_log.txt" in the dirs.cal directory
%   12/09/2016 Salinity and temperature choices have been added
%   03/01/2017 Auto slope & offset calculation added
%   03/10/2017 Re-process can now add non exisiting floats to FloatQCList
%   03/22/2017 Re-process now re-sets handles.info.qc_flag = 1 if QC O2
%              exists ~ line 1630
%   04/26/2017 Implement LIAR V2.2, LPHR and LINR, add button and code to
%              display credits / awknowledgements 
%   08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
%   08/15/2017 Added code for user dialogue box if multiple UW ID choices for
%               bottle data (ie 8514)
%   08/17/2017 - TM add map enhancements (global)
%   08/17/2017 Changed the bottle lookup table format (added MBARI_ID as
%              1st col) and changed code in SAGE to accomodate. MBARI_ID is
%              a better lookup paramater when dealing with duplicate UW_ID
%              #'s. This makes the 8/15 changes obsolete so they are
%              commented out
%   10/30/2017 Several changes to enable the use of text files created from Mprof.nc
%              picking a file can now occur from any directory and  shouldn't barf during
%              reprocessing. Flags to check for sensors are now operational vs false positive
%              which occured after switch to new text file format a while back.
%
%   11/14/17  MIT license added.
%   
%   01/10/18  Updated display, switched from MATLAB's guide, to GUI Layout
%             Toolbox (old version is archived as sage_version1.m)
% NOTES: 
%
% ************************************************************************
%
% ************************************************************************

 
% Declare shared variables
gui = createInterface();
handles = [];
% Define paths
dirs.user_dir = [getenv('USERPROFILE'),'\Documents\MATLAB\']; 
dirs.mfiles    = [dirs.user_dir,'ARGO_PROCESSING\MFILES\'];
dirs.mat       = [dirs.user_dir,'ARGO_PROCESSING\DATA\FLOATS\'];
dirs.cal       = [dirs.user_dir,'ARGO_PROCESSING\DATA\CAL\'];
dirs.FVlocal   = [dirs.user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ\'];
dirs.FV        = [dirs.user_dir,'ARGO_PROCESSING\DATA\FLOATVIZ\'];
%dirs.FV        = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
dirs.QCadj     = [dirs.user_dir,'ARGO_PROCESSING\DATA\CAL\QC_LISTS\'];
dirs.bottle    = [dirs.user_dir,'ARGO_PROCESSING\DATA\SHIPBOARD\'];
dirs.QC_images = [dirs.user_dir,'ARGO_PROCESSING\DATA\QC_images\'];
dirs.ERA       = [dirs.user_dir,'ARGO_PROCESSING\DATA\ERA_INT\'];
dirs.temp      = 'C:\temp\';
dirs.msg       = '\\atlas\ChemWebData\floats\';
dirs.alt       = '\\atlas\ChemWebData\floats\alternate\'; %alternate msg file directory (MBARI).  Comment out if not used.
dirs.msg_comb  = '\\atlas\ChemWebData\floats\combined\';
% Initialize variables
DATA = [];
inputs = [];


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


% SET UP INFO STRUCTURE TO STORE SENSOR EXIST FLAGS
handles.info.NO3_sensor = 0;
handles.info.O2_sensor  = 0;
handles.info.PH_sensor  = 0;
handles.info.S_sensor  = 0;
handles.info.T_sensor  = 0;

% LOAD FLOAT ID LIST
%MBARI name, UW_ID, WMO#, type
handles.float_IDs = MBARI_float_list(dirs); 
%-------------------------------------------------------------------------%
function gui = createInterface( ~ )
        % Create the user interface for the application and return a
        % structure of handles for global use.
        gui = struct();
        % Open a window and add some menus
        scrsz = get(0,'ScreenSize');
        b=0.1;
        l=0.1;
        w=0.8;
        h=0.8;
        gui.Window = figure( ...
            'Name', 'SAGE (SOCCOM Assessment and Graphical Evaluation)', ...
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
                        'String','Surface','tag','1','Value',0,'Callback',@plottype_onClicked ); 
                    rb2(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Deep','tag','2','Value',1,'Callback',@plottype_onClicked ); 
                    rb2(3) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Profile','tag','3','Value',0,'Callback',@plottype_onClicked ); 
                gui.rb2 = rb2;
                            
           % Control Box 3 (Reference Data)
           VB3 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);         
           VP3 = uix.Panel('Parent',VB3,'Padding',5,'title','Reference Data');        
                bbox = uix.VButtonBox('Parent',VP3,'Padding',5,...
                'ButtonSize',[200,20],'HorizontalAlignment','left');
                rb3(1) = uicontrol('Parent',bbox,'Style','radiobutton',...
                    'String','LIR','tag','LIR','Value',1,'Callback',@ref_onClicked ); 
                rb3(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
                    'String','CANYON','tag','CANYON','Value',0,'Callback',@ref_onClicked ); 
                rb3(3) = uicontrol('Parent',bbox,'Style','radiobutton',...
                    'String','WOA2013','tag','WOA','Value',0,'Callback',@ref_onClicked ); 
                rb3(4) = uicontrol('Parent',bbox,'Style','radiobutton',...
                    'String','Williams_50Sto80S','tag','MLR W50to80','Value',0,'Callback',@ref_onClicked ); 
                rb3(5) = uicontrol('Parent',bbox,'Style','radiobutton',...
                    'String','Williams_30Sto50S','tag','MLR W30to50','Value',0,'Callback',@ref_onClicked ); 
                gui.rb3 = rb3;
                
            % Control Box 5 (Parameter to plot) - This is really control box 4
            % as you go down the layout, but added last, so
            % calling it 5 
           VB5 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);         
           VP5 = uix.Panel('Parent',VB5,'Padding',5,'title','Parameter to Plot');        
                bbox = uix.VButtonBox('Parent',VP5,'Padding',5,...
                    'ButtonSize',[100,20],'HorizontalAlignment','left');
                    rb5(1) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Nitrate','tag','NO3','Value',1,'Callback',@plotparam_onClicked ); 
                    rb5(2) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Oxygen','tag','O2','Value',0,'Enable','on','Callback',@plotparam_onClicked ); 
                    rb5(3) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','pH','tag','PH','Value',0,'Enable','on','Callback',@plotparam_onClicked ); 
                    rb5(4) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Salinity','tag','S','Value',0,'Enable','on','Callback',@plotparam_onClicked ); 
                    rb5(5) = uicontrol('Parent',bbox,'Style','radiobutton',...
                        'String','Temperature','tag','T','Value',0,'Enable','on','Callback',@plotparam_onClicked ); 
                gui.rb5 = rb5;
                
            % Control Box 4 (QC Adjustments)
            VP4 = uix.VBox('Parent',controlLayout,'Padding',5,'BackgroundColor',BGC);     
                P4 = uix.BoxPanel('Parent',VP4,'title','Data Adjustments:');
%                 set( gui.P4{1}, 'DockFcn', {@nDock,1} );
                    P4v = uix.VBox('Parent',P4,'Padding',3,'Spacing',3);
                        gui.calcadjs = uicontrol('Parent',P4v,'Style','pushbutton','string',...
                            'CALC ADJUSTMENTS','Enable','on','Callback',@on_calcadj);
                        P4h = uix.HBox('Parent',P4v,'Padding',5,'Spacing',5);
                            gui.addrow = uicontrol('Parent',P4h,'Style','pushbutton','string',...
                                'ADD ROW','Enable','on','Callback',@on_addrow);
                            gui.removerow = uicontrol('Parent',P4h,'Style','pushbutton','string',...
                                'REMOVE ROW','Enable','on','Callback',@on_removerow);
                        gui.tbl=uitable('Parent',P4v,'Data',[1 1 0],'Enable','on',...
                            'CellEditCallback',@on_celledit);
                        gui.tbl.ColumnName = {'Cycle','Gain','Offset','Drift'};
                        gui.tbl.ColumnEditable = [true true true true];
                        gui.tbl.ColumnWidth = {50 50 50 50};
                        P4b = uicontrol('parent',P4v,'style','pushbutton','string',...
                            'RELOAD FloatQCList.txt','Callback',@on_reloadQC,'Enable','on');
                        P4c = uicontrol('parent',P4v,'style','pushbutton','backgroundcolor','red','string',...
                            'REPROCESS','fontsize',18,'Callback',@on_reprocess,'Enable','on');           

        set( controlLayout, 'Heights', [-4 -1.5 -2 -2 -5] ); %Main control vbox heights
        set( vb, 'Heights', [-1 -1 -1.75 -1.75 -1.75] ); %Float spec box heights
        set( P4v, 'Heights', [-.75 -1 -3 -.75 -.75] ); %QC adj box heights
        
        % + Create the view
        p = gui.ViewContainer;
        t = uiextras.TabPanel('Parent',p,'BackgroundColor','w');
                t1 = uix.VBox('Parent',t,'padding',10,'spacing',0.5);
                    t1a = uix.HBox('Parent',t1,'padding',0,'spacing',0);
                        t1a1 = uicontainer('Parent',t1a);
                            t1axes(1) = axes('Parent',t1a1,'Visible','off');
                        t1a2 = uicontainer('Parent',t1a);
                            t1axes(2) = axes('Parent',t1a2,'Visible','off');
                        set(t1a,'widths',[-0.85 -1])
                     t1b = uix.HBox('Parent',t1,'padding',0,'spacing',0.5);
                        t1b1 = uicontainer('Parent',t1b);
                            t1axes(3) = axes('Parent',t1b1,'Visible','off');
                        t1b2 = uicontainer('Parent',t1b);
                            t1axes(4) = axes('Parent',t1b2,'Visible','off');   
                        set(t1b,'widths',[-0.85 -1])

                t2 = uix.VBox('Parent',t,'padding',10,'spacing',0.5);
                    t2a = uix.HBox('Parent',t2,'padding',0,'spacing',0);
                        t2a1 = uicontainer('Parent',t2a);
                            t2axes(1) = axes('Parent',t2a1,'Visible','off');
                        t2a2 = uicontainer('Parent',t2a);
                            t2axes(2) = axes('Parent',t2a2,'Visible','off');
                            set(t2a,'widths',[-0.85 -1])
                     t2b = uix.HBox('Parent',t2,'padding',0,'spacing',0.5);
                        t2b1 = uicontainer('Parent',t2b);
                            t2axes(3) = axes('Parent',t2b1,'Visible','off');
                        t2b2 = uicontainer('Parent',t2b);
                            t2axes(4) = axes('Parent',t2b2,'Visible','off');
                            set(t2b,'widths',[-0.85 -1])
            gui.t1=t1;
            gui.t2=t2;
            t.TabNames = {'Raw','QC'};
            t.TabSize = 75;
            t.FontSize = 16;
            t.SelectedChild = 1;
            gui.t = t;
            gui.t1ax = t1axes;
            gui.t2ax = t2axes;
        % POSITION AXES, MOVE PLOTS LEFT SO ROOM ON RIGHT FOR STATS
        % RAW DATA TAB POSITIONING: [l b w h]
        pos1 = get(gui.t1ax(1),'Position');
        pos1(2) = 0.9*pos1(2);
        pos1(4) = 0.9*pos1(4);
        pos1(1) = 0.9*pos1(1);
        set(gui.t1ax(1),'Position',pos1);
        pos2 = get(gui.t1ax(2),'Position');
        pos2(3) = 0.85*pos2(3);
        pos2(1) = 0.6*pos2(1);
        pos2(2) = 0.9*pos2(2);
        pos2(4) = 0.9*pos2(4);
        set(gui.t1ax(2),'Position',pos2);
        pos3 = get(gui.t1ax(3),'Position');
        pos3(1) = 0.9*pos3(1);
        pos3(2) = 0.9*pos3(2);
        pos3(4) = 0.9*pos3(4);
        set(gui.t1ax(3),'Position',pos3);
        pos4 = get(gui.t1ax(4),'Position');
        pos4(3) = 0.85*pos4(3);
        pos4(1) = 0.6*pos4(1);
        pos4(2) = 0.91*pos4(2);
        pos4(4) = 0.9*pos4(4);
        set(gui.t1ax(4),'Position',pos4);
%         % QC DATA TAB POSITIONING
        pos1 = get(gui.t2ax(1),'Position');
        pos1(2) = 0.9*pos1(2);
        pos1(4) = 0.9*pos1(4);
        pos1(1) = 0.9*pos1(1);
        set(gui.t2ax(1),'Position',pos1);
        pos2 = get(gui.t2ax(2),'Position');
        pos2(3) = 0.85*pos2(3);
        pos2(1) = 0.6*pos2(1);
        pos2(2) = 0.9*pos2(2);
        pos2(4) = 0.9*pos2(4);
        set(gui.t2ax(2),'Position',pos2);
        pos3 = get(gui.t2ax(3),'Position');
        pos3(1) = 0.9*pos3(1);
        pos3(2) = 0.9*pos3(2);
        pos3(4) = 0.9*pos3(4);
        set(gui.t2ax(3),'Position',pos3);
        pos4 = get(gui.t2ax(4),'Position');
        pos4(3) = 0.85*pos4(3);
        pos4(1) = 0.6*pos4(1);
        pos4(2) = 0.91*pos4(2);
        pos4(4) = 0.9*pos4(4);
        set(gui.t2ax(4),'Position',pos4);   
    end % createInterface
%-------------------------------------------------------------------------%
    function updateInterface()
        % Update various parts of the interface in response to the inputs
        % being changed.
        
       if inputs.rorq == 1
           DATA.datatype =  handles.raw_data; 
            gui.whichAX = gui.t1ax;
       elseif inputs.rorq == 2 && isfield(handles,'new_qc_data') == 1
           DATA.datatype = handles.new_qc_data;
            gui.whichAX = gui.t2ax;
       elseif inputs.rorq == 2 && isfield(handles,'new_qc_data') == 0
           DATA.datatype = handles.qc_data;
            gui.whichAX = gui.t2ax;
       end
        
        set(gui.profmin,'String',num2str(inputs.profedit(1)));
        if max(DATA.datatype.data(:,2)) < inputs.profedit(2);
            set(gui.profmax,'String',num2str(max(DATA.datatype.data(:,2))));
        else
            set(gui.profmax,'String',num2str(inputs.profedit(2))); %last cast #
        end
        set(gui.Pmin,'String',num2str(inputs.depthedit(1))); 
        set(gui.Pmax,'String',num2str(inputs.depthedit(2))); 
%         if max(DATA.datatype.data(:,6)) > inputs.depthedit(2);
%             set(gui.Pmax,'String',num2str(ceil(max(DATA.datatype.data(:,6)))));
%         end
%         if max(DATA.datatype.data(:,6)) < inputs.depthedit(1) && inputs.SFLT == 1
%             hmsg = msgbox({'No Data within assigned depth limits','Expanding lower end of range.'},'Depth Range Adjustment');
%             set(gui.Pmin,'String',num2str(floor(max(DATA.datatype.data(:,6)))-100));
%         end

        set(gui.GLDPkm,'String',num2str(inputs.GLDPkm));
        depthmin = str2double(get(gui.Pmin,'String'));
        depthmax = str2double(get(gui.Pmax,'String'));
        PROFmin = str2double(get(gui.profmin,'String'));
        PROFmax = str2double(get(gui.profmax,'String'));
        % SET X LIMITS AND TICK LOCATIONS FOR PLOTTING (AS DATE AND CYCLE#)
        if PROFmin ~= PROFmax
            DATA.xlims{1} = [DATA.track(DATA.track(:,2)==PROFmin,1) DATA.track(DATA.track(:,2)==PROFmax,1)];
            DATA.xlims{2} = [DATA.track(DATA.track(:,2)==PROFmin,2) DATA.track(DATA.track(:,2)==PROFmax,2)];
        else
            DATA.xlims{1} = [DATA.track(DATA.track(:,2)==PROFmin,1)-1 DATA.track(DATA.track(:,2)==PROFmax,1)+1];
            DATA.xlims{2} = [DATA.track(DATA.track(:,2)==PROFmin,2)-1 DATA.track(DATA.track(:,2)==PROFmax,2)+1];
        end
        
        if nanmax(DATA.xlims{2}(2)) >= 6 % at least 6 cycles
            xtckfac{2} = floor((DATA.xlims{2}(2)-DATA.xlims{2}(1))/5);
        else
            xtckfac{2} = 1; %else increment x-axis by 1 cycle
        end
        DATA.xticks{2} = DATA.xlims{2}(1):xtckfac{2}:DATA.xlims{2}(2);
        [~,xti,~] = intersect(DATA.track(:,2),DATA.xticks{2});
        DATA.xticks{1} = DATA.track(xti,1);

        % GET FLOAT DATA, SET MISSING VALUES TO NaN, GET INDICES
        qc_flag = handles.info.qc_flag;
% % %         dRAW = handles.raw_data;
        DATA.datatype.data(DATA.datatype.data == -1e10) = NaN; % missing values

        iP    = find(strcmp('Pressure[dbar]',DATA.datatype.hdr)   == 1);

        % CHECK IF QC DATA EXISTS (O,N,PH) IF NOT IT HAS BEEN SET TO RAW
        inputs.data_str  = '';
        if qc_flag == 0
            inputs.data_str = ' !!! NO QC DATA - PLOTTING RAW !!!';
        end
       

        % CACULATE MLR OR GET WOA(USE QC DATA FOR THIS)
%         tNaN_QCO2 = isnan(handles.qc_data.data(:,iO)); %Any NaN's in QC O2 (for MLR)
%         tNaN_QCS = isnan(handles.qc_data.data(:,iS)); %Any NaN's in QC Salinity (for MLR)
        if strcmp(DATA.paramtag,'O2') == 1 || strcmp(DATA.paramtag,'S') == 1 || strcmp(DATA.paramtag,'T') == 1 || strcmp(DATA.reftag,'NOQC') == 1 
            DATA.refdata = handles.qc_data.data(:,iP) * NaN; %replace with nans for sal, temp, oxygen
        else
            switch DATA.reftag
                case 'WOA'
                    DATA.refdata = DATA.reftemp;
                case 'CANYON'
                    DATA.refdata = DATA.reftemp(:,DATA.CIND);
                case 'LIR'
                    DATA.refdata = DATA.reftemp(:,DATA.LIND);
                case {'MLR W50to80','MLR W30to50'}
                    DATA.refdata = DATA.MLRdata.(DATA.paramtag).(DATA.refs);
%                     MLR = DATA.refs.(DATA.paramtag);
%                     if ~isempty(MLR) %will be empty for salinity, temp, oxygen
%                         tMLR = isnan(handles.qc_data.data(:,iO)) | isnan(handles.qc_data.data(:,iS));
%                         potT = theta(handles.qc_data.data(:,iP), handles.qc_data.data(:,iT), handles.qc_data.data(:,iS),0);
%                         sig_theta  = density(handles.qc_data.data(:,iS), potT)-1000; %density at p =0 t= pot temp
%                         DATA.refdata = MLR.cC + handles.qc_data.data(:,iO)*MLR.cO + handles.qc_data.data(:,iS)*MLR.cS + ...
%                             handles.qc_data.data(:,iT)*MLR.cT + sig_theta*MLR.cST + handles.qc_data.data(:,iP)*MLR.cP;
%                         DATA.refdata(tMLR) = NaN;
%                     else
%                         DATA.refdata = handles.qc_data.data(:,iP) * NaN; %replace with nans
%                     end
            end
        end

        if ~isempty(DATA.IND) && ~isempty(DATA.refdata)
            inputs.y_label = DATA.datatype.hdr{DATA.IND};
            DATA.DIFF_X = DATA.datatype.data(:,DATA.IND) - DATA.refdata;
        else
            inputs.y_label = DATA.datatype.hdr{DATA.IND};
            DATA.DIFF_X = DATA.refdata .* NaN;
        end

        
                % SUBSET DATA SETS WITHIN DEPTH WINDOW (profile subset will
                % happen in PlotGuiData.  Do this because for calculating 
%         DATA.datasub = DATA.datatype.data(DATA.datatype.data(:,iP) >= depthmin & ...
%             DATA.datatype.data(:,iP) <= depthmax & DATA.datatype.data(:,2) >= PROFmin &...
%             DATA.datatype.data(:,2) <= PROFmax,:);
%         DATA.refsub = DATA.refdata(DATA.datatype.data(:,iP) >= depthmin & ...
%             DATA.datatype.data(:,iP) <= depthmax & DATA.datatype.data(:,2) >= PROFmin &...
%             DATA.datatype.data(:,2) <= PROFmax,:);
%         DATA.diffsub = DATA.DIFF_X(DATA.datatype.data(:,iP) >= depthmin & ...
%             DATA.datatype.data(:,iP) <= depthmax & DATA.datatype.data(:,2) >= PROFmin &...
%             DATA.datatype.data(:,2) <= PROFmax,:); 
        DATA.rawsub = handles.raw_data.data(handles.raw_data.data(:,iP)>=depthmin & ...
            handles.raw_data.data(:,iP)<=depthmax,:);
        DATA.qcsub = handles.qc_data.data(handles.qc_data.data(:,iP)>=depthmin & ...
            handles.qc_data.data(:,iP)<=depthmax,:); 
        
       DATA.datasub = DATA.datatype.data(DATA.datatype.data(:,iP) >= depthmin & ...
            DATA.datatype.data(:,iP) <= depthmax,:);
        DATA.refsub = DATA.refdata(DATA.datatype.data(:,iP) >= depthmin & ...
            DATA.datatype.data(:,iP) <= depthmax,:);
        DATA.diffsub = DATA.DIFF_X(DATA.datatype.data(:,iP) >= depthmin & ...
            DATA.datatype.data(:,iP) <= depthmax,:); 
        DATA.glosub = DATA.G.data(DATA.G.data(:,DATA.iGP)>=depthmin & DATA.G.data(:,DATA.iGP)<=depthmax,:);

        
                     
    end % updateInterface

%-------------------------------------------------------------------------%
    function onhelp( ~, ~ )
        % User has asked for the documentation
        open([dirs.mfiles,'GUIS\SAGE\README_sage.txt'])
    end % onHelp

%-------------------------------------------------------------------------%
    function onack( ~, ~ )
        % User has asked for the documentation
        open([dirs.mfiles,'GUIS\SAGE\SAGE_Awknowledgements.txt'])
    end % onHelp

%-------------------------------------------------------------------------%
 %%
    function on_selectfloat( ~, ~ )
        % select float data dir from dialog box
        % first set some limits and defaults
%         handles=[];
        if isfield(handles,'CGOD')
            handles = rmfield(handles,'CGOD');
        end
        if isfield(handles,'new_qc_data')
            handles = rmfield(handles,'new_qc_data');
        end
        inputs=[];
        inputs.isprof = 0;
        gui.t.SelectedChild = 1; %default to Raw upon float selection
        gui.t.SelectionChangedFcn = @on_RAWorQC;
        inputs.rorq = 1;
        inputs.depthedit = [1480 1520]; %default pressure range (deep)
        inputs.GLDPkm = 30;
        set(gui.Mbutton,'Enable','on');
        set(gui.rb3(1),'Enable','on')
        set(gui.rb3(2),'Enable','on')
        set(gui.rb3(3),'Enable','on')
        set(gui.rb3(4),'Enable','on')
        set(gui.rb3(5),'Enable','on')
        % CHOOSE FILE
        [fn,pn] = uigetfile([dirs.FVlocal,'*.txt'],'SELECT FILE');
        if ~isequal(fn, 0)  
            str = [pn,fn];
            set( gui.Fbutton,'String',fn);
            handles.info.file_name  = fn;
            handles.info.file_path  = pn
            handles.info.float_name = regexpi(fn,'\w+(?=\QC.txt)|\w+(?=\.txt)', ...
                'match', 'once');
            handles.info.UW_ID  = regexp(handles.info.float_name, ...
                '^\d{3}\d+(?=\w+)','match', 'once'); %#'s but chars follow
            handles.info.QCadj_file   = [handles.info.float_name,'_FloatQCList.txt']; 

            % Set flag for ODV file created from Mprof netcdf vs msg files
            % based on file name
            handles.info.Mprof = 0; 
            if regexp(handles.info.file_name,'^ODV', 'once')
                handles.info.Mprof = 1;
            end

            % DEFINE DATA QUALITY
            if regexpi(handles.info.file_name, 'QC.TXT', 'once')
                handles.info.data_quality = 'QC';
                handles.RAWorQC.String    = 'QC DATA';
            else
                handles.info.data_quality = 'RAW';
                handles.RAWorQC.String    = 'RAW DATA';
            end

            % GET DATA (RAW AND QC)
            % EXCLUDE FLBB & CARBONATE SYSTEM VARIABLES - ONLY CHECKING O, N, pH

            % RAW DATA FIRST
            %fv_path = [handles.dirs.FVlocal, handles.info.float_name,'.TXT'];
            fv_path = [pn, handles.info.float_name,'.TXT'];
            d = get_FloatViz_data(fv_path);

            % GET SOME RAW INDICES
            DATA.iP    = find(strcmp('Pressure[dbar]', d.hdr)  == 1);
            DATA.iT    = find(strcmp('Temperature[°C]', d.hdr)  == 1);
            DATA.iS    = find(strcmp('Salinity[pss]', d.hdr)  == 1);
            DATA.iZ    = find(strcmp('Depth[m]', d.hdr)  == 1);
            DATA.iO    = find(strcmp('Oxygen[µmol/kg]', d.hdr)  == 1);
            DATA.iOsat = find(strcmp('OxygenSat[%]', d.hdr)     == 1);
            DATA.iN    = find(strcmp('Nitrate[µmol/kg]', d.hdr) == 1);
            DATA.iPH   = find(strcmp('pHinsitu[Total]', d.hdr)  == 1);
            DATA.iCHL  = find(strcmp('Chl_a[mg/m^3]', d.hdr)  == 1);

            handles.info.CHL_sensor = 0;
            if ~isempty(DATA.iCHL) % TEST FOR CHL DATA USED WITH MPROF ADJUSTMENT FILE
                t1 = d.data(:,DATA.iCHL) ~= -1e10 & ~isnan(d.data(:,DATA.iCHL));
                if sum(t1) > 0
                    handles.info.CHL_sensor = 1;
                end
                clear t1
            end

            % CONDENSE THE DATA empty & empty+1 will be ignored
            handles.raw_data.hdr  = d.hdr([1:DATA.iZ+1,DATA.iO:DATA.iOsat+1,DATA.iN,DATA.iN+1,DATA.iPH,DATA.iPH+1]);
            handles.raw_data.data = d.data(:,[1:DATA.iZ+1,DATA.iO:DATA.iOsat+1,DATA.iN,DATA.iN+1,DATA.iPH,DATA.iPH+1]);
            handles.raw_data.data(handles.raw_data.data == -1e10) = NaN;
            clear d     

            % TRY TO GET QC DATA NEXT
            %fv_path = [handles.dirs.FVlocal,'QC\' handles.info.float_name,'QC.TXT'];
            fv_path = [pn,'QC\' handles.info.float_name,'QC.TXT'];
            d = get_FloatViz_data(fv_path);
            if ~isempty(d)
                handles.info.qc_flag = 1; % 1 = QC data exists
                DATA.iPH   = find(strcmp('pHinsitu[Total]', d.hdr)  == 1);
                handles.qc_data.hdr  = d.hdr([1:DATA.iZ+1,DATA.iO:DATA.iOsat+1,DATA.iN,DATA.iN+1,DATA.iPH,DATA.iPH+1]);
                handles.qc_data.data = d.data(:,[1:DATA.iZ+1,DATA.iO:DATA.iOsat+1,DATA.iN,DATA.iN+1,DATA.iPH,DATA.iPH+1]);
                handles.qc_data.data(handles.qc_data.data == -1e10) = NaN;
            else
                handles.info.qc_flag = 0; % NO QC DATA = 0
                handles.qc_data.hdr  = handles.raw_data.hdr;
                handles.qc_data.data = handles.raw_data.data;        
            end
            clear d

            % REDO AFTER SUBSETTING
            DATA.iO    = find(strcmp('Oxygen[µmol/kg]', handles.raw_data.hdr)  == 1);
            DATA.iOsat = find(strcmp('OxygenSat[%]', handles.raw_data.hdr)     == 1);
            DATA.iN    = find(strcmp('Nitrate[µmol/kg]', handles.raw_data.hdr) == 1);
            DATA.iPH   = find(strcmp('pHinsitu[Total]', handles.raw_data.hdr)  == 1);   

            % ONLY WANT GOOD DATA FOR QC PURPOSES & SET SENSOR EXIST FLAGS
            if ~isempty(DATA.iS)
                tbad = handles.raw_data.data(:,DATA.iS+1) == 8; % SET bad S to NaN
                handles.raw_data.data(tbad,DATA.iS) = NaN;

                tbad = handles.qc_data.data(:,DATA.iS+1) == 8;
                handles.qc_data.data(tbad,DATA.iS) = NaN;
            end

            if ~isempty(DATA.iO)
                if all(isnan(handles.raw_data.data(:,DATA.iO))); % NO DATA!!
                    handles.info.O2_sensor = 0;
                else
                    handles.info.O2_sensor = 1;
                    tbad = handles.raw_data.data(:,DATA.iO+1) == 8;
                    handles.raw_data.data(tbad,DATA.iO) = NaN;
                    handles.raw_data.data(tbad,DATA.iOsat) = NaN;

                    tbad = handles.qc_data.data(:,DATA.iO+1) == 8;
                    handles.qc_data.data(tbad,DATA.iO) = NaN;
                    handles.qc_data.data(tbad,DATA.iOsat) = NaN;
                end
            end

            if ~isempty(DATA.iN)
                if all(isnan(handles.raw_data.data(:,DATA.iN))); % NO DATA!!
                    handles.info.NO3_sensor = 0;
                else
                    handles.info.NO3_sensor = 1;
                    tbad = handles.raw_data.data(:,DATA.iN+1) == 8;
                    handles.raw_data.data(tbad,DATA.iN) = NaN;

                    tbad = handles.qc_data.data(:,DATA.iN+1) == 8;
                    handles.qc_data.data(tbad,DATA.iN) = NaN;
                end
            end

            if ~isempty(DATA.iPH)
                if all(isnan(handles.raw_data.data(:,DATA.iPH))); % NO DATA!!
                    disp('no data')
                    handles.info.PH_sensor = 0;
                else
                    handles.info.PH_sensor = 1;
                    disp('ph data')
                    tbad = handles.raw_data.data(:,DATA.iPH+1) == 8;
                    handles.raw_data.data(tbad,DATA.iPH) = NaN;

                    tbad = handles.qc_data.data(:,DATA.iPH+1) == 8;
                    handles.qc_data.data(tbad,DATA.iPH) = NaN;
                end
            end

            clear tbad 

            % LOAD CALIBRATION DATA - IT WILL BE USED LATER FOR SENSOR CHECKS
            % WHEN UPDATING FloatQCList AND MAYBE MORE STUFF
            if exist([dirs.cal,'cal',handles.info.float_name,'.mat'],'file')
                handles.info.cal = load([dirs.cal,'cal',handles.info.float_name,'.mat']);
%                 handles.info.cal = cal;
            end

            % SAVE FLOAT TRACK & GET WOA 2013 NITRATE DATA FOR TRACK
            d = handles.raw_data.data; % Get raw data
            [~,ia,~] = unique(d(:,2));
            DATA.track = d(ia,1:4); % sdn cycle lon lat
            inputs.profedit = [1 DATA.track(end,2)]; 
            DATA.iN    = find(strcmp('Nitrate[µmol/kg]',handles.raw_data.hdr) == 1);

%             set(handles.recumpute_text,'Visible','on')
%             set(handles.recumpute_text, ...
%                 'String','LOADING WOA 2013 NITRATE DATA ....')

            WOA_NO3 = get_WOA2013_local(DATA.track(:,[1,4,3]), [0 2000], 'NO3');
            % NOW MATCH WOA DATA TO RAW PROFILE DATA, sample by sample

            WNO3 = ones(size(handles.raw_data.data(:,1))) * NaN; % predim
            Z = WOA_NO3.Z; % WOA depth grid
            N = WOA_NO3.d; % WOA nitrate, µM / L

            for cast_ct = 1:size(DATA.track,1) % step through profiles
                % INTRPOLATE WOA2013 ON TO FLOAT PRESSURE PROFILE
                t1  = d(:,2) == DATA.track(cast_ct,2); % get profile
                %t1  = d(:,2) == track(cast_ct,2) & ~isnan(d(:,iN)); % get profile
                tmpZ = handles.raw_data.data(t1,6); % Get float pressure profile for cast
                WNO3(t1) = interp1(Z, N(:,cast_ct),tmpZ);
            end


            % NOW CONVERT TO µmol/kg
            potT  = theta(handles.raw_data.data(:,6), ...
                handles.raw_data.data(:,8),handles.raw_data.data(:,10),0); % P,T,S,P0
            den  = density(handles.raw_data.data(:,10), potT);
            DATA.WOA_NO3 = WNO3./den*1000; %µmol/kg
            clear WOA_NO3 d WNO3 t1 t2 tmp WOAtmp Z N ia potT den

            % ********************************************************************
            % GET ANY GLODAP DATA THAT IS NEAR ANY FLOAT TRACK POINTS
%             set(handles.recumpute_text,'Visible','on')
%             set(handles.recumpute_text, ...
%                 'String','LOADING GLODAPv2 DATA ....')
%             drawnow
%             track = DATA.track;
            d = get_GLODAPv2_local(DATA.track (:,[1,2,4,3]), ...
                inputs.GLDPkm, [0 2000]);
            handles.GLODAP = d;
            clear track d cycles i t1
            % GET GLODAPv2 DATA & SET INDICES
            DATA.G = handles.GLODAP; % Could be empty if no crossover data

            DATA.iGcyc = find(strcmp('float cycle',  DATA.G.hdr) == 1);
            DATA.iGSDN = find(strcmp('Date',  DATA.G.hdr)        == 1);
            DATA.iGP   = find(strcmp('G2pressure',DATA.G.hdr)    == 1);
            DATA.iGT   = find(strcmp('G2temperature',DATA.G.hdr) == 1);
            DATA.iGS   = find(strcmp('G2salinity',DATA.G.hdr)    == 1);
            DATA.iGO   = find(strcmp('G2oxygen',DATA.G.hdr)      == 1);
            DATA.iGN   = find(strcmp('G2nitrate',DATA.G.hdr)     == 1);
            DATA.iGPH  = find(strcmp('ph_insitu',DATA.G.hdr)     == 1);
 
            
            [handles, DATA] = get_LIR_CAN_MLR(handles,DATA);
% % % %         
% % % %             % ********************************************************************
% % % %             % GET CANYON NEURAL NETWORK, LINR & LIPHER APROXIMATIONs FOR NO3 AND PH
% % % %             % NEED QC OXYGEN FOR THIS
% % % %             % out = CANYON_jp(gtime,lat,lon,pres,temp,psal,doxy,param)
% % % %             if handles.info.qc_flag == 1
% % % %                 d = handles.qc_data;
% % % %             else
% % % %                 d.hdr = {};
% % % %             end
% % % % 
% % % %             % CANYON
% % % %             if handles.info.qc_flag == 1 && ~isempty(DATA.iO)
% % % % %                 set(handles.recumpute_text,'Visible','on')
% % % % %                 set(handles.recumpute_text, ...
% % % % %                     'String','LOADING CANYON NEURAL NETWORK NO3  & PH ESTIMATES ....')
% % % % %                 drawnow
% % % %                 canyon_no3 = CANYON_jp(d.data(:,1),d.data(:,4),d.data(:,3), ...
% % % %                     d.data(:,6),d.data(:,8),d.data(:,10),d.data(:,DATA.iO),'NO3');
% % % %                 canyon_ph = CANYON_jp(d.data(:,1),d.data(:,4),d.data(:,3), ...
% % % %                     d.data(:,6),d.data(:,8),d.data(:,10),d.data(:,DATA.iO),'PH');     
% % % %                 handles.canyon.hdr  = [d.hdr([1,2,6]),'canyon_no3','canyon_ph'];
% % % %                 handles.canyon.data = [d.data(:,[1,2,6]), canyon_no3, canyon_ph];
% % % %             else
% % % % %                 set(handles.recumpute_text,'Visible','on')
% % % % %                 set(handles.recumpute_text, ...
% % % % %                     'String','NO CANYON NEURAL NETWORK NO3  OR PH ESTIMATES!!')
% % % % %                 drawnow
% % % %                 handles.canyon =[];
% % % %             end
% % % %             
% % % %             % GET CANYON NEURAL NETWORK ESTIMATES 
% % % %             DATA.C = handles.canyon; 
% % % %             if ~isempty(DATA.C)
% % % %                 DATA.iCN   = find(strcmp('canyon_no3',DATA.C.hdr)     == 1);
% % % %                 DATA.iCPH  = find(strcmp('canyon_ph',DATA.C.hdr)      == 1);
% % % %             else
% % % %                 DATA.iCN   = [];
% % % %                 DATA.iCPH  = [];    
% % % %             end
% % % % 
% % % %             % LINR & LIPHR ESTIMATES FOR NITRATE AND PH 
% % % %             if handles.info.qc_flag == 1 && ~isempty(DATA.iO)
% % % % %                 set(handles.recumpute_text,'Visible','on')
% % % % %                 set(handles.recumpute_text, ...
% % % % %                     'String','LOADING LINR & LIPHR  NO3  & PH ESTIMATES ....')
% % % % %                 drawnow
% % % % 
% % % %                 LXXX_pos = [d.data(:,3), d.data(:,4), d.data(:,DATA.iZ)]; % lon,lat,Z
% % % %                 ptemp    = theta(d.data(:,6), d.data(:,8) ,d.data(:,10),0);
% % % % 
% % % %         %         MeasIDVec    = [1 2 6]; % PSAL, Pot_TEMP, DOXY_ADJ,
% % % %         %         Measurements = [d.data(:,10), ptemp, d.data(:,iO)];   
% % % % 
% % % %                 MeasIDVec    = [1 6 7]; % PSAL, DOXY_ADJ, TEMP, ,
% % % %                 Measurements = [d.data(:,10), d.data(:,DATA.iO), d.data(:,8)];  
% % % % 
% % % %                 Equations    = 7; % S, Theta, AOU
% % % % 
% % % %         %         [NO3_Est, Uncert_Est, MinUncert_Equ] = LINR(LXXX_pos, ...
% % % %         %             Measurements, MeasIDVec, Equations, [], 1); % update 04/25/17
% % % %         %         
% % % %         %         [PH_Est, Uncert_Est, MinUncert_Equ] = LIPHR(LXXX_pos, ...
% % % %         %             Measurements, MeasIDVec, Equations, [], 1); % update 04/25/17 
% % % % 
% % % %                 [NO3_Est, Uncert_Est, MinUncert_Equ] = LINR(LXXX_pos, ...
% % % %                     Measurements, MeasIDVec,'Equations', Equations); % update 08/11/17
% % % % 
% % % %                 [PH_Est, Uncert_Est, MinUncert_Equ] = LIPHR(LXXX_pos, ...
% % % %                     Measurements, MeasIDVec, ...
% % % %                     'Equations', Equations,'OAAdjustTF',false); % update 08/11/17        
% % % % 
% % % %                 %
% % % % 
% % % %                 handles.LIR.hdr  = [d.hdr([1,2,6]),'LIR_no3','LIR_ph'];
% % % %                 handles.LIR.data = [d.data(:,[1,2,6]), NO3_Est, PH_Est];
% % % %             else
% % % % %                 set(handles.recumpute_text,'Visible','on')
% % % % %                 set(handles.recumpute_text, ...
% % % % %                     'String','NO LINR or LIPHER NO3  OR PH ESTIMATES!!')
% % % %                 handles.LIR =[];
% % % %             end
% % % % 
% % % %             clear canyon_no3 canyon_ph LXXX_pos  ptemp  MeasIDVec
% % % %             clear Measurements Equations NO3_Est PH_Est Uncert_Est MinUncert_Equ
% % % % 
% % % %         % GET LINR & LIPHER NO3 & PH ESTIMATES 
% % % %         DATA.L = handles.LIR ;
% % % %         if ~isempty(DATA.L)
% % % %             DATA.iLN   = find(strcmp('LIR_no3',DATA.L.hdr)     == 1);
% % % %             DATA.iLPH  = find(strcmp('LIR_ph',DATA.L.hdr)      == 1);
% % % %             DATA.reftemp = DATA.L.data;
% % % %             DATA.reftag = 'LIR'; %will be the default on SELECT FLOAT
% % % %         else
% % % %             DATA.iLN   = [];
% % % %             DATA.iLPH  = [];    
% % % %             DATA.reftag = 'NOQC'; %No QC has been done.
% % % %         end

            % GET CALIBRATION BOTTLE DATA IF IT EXISTS
            % LOOKUP TABLE  HEADER = [MBARI_ID UW_ID  WMO   CRUISE   STN   CAST   Data file]
            % LOAD BOTTLE DATA LOOK UP TABLE & STORE IN handles STRUCTURE
            % LOOK UP TABLE  HEADER = [MBARI_ID UW_ID  WMO   CRUISE   STN   CAST   Data file]
            fid = fopen([dirs.bottle,'BottleData_lookup_table.txt' ]);
            d   = textscan(fid, '%s %s %s %s %f %f %s','Delimiter', '\t','HeaderLines',1);
            fclose(fid);
            handles.bottle_lookup = d;
            clear d fid

            ind = strcmpi(handles.info.float_name, handles.bottle_lookup{1,1});

            if sum(ind) > 0; % float(s) exists in lookup table
                bottle_fname = handles.bottle_lookup{1,7}{ind};
                stn = handles.bottle_lookup{1,5}(ind);

                cst = handles.bottle_lookup{1,6}(ind);

                % data file & data exist for float
                if ~isempty(bottle_fname) && stn ~= -999 && cst ~= -999
                    d = get_shipboard_data([dirs.bottle,bottle_fname]);
                    DATA.iStn  = find(strcmp(d.hdr,'STNNBR') == 1);
                    DATA.iCast = find(strcmp(d.hdr,'CASTNO') == 1);
                    DATA.tStn  = d.data(:,DATA.iStn)  == stn;
                    DATA.tCast = d.data(:,DATA.iCast) == cst;
                    handles.bdata.cruise  = d.cruise;
                    handles.bdata.hdr     = d.hdr;
                    %handles.bdata.units   = d.units;
                    d.data(d.data == -999) = NaN;
                    handles.bdata.data    = d.data(DATA.tStn&DATA.tCast,:);
                else
                    handles.bdata.Cruise = '';
                    handles.bdata.hdr    = '';
                    handles.bdata.data   = [];
                end
                clear bottle_fname stn cst d IStn iCast tStn tCast
            else
                handles.bdata.Cruise = '';
                handles.bdata.hdr    = '';
                handles.bdata.data   = [];
            end

            b = handles.bdata; % Could be empty if no calibration data
            DATA.b = b;

            DATA.ibSDN = find(strcmp('DATE',  b.hdr) == 1);
            DATA.ibZ   = find(strcmp('DEPTH', b.hdr) == 1);
            DATA.ibP   = find(strcmp('CTDPRS',b.hdr) == 1);
            DATA.ibT   = find(strcmp('CTDTMP',b.hdr) == 1);
            DATA.ibS   = find(strcmp('CTDSAL',b.hdr) == 1);
            DATA.ibO   = find(strcmp('OXYGEN',b.hdr) == 1);
            DATA.ibN   = find(strcmp('NITRAT',b.hdr) == 1);
            DATA.ibPH1  = find(strcmp('PH_TOT_INSITU',b.hdr) == 1);
            DATA.ibPH2  = find(strcmp('PH_TOT_INSITU_ALKDIC',b.hdr) == 1);
        
           % MLR ESTIMATES FOR NITRATE AND PH 
%            DATA.MLR = LoadGuiMLR_GLT;

           % SET SOME MORE DEFAULTS UPON LOADING FLOAT
           set(gui.rb3(1),'Value',1)
           set(gui.rb3(gui.rb3~=gui.rb3(1)),'Value',0)
           set(gui.rb5(1),'Value',1)
           set(gui.rb5(gui.rb5~=gui.rb5(1)),'Value',0)
%            DATA.reftemp = DATA.L.data;
%            DATA.reftag = 'LIR';
           DATA.paramtag = 'NO3';
           DATA.IND  = DATA.iN;
           DATA.bIND = DATA.ibN;
           DATA.GIND = DATA.iGN;
           DATA.CIND = DATA.iCN;
           DATA.LIND = DATA.iLN;
%                    if max(handles.raw_data.data(:,6)) < inputs.depthedit(1) 
%             hmsg = msgbox({'No Data within assigned depth limits','Expanding lower end of range.'},'Depth Range Adjustment');
%             inputs.depthedit(1) = floor(max(DATA.datatype.data(:,6)))-100;
%         end
           
            % GET & DISPLAY QC ADJUSTMENTS IN TABLE
            handles.add_row.Enable = 'on'; % Don't activate until data loaded
            handles.remove_row.Enable = 'on';
            DATA.QCA_path = [dirs.QCadj,handles.info.QCadj_file];
            QCA      = get_QCA(DATA.QCA_path,handles.info.float_name);
            handles.QCA = QCA;
            if ~isempty(handles.QCA.(DATA.paramtag))
                DATA.tableDATA=handles.QCA.(DATA.paramtag);
            else
                DATA.tableDATA=[1 1 0 0]; % NO QC so add start row
            end
            set(gui.tbl,'Data',DATA.tableDATA)

           inputs.rorq = 1; %start with raw data
           DATA.paramrefnum = 1; %start with LIR
           % PLOT THE DATA
           inputs.SFLT = 1;
           updateInterface()
           inputs.SFLT = 2;
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
%             plot
        end
    end %end selectfloat
%%
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
        title(['FLOAT ',handles.info.float_name],'FontSize', 16)
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
       updateInterface()
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
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
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
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
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
       end
   end

%-------------------------------------------------------------------------%        
   function on_GLODAP( source, ~ ) 
       GDkm = get(source,'String');
       inputs.GLDPkm = str2double(GDkm);
       d = get_GLODAPv2_local(DATA.track (:,[1,2,4,3]), ...
        inputs.GLDPkm, [0 2000]);
        handles.GLODAP = d;
        % GET GLODAPv2 DATA & SET INDICES
        DATA.G = handles.GLODAP; % Could be empty if no crossover data

        DATA.iGcyc = find(strcmp('float cycle',  DATA.G.hdr) == 1);
        DATA.iGSDN = find(strcmp('Date',  DATA.G.hdr)        == 1);
        DATA.iGP   = find(strcmp('G2pressure',DATA.G.hdr)    == 1);
        DATA.iGT   = find(strcmp('G2temperature',DATA.G.hdr) == 1);
        DATA.iGS   = find(strcmp('G2salinity',DATA.G.hdr)    == 1);
        DATA.iGO   = find(strcmp('G2oxygen',DATA.G.hdr)      == 1);
        DATA.iGN   = find(strcmp('G2nitrate',DATA.G.hdr)     == 1);
        DATA.iGPH  = find(strcmp('ph_insitu',DATA.G.hdr)     == 1);
       updateInterface()
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
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
            PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
        else
            inputs.isprof = 0;
            updateInterface()
            PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)  
        end
   end

%-------------------------------------------------------------------------%
    function plotparam_onClicked( source, ~ ) 
        r5=gui.rb5;
        source.Value = 1; % select this 
        set( r5(r5~=source), 'Value', 0 ) % unselect others
        DATA.paramtag = get(source,'tag');
        % CHECK DATA TYPE
        switch DATA.paramtag
            case 'NO3'
                DATA.IND  = DATA.iN;
                DATA.bIND = DATA.ibN;
                DATA.GIND = DATA.iGN;
                DATA.CIND = DATA.iCN;
                DATA.LIND = DATA.iLN;
                set(gui.rb3(1),'Value',0,'Enable','on')
                set(gui.rb3(2),'Value',0,'Enable','on')
                set(gui.rb3(3),'Value',0,'Enable','on')
                set(gui.rb3(4),'Value',0,'Enable','on')
                set(gui.rb3(5),'Value',0,'Enable','on')
            case 'PH'
                DATA.IND  = DATA.iPH; 
                DATA.bIND = DATA.ibPH1;
                DATA.bIND2 = DATA.ibPH2;
                DATA.GIND = DATA.iGPH;
                DATA.CIND = DATA.iCPH;
                DATA.LIND = DATA.iLPH;
                set(gui.rb3(1),'Value',0,'Enable','on')
                set(gui.rb3(2),'Value',0,'Enable','on')
                set(gui.rb3(3),'Value',0,'Enable','off')
                set(gui.rb3(4),'Value',0,'Enable','on')
                set(gui.rb3(5),'Value',0,'Enable','on')
            case 'O2'
                DATA.IND  = DATA.iO;
                DATA.bIND = DATA.ibO;
                DATA.GIND = DATA.iGO;
                set(gui.rb3(1),'Value',0,'Enable','off')
                set(gui.rb3(2),'Value',0,'Enable','off')
                set(gui.rb3(3),'Value',0,'Enable','off')
                set(gui.rb3(4),'Value',0,'Enable','off')
                set(gui.rb3(5),'Value',0,'Enable','off')
             case 'S'
                DATA.IND  = DATA.iS;
                DATA.bIND = DATA.ibS;
                DATA.GIND = DATA.iGS;
                set(gui.rb3(1),'Value',0,'Enable','off')
                set(gui.rb3(2),'Value',0,'Enable','off')
                set(gui.rb3(3),'Value',0,'Enable','off')
                set(gui.rb3(4),'Value',0,'Enable','off')
                set(gui.rb3(5),'Value',0,'Enable','off')
            case 'T'
                DATA.IND  = DATA.iT;
                DATA.bIND = DATA.ibT;
                DATA.GIND = DATA.iGT;
                set(gui.rb3(1),'Value',0,'Enable','off')
                set(gui.rb3(2),'Value',0,'Enable','off')
                set(gui.rb3(3),'Value',0,'Enable','off')
                set(gui.rb3(4),'Value',0,'Enable','off')
                set(gui.rb3(5),'Value',0,'Enable','off')
        end
        set(gui.rb3(DATA.paramrefnum),'Value',1);
        % Set appropriate Adjustments:
        if isfield(handles.QCA,DATA.paramtag)
            set(gui.tbl,'Enable','on')
            if ~isempty(handles.QCA.(DATA.paramtag))
                DATA.tableDATA=handles.QCA.(DATA.paramtag);
            else
                DATA.tableDATA=[1 1 0 0]; % NO QC so add start row
            end
        else
            DATA.tableDATA=[1 1 0 0]; % NO QC so add start row
            set(gui.tbl,'Data',DATA.tableDATA)
            set(gui.tbl,'Enable','off')
        end
        set(gui.tbl,'Data',DATA.tableDATA)
        updateInterface()
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
       end
    end

%-------------------------------------------------------------------------%
    function ref_onClicked( source, ~ ) 
        r3=gui.rb3;
        source.Value = 1; % select this 
        set( r3(r3~=source), 'Value', 0 ) % unselect others 
        % get track data
        reftag = get(source,'tag');
        if (strcmp(reftag,'MLR W50to80')) == 1
            DATA.refs = 'Williams_50Sto80S';
            DATA.paramrefnum = 4;
        elseif (strcmp(reftag,'MLR W30to50')) == 1
            DATA.refs = 'Williams_30Sto50S';
            DATA.paramrefnum = 5;
        elseif (strcmp(reftag,'WOA')) == 1
            DATA.reftemp = DATA.WOA_NO3;
            DATA.paramrefnum = 3;
        elseif (strcmp(reftag,'CANYON')) == 1
            DATA.reftemp = DATA.C.data;
            DATA.paramrefnum = 2;
        elseif (strcmp(reftag,'LIR')) == 1
            DATA.reftemp = DATA.L.data;
            DATA.paramrefnum = 1;
        end
        DATA.reftag = reftag;
        updateInterface()
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
       end
    end
    
%-------------------------------------------------------------------------%
    function on_addrow( source, ~ ) 
        celldata = gui.tbl.Data;
        newstrt = celldata(end,1)+1;
        newrow = [newstrt 1 0 0];
        new_celldata = [celldata;newrow];
%         new_ends = [new_celldata(2:end,1);inputs.cyEND];
        DATA.tableDATA=new_celldata;
        set(gui.tbl,'Data',DATA.tableDATA);
    end 

%-------------------------------------------------------------------------%
    function on_removerow( source, ~ ) 
        celldata = gui.tbl.Data;
        if size(celldata,1) == 1 %do not clear entire table
            new_celldata = [1 1 0 0];
        else
            new_celldata = celldata(1:end-1,:);
        end
        set(gui.tbl,'Data',new_celldata);
    end

%-------------------------------------------------------------------------%
    function on_celledit( source, callbackdata ) 
        set(gui.tbl,'Data',source.Data)
        DATA.tableDATA = source.Data;
        handles.QCA.(DATA.paramtag) = DATA.tableDATA;
        handles.new_qc_data = apply_GUIQC_corr_GLT(handles,DATA);
        if strcmp(DATA.paramtag,'O2') == 1 %O2 gain value was modified
            Omsg = figure('Name','UPDATING O2 AND RECALCULATING REFERENCE FIELDS...','NumberTitle','off','units','pixels','position',[500 500 200 50],'windowstyle','modal');
            uicontrol('style','text','string','PLEASE WAIT.','units','pixels','position',[75 10 50 30]);    
            [handles, DATA] = get_LIR_CAN_MLR(handles,DATA);
        end
        updateInterface()
		if strcmp(DATA.paramtag,'O2') == 1 %O2 gain value was modified
			close(Omsg)
		end
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
       end
    end

%-------------------------------------------------------------------------%
    function on_calcadj( source, ~)
        handles.CGOD  = getGUIQC_M_B_GLT(DATA,handles);
        set(gui.tbl,'Data',handles.CGOD)
        handles.QCA.(DATA.paramtag) = handles.CGOD;
        DATA.tableDATA = handles.QCA.(DATA.paramtag);
        handles.new_qc_data = apply_GUIQC_corr_GLT(handles,DATA);
        updateInterface()
       if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
       end
    end

%-------------------------------------------------------------------------%
    function on_reloadQC( source, ~ ) 
%         DATA.QCA = get_QCA(inputs.qc_path,DATA.floatNAME);
        handles.QCA      = get_QCA(DATA.QCA_path,handles.info.float_name);
        if ~isempty(handles.QCA.(DATA.paramtag))
            DATA.tableDATA=handles.QCA.(DATA.paramtag);
        else
            DATA.tableDATA=[1 1 0 0]; % NO QC so add start row
            msgbox('ERROR: FLOAT HAS NOT YET BEEN QUALITY CONTROLLED');
        end
        set(gui.tbl,'Data',DATA.tableDATA)
        handles.new_qc_data = apply_GUIQC_corr_GLT(handles,DATA);
        updateInterface()
        if inputs.isprof == 1 %profile selected?
           PlotGuiData_profile_GLT(dirs,gui,DATA,inputs,handles)
       else
           PlotGuiData_GLT(dirs,gui,DATA,inputs,handles)
       end 
    end

%-------------------------------------------------------------------------%
    function on_reprocess( source, ~ )
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
            mymsg = figure('Name','UPDATING QC AND REPROCESSING...','NumberTitle','off','units','pixels','position',[500 500 200 50],'windowstyle','modal');
            uicontrol('style','text','string','PLEASE WAIT.','units','pixels','position',[75 10 50 30]);       
            %             mymsg = msgbox('UPDATING QC AND REPROCESSING....');
            handles.info.QCadj_log    = 'FloatQCList_log.txt';
            fid = fopen([dirs.cal, handles.info.QCadj_log],'a');
            fprintf(fid,'%s\t%s\t%s\t%s\r\n',handles.info.UW_ID,sdn,USER, ...
            user_input{1});
            fclose(fid);
            clear USER sdn title_txt user_input fid
            
        % MAKE NEW QC ADJUSTMENT LIST
            tf = NewFloatQCList_GLT(handles,dirs); % Make New QC list
            
        % REPROCESS FLOATVIZ QC DATAFILE
            if handles.info.Mprof == 0
                tf = Process_GUI_float_GLT(handles,dirs);
            elseif handles.info.Mprof == 1
               tf = Make_Mprof_ODVQC(handles);     
            end
        
        end
        close(mymsg)
    end
%-------------------------------------------------------------------------%
% %     function nDock( eventSource, eventData, whichpanel )
% %         gui.P4{whichpanel}.IsDocked = ~gui.P4{whichpanel}.IsDocked;
% %          if gui.P4{whichpanel}.IsDocked
% %             % Put it back into the layout
% %             newfig = get( gui.P4{whichpanel}, 'Parent' );
% %             set( gui.P4{whichpanel}, 'Parent', vbox );
% %             delete( newfig );
% %         else 
% %             % Take it out of the layout
% %             pos = getpixelposition( gui.P4{whichpanel} );
% %             newfig = figure( ...
% %                 'Name', get( gui.P4{whichpanel}, 'Title' ), ...
% %                 'NumberTitle', 'off', ...
% %                 'MenuBar', 'none', ...
% %                 'Toolbar', 'none', ...
% %                 'CloseRequestFcn', {@nDock, whichpanel} );
% %             figpos = get( newfig, 'Position' );
% %             set( newfig, 'Position', [figpos(1,1:2), pos(1,3:4)] );
% %             set( gui.P4{whichpanel}, 'Parent', newfig, ...
% %                 'Units', 'Normalized', ...
% %                 'Position', [0 0 1 1] );
% %         end 
% %     end % nDock
%-------------------------------------------------------------------------%

end % end sageGLT.m
