function d = Merge_dura_msgs(MBARI_ID_str, dirs)
% ************************************************************************
% PURPOSE:
%    This function processes raw dura message files for a given APEX float
%    (*.dura), calculates useful values from the raw signals and then
%    merges data for all profiles.
%
% USAGE:
%	d = Merge_dura_msgs(MBARI_ID_str, dirs)
%
% INPUTS:
%   UW_ID_str  = UW/MBARI Float ID # as a string
%
%
%   dirs       = Either an empty variable ora structure with directory
%                strings where files are located. It must contain these
%                fields at a minimum or be empty. If dirs is empty default
%                paths will be used.
%
% OUTPUTS:
%
% REVISONS:
% 12/13/20 - jp - Forced all fopen writes to UTF-8, because that is the
%     new default for Matlab 2020 and better cross platform sharing
% 12/21/20 - JP, Added header line to txt files to alert ODV that format is UTF-8, //<Encoding>UTF-8</Encoding>
% 03/202021 - JP - updated code for new GOBGC naming conventions
plot_flag = 0;
print_flag = 1;
d = [];
% ************************************************************************
% FOR TESTING

% MBARI_ID_str  = 'ua19719';
% dirs =[];

% ************************************************************************
% SET FORMATS, DEFAULT DIRS, PREDIMENSION STRUCTURE, BOOKKEEPING
% ************************************************************************
RCR.S    = [26 38]; % from argo parameter list
RCR.T    = [-2.5 40]; % from argo parameter list
RCR.PH   = [7.0 8.8];
RCR.PHV  = [-1.2 -0.7]; % Range check on pH volts
fclose all; % CLOSE ANY OPEN FILES
tf_float = 0; % Default flag for float processing success 0 = no good
MVI_str = '-1e10';

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\'];
    dirs.msg   = '\\seaecho.shore.mbari.org\floats\';
    %dirs.mat   = 'C:\Users\jplant\Documents\MATLAB\ARGO\DURA\';
    dirs.cal   = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
    dirs.temp  = 'C:\temp\';
%     dirs.FV    = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
%                   'DATA\FLOATVIZ\'];
    dirs.FV    = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\';
                             
    dirs.save = [user_dir,'ARGO_PROCESSING\', ...
                  'DATA\PH_DIAGNOSTICS\'];
    load('\\atlas\Chem\ARGO_PROCESSING\DATA\PH_DIAGNOSTICS\MBARI_pHsensor_versions.mat')
    dirs.ph_ver = d;
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end
%INFO = dirs.INFO;

     
% ************************************************************************
%LOAD CAL FILE TO CALC pH with later
if exist([dirs.cal,'cal',MBARI_ID_str,'.mat'],'file')
    load([dirs.cal,'cal',MBARI_ID_str,'.mat'])
else
    disp(['Could not find cal file for ',MBARI_ID_str])
    return
end

% CHECK FOR pH CAL
if ~isfield(cal,'pH')
    disp(['No pH calibration found for ',MBARI_ID_str,': exiting'])
    return
end
    
INFO.INST_ID = cal.info.INST_ID; % a string
INFO.WMO     = cal.info.WMO_ID; % a string
INFO.Program = cal.info.Program; % a string
INFO.Region  = cal.info.Region; % a string
INFO.msg_dir = cal.info.msg_dir; % a string
% use regexp to strip any letters out (Seabird SN sometimes)
INFO.DF_num = str2double(regexp(cal.pH.SN,'\d+','once','match'));
INFO.k0     = cal.pH.k0;
INFO.k2     = cal.pH.k2;
INFO.pcoef  = cal.pH.pcoefs;
INFO.VER    = NaN;

%CHECK FOR PH BUILD VERSION TABLE
% ph_ver = get_ph_version; % FOR TESTING, REMOVE WHEN ALL GOOD
% if 1 ==1 % TESTING
if isfield(dirs,'ph_ver')
   ph_ver = dirs.ph_ver;
    %{'DF ID','UW ID','MSC','K2','K0_pON','sensor version','verion notes','Build notes'}
    vUW  = find(strcmp(ph_ver.hdr,'UW ID') == 1);
    VUW  = cell2mat(ph_ver.data(:,vUW));
    vVER = find(strcmp(ph_ver.hdr,'sensor version') == 1);
    t1   = VUW == str2double(INFO.INST_ID); % dup float will flag twice
    if sum(t1) == 1
        INFO.VER = cell2mat(ph_ver.data(t1,vVER));
    elseif sum(t1) > 1
        INFO.VER = cell2mat(ph_ver.data(t1,vVER));
        INFO.VER = INFO.VER(1);
    end
end

% INFO.DF_num = cal.pH.SN;
% DF_info = [INFO.DF_ver, INFO.K2, INFO.KO];

% ************************************************************************
% BUILD *.DURA LIST
% ************************************************************************
dlist = get_msg_list(cal.info,'dura'); % dura file list

% CHECK FOR ODD FILE NEMES & REMOVE (i.e. 19719 has "._19719.001.dura"
t1 = cellfun(@isempty,regexp(dlist.list(:,1),'^\d','once'));
if sum(t1) > 0 % non standard dura files
    dlist.list = dlist.list(~t1,:);
    disp(['Odd dura file removed from list for ',MBARI_ID_str])
end

% ************************************************************************
% COPY FLOAT MESSAGE FILES TO LOCAL TEMPORARY DIRECTORY
% Three file groups *.msg, *.isus, *.dura
% ************************************************************************
if isempty(dlist.list)
    disp(['No float dura files found to process for float ',MBARI_ID_str]);
    return
    end
    
% COPY *.DURA files
    if ~isempty(ls([dirs.temp,'*.dura']))
        delete([dirs.temp,'*.dura']) % Clear any message files in temp dir
    end

for i = 1:size(dlist.list,1)
    fp = fullfile(dlist.list{i,2}, dlist.list{i,1});
    copyfile(fp, dirs.temp);
end

% ************************************************************************
% PROCESS APEX MESSAGE FILES FOR GIVEN FLOAT
% ************************************************************************
msg_list   = ls([dirs.temp,'*.dura']); % get list of file names to process

disp(['Processing ARGO float ',MBARI_ID_str, '.........'])

merge_pH =[];
pH_hdr_chk = 0;
for msg_ct = 1:size(msg_list,1)
    pH_file = strtrim(msg_list(msg_ct,:));
    % find block of numbers then look ahead to see if '.dura' follows
    cast_num = regexp(pH_file,'\d+(?=\.dura)','once','match');
    cast_num = str2double(cast_num);
    
    % ****************************************************************
    % CALCULATE pH (µmol / kg scale)
    % ****************************************************************
    % GET  DATA FROM *.dura FOR DIAGNOSTICS
    dura = parse_pHmsg([dirs.temp,pH_file]);
    
    if isempty(dura.data)
        disp(['No pH data returned for ',pH_file])
        delete([dirs.temp,pH_file]);
        continue
    end
    
    if pH_hdr_chk == 0
        pH_hdr  = dura.hdr;
        pH_hdr  = ['Cast #', pH_hdr, 'pH_tot'];
        MSC_ver = dura.MSC_FW_ver;
        pH_hdr_chk = 1;
    end
    pH_data = dura.data;
    
   
    
    % CHECK FOR NAN's IN TIME & REMOVE - INCOMPLETE FILE ISSUE with pH
    % msg parser. I SHould really fix there eventually
    tnan = isnan(pH_data(:,1));
    pH_data(tnan,:) =[];
    
    %         {'Date/Time','CTD Pres','CTD Temp','CTD Sal','Internal Temp',
    %             'Internal Humidity','System Power, V','System Power, I',
    %             'Backup Batt, V','Vrs','std Vrs','Vk','std Vk','Ik','Ib'}
    
    % FIX GPS TIME STAMP BUG IN 9634
    if strcmp(MBARI_ID_str,'9634SOOCN') && cast_num > 137
        fprintf('%s Cycle %0.0f GPS time bug has been corrected\n', ...
            MBARI_ID_str,  cast_num)
        pH_data(:,1) = pH_data(:,1) + 1024*7;
    end
    
    % LOOK FOR REALLY BAD SALINITY VALUES & EXCLUDE FROM PH CALC
    % otherwise imaginary number land
    tS = pH_data(:,4) > 0;
    if sum(~tS) > 0
        disp(['Negative salinity values detected & excluded from ph ', ...
            'calculation for cycle ', num2str(cast_num)]);
    end
    
    phtot = ones(size(pH_data,1),1) * NaN;
    
    [phfree,phtot_tmp] = phcalc(pH_data(tS,10), pH_data(tS,2), ...
        pH_data(tS,3), pH_data(tS,4), cal.pH.k0, cal.pH.k2, cal.pH.pcoefs);
    phtot(tS) = phtot_tmp;
    
    realPH = isreal(phtot);
    badPH  = phtot < RCR.PH(1) & phtot > RCR.PH(2);
    if ~realPH
        phtot = pH_data(:,1)*NaN;
        disp(['pH calc contains imaginary numbers - sensor failure ', ...
            'setting ph values = NaN'])
    end
    
    if sum(badPH > 0)
        phtot = pH_data(:,1)*NaN;
        disp(['pH out of range - setting =NAN']);
    end
    
    pH_data = [ones(size(pH_data(:,1)))*cast_num,pH_data,phtot];
    
    delete([dirs.temp,pH_file]);
    
    merge_pH = [merge_pH;pH_data];
    clear dura pH_data ph_tmp
    
end
d.hdr = pH_hdr;
d.data = merge_pH;

if plot_flag == 1 && size(unique(merge_pH(:,1)),1) > 1 % need >1 profiles for contour

    iSDN = find(strcmp(pH_hdr,'Date/Time') == 1);
    iC   = find(strcmp(pH_hdr,'Cast #') == 1);
    iP   = find(strcmp(pH_hdr,'CTD Pres') == 1);
    iH   = find(strcmp(pH_hdr,'Internal Humidity') == 1);
    iBB  = find(strcmp(pH_hdr,'Backup Batt, V') == 1);
    iVRS = find(strcmp(pH_hdr,'Vrs') == 1);
    iSVRS = find(strcmp(pH_hdr,'std Vrs') == 1);
    iVK  = find(strcmp(pH_hdr,'Vk') == 1);
    iSVK = find(strcmp(pH_hdr,'std Vk') == 1);
    iIK  = find(strcmp(pH_hdr,'Ik') == 1);
    iIB  = find(strcmp(pH_hdr,'Ib') == 1);

% ****************************************************************

    % MAKE SOME CONTOUR PLOT
    % {'Cast #','Date/Time','CTD Pres','CTD Temp','CTD Sal','Internal Temp',
    %     'Internal Humidity','System Power, V','System Power, I','Backup Batt, V',
    %     'Vrs','std Vrs','Vk','std Vk','Ik','Ib'}
    
    P          = [0:5:100,110:10:400,450:50:1000,1100:100:2000]';
    rP         = size(P,1);
    cycles     = unique(merge_pH(:,iC));
    r_cycles   = size(cycles,1);
    interp_ph  = ones(rP*r_cycles, size(pH_hdr,2))*NaN;
    line_ct    = 1;
    for cyc_ct = 1:r_cycles
        t1       = merge_pH(:,iC) == cycles(cyc_ct);
        tmp      = merge_pH(t1,:);
        [~,ia,~] = unique(tmp(:,iP));
        tmp      = tmp(ia,:); % sorted and unique
        rtmp     = size(tmp,1);
        tg       = ~isnan(tmp(:,iP));
        
        % CHECK DATA SIZE FOR INTERP
        if size(tmp(tg,iP),1) < 2
            disp(['Less than 2 cycles for interp - skipping cycle ', ...
                num2str(cycles(cyc_ct))]);
            continue
        end

        interp_ph(line_ct:line_ct+rP-1,:) = interp1(tmp(tg,iP),tmp(tg,:),P);
        line_ct    = line_ct+rP;
    end
    
    % ****************************************************************
    % NOW MAKE SOME DAIGNOSTIC DATA PLOTS
    F10 = figure(10);
    F10.Units = 'Normalized';
    %F10.Position = [0.0917 0.0713 0.8000 0.8000];
    %F10.Position =[0.0917 0.3278 0.4515 0.5435];
    F10.Position =[0.0917 0.2972 0.4630 0.5741];
    %F10.Position =[0.0917 0.0972 0.6057 0.7741];
    ax(1) = subplot(3,2,1);
    ax(2) = subplot(3,2,2);
    ax(3) = subplot(3,2,3);
    ax(4) = subplot(3,2,4);
    ax(5) = subplot(3,2,5);
    ax(6) = subplot(3,2,6);
    
    ax(1).Position = ax(1).Position + [0 -0.05 0 0.05];
    ax(2).Position = ax(2).Position + [0 -0.05 0 0.05];
    ax(3).Position = ax(3).Position + [0 -0.1 0 0.05];
    ax(4).Position = ax(4).Position + [0 -0.1 0 .05];
    ax(5).Position = ax(5).Position + [0 -0.05 0 -0.025];
    ax(6).Position = ax(6).Position + [0 -0.05 0 -0.035];
    
    y_lim = [0 max(interp_ph(:,iP))];
        
    
    % ax(1) = subplot(3,2,1);
    % ax(1).Position = ax(1).Position + [0 -0.05 0 0.05];
    IX    = iVRS;
    Z     = reshape(interp_ph(:,IX), rP, r_cycles);
    t_str = [MBARI_ID_str,'    ',pH_hdr{IX}];
    if IX == iIB || IX == iIK
        Z = Z * 1e9;
        Z(abs(Z) > 1000) = NaN;
        t_str =[t_str,' * 1e9'];
    elseif IX == iVRS || IX == iVK
        Z(Z > -0.4 | Z <-1.6) = NaN;
    end
    contourf(ax(1), cycles', P, Z, 30,'LineColor','none')
    set(ax(1),'YDir','Reverse')
    ylim(ax(1),y_lim);
    title(ax(1), t_str);
    cb = colorbar(ax(1));
    
    % ax(2) = subplot(3,2,2);
    % ax(2).Position = ax(2).Position + [0 -0.05 0 0.05];
    IX    = iVK;
    Z     = reshape(interp_ph(:,IX), rP, r_cycles);
    t_str = [MBARI_ID_str,'    ',pH_hdr{IX}];
    if IX == iIB || IX == iIK
        Z = Z * 1e9;
        Z(abs(Z) > 1000) = NaN;
        t_str =[t_str,' * 1e9'];
    elseif IX == iVRS || IX == iVK
        Z(Z > -0.4 | Z <-1.6) = NaN;
    end
    contourf(ax(2), cycles', P, Z, 30,'LineColor','none')
    set(ax(2),'YDir','Reverse')
    ylim(ax(2),y_lim);
    title(ax(2), t_str);
    cb = colorbar(ax(2));
    
    % ax(3) = subplot(3,2,3);
    % ax(3).Position = ax(3).Position + [0 -0.1 0 0.05];
    IX    = iIK;
    IX    = iSVRS;
    Z     = reshape(interp_ph(:,IX), rP, r_cycles);
    t_str = [MBARI_ID_str,'    ',pH_hdr{IX}];
    if IX == iIB || IX == iIK
        Z = Z * 1e9;
        Z(abs(Z) > 1000) = NaN;
        t_str =[t_str,' * 1e9'];
    elseif IX == iVRS || IX == iVK
        Z(Z > -0.4 | Z <-1.6) = NaN;
    end
    contourf(ax(3), cycles', P, Z, 30,'LineColor','none')
    set(ax(3),'YDir','Reverse')
    ylim(ax(3),y_lim);
    title(ax(3), t_str);
    cb = colorbar(ax(3));
    
    % ax(4) = subplot(3,2,4);
    % ax(4).Position = ax(4).Position + [0 -0.1 0 0.05];
    IX    = iIB;
    IX    = iSVK;
    Z     = reshape(interp_ph(:,IX), rP, r_cycles);
    t_str = [MBARI_ID_str,'    ',pH_hdr{IX}];
    if IX == iIB || IX == iIK
        Z = Z * 1e9;
        Z(abs(Z) > 1000) = NaN;
        t_str =[t_str,' * 1e9'];
    elseif IX == iVRS || IX == iVK
        Z(Z > -0.4 | Z <-1.6) = NaN;
    end
    contourf(ax(4), cycles', P, Z, 30,'LineColor','none')
    set(ax(4),'YDir','Reverse')
    ylim(ax(4),y_lim);
    title(ax(4),t_str);
    cb = colorbar(ax(4));
    
    
    
    IX    = iIK;
    plot(ax(5),interp_ph(:,iC),interp_ph(:,IX),'ko-','MarkerSize', 4, ...
        'MarkerFaceColor','r')
    ylabel(ax(5),pH_hdr{IX})
    ax(5).Position(3) = ax(3).Position(3);
    
    
    IX    = iIB;
    plot(ax(6),interp_ph(:,iC),interp_ph(:,IX),'ko-','MarkerSize',4, ...
        'MarkerFaceColor','r')
    ylabel(ax(6), pH_hdr{IX})
    ax(6).Position(3) = ax(4).Position(3);
    
    % ax(5) = subplot(3,2,5);
    % ax(5).Position = ax(5).Position + [0 -0.1 0 -0.05];
    %
    % ax(6) = subplot(3,2,6);
    % ax(6).Position = ax(6).Position + [0 -0.1 0 -0.05];
    
    %linkaxes(ax,'x');
    % disp('I put a "return" in at line 304')
    % return
    
    print(gcf,[dirs.save,'PLOTS\',MBARI_ID_str],'-dpng');
    
end

% *********************************************************************
% CLEAN UP
if ~isempty(ls([dirs.temp,'*.dura']))
    delete([dirs.temp,'*.dura']);
end


% ************************************************************************
% PRINT MERGED DURA FILES TO TEXT
% ************************************************************************
if print_flag ==1
    % ********************************************************************
    % SET SOME QF RANGE CHECKS & ASSIGN QF's
    RC.S8   = [26 38]; % from argo parameter list
    RC.T8   = [-2.5 40]; % from argo parameter list
    RC.VK4  = [-1.18 -0.70];    % VK
    RC.VK8  = [-1.30 -0.58];
    RC.SVK4 = 1e-4;     % STD VK 100 uA
    RC.SVK8 = 3e-4; 
    RC.VRS4 = [-1.016 -0.824]; % VRS
    RC.VRS8 = [-1.05 -0.75];
    RC.IK4  = 2e-8;           % IK
    RC.IK8  = 1e-7;
    RC.IB4  = 5e-9;           % IB
    RC.IB8  = 1e-7;
    RC.VSYS4  = 13;    % VK
    RC.VSYS8  = 10.8;
    
    [rr,cc] = size(merge_pH); %predim
    pdata = ones(rr,cc*2)*NaN;
    pdata(:,1:2:2*cc-1) = merge_pH;
    pdata(:,2:2:2*cc) = 0;
    
    phdr  = cell(1,2*cc);
    phdr(1:2:2*cc-1) = pH_hdr;
    phdr(2:2:2*cc) = {'QF'};
    
    tQF = strcmp(phdr,'QF');
    pvars = phdr(~tQF);

    clear rr cc
    
    for i = 1:size(pvars,2)
        X  = pvars{1,i};
        iX = find(strcmp(phdr,X) == 1);
        tnan = isnan(pdata(:,iX));
        if ~strcmp(X,'QF') % reoces NaN's with fill value and set QF=1
            %pdata(tnan,iX) = -1e10; % will put fill value in layer as str
            pdata(tnan,iX+1) = 1;
        end
            
        if strcmp(X,'QF')
            continue
        elseif strcmp(X,'Vk')     
            t4 = pdata(:,iX) < RC.VK4(1) | pdata(:,iX) > RC.VK4(2);
            t8 = pdata(:,iX) < RC.VK8(1) | pdata(:,iX) > RC.VK8(2);
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'Vrs')
            t4 = pdata(:,iX) < RC.VRS4(1) | pdata(:,iX) > RC.VRS4(2);
            t8 = pdata(:,iX) < RC.VRS8(1) | pdata(:,iX) > RC.VRS8(2);
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'std Vk') % should never be neg but 9099 cycle 119 some -3's
            t4 = abs(pdata(:,iX)) > RC.SVK4(1); 
            t8 = abs(pdata(:,iX)) > RC.SVK8(1);
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'Ik')
            t4 = abs(pdata(:,iX)) > RC.IK4(1);
            t8 = abs(pdata(:,iX)) > RC.IK8(1);
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;   
        elseif strcmp(X,'IB')
            t4 = abs(pdata(:,iX)) > RC.IB4(1);
            t8 = abs(pdata(:,iX)) > RC.IB8(1);
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;  
        elseif strcmp(X,'Salinity[pss]')
            t8 = pdata(:,iX) < RC.S8(1) | pdata(:,iX) > RC.S8(2);
            pdata(t8&~tnan,iX+1) = 8;  
        elseif strcmp(X,'Temperature[deg_C]')
            t8 = pdata(:,iX) < RC.T8(1) | pdata(:,iX) > RC.T8(2);
            pdata(t8&~tnan,iX+1) = 8;  
        elseif strcmp(X,'System_Power_V')
            t4 = abs(pdata(:,iX)) < RC.VSYS4(1);
            t8 = abs(pdata(:,iX)) < RC.VSYS8(1);
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;  
        end
    end

    iSDN = find(strcmp(phdr,'Date/Time') == 1);
    iC   = find(strcmp(phdr,'Cast #') == 1);
    iP   = find(strcmp(phdr,'CTD Pres') == 1);
    iH   = find(strcmp(phdr,'Internal Humidity') == 1);
    iBB  = find(strcmp(phdr,'Backup Batt, V') == 1);
    iVRS = find(strcmp(phdr,'Vrs') == 1);
    iSVRS = find(strcmp(phdr,'std Vrs') == 1);
    iVK  = find(strcmp(phdr,'Vk') == 1);
    iSVK = find(strcmp(phdr,'std Vk') == 1);
    iIK  = find(strcmp(phdr,'Ik') == 1);
    iIB  = find(strcmp(phdr,'Ib') == 1);
    
    %dstr  = datestr(pdata(:,iSDN),'yyyy-mm-dd HH:MM');
    dstr  = datestr(pdata(:,iSDN),'mm/dd/yyyy');
    hstr  = datestr(pdata(:,iSDN),'HH:MM');
    
    % GET LAT AND LON FROM FLOATVIZ FILES
    d = get_FloatViz_data([dirs.FV, INFO.WMO,'.TXT']);
    uC = unique(pdata(:,iC)); % merge file cycle #'s
    LL = ones(size(pdata,1),3)*NaN; % predim lon lat QF
    
    for i = 1: size(uC,1)
        t1 = d.data(:,2) == uC(i);
        tmp = d.data(t1,:);
        if isempty(tmp)
            disp(['Cycle ',num2str(uC(i)), ' exists for dura msg file ', ...
                'but not for FloatViz file - settin lon & lat to NaN'])
            disp('Are your Float Viz files up to date on you local computer?')
            LL(t2,1:2) = NaN;
            LL(t2,3)   = 1; % QF missing
        else
            t2 = pdata(:,iC) == uC(i);
            LL(t2,1) = tmp(1,3);
            LL(t2,2) = tmp(1,4);
            LL(t2,3) = 0; %QF good
        end
    end

    % PRINT HEADER FIRST
    
    
%     print_hdr ={'Cruise' 'Station' 'Type' 'yyyy-mm-dd hh:mm' ...
%         'Longitude [degrees_east]' 'Latitude [degrees_north]' 'QF'...
%         'Pressure[dbar]' 'QF' 'Temperature[°C]' 'QF' 'Salinity[pss]' 'QF'};
    
%     print_hdr ={'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm'...
%         'Lon[deg E]' 'Lat[deg N]' 'QF'...
%         'Pressure[dbar]' 'QF' 'Temperature[deg_C]' 'QF' 'Salinity[pss]' 'QF'};
    
    print_hdr ={'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm'...
        'Lon [°E]' 'Lat [°N]' 'QF'...
        'Pressure[dbar]' 'QF' 'Temperature[deg_C]' 'QF' 'Salinity[pss]' 'QF'};
    
    hdr_tmp = phdr(11:end);
    hdr_tmp = regexprep(hdr_tmp,',',''); % remove commas
    hdr_tmp = regexprep(hdr_tmp,' ','_'); % replace spaces with underscore
    
    %print_hdr = [print_hdr, phdr(6:end)];
    print_hdr = [print_hdr, hdr_tmp];
    %print_hdr = [print_hdr, 'DF_SN' 'QF' 'Stem_version' 'QF' 'K2' 'QF' 'K0' 'QF'];
    print_hdr = [print_hdr, 'DF_SN' 'QF' 'Stem_version' 'QF' 'K2' 'QF', ...
        'K0' 'QF' 'MSC_Version' 'QF'];
    
    % MODIFY HEADER - ADD UNITS TO COLUMNS
    for i = 1:size(print_hdr,2)
        X  = print_hdr{i};
        if strcmp(X,'QF') % reoces NaN's with fill value and set QF=1
            continue
        elseif regexp(X,'Power\_V$','once')
            print_hdr{i} = regexprep(X,'Power_V','Volts[V]');
        elseif regexp(X,'Batt\_V$','once')
            print_hdr{i} = regexprep(X,'_V','[V]');       
        elseif regexp(X,'Power\_I$','once')
            print_hdr{i} = regexprep(X,'Power_I','Amps[Amps]');
        elseif regexp(X,'Ik|Ib','once')
            print_hdr{i} =[X,'[Amps]'];
        elseif regexp(X,'Vk|Vrs','once')
            print_hdr{i} =[X,'[V]'];
        elseif regexp(X,'Temp$','once')
            print_hdr{i} =[X,'[deg_C]'];
        elseif regexp(X,'Humidity$','once')
            print_hdr{i} =[X,'[%]'];
        end
    end
            
    out_name = [INFO.WMO,'_DURA.TXT'];
    save_fp  = [dirs.save,'DATA\', out_name];
    disp(['Printing ph diagnostic data data to: ', save_fp]);
    fid  = fopen([dirs.save,'DATA\', out_name],'W','n','UTF-8');
    
    fprintf(fid,'//0\r\n');
    fprintf(fid,'//<Encoding>UTF-8</Encoding>\r\n');
    fprintf(fid,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
        '\r\n']);
    
    fprintf(fid,['//WMO ID: ',INFO.WMO,'\r\n']);
    fprintf(fid,['//Institution ID: ',INFO.INST_ID,'\r\n']);
    fprintf(fid,['//MBARI ID: ',MBARI_ID_str,'\r\n']);
    fprintf(fid,['//Project Name: ',INFO.Program,'\r\n']);
    fprintf(fid,['//Region: ',INFO.Region,'\r\n//\r\n']);
    
    fprintf(fid,'//STEM VERSION DEFINITIONS:\r\n');
    fprintf(fid,'//1.0 = Ag wire, thinner ISFET cover (SN 1-127)\r\n');
    fprintf(fid,'//1.1 = Ag wire, thicker ISFET cover (SN 128-163)\r\n');
    fprintf(fid,'//1.2 = Pt wire, thicker ISFET cover (SN 164 -174)\r\n');
    fprintf(fid,'//2.0 = Double O ring, Pt wire, thicker ISFET cover (SN 175-197)\r\n');
    fprintf(fid,['//3.0 = Smaller PT feed holes , double O ring, Pt wire, ', ...
        'thicker ISFET cover (SN 198-250)\r\n']);
    fprintf(fid,['//4.0 = roll pin, Smaller PT feed holes , double O ring, ', ...
        'Pt wire, thicker ISFET cover(SN 250-)\r\n']);
    fprintf(fid,'//0.0 = SBE stem\r\n//\r\n');
    
    fprintf(fid,'//QUALITY FLAG RANGE LIMITS:\r\n');
    fprintf(fid,['//Vk QF=4: [%0.3f< Vk >%0.3f], Vk QF=8: ', ...
       '[%0.3f> Vk >%0.3f]\r\n'],RC.VK4, RC.VK8);
    fprintf(fid,['//Vrs QF=4: [%0.3f< Vrs >%0.3f], Vrs QF=8: ', ...
       '[%0.3f> Vrs >%0.3f]\r\n'],RC.VRS4, RC.VRS8); 
    fprintf(fid,['//std Vk QF=4: std Vk >%0.1e, std Vk QF=8: ', ...
       'std Vk >%0.1e\r\n'],RC.SVK4, RC.SVK8); 
    fprintf(fid,['//Ik QF=4: abs(Ik) >%0.1e, Ik QF=8: ', ...
        'abs(Ik) >%0.1e\r\n'],RC.IK4, RC.IK8); 
    fprintf(fid,['//Ib QF=4: abs(Ib) >%0.1e, Ib QF=8: ', ...
        'abs(Ib) >%0.1e\r\n'],RC.IB4, RC.IB8); 
    fprintf(fid,['//System Voltage QF=4: Vsys <%0.1f, System Voltage QF=8: ', ...
        'Vsys <%0.1f\r\n//\r\n'],RC.VSYS4, RC.VSYS8); 
    
    
    fprintf(fid,['//Data quality flags: 0=Good, 4=Questionable, ', ...
        '8=Bad, 1=Missing value\r\n']);
    fprintf(fid,['//Missing data value = ',MVI_str,'\r\n//\r\n']);
    
    hdr_size  = size(print_hdr,2);
    data_rows = size(pdata,1);
    
    % PRINT HEADER & BUILD FORMAT STRING
    fstr = '';
    for i = 1:size(print_hdr,2)-1
        if regexp(print_hdr{i},'Station|QF|DF\_SN|MSC_Ver','once')
            fstr =[fstr,'%d\t'];
        elseif regexp(print_hdr{i},'Cruise|Type|^mon|^hh','once')
            fstr =[fstr,'%s\t'];
        elseif regexp(print_hdr{i},'^Vrs|^Vk','once')
            fstr =[fstr,'%0.6f\t'];
        elseif regexp(print_hdr{i},'^std|Ik|Ib|K2','once')
            fstr =[fstr,'%0.4e\t'];    
        elseif regexp(print_hdr{i},'^Pres|Humid|System_Volts','once')
            fstr =[fstr,'%0.2f\t'];
        else
            fstr =[fstr,'%0.4f\t'];
        end
        fprintf(fid,'%s\t',print_hdr{i});
    end
    fprintf(fid,'%s\r\n',print_hdr{i+1});
    %fstr =[fstr,'%d\r\n'];
    fstr =[fstr,'%d'];
    
    ifstr    = regexp(fstr,'s.+s.+s.t','once','end');
    fstr_txt = fstr(1:ifstr);
    fstr_num = fstr(ifstr+1:end);
    
    % SET DF_num to string because of SBE DF_SN can have letters
%     if ~ischar(INFO.DF_num)
%         INFO.DF_str = num2str(INFO.DF_num,'%0.0f');
%     else
%         INFO.DF_str = INFO.DF_num;
%     end
%     DF_info = [INFO.DF_ver, INFO.K2, INFO.KO];
%     DF_info(isnan(DF_info)) = -1e10; % nan's to mvi

%DF_info = [INFO.DF_num, 0, INFO.VER, 0, INFO.k2, 0, INFO.k0, 0];
DF_info = [INFO.DF_num, 0, INFO.VER, 0, INFO.k2, 0, INFO.k0, 0, MSC_ver, 0];

tnan = isnan(DF_info); % check for missing cal info
if sum(tnan,2)> 0
    %DF_info(tnan) = -1e10;  %MVI % do this later as str replacement
    DF_info(logical([0,tnan(1:end-1)])) = 1; % MVI QF
end
    
    for ct = 1 : data_rows  
        
        pos_fix = LL(ct,:);
        tnan = isnan(pos_fix); % check for missing Lat Lon
        if sum(tnan,2) > 0
            pos_fix= [-1e10, -1e10, 1]; % MVI & MVI QF
        end
        
%         fprintf(fid, fstr, MBARI_ID_str, pdata(ct,iC), 'B', ...
%             dstr(ct,:), hstr(ct,:), pos_fix, pdata(ct,iP:end), ...
%             DF_info);
        
        fprintf(fid, fstr_txt, INFO.WMO, pdata(ct,iC), 'B', ...
            dstr(ct,:), hstr(ct,:));
        
        str = sprintf(fstr_num, pos_fix, pdata(ct,iP:end), DF_info);
        fprintf(fid,'%s\r\n', regexprep(str,'NaN','-1e10'));
        

    end
    
    fclose(fid);
    
    % MAKE CONFIG FILE
    %fid  = fopen([dirs.save,'DATA\', regexprep(out_name,'TXT','CFG')],'W');
    fid  = fopen([dirs.save,'DATA\', regexprep(out_name,'TXT','CFG')],...
        'W','n','UTF-8');
    
    fprintf(fid,'//%0.0f\r\n',data_rows);
    fclose(fid);
    
end














