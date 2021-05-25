function d = Merge_isus_msgs(MBARI_ID_str, dirs)
% ************************************************************************
% PURPOSE:
%    This function processes raw isus files for a given APEX float
%    (.isus), calculates useful values from the raw signals
%    (P, T, S, O2, NO3, FIT DIAGNOSTICS) and then merges these results
%    for a given profile.
%
%
% USAGE:
%	tf_float = Process_APEX_float(UW_ID_str, dirs, update_str)
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
%   tf_float =  1 if processing was a success, 0 otherwise
%
% REVISONS:
% 12/13/20 - jp - Forced all fopen writes to UTF-8, because that is the
%     new default for Matlab 2020 and better cross platform sharing
% 12/21/20 - JP, Added header line to txt files to alert ODV that format is UTF-8, //<Encoding>UTF-8</Encoding>
% 03/23/21 - JP/TM modifications to bring in line with new WMO-based file naming system (GOBGC!)
print_flag = 1;
d = [];
% ************************************************************************
% FOR TESTING
%MBARI_ID_str  = 'ua6966';
% MBARI_ID_str  = 'ua7622';
%MBARI_ID_str  = 'ua19719';
% dirs =[];

% ************************************************************************
% SET FORMATS, DEFAULT DIRS, PREDIMENSION STRUCTURE, BOOKKEEPING
% ************************************************************************
% RCR.S    = [26 38]; % from argo parameter list
% RCR.T    = [-2.5 40]; % from argo parameter list
% RCR.NO3  = [-15 65];

fclose all; % CLOSE ANY OPEN FILES
tf_float = 0; % Default flag for float processing success 0 = no good
MVI_str = '-1e10';

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
    user_dir = [user_dir, '\Documents\MATLAB\'];
    dirs.msg   = '\\seaecho.shore.mbari.org\floats\';
    %dirs.cal   = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\';
    %dirs.cal   = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL';
    dirs.cal   = [user_dir,'ARGO_PROCESSING\DATA\CAL\'];
    
    dirs.temp  = 'C:\temp\';
    %     dirs.FV    = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
    %                   'DATA\FLOATVIZ\'];
    dirs.FV    = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\';
    
    %    dirs.save = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\', ...
    %                  'DATA\SENSOR_STATS\NO3 DIAGNOSTICS\'];
    dirs.save = [user_dir,'ARGO_PROCESSING\', ...
        'DATA\NO3_DIAGNOSTICS\'];
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end


% ************************************************************************
%LOAD CAL FILE TO CALC NO3 with later
if exist([dirs.cal,'cal',MBARI_ID_str,'.mat'],'file')
    load([dirs.cal,'cal',MBARI_ID_str,'.mat'])
else
    disp(['Could not find cal file for ',MBARI_ID_str])
    return
end

if ~isfield(cal,'N')
    disp(['No nitrate calibration file found for ',MBARI_ID_str])
    return
end

ncal            = cal.N;

INFO.INST_ID = cal.info.INST_ID; % a string
INFO.WMO     = cal.info.WMO_ID; % a string
INFO.Program = cal.info.Program; % a string
INFO.Region  = cal.info.Region; % a string
INFO.msg_dir = cal.info.msg_dir; % a string

% ************************************************************************
% BUILD *.ISUS LIST
% ************************************************************************
nlist = get_msg_list(cal.info,'isus'); % isus file list

% CHECK FOR ODD FILE NEMES & REMOVE (i.e. 19719 has "._19719.001.dura")
t1 = cellfun(@isempty,regexp(nlist.list(:,1),'^\d','once'));
if sum(t1) > 0 % non standard dura files
    nlist.list = nlist.list(~t1,:);
    disp(['Odd dura file removed from list for ',MBARI_ID_str])
end

% ************************************************************************
% COPY FLOAT MESSAGE FILES TO LOCAL TEMPORARY DIRECTORY
% Three file groups *.msg, *.isus, *.dura
% ************************************************************************
if isempty(nlist.list)
    disp(['No float dura files found to process for float ',MBARI_ID_str]);
    return
    end
    
% COPY *.ISUS files
    if ~isempty(ls([dirs.temp,'*.isus']))
        delete([dirs.temp,'*.isus']) % Clear any message files in temp dir
    end

for i = 1:size(nlist.list,1)
    fp = fullfile(nlist.list{i,2}, nlist.list{i,1});
    copyfile(fp, dirs.temp);
end

% ************************************************************************
% PROCESS MESSAGE FILES FOR GIVEN FLOAT
% ************************************************************************
msg_list   = ls([dirs.temp,'*.isus']); % get list of file names to process
disp(['Processing float ',MBARI_ID_str, '.........'])

merge_NO3_hdr = {'WMO ID' 'Cycle' 'SDN' 'Dark_Current' 'Pres' 'Temp', ...
    'Sal' 'NO3' 'BL_intercept' 'BL_slope' 'RMS_ERROR' 'WL240', ...
    'ABS_240' 'INTENSITY_240' 'MSC_Version'};
merge_NO3 = ones(300*70,size(merge_NO3_hdr,2))*NaN;
line_ct = 1;
INST_ID_num = str2double(INFO.INST_ID);

for msg_ct = 1:size(msg_list,1)
    NO3_file = strtrim(msg_list(msg_ct,:));
    % find block of numbers then look ahead to see if '.isus' follows
    cast_num = regexp(NO3_file,'\d+(?=\.isus)','once','match');
    cast_num = str2double(cast_num);
    
    % ****************************************************************
    % CALCULATE NO3 (µmol / kg scale)
    % ****************************************************************
    spec = parse_NO3msg([dirs.temp,NO3_file]); % return struct
    if ~isempty(spec.UV_INTEN)
        
        % FIX GPS TIME STAMP BUG IN 9634
        if strcmp(MBARI_ID_str,'9634SOOCN') && cast_num > 137
            fprintf('%s Cycle %0.0f GPS time bug has been corrected\n', ...
                MBARI_ID_str,  cast_num)
            spec.SDN = spec.SDN + 1024*7;
        end
        
        if strcmp(cal.info.float_type, 'APEX')
            if size(spec.pix_fit_win,2) == 2
                iINT = spec.pix_fit_win(2) - spec.spectra_pix_range(1); % intensity index ~240
            else
                str = sprintf(['%s is an APEX float but no pix_fit_win ', ...
                    'field!! Setting pix_fit_win to [217 240]'],MBARI_ID_str);
                disp(str)
                spec.pix_fit_win    = [NaN NaN];
                spec.pix_fit_win(1) = find(ncal.WL >= 217,1,'first');
                spec.pix_fit_win(2) = find(ncal.WL <= 240 ,1,'last');
                
                %                 pfit_low = find(cal.N.WL >= spec.WL_fit_win(1),1,'first');
                %                 pfit_hi  = find(cal.N.WL <= spec.WL_fit_win(2),1,'last');
                %                 spec.pix_fit_win = [pfit_low pfit_hi];
                iINT = spec.pix_fit_win(2) - spec.spectra_pix_range(1); % intensity index ~240
                clear pfit_low pfit_hi
            end
        elseif strcmp(cal.info.float_type, 'NAVIS')
            pfit_low = find(cal.N.WL >= spec.WL_fit_win(1),1,'first');
            pfit_hi  = find(cal.N.WL <= spec.WL_fit_win(2),1,'last');
            spec.pix_fit_win = [pfit_low pfit_hi];
            iINT = spec.pix_fit_win(2) - spec.spectra_pix_range(1); % intensity index ~240
            clear pfit_low pfit_hi
        end
        
        UV_INTEN = spec.UV_INTEN;
        %[SDN, DarkCur,Pres,Temp,Sal,NO3,BL_B,BL_M,RMS ERROR,WL~240,ABS~240]
        NO3 = calc_FLOAT_NO3(spec, ncal, 1); % ESW P corr, umol/L
        rN = size(NO3,1);
        merge_NO3(line_ct: rN+line_ct-1,:) = [ones(rN,1)*[INST_ID_num, cast_num], ...
            NO3, UV_INTEN(:,iINT), ones(rN,1)*spec.MSC_FW_ver];
        line_ct = line_ct + rN;
    end
end
merge_NO3 = merge_NO3(1:line_ct-1,:);
d.hdr     = merge_NO3_hdr;
d.data    = merge_NO3;

% ************************************************************************
% PRINT MERGED ISUS FILES TO TEXT
% ************************************************************************
if print_flag == 1
    % ********************************************************************
    % SET SOME QF RANGE CHECKS & ASSIGN QF's
    RC.S8   = [26 38]; % from argo parameter list
    RC.T8   = [-2.5 40]; % from argo parameter list
    RC.N8   = [-15 65]; %NO3
    RC.ABS4 = 0.8; % ABS 240nm
    RC.ABS8 = 1.1; % ABS 240nm
    RC.FIT8 = 0.003; % RMS FIT ERROR
    
    [rr,cc] = size(d.data); %predim
    pdata = ones(rr,cc*2)*NaN;
    pdata(:,1:2:2*cc-1) = d.data;
    pdata(:,2:2:2*cc) = 0;
    
    Nhdr  = cell(1,2*cc);
    Nhdr(1:2:2*cc-1) = d.hdr;
    Nhdr(2:2:2*cc) = {'QF'};
    
    tQF = strcmp(Nhdr,'QF');
    pvars = Nhdr(~tQF);
    
    clear rr cc
    
    % GET INDICES
    iSDN  = find(strcmp(Nhdr,'SDN') == 1);
    iC    = find(strcmp(Nhdr,'Cycle') == 1);
    iP    = find(strcmp(Nhdr,'Pres') == 1);
    iT    = find(strcmp(Nhdr,'Temp') == 1);
    iS    = find(strcmp(Nhdr,'Sal') == 1);
    iN    = find(strcmp(Nhdr,'NO3') == 1);
    iA240 = find(strcmp(Nhdr,'ABS_240') == 1);
    iFIT  = find(strcmp(Nhdr,'RMS_ERROR') == 1);
    iDC   = find(strcmp(Nhdr,'Dark_Current') == 1);
    
    % merge_NO3_hdr = {'UW ID' 'Cycle' 'SDN' 'Dark_Current' 'Pres' 'Temp', ...
    %     'Sal' 'NO3' 'BL_intercept' 'BL_slope' 'RMS_ERROR' 'WL240', ...
    %     'ABS_240' 'INTENSITY_240'};
    
    % SET SAME QUALITY FLAGS
    for i = 1:size(pvars,2)
        X  = pvars{1,i};
        iX = find(strcmp(Nhdr,X) == 1);
        tnan = isnan(pdata(:,iX));
        if ~strcmp(X,'QF') % replaces NaN's with fill value and set QF=1
            %pdata(tnan,iX) = -1e10; % will put fill value in later as str
            pdata(tnan,iX+1) = 1;
        end
        
        if strcmp(X,'QF')
            continue
        elseif strcmp(X,'Temp')
            t8 = pdata(:,iX) < RC.T8(1) | pdata(:,iX) > RC.T8(2);
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'Sal')
            t8 = pdata(:,iX) < RC.S8(1) | pdata(:,iX) > RC.S8(2);
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'NO3')
            %t4 = pdata(:,iA240) > RC.ABS4;
            t8 = pdata(:,iX) < RC.N8(1) | pdata(:,iX) > RC.N8(2) | ...
                pdata(:,iS+1) == 4 | pdata(:,iA240) > RC.ABS8 | ...
                pdata(:,iFIT) > RC.FIT8;
            pdata(~tnan,iX+1) = 4; % set all NO3 as questionable
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'ABS_240')
            t4 = pdata(:,iX) > RC.ABS4;
            t8 = pdata(:,iX) > RC.ABS8;
            pdata(t4&~tnan,iX+1) = 4;
            pdata(t8&~tnan,iX+1) = 8;
        elseif strcmp(X,'RMS_ERROR')
            t8 = pdata(:,iX) > RC.FIT8;
            pdata(t8&~tnan,iX+1) = 8;
        end
    end
    
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
    
    % merge_NO3_hdr = {'UW ID' 'Cycle' 'SDN' 'Dark_Current' 'Pres' 'Temp', ...
    %     'Sal' 'NO3' 'BL_intercept' 'BL_slope' 'RMS_ERROR' 'WL240', ...
    %     'ABS_240' 'INTENSITY_240'};
    
    % PRINT HEADER FIRST
    out_name = [INFO.WMO,'_NO3.TXT'];
    
    
    print_hdr ={'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm'...
        'Lon [°E]' 'Lat [°N]' 'QF'...
        'Pressure[dbar]' 'QF' 'Temperature[deg_C]' 'QF' 'Salinity[pss]' 'QF'};
    
    hdr_tmp = Nhdr;
    hdr_tmp = regexprep(hdr_tmp,' ','_'); % replace spaces with underscore
    print_hdr = [print_hdr, hdr_tmp(iDC:iDC+1),hdr_tmp(iN:end)];
    
    % MODIFY HEADER - ADD UNITS TO COLUMNS
    for i = 1:size(print_hdr,2)
        X  = print_hdr{i};
        if strcmp(X,'QF')
            continue
        elseif regexp(X,'NO3','once')
            print_hdr{i} = regexprep(X,'NO3','Nitrate[µmol/L]');
        end
    end
    
    
    disp(['Printing nitrate diagnostic data data to: ',dirs.save,out_name]);
    %fid  = fopen([dirs.save,'DATA\', out_name],'W');
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
    fprintf(fid,'//WARNING:\r\n');
    fprintf(fid,'//  NITRATE CONCENTRATION IN THIS FILE IS NOT CORRECTED!\r\n');
    fprintf(fid,['//  Users of this data file should have a good ',...
        'understanding how the\r\n//  raw seawater absorption spectra ',...
        'are used to calculate nitrate\r\n']);
    fprintf(fid,'//\r\n');
    fprintf(fid,'//QUALITY FLAG RANGE LIMITS:\r\n');
    fprintf(fid,'//  Temperature QF=8: %0.0f< T(C) >%0.0f\r\n', RC.T8);
    fprintf(fid,'//  Salinity QF=8: %0.0f< T(C) >%0.0f\r\n', RC.S8);
    fprintf(fid,'//  Nitrate QF=4: %0.1f< ABS240 >%0.1f\r\n',RC.ABS4, RC.ABS8);
    fprintf(fid,['//  Nitrate QF=8: [%0.0f< NO3 >%0.0f] | SALINITY QF = 8 | ', ...
        'ABS240 >%0.1f | Fit error >%0.3f  \r\n'],RC.N8, RC.ABS8,RC.FIT8);
    fprintf(fid,['//  ABS240 QF=4: [%0.1f< ABS240 >%0.1f], ABS240 QF=8: ', ...
        '[ABS240 >%0.1f]\r\n'],RC.ABS4, RC.ABS8,RC.ABS8);
    fprintf(fid,'//  RMS Fit error QF=8: RMS ERROR >%0.3f\r\n', RC.FIT8);
    
    fprintf(fid,['//Data quality flags: 0=Good, 4=Questionable, ', ...
        '8=Bad, 1=Missing value\r\n']);
    fprintf(fid,['//Missing data value = ',MVI_str,'\r\n//\r\n']);
    
    hdr_size  = size(print_hdr,2);
    data_rows = size(pdata,1);
    
    % PRINT HEADER & BUILD FORMAT STRING
    fstr = '';
    for i = 1:size(print_hdr,2)-1
        if regexp(print_hdr{i},'Station|QF|^INTENSITY|MSC_Ver','once')
            fstr =[fstr,'%d\t'];
        elseif regexp(print_hdr{i},'Cruise|Type|^mon|^hh','once')
            fstr =[fstr,'%s\t'];
        elseif regexp(print_hdr{i},'^Lon|^Lat|','once')
            fstr =[fstr,'%0.4f\t'];
        elseif regexp(print_hdr{i},'^Pressure|^Nitrate','once')
            fstr =[fstr,'%0.2f\t'];
        elseif regexp(print_hdr{i},'^Dark','once')
            fstr =[fstr,'%0.0f\t'];
        else
            fstr =[fstr,'%0.4g\t'];
        end
        fprintf(fid,'%s\t',print_hdr{i});
    end
    fprintf(fid,'%s\r\n',print_hdr{i+1});
    %fstr =[fstr,'%d\r\n'];
    fstr =[fstr,'%d'];
    
    % Find start of number only formating
    ifstr    = regexp(fstr,'s.+s.+s.t','once','end');
    fstr_txt = fstr(1:ifstr);
    fstr_num = fstr(ifstr+1:end);
    
    % PRINT DATA LINES
    for ct = 1 : data_rows
        pos_fix = LL(ct,:);
        tnan = isnan(pos_fix); % check for missing Lat Lon
        if sum(tnan,2) > 0
            pos_fix= [-1e10, -1e10, 1]; % MVI & MVI QF
        end
        
        % ADD SOME META DATA
        fprintf(fid, fstr_txt, INFO.WMO, pdata(ct,iC), 'B', ...
            dstr(ct,:), hstr(ct,:));
        
        % BUILD DATA LINE AS STRING
        str = sprintf(fstr_num, pos_fix, pdata(ct,iP:iP+5), ... %position & PTS
            pdata(ct,iDC:iDC+1),pdata(ct,iN:end));
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














