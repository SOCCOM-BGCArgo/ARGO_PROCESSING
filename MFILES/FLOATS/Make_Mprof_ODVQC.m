function tf = Make_Mprof_ODVQC(handles)
% FOR USE IN SAGE AND SAGE-O2Argo
% CREATE AN ADJUSTED ODV STYLE TEXT FILE FROM RAW ODV FILE GIVEN A
% CORRECTIONS FILE. THIS IS DESIGNED ONLY FOR MPROF FILES WHERE NO MESSAGE
% FILES EXIST & NO ARGO FILES.

fp = filesep; % File separator for current platform
tf = 0;
% CHECK TO MAKE SURE THIS IS AN MPROF FILE
if handles.info.Mprof == 0
    disp('File does not appear to originated from a Mprof NetCDF')
    return
end

% DEFINE RANGE CHECK TABLE
RC.S     = [26 38]; % from argo parameter list
RC.T     = [-2.5 40]; % from argo parameter list
RC.O     = [-5 550]; % from argo parameter list
RC.CHL   = [-0.1 150]; % argoBGC QC manual 09July2016
RC.BB700 = [-0.000025 0.1]; % argoBGC QC manual 09July2016
RC.BB532 = [-0.000005 0.1]; % argoBGC QC manual 09July2016
RC.CDOM  = [-1000 1000]; % Place Holder
RC.NO3   = [-15 65];
RCR.PH   = [7.0 8.8];

% GET THE RAW DATA
fv_path  = [handles.info.file_path, handles.info.float_name,'.TXT'];
d        = get_FloatViz_data(fv_path);
raw_hdr      = d.hdr;
raw_data = d.data;
qc_data  = d.data; % start with raw as QC template
qc_hdr = raw_hdr;
clear d

iSDN  = find(strcmp('SDN',qc_hdr)                == 1);
iSTN  = find(strcmp('Station',qc_hdr)            == 1);
iP    = find(strcmp('Pressure[dbar]',qc_hdr)     == 1);
iT    = find(strcmp('Temperature[°C]',qc_hdr)    == 1);
iS    = find(strcmp('Salinity[pss]',qc_hdr)      == 1);
iZ    = find(strcmp('Depth[m]',qc_hdr)           == 1);
iO    = find(strcmp('Oxygen[µmol/kg]',qc_hdr)    == 1);
iN    = find(strcmp('Nitrate[µmol/kg]',qc_hdr)   == 1);
iPH   = find(strcmp('pHinsitu[Total]',qc_hdr)    == 1);
iCHL  = find(strcmp('Chl_a[mg/m^3]',qc_hdr)      == 1);
i380  = find(strcmp('D_IRRAD380[W/m^2/nm]',qc_hdr) == 1);
i412  = find(strcmp('D_IRRAD412[W/m^2/nm]',qc_hdr) == 1);
i490  = find(strcmp('D_IRRAD490[W/m^2/nm]',qc_hdr) == 1);
iPAR  = find(strcmp('D_PAR[W/m^2/nm]',qc_hdr)    == 1);
iHS   = find(strcmp('Bisulfide[µmol/kg]',qc_hdr) == 1);


% SET MASTER LIST OF POSSIBLE CORRECTABLE PARAMETERS & TEST FOR DATA
qc_params = {'Oxygen[µmol/kg]', 'Nitrate[µmol/kg]','pHinsitu[Total]', ...
             'Chl_a[mg/m^3]'};
    
% GET CORRECTIONS FROM GUI STRUCTURE
QCA       = handles.QCA;
% qca_params = fieldnames(QCA);
for qc_ct = 1: size(qc_params,2)
    switch(qc_params{qc_ct})
        
        case 'Oxygen[µmol/kg]'
            if ~isfield(QCA,'O2') 
                disp('No O2 gain correction found.')
                if isfield(QCA,'O2')
                    QCA = rmfield(QCA,'O2');
                end
                continue
            end
            if isempty(iO) || isempty(QCA.O2)
                disp('No O2 index found...skipping O2 QC application.')
                if isfield(QCA,'O2')
                    QCA = rmfield(QCA,'O2');
                end
                continue
            end
            tdata = ~isnan(qc_data(:,iO)) & qc_data(:,iO) ~= -1e10;
            if sum(tdata) > 0 % data to correct
                qca = QCA.O2;
                Guse = [];
                for gi = 1:size(qca,1) %perform operation for each row of O2 gains (usually only one single drift, so up to 2 rows, but allow for more complex)          
                    data_tmp = qc_data(tdata,:);
                    if gi == size(qca,1) %either single gain, or last row in qca gain matrix
                        data_tmp = data_tmp(data_tmp(:,iSTN)>= qca(gi,1) & data_tmp(:,iSTN)<=nanmax(data_tmp(:,iSTN)),:);
                    else
                        data_tmp = data_tmp(data_tmp(:,iSTN)>= qca(gi,1) & data_tmp(:,iSTN)<qca(gi+1,1),:);
                    end
                    Tdiff = data_tmp(:,iSDN) - nanmin(data_tmp(:,iSDN));
                    Gtmp = qca(gi,2) + qca(gi,end)./365.*Tdiff;
                    Guse = [Guse;Gtmp]; %will result in a gain value to apply per sample.
                end
%                 qc_data(tdata,iO) = qc_data(tdata,iO) * qca(1,2); %if only static mean gain ever applied.
                qc_data(tdata,iO) = qc_data(tdata,iO) .* Guse;
                qc_data(tdata,iO+2) = qc_data(tdata,iO) ./ ...
                   oxy_sol(qc_data(tdata,iT), qc_data(tdata,iS),0) * 100;
                tRC = qc_data(:,iO) >= RC.O(1) & qc_data(:,iO) <= RC.O(2);
                tbad = qc_data(:,iO+1) == 8;
                qc_data(tdata & tRC & ~tbad, iO+1) = 0; % GOOD DATA!
                qc_data((tdata & ~tRC) | tbad, iO+1) = 8; % BAD - OUT OF RANGE DATA!
                qc_data(:,iO+3) = qc_data(:,iO+1); % o2 sat qc = O2 qc
            end
            clear tdata qca tRC
            
        case 'Nitrate[µmol/kg]'
            if ~isfield(QCA,'NO3')
                disp('No Nitrate QC corrections found.')
                if isfield(QCA,'NO3')
                    QCA = rmfield(QCA,'NO3');
                end                
                continue
            end
            if isempty(iN) || isempty(QCA.NO3)
                disp('No Nitrate index exists...skipping Nitrate QC application.')
                if isfield(QCA,'NO3')
                    QCA = rmfield(QCA,'NO3');
                end                
                continue
            end
            tdata = ~isnan(qc_data(:,iN)) & qc_data(:,iN) ~= -1e10;
            if sum(tdata) > 0 % data to correct
                qc_tmp = qc_data(tdata,:);
                cor = ones(size(qc_tmp(:,1)))*0; % predim w/ zeros
                qca = QCA.NO3; % corections
                rr  = size(qca,1); % # adjustment rows
                last_cor = 0;
                for i = 1:rr
                    if i < rr
                        t1 = qc_tmp(:,2) >= qca(i,1) & ...
                             qc_tmp(:,2) <= qca(i+1,1); %get data block
                    else
                        t1 = qc_tmp(:,2) >= qca(i,1); %get data block
                    end
                    % elapsed years in block
                    dt = (qc_tmp(t1,1) - min(qc_tmp(t1,1))) ./ 365;
                    tmp_cor  = last_cor + qca(i,3) + qca(i,4) .* dt; %offset + drift
                    last_cor = tmp_cor(end);
                    cor(t1)  = tmp_cor;
                end
                qc_tmp(:,iN)   = (qc_tmp(:,iN) - cor) ./ qca(i,2);
                
                % SET BAD QC DATA TO NaN
                tRC  = qc_tmp(:,iN) >= RC.NO3(1) & ...
                       qc_tmp(:,iN) <= RC.NO3(2);
                tbad = qc_tmp(:,iN+1) == 8; % carry over from raw
                qc_tmp(tRC & ~tbad, iN+1) = 0; % GOOD DATA!
                qc_tmp((~tRC) | tbad, iN+1) = 8; % BAD!
                qc_data(tdata,:) = qc_tmp;
            end
            clear tdata cor qca rr last_cor i t1 dt tmp_cor  tRC tbad qc_tmp
            
            
        case 'pHinsitu[Total]'
            if ~isfield(QCA,'PH')
                disp('No pH QC correctionc found.')
                if isfield(QCA,'PH')
                    QCA = rmfield(QCA,'PH');
                end
                continue
            end
            if isempty(iPH)  || isempty(QCA.PH)
                disp('No pH QC index exists...skipping pH QC application.')
                if isfield(QCA,'PH')
                    QCA = rmfield(QCA,'PH');
                end
                continue
            end
            tdata = ~isnan(qc_data(:,iPH)) & qc_data(:,iPH) ~= -1e10;
            if sum(tdata) > 0 % data to correct
                qc_tmp = qc_data(tdata,:);
                cor = ones(size(qc_tmp(:,1)))*0; % predim w/ zeros
                qca = QCA.PH; % corections
                rr  = size(qca,1); % # adjustment rows
                
                if isempty(QCA.PH_OFFSET)
                    pump_offset = 0;
                else
                    pump_offset = QCA.PH_OFFSET;
                end
                tP          = qc_tmp(:,iP) < 980; % logical flag for pump on calc
                pump_offset = pump_offset * tP; % array of 0's and offset
                TCOR        = (2 + 273.15)./(qc_tmp(:,iT) + 273.15); % drift + offset is f(T) too????
                
                last_cor = 0;
                for i = 1:rr
                    if i < rr
                        t1 = qc_tmp(:,2) >= qca(i,1) &  ...
                             qc_tmp(:,2) <= qca(i+1,1); %get block
                    else
                        t1 = qc_tmp(:,2) >= qca(i,1); %get block
                    end
                    dt = (qc_tmp(t1,1) - min(qc_tmp(t1,1))) ./ 365; % elapsed years in block
                    tmp_cor  = last_cor + qca(i,3) + qca(i,4) .* dt; %offset + drift
                    last_cor = tmp_cor(end);
                    cor(t1)  = tmp_cor;
                end
                qc_tmp(:,iPH) = qc_tmp(:,iPH) +(pump_offset - cor) ...
                                 .* TCOR;
                
                % SET BAD QC DATA TO NaN
                tRC  = qc_tmp(:,iPH) >= RC.PH(1) & ...
                       qc_tmp(:,iPH) <= RC.PH(2);
                tbad = qc_tmp(:,iPH+1) == 8;
                qc_tmp(tRC & ~tbad, iPH+1) = 0; % GOOD DATA!
                qc_tmp((~tRC) | tbad, iPH+1) = 8; % BAD!
                qc_data(tdata,:) = qc_tmp;
            end
            clear tdata cor qca rr last_cor i t1 dt tmp_cor  tRC tbad qc_tmp
            
        case 'Chl_a[mg/m^3]' % GAIN CORRECTION ONLY
            if isempty(iCHL) || ~isfield(QCA,'CHL')
                disp('No CHL index or QC corrections found')
                if isfield(QCA,'CHL')
                    QCA = rmfield(QCA,'CHL');
                end                
                continue
            end
            tdata = ~isnan(qc_data(:,iCHL)) & qc_data(:,iCHL) ~= -1e10;
            if sum(tdata) > 0 % data to correct
                qca = QCA.NO3; % corections
                qc_data(tdata,iCHL) = qc_data(tdata,iCHL) * qca(1,2);
                tRC = qc_data(:,iCHL) >= RC.CHL(1) & ...
                    qc_data(:,iCHL) <= RC.CHL(2);
                tbad = qc_data(:,iCHL+1) == 8; % carry over from raw
                qc_data(tdata & tRC & ~tbad, iCHL+1) = 0; % GOOD DATA!
                qc_data((tdata & ~tRC) | tbad, iCHL+1) = 8; % BAD - OUT OF RANGE DATA!
            end
            clear tdata qca tRC tbad
    end
end

[qc_r, qc_c] = size(qc_data);
% ************************************************************************
% ************************************************************************
%                QC ADJUSTMENTS MADE - NOW PRINT TO FILE
% ************************************************************************
% ************************************************************************
MVI_str = '-1e10'; % MISSING VALUE INDICATOR FOR ODV

% ************************************************************************
% BUILD LOOK UP CELL ARRAY to match variables and set format string
% ************************************************************************
%RAW ODV FILE
ODV_adj(1,:)  = {'Pressure[dbar]'        '%0.2f' 'PRES' '' '' ''}; % ?
ODV_adj(2,:)  = {'Temperature[°C]'       '%0.4f' 'TEMP' '' '' ''};   
ODV_adj(3,:)  = {'Salinity[pss]'         '%0.4f' 'PSAL' '' '' ''};   
ODV_adj(4,:)  = {'Sigma_theta[kg/m^3]'   '%0.3f' 'SIGMA_THETA' '' '' ''};
ODV_adj(5,:)  = {'Depth[m]'              '%0.3f' 'DEPTH' '' '' ''};
ODV_adj(6,:)  = {'Oxygen[µmol/kg]'       '%0.3f' 'DOXY' '' '' ''};   
ODV_adj(7,:)  = {'OxygenSat[%]'          '%0.3f' 'DOXY_%SAT' '' '' ''};
ODV_adj(8,:)  = {'Nitrate[µmol/kg]'      '%0.2f' 'NITRATE' '' '' ''}; 
ODV_adj(9,:)  = {'Chl_a[mg/m^3]'         '%0.4f' 'CHLA' '' '' ''};   
ODV_adj(10,:) = {'b_bp700[1/m]'          '%0.6f' 'BBP700' '' '' ''};
ODV_adj(11,:) = {'CDOM[ppb]'             '%0.2f' 'CDOM' '' '' ''}; 
ODV_adj(12,:) = {'CP660[1/m]'            '%0.4f' 'CP660' '' '' ''}; 

% ADD THESE FOR ODV FLAVOR #2 -PROVOR
ODV_adj(13,:) = {'D_IRRAD380[W/m^2/nm]'  '%4.4f' 'DOWN_IRRADIANCE380' '' '' ''}; 
ODV_adj(14,:) = {'D_IRRAD412[W/m^2/nm]'  '%4.4f' 'DOWN_IRRADIANCE412' '' '' ''}; 
ODV_adj(15,:) = {'D_IRRAD490[W/m^2/nm]'  '%4.4f' 'DOWN_IRRADIANCE490' '' '' ''}; 
ODV_adj(16,:) = {'D_PAR[W/m^2/nm]'       '%4.4f' 'DOWNWELLING_PAR' '' '' ''}; 
ODV_adj(17,:) = {'Bisulfide[µmol/kg]'    '%4.4f' 'BISULFIDE' '' '' ''}; 

% ************************************************************************

if isempty(i380) && isempty(i412) && isempty(i490) && isempty(iPAR)
    ind1 = find(strcmp('D_IRRAD380[W/m^2/nm]',ODV_adj(:,1))   == 1);
    ind2 = find(strcmp('D_IRRAD412[W/m^2/nm]',ODV_adj(:,1))   == 1);
    ind3 = find(strcmp('D_IRRAD490[W/m^2/nm]',ODV_adj(:,1))   == 1);
    ind4 = find(strcmp('D_PAR[W/m^2/nm]',ODV_adj(:,1))   == 1);
    ODV_adj([ind1,ind2,ind3,ind4],:) = [];
end

if isempty(iHS)
    ind1 = find(strcmp('Bisulfide[µmol/kg]',ODV_adj(:,1))   == 1);
    ODV_adj(ind1,:) = [];
end

qc_var_ct = size(ODV_adj,1);

if ~isdir([handles.info.file_path, 'QC',fp])
    mkdir([handles.info.file_path, 'QC',fp]);
end
    
fvqc_path  = [handles.info.file_path, 'QC',fp,handles.info.float_name, ...
              'QC.TXT'];
fid_qc  = fopen(fvqc_path, 'W');
fid_raw = fopen(fv_path);

disp(['Printing adjusted data to: ',fvqc_path]);

fprintf(fid_qc,'//0\r\n');
fprintf(fid_qc,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
       '\r\n']);
 
 tline    = ' ';  
 print_it = 0;
 while ischar(tline)
     if regexp(tline,'^//File','once')
         print_it = 1;
     elseif regexp(tline,'^//Missing','once')
         break
     elseif print_it == 1
         fprintf(fid_qc, '%s\r\n', tline);
         if regexp(tline,'^//WMO', 'once')
             WMO = regexp(tline,'\d+', 'match','once');
         end
     end
     tline = fgetl(fid_raw);
 end
 fclose(fid_raw); % Don't need raw file any more
 
 % PRINT OUT FLOAT VARIABLE QC CORRECTION INFO
 fprintf(fid_qc,'//QUALITY CONTROLLED DATA CORRECTIONS:\r\n');
 fprintf(fid_qc,'//Measurement\tStation\tGain\tOffset\tDrift\r\n');
 possible_fields = {'O2' 'NO3' 'PH' 'CHL' 'BB' 'CDOM'};
 possible_hdr_names = {'Oxygen' 'Nitrate' 'pH' 'Chl' 'BB' 'CDOM'};
 for i = 1 : size(possible_fields,2)
     if isfield(QCA, possible_fields{i}) %
         qca = QCA.(possible_fields{i});
         rr  = size(qca,1); % # adjustment rows
         for j = 1: rr
             fprintf(fid_qc, ['//',possible_hdr_names{i},'\t', ...
                 '%0.0f\t%0.4f\t%0.4f\t%0.4f\r\n'], qca(j,1:4));
         end
     end
 end
 fprintf(fid_qc,'//\r\n');
 
 
 fprintf(fid_qc,['//Missing data value = ',MVI_str,'\r\n']);
 fprintf(fid_qc,['//Data quality flags: 0=Good, 4=Questionable, 8=Bad, '...
     '1=Missing or not inspected \r\n']);
 
% NOW PRINT THE QC DATA HEADER
std_ODV_vars   = {'Cruise' 'Station' 'Type' 'mon/day/yr' 'hh:mm' ...
                  'Lon [°E]' 'Lat [°N]' 'QF'}; % SIZE = 8
std_size = size(std_ODV_vars,2);

for i = 1:std_size % PRINT STANDARD HEADER VARS
    fprintf(fid_qc,'%s\t',std_ODV_vars{1,i}); % std vars
end

for i = 1:qc_var_ct % PRINT FLOAT SPECIFIC HEADER VARS
    if i < qc_var_ct
        fprintf(fid_qc,'%s\t%s\t',ODV_adj{i,1},'QF'); % std vars
    else
        fprintf(fid_qc,'%s\t%s\r\n',ODV_adj{i,1},'QF'); % std vars
    end
end
 
% BUILD DATA FORMAT STRING & CREATE DUMMY MATRIX TO DUMP DATA TO
% THIS CATCHES FLOATS WHICH DON'T HAVE SENSORS FOR GIVEN VARIABLES
% THESE COLUMNS WILL GET FILLED WITH -1e10 and QC FLAG = 1
dummy_out  = ones(qc_r, qc_var_ct) * NaN;
fill_MVI = ones(qc_r, 1) * -1e10;
fill_QC  = ones(qc_r, 1);
ODV_std_f = '%s\t%0.0f\t%s\t%s\t%s\t%0.3f\t%0.3f\t%0.0f\t'; %std_vars
ODV_qc_f =''; ODV_space_f = '';
for i = 0:qc_var_ct-1
    c_ct = i*2+1; % Need to add data and QC col
    ind = find(strcmp(ODV_adj{i+1,1}, qc_hdr) == 1);
    
    if ~isempty(ind)
        if i < qc_var_ct-1
            ODV_qc_f   = [ODV_qc_f, ODV_adj{i+1,2},'\t%0.0f\t'];
            ODV_space_f = [ODV_space_f,'\t\t'];
        else
            ODV_qc_f   = [ODV_qc_f, ODV_adj{i+1,2},'\t%0.0f\r\n'];
            ODV_space_f = [ODV_space_f,'\t\r\n'];
        end
        dummy_out(:,c_ct:c_ct+1)  = qc_data(:,ind:ind+1);
    else
        if i < qc_var_ct-1
            ODV_qc_f   = [ODV_qc_f, '%1.0E\t%0.0f\t'];
            ODV_space_f = [ODV_space_f,'\t\t'];
        else
            ODV_qc_f   = [ODV_qc_f, '%1.0E\t%0.0f\r\n'];
            ODV_space_f = [ODV_space_f,'\t\r\n'];
        end
        dummy_out(:,c_ct:c_ct+1)  = [fill_MVI,fill_QC]; % No data, but need cols
    end
end 

% NOW PRINT DATA LINES TO FILE
cast_num  = 0; %initalize
prof_num  = 0;
line_ct   = 0;
for sample_ct = 1 : qc_r
    if qc_data(sample_ct,2) - cast_num > 0 % build standard cast part
        if sample_ct > 1
            fprintf(fid_qc,[std_str,ODV_space_f]); % add profile spacer line
            line_ct = line_ct+1;
        end
        cast_num = qc_data(sample_ct,2);
        date_str = datestr(qc_data(sample_ct,1),'mm/dd/yyyy');
        time_str = datestr(qc_data(sample_ct,1),'HH:MM');
        std_str  = sprintf(ODV_std_f, WMO, cast_num, 'C', ...
            date_str, time_str, qc_data(sample_ct,3:5));
        std_str = regexprep(std_str,'NaN',MVI_str);
    end   
    
    Dout = dummy_out(sample_ct,:);
    Dout(Dout<-100000000) = nan;
    data_str = sprintf(ODV_qc_f, Dout);
    out_str = [std_str,data_str];

    % replace NaN' w/ missing value indicator
%     out_str = regexprep(out_str,'NaN',MVI_str);  % replace NaN' w/ missing value indicator
    out_str = regexprep(out_str,'NaN',MVI_str);  % replace -10000000000 w/ missing value indicator
    fprintf(fid_qc, '%s', out_str);
    line_ct = line_ct+1;
    clear out_str data_str date_str time_str
    if sample_ct == qc_r
        fprintf(fid_qc,[std_str,ODV_space_f]); % add profile spacer line
        line_ct = line_ct+1;
    end    
end
fclose(fid_qc);
clear fid_raw cast_num sample_ct

% MAKE CONFIG FILE FOR FLOATVIZ
fvqc_path  = [handles.info.file_path, 'QC',fp,handles.info.float_name, ...
              'QC.CFG'];
fid_qc  = fopen(fvqc_path,'w');
fprintf(fid_qc,'//%0.0f\r\n',line_ct);
fclose(fid_qc);
disp(['DONE printing adjusted data to: ',fvqc_path]);
clear fid_raw line_ct dummy_out


            


