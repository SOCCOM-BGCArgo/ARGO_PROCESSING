function tf = NewFloatQCList_GLT(handles,dirs)
% This function updates the FloatQCList. If float doesn't exist on list it
% is added. If sensor exists, but there is no QC line for the sensor, a 
% default line is added that doesn't change the QC data value for that
% sensor
%
%UPDATES:
% 08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
% 10/30/2017 - Addapted to use text files generated from Mprof.nc files & not to print QC line if sensor 
%              doesn't exist     


tf = 1;
% chl_flag = 0;
% new_float = 0;
% Make new  FloatQcList.txt locally from float_qc GUI data table changes
fname        = [dirs.QCadj,handles.info.QCadj_file]; % Existing adjustment file
new_fname    = [dirs.temp,filesep,'new_FloatQCList.txt'];  %file that will be updated
float_name   = handles.info.float_name;

if exist(fname) == 0 %does not exist (ie new float), create file
    fid = fopen(fname,'w+');
    wasempty = 1;
else
    fid        = fopen(fname); % Existing list
    wasempty = 0;
end

fid_new    = fopen(new_fname,'a'); % new list

QCA = handles.QCA;
field_names = fieldnames(QCA);
%cal = handles.info.cal; %used to check sensors on float

fprintf(fid_new,['%s\t',datestr(now,'mm/dd/yy HH:MM'),'\r\n'],...
                 float_name); % Update date on existing float line
if isfield(QCA,'O2') && handles.info.O2_sensor == 1
    rr = size(QCA.O2,1);
    if rr == 0 % no adjustments yet
        fprintf(fid_new,'Oxygen,  %0.0f,%1.4f,%1.4f,%1.4f\r\n',[1, 1, 0, 0]);
    else
        for i = 1 : rr
            % cycle gain offset drift (offset will be zero)
            fprintf(fid_new,'Oxygen, %0.0f,%1.4f,%1.4f,%1.4f\r\n',QCA.O2(i,:));
        end
    end
end

if isfield(QCA,'NO3') && handles.info.NO3_sensor == 1
    rr = size(QCA.NO3,1);
    if rr == 0 % no adjustments yet
        fprintf(fid_new,'Nitrate, %0.0f,%1.4f,%1.4f,%1.4f\r\n',1,1,0,0);
    else
        for i = 1 : rr
%             if i == 1
%                 fprintf(fid_new,'Nitrate, %0.0f,%1.3f,%1.3f,%1.3f,%s\r\n',...
%                     QCA.NO3(i,:),handles.info.comp_data);
%             else
            fprintf(fid_new,'Nitrate, %0.0f,%1.4f,%1.4f,%1.4f\r\n', ...
                QCA.NO3(i,:));         
        end
    end
end

if isfield(QCA,'PH_OFFSET') && handles.info.PH_sensor == 1
    rr = size(QCA.PH_OFFSET,1);
    if rr == 0
        fprintf(fid_new,'pH, offset, %1.4f\r\n',0);
    else
        for i = 1 : rr
            fprintf(fid_new,'pH, offset, %1.4f\r\n',QCA.PH_OFFSET(i,1));
        end
    end
end

if isfield(QCA,'PH') && handles.info.PH_sensor == 1
    rr = size(QCA.PH,1);
    if rr == 0 % no adjustments yet
        fprintf(fid_new,'pH, %0.0f,%1.4f,%1.4f\r\n',1,0,0);
    else
        for i = 1 : rr
            % CYCLE OFFSET DRIFT
            fprintf(fid_new,'pH, %0.0f,%1.4f,%1.4f\r\n',QCA.PH(i,[1,3,4]));
        end
    end
end

% BIO-OPTICS NOT STORED IN QCA, ONLY O, N, and pH
            % no QC CHL line, but CHL exists
if handles.info.Mprof == 0 %&& isfield(handles.info.cal,'CHL') %TM 12/10/20; don't populate QClist with CHL anymore.  May change in future?
%    fprintf(fid_new,'%s\r\n', 'CHL, 0.5, 0'); %TM 12/10/20; don't populate QClist with CHL anymore.  May change in future?
% MPROF FILE - NO ASSOCIATED CAL FILE TO CHECK    
elseif handles.info.Mprof == 1 && handles.info.CHL_sensor == 1
    fprintf(fid_new,'%s\r\n', 'CHL, 0.5, 0');
end

if wasempty == 0 % not a new float; append previous
    fprintf(fid_new,'\r\n')
    fprintf(fid_new,'PREVIOUS\r\n')

    % Now copy all previous QC corrections to new file, latest correction at
    % top
    if wasempty == 1 %new float, floatQC file was empty 
        %do nothing
    else
        while 1
            mytline = fgetl(fid);
            if ~ischar(mytline), break, end
            fprintf(fid_new, '%s\r\n',mytline);
        end
    end
else
    % new float; do nothing more
end

fclose(fid);
fclose(fid_new);


% NOW SAVE OLD AND RENAME NEW to default name
% copyfile([handles.dirs.QCadj,fname],[handles.dirs.QCadj,'old_',fname]);
copyfile(new_fname,fname);
delete(new_fname);

