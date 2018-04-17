function tf_float = Merge_dura_msgs(UW_ID_str, dirs)
% ************************************************************************
% PURPOSE:
%    This function processes raw message files for a given APEX float
%    (.msg, .isus, .dura), calculates useful values from the raw signals
%    (P, T, S, O2, NO3, pH, CHl, Bbp, CDOM) and then merges these results
%    for a given profile.
%
%    Each profile is saved as a *.mat file with the following format:
%    WMO_ID.PROFILE_#.mat (Ex 5904683.006.mat) in a directory using the
%    WMO_ID for its name. Each file contains 3 structures with variables
%    inside named using ARGO definitions:
%       LR   - Low resolution data
%       HR   - High resolution data
%       INFO - Some information about the profile and float
%
% USAGE:
%	tf_float = Process_APEX_float(UW_ID_str, dirs, update_str)
%
% INPUTS:
%   UW_ID_str  = UW/MBARI Float ID # as a string
%
%   update_str = "all" or "update"
%                 all    - to process all available msg files
%                 update - only process new msg files
%
%   dirs       = Either an empty variable ora structure with directory
%                strings where files are located. It must contain these
%                fields at a minimum or be empty. If dirs is empty default
%                paths will be used.
%
% OUTPUTS:
%   tf_float =  1 if processing was a success, 0 otherwise
%

% ************************************************************************
% FOR TESTING
% UW_ID_str  = '9095'; % very low intensites
% UW_ID_str  = '9125';
% dirs =[];

% ************************************************************************
% SET FORMATS, DEFAULT DIRS, PREDIMENSION STRUCTURE, BOOKKEEPING
% ************************************************************************
fclose all; % CLOSE ANY OPEN FILES
tf_float = 0; % Default flag for float processing success 0 = no good
MVI_str = '-1e10';

% **** DEFAULT STRUCTURE FOR DIRECTORY PATHS ****
if isempty(dirs)
    dirs.msg       = '\\atlas\ChemWebData\floats\';
    dirs.mat       = 'C:\Users\jplant\Documents\MATLAB\ARGO\DURA\';
    dirs.temp      = 'C:\temp\';
elseif ~isstruct(dirs)
    disp('Check "dirs" input. Must be an empty variable or a structure')
    return
end

% ************************************************************************
% BUILD *.DURA LIST
% ************************************************************************
dlist = get_msg_list(UW_ID_str,dirs,'dura'); % isus file list

% ************************************************************************
% COPY FLOAT MESSAGE FILES TO LOCAL TEMPORARY DIRECTORY
% Three file groups *.msg, *.isus, *.dura
% ************************************************************************
if isempty(dlist)
    disp(['No float dura files found to process for float ',UW_ID_str]);
    return
elseif isempty(dlist.reg_list) && isempty(dlist.alt_list)
    disp(['No float dura files found to process for float ',UW_ID_str]);
    return
else
    disp(['Copying dura message files to ', dirs.temp, '  .......'])
    % CHECK FLOAT TYPE
    if regexp(dlist.reg_dir,'floats\\f','once') % APEX
        float_type = 'APEX';
        disp([UW_ID_str,' is an APEX float'])
    end
    
    % COPY *.dura files
    if ~isempty(ls([dirs.temp,'*.dura']))
        delete([dirs.temp,'*.dura']) % Clear any message files in temp dir
    end
    if ~isempty(dlist)
        if ~isempty(dlist.reg_list)
            for i = 1:length(dlist.reg_list)
                copyfile([dlist.reg_dir, dlist.reg_list{1,i}],dirs.temp);
            end
        end
        if ~isempty(dlist.alt_list)
            for i = 1:length(dlist.alt_list)
                copyfile([dlist.alt_dir, dlist.alt_list{1,i}],dirs.temp);
            end
        end
        clear dlist i
    end
end

% ************************************************************************
% PROCESS APEX MESSAGE FILES FOR GIVEN FLOAT
% ************************************************************************
msg_list   = ls([dirs.temp,'*.dura']); % get list of file names to process

disp(['Processing ARGO float ',UW_ID_str, '.........'])

merge_pH =[];
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
        continue
    else
        if msg_ct == 1
            pH_hdr = dura.hdr;
            pH_hdr = ['Cast #', pH_hdr];
        end
        pH_data = dura.data;
        pH_data = [ones(size(pH_data(:,1)))*cast_num,pH_data];
    end
    merge_pH = [merge_pH;pH_data];
    clear dura pH_data
end

merge_pH(:,2) = excelsdn(merge_pH(:,2)); %Excell date num for Hans
tNaN = isnan(merge_pH);
merge_pH(tNaN) = -1e10;


%                 ph_rows = size(pH_data,1);
%                 IX      = (ph_rows:-1:1)';
%                 %[B,IX]  = sort(pH_data(:,ipH_p));
%                 pH_data = pH_data(IX,:);
%                 clear dura B IX ph_rows


% *********************************************************************
% CLEAN UP

if ~isempty(ls([dirs.temp,'*.dura']))
    delete([dirs.temp,'*.dura']);
end
%end

% ************************************************************************
% PRINT MERGED DURA FILES TO TEXT
% ************************************************************************
% PRINT HEADER FIRST
out_name = [UW_ID_str,'_ALL_DURA.txt'];
disp(['Printing raw data to: ',dirs.mat,out_name]);
fid_raw  = fopen([dirs.mat, out_name],'W');
fprintf(fid_raw,'//0\r\n');
fprintf(fid_raw,['//File updated on ',datestr(now,'mm/dd/yyyy HH:MM'), ...
       '\r\n']);
fprintf(fid_raw,['//Univ. of Washington ID: ',UW_ID_str,'\r\n']);
fprintf(fid_raw,['//Missing data value = ',MVI_str,'\r\n']);


hdr_size  = size(pH_hdr,2);
data_size = size(merge_pH,1);

% PRINT HEADER
for i = 1:hdr_size % PRINT STANDARD HEADER VARS
    if i < hdr_size
        fprintf(fid_raw,'%s\t',pH_hdr{1,i}); % std vars
    else
        fprintf(fid_raw,'%s\r\n',pH_hdr{1,i}); % std vars
    end
end

% BUILD DATA FORMAT STRING NEXT
data_f ='';
for i = 1:hdr_size
    if i == hdr_size
        data_f   = [data_f,'%e\r\n'];
    elseif regexp(pH_hdr{i},'Ik|Ib|std', 'once')
        data_f   = [data_f,'%e\t'];
    else
        data_f   = [data_f,'%f\t'];
    end
end

% NOW PRINT DATA LINES TO FILE
for sample_ct = 1 : data_size
    fprintf(fid_raw, data_f, merge_pH(sample_ct,:));
end

fclose(fid_raw);
















