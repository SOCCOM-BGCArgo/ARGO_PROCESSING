function d = get_ph_version()

% SCRIPT TO BUILD LIST OF PH SENSOR SN & VERSION NUMBER
d =[];



% ************************************************************************
% SET DIR STRUCTURE
% ************************************************************************
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir = [user_dir, '\Documents\MATLAB\ARGO_PROCESSING\DATA\'];

dirs.temp      = 'C:\temp\';
dirs.ph      = '\\atlas\Chem\DuraFET\APEX+pH Calibrations\pHlogFiles\';
% ************************************************************************
% LOAD AND PARSE pH_log file
% This is a messy excel sheet, condensing down to useful info for APEX
% floats and then adding version defination and descriptions
% ************************************************************************

disp('Trying to extract info from pHlog.xlsx calibration file ...')
[status,msg] = copyfile([dirs.ph, 'pHlog.xlsx'], dirs.temp); % copy to local

[~,~,PH_raw] = xlsread([dirs.temp, 'pHlog.xlsx'], ...
     'CALIBRATION SUMMARY','A5:AH500'); % AH500 # may need to be extended
ph_ver_hdr = {'DF ID' 'UW ID' 'MSC' 'K2' 'K0_pON' 'Build notes'};
ph_ver_data = PH_raw(:,[3,4,6,16,34,9]);

% NOW look for all NaN rows & remove. A little tricky beacuse of mixed data
% types in cell array
tgood = ones(size(ph_ver_data,1),1);
for i = 1:size(ph_ver_data,1)
    try
        if ischar(ph_ver_data{i,2}) % char string in UW ID col (no UW ID #)
            tgood(i) = 0;
            %continue 
        elseif sum(~isnan(cell2mat(ph_ver_data(i,2))),2) == 0 % No  UW ID
            tgood(i) = 0;
            %continue
        elseif ischar(ph_ver_data{i,2}) % char string in UW ID col
            tgood(i) = 0;
            %continue
        elseif sum(~isnan(cell2mat(ph_ver_data(i,:))),2) == 0 % rows all NaN
            tgood(i) = 0;
        elseif cell2mat(ph_ver_data(i,2)) == 8514 && isnan(cell2mat(ph_ver_data(i,6)))
            tgood(i) = 0;
        end       
    catch
        continue
    end
end
ph_ver_data = ph_ver_data(logical(tgood),:);

% NOW DO A SECOND LOOP AND ASSIGN Version type
build_type = cell(size(ph_ver_data,1),2);
iDF = find(strcmp(ph_ver_hdr,'DF ID'));
for i = 1:size(ph_ver_data,1)
    DF = ph_ver_data{i, iDF};
    if isnan(DF)
        continue
    elseif ischar(DF) || DF > 10000
        build_type{i,1} = 0.0;
        build_type{i,2} = 'SBE stem';
    elseif DF < 128
        build_type{i,1} = 1.0;
        build_type{i,2} = 'Ag wire, thinner ISFET cover';
    elseif DF < 164
        build_type{i,1} = 1.1;
        build_type{i,2} = 'Ag wire, thicker ISFET cover';
    elseif DF < 175
        build_type{i,1} = 1.2;
        build_type{i,2} = 'Pt wire, thicker ISFET cover';
    elseif DF < 198
        build_type{i,1} = 2.0;
        build_type{i,2} = 'Double ring, Pt wire, thicker ISFET cover';
    elseif DF < 251 % not certain about this divission
        build_type{i,1} = 3.0;
        build_type{i,2} = 'Smaller PT feed holes + V2.0';     
     elseif DF >= 251 % not certain about this divission
        build_type{i,1} = 4.0;
        build_type{i,2} = 'Roll pin, Smaller PT feed holes + V2.0';
    end
end

%FINALIZE
ph_ver_hdr = [ph_ver_hdr(1:5), 'sensor version', ...
    'verion notes', ph_ver_hdr(6)];
ph_ver_data = [ph_ver_data(:,1:5), build_type, ph_ver_data(:,6)];

d.hdr = ph_ver_hdr;
d.data = ph_ver_data;

save([user_dir,'\PH_DIAGNOSTICS\MBARI_pHsensor_versions.mat'],'d')
clearvars -except d

