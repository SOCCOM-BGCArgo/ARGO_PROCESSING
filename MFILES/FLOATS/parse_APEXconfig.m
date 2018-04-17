% parse_APEXconfig.m

% step through all the APEX floats (UW-MBARI floats) directories and and
% parse each floatconfig.txt file. Extract all calibration coefficients
% and store the coefficients in a *.mat variable.

% ************************************************************************
% PATHS AND NAMES
% ************************************************************************
data_dir    = '\\atlas\chem\ISUS\Argo\'; % config (txt calfile to parse)
local_dir   = 'C:\temp\';
config_file = 'FloatConfig.txt';

O2_3830_vars = {'C0', 'C1', 'C2', 'C3', 'C4', 'P'};
% O2_4330_vars = {'PhaseCoef', 'FoilCoefA', 'FoilCoefB',...
%                 'FoilPolyDegT', 'FoilPolyDegO'};

dc_vars     = {'ChlDC' 'BetabDC'};
scale_vars  = {'ChlScale' 'BetabScale'};
pH_vars     = {'k0' 'k2' 'k3xP' 'k4xP^2' 'k5xP^3' 'k6xP^4' 'k7xP^5' 'k8xP^6'};
pH_pattern  = {'k0' 'k2' 'k3xP' 'k4xP' 'k5xP' 'k6xP' 'k7xP' 'k8xP'};

% ************************************************************************
% GET FLOAT LIST FROM DIRECTORY NAMES
% THIS ASSUMES THAT THE FLOAT NAME STARTS WITH FORM: ^\d{4}[a-z]
floats = dir(data_dir); % get structures
floats = floats(find(~cellfun(@isempty,{floats(:).isdir}))); % only dirs
floats ={floats.name}; % dir names
t1 = regexpi(floats,'^\d{4}[a-z]'); % only dir names w 4#'s followed by letters
floats = floats(find(~cellfun(@isempty,t1)));
floats = floats'; % because easier for me to look at columns

% floats ={'7620SoOcn'}; % for testing
% floats ={'7601STNP'};
% floats ={'9094SoOcn'}
% floats ={'6963HOTpH'}
% ************************************************************************
% LOOP THROUGH EACH FLOAT DIR AND PARSE CONFIG FILE
% CAL ORDER = O2,pH, chl, backscatter
% ************************************************************************
no_config ={}; ct=1;
for i = 1:length(floats)
    sf   = ones(size(scale_vars))*NaN; % predim
    dc   = ones(size(dc_vars))*NaN;
    phcf = ones(size(pH_vars))*NaN;
    O2cf = struct;
    
    % CHECK EXISTANCE AND THEN COPY CONFIG FILE TO LOCAL
    config_path = [data_dir, floats{i}, '\', config_file];
    local_path  = [local_dir, config_file];
    if ~exist(config_path,'file')   
        disp(['NO FILE FOUND AT: ', config_path])
        no_config(ct,1) = {['NO FILE FOUND AT: ', config_path]};
        ct = ct+1;
        continue   
    else
        status    = copyfile(config_path, local_dir);
        if status == 1
            disp([floats{i},' config file copied to local(', local_dir, ')'])
        else
            disp(['Could not copy ', config_path, 'to local'])
            no_config(ct,1) = {['Could not copy ', config_path, 'to local']};
            ct = ct+1;
            continue
        end
    end
    
    % ***********************
    % DO OXYGEN FIRST
    % ***********************
    
    % ********************************************************************
    % PARSE CONFIG FILE TO GET OXYGEN CALIBRATION COEFFICIENTS
    % Examples of OptodeSn output:
    % 'OptodeSn = 737'        (old style)
    % 'OptodeSn = 4330 1168'  (new style 4330)
    % ********************************************************************
    fid   = fopen(local_path); % open config file
    tline = ' ';% initialize 

    while isempty(regexp(tline,'OptodeSn =','once'))
        tline = fgetl(fid);
    end
    % parse any #'s from SN line
    opt_type = sscanf(tline,'%*s %*s %f %f');
    
    if length(opt_type) == 1 && opt_type(1) ~= 4330
        disp(['    model 3830 optode: ', tline])
        % 22 coefficients to extract
        coef = textscan(fid,'%s %*s %f',22); %{1}= ID, {2}= value
        pointer = ftell(fid);
        
        %BUILD COEFF ARRAYS for POLYVAL
        for j = 1:length(O2_3830_vars)
            % Aanderaa = ascending powers
            ind1 = find(strncmp(O2_3830_vars{j},coef{1},...
                   length(O2_3830_vars{j})));
            ind1 = sort(ind1,'descend')'; % Matlab needs descending powers
            eval(['O2cf.p',O2_3830_vars{j},'=coef{2}(ind1);']);
        end;
        O2cf.type = '3830';
    elseif length(opt_type) == 2 && opt_type(1) == 4330
        disp(['    model 4330 optode: ', tline])                 
        while ischar(tline); % EXTRACT COEFFICIENTS
            if ~isempty(regexp(tline,'PhaseCoef','once'))
                ind = regexp(tline,'PhaseCoef','once')+ 19;
                O2cf.PCoef = cell2mat((textscan(tline(ind:end),'%f',4)));
            end
            if ~isempty(regexp(tline,'FoilCoefA','once'))
                ind = regexp(tline,'FoilCoefA','once')+ 19;
                FCoefA = cell2mat((textscan(tline(ind:end),'%f',14)));
            end
            if ~isempty(regexp(tline,'FoilCoefB','once'))
                ind = regexp(tline,'FoilCoefB','once')+ 19;
                FCoefB = cell2mat((textscan(tline(ind:end),'%f',14)));
            end
            if ~isempty(regexp(tline,'FoilPolyDegT','once'))
                ind = regexp(tline,'FoilPolyDegT','once')+ 22;
                O2cf.PolyDegT = cell2mat((textscan(tline(ind:end),'%f',28)));
            end
            if ~isempty(regexp(tline,'FoilPolyDegO','once'))
                ind = regexp(tline,'FoilPolyDegO','once')+ 22;
                O2cf.PolyDegO = cell2mat((textscan(tline(ind:end),'%f',28)));
                pointer = ftell(fid);
                break % got all the coefficients so move on
            end
            tline = fgetl(fid);
        end
        O2cf.FCoef =[FCoefA;FCoefB];
        t1 = O2cf.FCoef == 0; % Look for zero, remove to shorten calc, 0 * x = 0
        O2cf.FCoef(t1)    = [];
        O2cf.PolyDegT(t1) = [];
        O2cf.PolyDegO(t1) = [];
        O2cf.type = '4330';
    else
        s1=['CAN NOT DETERMINE OPTODE TYPE FOR ',floats{i},...
            ' MOVING TO NEXT FLOAT'];
        disp(s1)
        no_config(ct,1) = {s1};
        ct =ct+1;
        fclose(fid);
        continue
    end
    clear ind ind1 FCoefA FCoefB

    % *****************************************
    % LOOK FOR pH, IF EXISTS PARSE COEFFICIENTS
    % *****************************************
    ph_test  = 0;
    chl_test = 0;
    while ischar(tline)
        if ph_test == 1
            for j = 1:length(pH_vars) % look for ph var names in tline
                ind1 = regexp(tline,pH_pattern{j},'once');
                if ~isempty(ind1)
                    tmp  = textscan(tline,'%f %s',1,'Delimiter',',');
                    phcf(strcmp(strtrim(tmp{1,2}),pH_vars))  = tmp{1,1};
                    continue % found a match move on to next tline
                end
            end
        elseif ~isempty(regexp(tline,'^CHLFLUOR','once'))
            disp(['    ',tline])
            chl_test = 1;
            break % NO pH BUT CHL EXISTS
        elseif ~isempty(regexp(tline,'^Durafet','once')) % pH exists so extract
            disp(['    ',floats{i},'   ',tline])
            ph_test =1;
        end
        tline = fgetl(fid);
    end

    if ~ischar(tline) % NO pH - reset pointer to end of O2 & look for CHL
        fseek(fid, pointer,'bof');
    end
    
    % *****************************************
    % LOOK FOR CHL, IF EXISTS PARSE COEFFICIENTS
    % *****************************************
    while  chl_test == 1 && ischar(tline)
        ind1 = regexp(tline,'\w+Scale','once');       % scale values
        ind2 = regexp(tline,'\w+DC','once');          % DC values
        if ~isempty(ind1) %SCALE VALUES
            tmp = textscan(tline,'%f %s',1); %use format once
            sf(strcmp(tmp{1,2},scale_vars))  = tmp{1,1};
        elseif ~isempty(ind2) %DC VALUES
            tmp = textscan(tline,'%f %s',1); %use format once
            dc(strcmp(tmp{1,2},dc_vars))  = tmp{1,1};
        end
        tline = fgetl(fid);
    end
    
    if ph_test == 0
        phcf =[]; %set to empty
    end
    
    if chl_test == 0; % NO CHL 
        sf =[]; %set to empty
        dc =[];
    end
    
    fclose(fid);
    save([local_dir,floats{i},'_cal.mat'],'O2cf','phcf', 'sf','dc');
    %clear fid
    delete(local_path) % Remove config file from local once done with it
    %pause
end
clear t1 status ct i fid
    