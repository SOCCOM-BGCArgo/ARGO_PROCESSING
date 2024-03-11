function synthfull = build_BSOLO_syn_prof(mat_file, var_list)
% This function merges all BGC sensors axes for a given cycle. It is based
% Henry Bittig's Sprof profile code ("ARGO_simplified_profile.m")

% CHANGES
% 11/06/2023 - JP - Added check block to make sure S.LR subfields contain data.
%    If all fields are empty remove LR from S. Fix for SS0002.102 which 
%    has no LR data
% 12/19/2023 - JP - Commented out previous fix which was due to a bug in
%                   Johns's code starting around ~Line 84
% 02/12/2024 - JP - minor fix to deal with repeated pressure values at
%                   surface for s4003 cycle 29 which was breaking some
%                   interploation code block (used unique instead of sort).

% *************************************************************************
% TESTING
% fd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\4903026\';
% fn = '4903026.009.mat';

% fd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\5906320\';
% fn = '5906320.009.mat';

% fd = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATS\5906765\';
% fd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\5906765\';
% fn = '5906765.102.mat'; % ss0002
% 

% fd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\2903886\';
% fn = '2903886.029.mat'; % ss0002
% mat_file = [fd,fn];
% var_list = [];
% *************************************************************************


% *************************************************************************
% DO SOME PREP WORK
% *************************************************************************
includeTSflag   = 0;
addTSeverywhere = 1;
addoffsetflag   = 1;
tf_APEX_OCR     = 0;

% This will always be the axis order, HR field contains primary P axis and will
% always be there. Other field may or may not exist
field_order   = {'HR' 'LR' 'BGC01' 'BGC02' 'BGC03' 'BGC04' 'BGC05' 'BGC06'};

if ~isempty(regexp(mat_file,'5906320|5906446','once')) % APEX OCR 19191 & 19314
    field_order = {'HR' 'LR' 'OCR'};
    tf_APEX_OCR = 1;
end

% SENSOR OFFSETS meters below CTD (+ = below), ECO strapped to float midsection
% ECO_offset  = 0.5;
% NO3_offset  = 0.8;

ECO_offset  = 0;
NO3_offset  = 0;

if isempty(var_list) % Default for solo if nothing defined
    %  BUILD LIST OF POSSIBLE BGC SOLO ARGO VARIABLES (UPDATE WHEN NEW VARIABLES ARISE)
    %      *_QC, *_ADJUSTED, & *_ADJUSTED_QC BUILT FROM RAW PARAMETER NAME
    raw_vars{1,1}  = 'PRES';
    raw_vars{2,1}  = 'TEMP';
    raw_vars{3,1}  = 'PSAL';
    raw_vars{4,1}  = 'DOXY';
    raw_vars{5,1}  = 'DOXY2';
    raw_vars{6,1}  = 'NITRATE';
    raw_vars{7,1}  = 'CHLA';
    raw_vars{8,1}  = 'BBP700';
    raw_vars{9,1}  = 'CDOM';
    raw_vars{10,1} = 'BBP532';
    raw_vars{11,1} = 'PH_IN_SITU_TOTAL';
    raw_vars{12,1} = 'DOWN_IRRADIANCE380';
    raw_vars{13,1} = 'DOWN_IRRADIANCE412';
    raw_vars{14,1} = 'DOWN_IRRADIANCE490';
    raw_vars{15,1} = 'DOWNWELLING_PAR';

    var_list       = raw_vars;
end
master_bgc_params = setdiff(var_list,{'PRES' 'TEMP' 'PSAL'},'stable');

% ************************************************************************
% ************************************************************************
%                     PROCESS *.MAT FILES PART 1:
%              ALIGN MULTIPLE PRES AXES & DATA ONTO 1 AXIS
% ************************************************************************
% ************************************************************************
S       = load(mat_file); % load cycle mat vars into structure

% ADD BLOCK HERE TO CHECK IF SECONDARY CTD AXES ("LR") HAS DATA. IF IT
% DOES NOT REMOVE FIELD. THIS CATCHES ss0002 cycle 102
% Used David Hill's approach to merge all data from structure fields into an array
% https://www.mathworks.com/matlabcentral/answers/852330-converting-struct-field-to-array

% JP 12/19/23 commented out since John's code now returning all axes.
% cycle 102 returned 1 less
% if isempty(cell2mat(struct2cell(S.LR)))
%     S = rmfield(S,'LR');
% end

% GATHER SOME CYCLE INFO
fields    = fieldnames(S);
cycle     = S.INFO.cast;
gps       = S.INFO.gps;

if ~tf_APEX_OCR
    EOP_sdn   = S.INFO.EOP_sdn;
    p_axes_ct = S.INFO.pres_axes_ct;
    sensors   = S.INFO.sensors;
else
    EOP_sdn   = S.INFO.sdn;
    p_axes_ct = 3;
    sensors   = {'SBE41CP' 'SBE41CP' 'OCR'};
end

% ALTER SENSORS STRINGS FOR AXES PRIORITY TESTING
% A BIT OF A FUDGE NOT SURE WHAT WILL HAPPEN WITH MULTI BGC ON ONE AXIS
% MAY NEED TO LOOP THROUGH & BUILD HENRY LIKE MERGED STRINGS
sens_priority = regexprep(sensors,'SBE41CP','CTD'); % CTD
sens_priority = regexprep(sens_priority,'ECO','BBP'); % ECOsens_priority = regexprep(sens_priority,'ECO','BBP'); % ECO
sens_priority = regexprep(sens_priority,'ALK','TRANSISTOR_PH'); % ECO
sens_priority = regexprep(sens_priority,'DOXY','OPTODE_DOXY'); % ECO
sens_priority = regexprep(sens_priority,'NO3','SPECTROPHTOMETER_NITRATE'); % ECO
sens_priority = regexprep(sens_priority,'OCR','RADIOMETER'); % ECO
[~,asort]     = sort(sens_priority(2:end));
asort         = [1 asort+1];


% RE-ORDER & REDUCE FIELD NAMES TO MEASURED AXES ORDER, THIS IS BASED
% ON PREDEFINED LIST ("field_order") DEFINED ABOVE
[~,Lib] = ismember(field_order',fields);
Lib(Lib == 0) =[];
pfields =  fields(Lib);  % fields now in sequential axis order
if size(pfields,1) ~= p_axes_ct
    str = fprintf(['\nWarning WMO %s pressure axes count & pressure ',...
        'axes field counts do not match for cycle %0.0f\n'],...
        S.INFO.WMO_ID,cycle);
    return
end


% FIND ALL BGC B PARAMETERS IN CYCLE. COMPARE TO USER DEFINED LIST OF
% BGC PARAMETERS TO LOOK FOR ("master_bgc_params")
bgcparams = cell(1,p_axes_ct); % cell(i) will be n x 1 if ax has multi params
bgcP_idx  = bgcparams;
bgcflag   = false(1,p_axes_ct);
for i = 1:p_axes_ct % step through pressure axes
    bgcparams{i} = intersect(master_bgc_params, fieldnames(S.(pfields{i})));
    if ~isempty(bgcparams{i})
        bgcflag(i) = 1;
    end
    bgcP_idx{i}  = ones(size(bgcparams{i}))*i; % axis index for param
end
% convert to N x 1 cell & numeric arrays, keeping axes order
ubgcparams = unique(cat(1, bgcparams{:}),'stable'); % List of bgc params only
ubgcP_idx  = cat(1, bgcP_idx{:}); %axis index array linking field to param

% ********************************************************************
% ********************************************************************
% IF SPECIAL CASE APEX OCR FLOAT, DEAL WITH NON UNIQUE & REPEATED
% SURFACE PRESSURE VALUES (UPPER 4 M)- Take MEDIAN OF REPEATED VALUES
% HENRY'S QUICK CHECK DOESN'T WORK FULLY

if tf_APEX_OCR
    S.INFO.file_name    = mat_file; % NOT IN OCR FILED INFO
    S.INFO.BOP_sdn      = NaN;
    S.INFO.EOP_sdn      = S.INFO.sdn;
    S.INFO.sensors      = sensors;
    S.INFO.pres_axes_ct = p_axes_ct;
    S.INFO.bgc_pres_axes_ct = 1;

    P = S.OCR.PRES;
    tP = P<4; % less than 4 dbar
    PT = P(tP);
    [uP,ia,ic] = unique(PT); %ic top 5m indices

    OCR_flds = fieldnames(S.OCR);
    tg = ~cellfun(@isempty, regexp(OCR_flds,'\d+$$|PAR$','once'));
    OCR_flds = OCR_flds(tg);

    for i = 1: size(OCR_flds,1)
        bgcp = OCR_flds{i}; % get apram name

        for j = 1:max(ic) % step through repeated surface P levels
            inds = find(ic == j); % repeated values at uP(j)
            S.OCR.(bgcp)(j) = median(S.OCR.(bgcp)(inds),1,'omitnan');
            S.OCR.([bgcp '_QC'])(j) = ...
                max(S.OCR.([bgcp '_QC'])(inds),[],1,'omitnan');

            if ~strncmp(bgcp,'RAW',3)
                S.OCR.([bgcp '_ADJUSTED'])(j) = ...
                    median(S.OCR.([bgcp '_ADJUSTED'])(inds),1,'omitnan');
                S.OCR.([bgcp '_ADJUSTED_QC'])(j) = ...
                    max(S.OCR.([bgcp '_ADJUSTED_QC'])(inds),[],1,'omitnan');
                S.OCR.([bgcp '_ADJUSTED_ERROR'])(j) = ...
                    max(S.OCR.([bgcp '_ADJUSTED_ERROR'])(inds),[],1,'omitnan');
            end

            if i == 1
                S.OCR.PRES_QC(j) = max(S.OCR.PRES_QC(inds),[],1,'omitnan');
                S.OCR.PRES_ADJUSTED_QC(j) =  ...
                    max(S.OCR.PRES_ADJUSTED_QC(inds),[],1,'omitnan');
            end
        end
        % MERGE MEDAIN VALUES WITH EXISTING DEEPER VALUES
        S.OCR.(bgcp) = [S.OCR.(bgcp)(1:j);S.OCR.(bgcp)(~tP)];
        S.OCR.([bgcp '_QC']) = ...
            [S.OCR.([bgcp '_QC'])(1:j);S.OCR.([bgcp '_QC'])(~tP)];

        if ~strncmp(bgcp,'RAW',3)
            S.OCR.([bgcp '_ADJUSTED']) = ...
                [S.OCR.([bgcp '_ADJUSTED'])(1:j);S.OCR.([bgcp '_ADJUSTED'])(~tP)];
            S.OCR.([bgcp '_ADJUSTED_QC']) = ...
                [S.OCR.([bgcp '_ADJUSTED_QC'])(1:j);S.OCR.([bgcp '_ADJUSTED_QC'])(~tP)];
            S.OCR.([bgcp '_ADJUSTED_ERROR']) = ...
                [S.OCR.([bgcp '_ADJUSTED_ERROR'])(1:j);S.OCR.([bgcp '_ADJUSTED_ERROR'])(~tP)];
        end

        if i == 1
            S.OCR.PRES_QC = [S.OCR.PRES_QC(1:j);S.OCR.PRES_QC(~tP)];
            S.OCR.PRES_ADJUSTED_QC = ...
                [S.OCR.PRES_ADJUSTED_QC(1:j);S.OCR.PRES_ADJUSTED_QC(~tP)];
        end
    end
    S.OCR.PRES = [uP; S.OCR.PRES(~tP)];
    S.OCR.PRES_ADJUSTED = S.OCR.PRES;
end
% ********************************************************************
% ********************************************************************

% STEP THROUGH ALL PRESSURE AXES & BUILD r x n matrices of PRES & PRES_QC
allP   = ones(1000, p_axes_ct)*NaN;
allPQC = allP; % sensor offset array

for i = 1:p_axes_ct
    P   = S.(pfields{i}).PRES;
    PQC = S.(pfields{i}).PRES_QC;
    rP  = size(P,1);

    % CHECK FOR ANY PRESSURE INVERSIONS. IF FOUND SET FLAG TO BAD
    ind             = ~isnan(P) & ismember(PQC,[0 1 2 3]);
    pinversion      = logical([diff(P(ind))<=0;0]);
    PQC(pinversion) = 4;

    % FILL IN MATRIX
    allP(1:rP,i)   = P;
    allPQC(1:rP,i) = PQC;
end

tnan       = sum(~isnan(allP),2) == 0;
allP       = allP(~tnan,:);
allPQC     = allPQC(~tnan,:);

tgP        = ismember(allPQC,[0 1 2 3]);
allP(~tgP) = NaN;
allP0      = allP; % Make a copy before any offset alterations
clear P PQC rP

% check where BGC samples are present (and not all NaN)
% build sensor offset array
bgcpresence = false(size(allP));
bgc_offsets = ones(1, size(ubgcparams,1)) * 0; % predim offsets array
for i=1:length(ubgcparams)
    axis_idx = ubgcP_idx(i); % get approriate structure
    bgc_data = S.(pfields{axis_idx}).(ubgcparams{i}); % param data
    tg       = ~isnan(bgc_data); % any data?
    rtg      = size(tg,1);
    bgcpresence(1:rtg, axis_idx) = bgcpresence(1:rtg,axis_idx) | tg;

    % ADD SENSOR OFFSETS TO ARRAY IF NEED BE (HARD WIRED)
    if regexp(ubgcparams{i},'^NITRATE','once')
        bgc_offsets(i) = NO3_offset;
    elseif regexp(ubgcparams{i},'^BBP|^CDOM|^CHL','once')
        bgc_offsets(i) = ECO_offset;
    end
end
clear i axis idx tg rtg


% INCLUDE T & S ONLY P AXES?
if ~includeTSflag
    allP(:,~bgcflag)   = NaN;
    allP(~bgcpresence) = NaN;
end

% split up pressure in each N_PROF for each parameter/sensor
% e.g. NITRATE and DOXY in same N_PROF can have different vertical offsets
if addoffsetflag
    linind = cellfun(@(x,y)repmat(x,length(y),1), ...
        num2cell(1:p_axes_ct,1),bgcparams,'uniform',0); % index to N_PROF in #repetitions like bgcparams
    linbgcparams = cat(1, bgcparams{:});
    linind       = cat(1, linind{:});
    allP         = allP(:,linind); % repeat pressure for each b-parameter
    allP         = allP + ones(size(allP(:,1)))* bgc_offsets; % add offsets
end % split up N_PROFs only if vertical offsets are added per sensor

% ********************************************************************
% BUILD SYNTHETIC PRESSURE AXIS. If apply offset is true it seems CTD
% only press axes are not used? Is this the intent?
uallP = unique(allP(~isnan(allP))); % unique bgc pres values

% *********************************************************************
% Check which pressure levels are present in which profile & get dP's
prespresent = false(length(uallP),size(allP,2));
for i = 1:size(allP,2)
    prespresent(:,i) = ismember(uallP, allP(:,i));
end

% get pressure differences between the levels in each N_PROF; use the
% minimum of preceeding/succeeding deltaPRES
valdp = ones(length(uallP), size(allP,2))*NaN;
for i=1:size(allP,2)
    pres = allP(~isnan(allP(:,i)),i);
    if ~isempty(pres)
        if length(pres)>1
            dpres = diff(pres);
            valdp(prespresent(:,i),i) = min(abs([dpres(1);dpres]),abs([dpres;dpres(end)]));
        else
            valdp(prespresent(:,i),i)=NaN;
        end
    end
    clear pres dpres
end

% *********************************************************************
%  NOW BUILD INDEX ARRAY FOR SYNTHETIC PRESSURE AXIS
% *********************************************************************
% Start from bottom & work towards the surface

% cycle through record from bottom
useind   = false(size(uallP)); %predim index array
i        = length(uallP); % start at bottom
niter    = 0;
nitermax = length(uallP)+1;

if all(isnan(valdp(:))) % each BGC obs has only one level - no way to jump/construct synthetic pressure axis
    presaxis=uallP(end); % take deepest pressure and align rest to this level
else % cycle pressure record
    while ~isempty(i)
        niter = niter+1;
        useind(i) = 1; % add current level to synthetic axis
        ind=find(uallP > uallP(i) - min(valdp(i,:)) & uallP <= uallP(i)); % get pressures that are within current level (included) and probably next level-min(dPRES) (excluded)
        % check if any of the intermittent levels has such a small dPRES, that
        % there will be a second observation within the current level-min(dPRES) "jump"
        obspresence = ~isnan(valdp(ind,:)); % make presence of non-FillValue a logical array
        if ~isempty(ind) && any(sum(obspresence,1)>1) % there are other levels in between and they have more than one observation, i.e., a denser sampling interval
            % go to deepest upres that features a second observation in the
            % same N_PROF: sum #obs from bottom in each N_PROF, get max in each
            % line (upres), and jump to (deepest) line that has >1
            i=ind(find(max(flipud(cumsum(flipud(obspresence),1)),[],2)>1,1,'last'));
            if isempty(i)
                disp(['S-PROF_WARNING: File ' bfilestr '.nc: Trouble',...
                    ' during creation of synthetic pressure axis. ',...
                    'Create synthetic profile only with available core data.'])
                useind=[]; % default to empty presaxis if failed
                break % loop of pressure levels from bottom
            end
        else % jump by at least current level+min(dPRES)
            if all(isnan(valdp(i,:))) % no min(dPRES) available for current level:
                % jump to next level that has a non-NaN valdp
                i=find(uallP<uallP(i) & any(~isnan(valdp),2),1,'last');
            else % jump by at least current level+min(dPRES)
                i=find(uallP<=uallP(i)-min(valdp(i,:)),1,'last');
            end
        end
        clear obspresence
        if niter>nitermax
            disp(['S-PROF_WARNING: File ' bfilestr '.nc: Exceeded maximum number of iterations in selection of synthetic pressure levels. Should not happen... Create synthetic profile only with available core data.'])
            useind = []; % default to empty presaxis if failed
            break % loop of pressure levels from bottom
        end
    end
    clear niter nitermax
    presaxis = uallP(useind);
end % pressure axis cycle from bottom possible

% *********************************************************************
%  SYNTHETIC PRESSURE AXIS BUILT - NOW PREP TO FILL / ALIGN BGC DATA
% *********************************************************************
coreparams = {'PRES';'TEMP';'PSAL'}; % define core parameters

% get non-overlapping pressure axis for core parameters
xpres   = allP0(:,asort); % keep original pressure
xpresqc = tgP(:,asort); % and flag to pressure that are not QC=4

% build argo like param matrix for P, T, S to fit Henry's flow
for i=1:length(coreparams)
    c_param    = coreparams{i};
    %disp(c_param)
    c_param_qc = [c_param ,'_QC'];
    y       = allP0 * NaN; % predim
    yqc     = y;
    for j   = 1:p_axes_ct % step through axes & look for c_param
        fnames = fieldnames(S.(pfields{j}));
        tf     = strcmp(fnames, c_param); % is param in structure?
        if sum(tf) == 0
            continue
        end
        x           = S.(pfields{j}).(c_param);
        xqc         = S.(pfields{j}).(c_param_qc);
        rx          = size(x,1);
        y(1:rx,j)   = x;
        yqc(1:rx,j) = xqc;
    end
    y   = y(:,asort);
    yqc = yqc(:,asort);
    clear j tf x xqg rx


    % extract data

    if ismember(coreparams{i},{'PRES'})
        overlap=true(size(y)); % pressures of different N_PROFs (almost) necessarily overlap
    else % not PRES
        % get max and min ranges per nprof
        xrange=ones(2,1+p_axes_ct)*NaN;
        overlap=false(size(y)); % and flag portions that don't overlap
        for k=1:p_axes_ct
            ind=~isnan(y(:,k)) & xpresqc(:,k);
            if any(ind) % data in current nprof and with proper pres
                xrange(:,1+k)=[min(xpres(ind,k)); max(xpres(ind,k))];
                % check if it should be added or not
                if isnan(xrange(1,1)) % first nprof with data
                    xrange(:,1)=xrange(:,1+k);
                    overlap(:,k)=ind;
                else % more than one nprof with data: keep only data outside existing range
                    overlap(:,k)=ind & (xpres(:,k)>xrange(2,1) | xpres(:,k)<xrange(1,1));
                    xrange(1,1)=min([xrange(1,1) xrange(1,1+k)]);
                    xrange(2,1)=max([xrange(2,1) xrange(2,1+k)]);
                end
            end
        end % cycle all nprofs
        clear ind xrange
    end % overlapping portion

    % do not use <PARAM>_QC of 8 or FillValue
    yflagqc   = ismember(yqc,[0 1 2 3 4 5]);
    yflagnoFV = ~isnan(y);

    % clean up: only data of current parameter without QC 8, only pflag
    x       = xpres(xpresqc & yflagqc & yflagnoFV);
    y       = y(xpresqc & yflagqc & yflagnoFV);
    overlap = overlap(xpresqc & yflagqc & yflagnoFV);
    clear yflagqc yflagnoFV yqc

    if ~isempty(y) % not only NaN/FV data
        % use only non-overlapping portion
        x=x(overlap);
        clear overlap
        % make monotonic for interpolation (looses nprof priority!)
        [~,ind]=sort(x); % monotonic sorting
        x=x(ind);
        full.(coreparams{i}).x=x;
    else % empty pressure axis
        full.(coreparams{i}).x = [];
    end % not only NaN/FV data
    clear x y
end

% get all non-overlapping (HR-)TEMP/PSAL levels and synthetic pressure axis levels
presmerge    = union(union(full.TEMP.x,full.PSAL.x),presaxis); % 'HR' T/S pressure axis and presaxis
% get indices for synthetic pressure axis
nosf         = length(presmerge);
[~,synthind] = intersect(presmerge,presaxis);
%clear full


% *************************************************************************
%                    FILL IN DATA ON SYNTHETIC AXIS
% *************************************************************************
ubgcparams = cat(1, {'PRES';'TEMP';'PSAL'},ubgcparams(:));
bgc_offsets = [0;0;0; bgc_offsets'];
isynth     = (1:length(presaxis))'; % index 1..length synthetic pressure axis
nos        = length(isynth);
% DOXY/<PARAM> can sit in more than one N_PROF, so need to be a bit more
% clunky than just simply interpolating a single N_PROF which contains
% <PARAM>

xpresqc          = tgP(:,asort); % and flag to pressure that are not QC=4
synth.PRES.value = presaxis;

if isempty(presmerge) % core and presaxis empty
    for i = 1:length(ubgcparams) % tap all with FillValue (avoid N_LEVEL dimension of 0)
        synth.(ubgcparams{i}).value                     = NaN;
        synth.([ubgcparams{i} '_QC']).value             = NaN;
        synth.([ubgcparams{i} '_ADJUSTED']).value       = NaN;
        synth.([ubgcparams{i} '_ADJUSTED_QC']).value    = NaN;
        synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value = NaN;
        synth.([ubgcparams{i} '_dPRES']).value          = NaN;
    end
    nosf     = 1;
    synthind = 1; % avoid N_LEVEL dimension of 0 for BGC parameters
else
    for i=1:length(ubgcparams)
        b_param = ubgcparams{i};
        if ismember(b_param,{'PRES';'TEMP';'PSAL'}) % check which pressure axis to use
            sind=1:length(presmerge); % synthetic levels and full-HR PTS data for core
        else % BGC
            sind=synthind; % synthetic levels only
        end
        spresaxis = presmerge(sind); % pressure axis for current parameter
        spresaxis = spresaxis(:); % make sure it's a column vector

        if addoffsetflag
            s_offset = 0; % default to 0
            % re-apply vertical offset from scratch for particular parameter
            % e.g., N_PROF with SUNA & CTD has large offset for NITRATE, but none
            % for TEMP & PSAL, and only one S.PRES.value / inpres for both
            %if ismember(ubgcparams{i},voff.linbgcparams)
            % check to which N_PROF the current parameter  fits
            ind   = cellfun(@(x)ismember(b_param,x),bgcparams); % for P axis
            o_ind = strcmp(ubgcparams, b_param); % for offset grab
            if sum(ind) == 1
                s_offset = bgc_offsets(o_ind);
            end

            xpres = allP0;

            % NaN+offset = NaN so no indexing needed
            xpres(:,ind)=xpres(:,ind) + ones(size(xpres(:,1)))*s_offset;
            %xpres(inpres0==FV)=FV; % keep FV
            xpres=xpres(:,asort);
            clear ind
            %else % no match found
            %xpres=inpres0(:,asort); % keep original pressure
            %end % found vertical offsets in config

            %             if strcmp(ubgcparams{i},'PRES')
            %                 xpres=inpres0(:,asort); % forget everything above, use all unadjusted pressures
            %                 %xpres=S.PRES.value(:,asort); % forget everything above, use all unadjusted pressures
            %             end
        else % no addoffsetflag
            xpres=allP0(:,asort); % keep original pressure
        end % addoffsetflag


        % NEED TO BUILD ARGO LIKE MATRICES FOR EACH BGC PARAM TO FIT
        % HENRY's FLOW
        y   = allP0 * NaN; yqc = y; yadj = y; yadjqc = y; yadjerr = y;% predim

        for j = 1:p_axes_ct % step through axes & look for c_param
            fnames = fieldnames(S.(pfields{j}));
            tf     = strcmp(fnames, b_param); % is param in structure?
            if sum(tf) == 0
                continue
            end
            r              = size(S.(pfields{j}).(b_param),1);
            y(1:r,j)       = S.(pfields{j}).(b_param);
            yqc(1:r,j)     = S.(pfields{j}).([b_param,'_QC']);
            yadj(1:r,j)    = S.(pfields{j}).([b_param,'_ADJUSTED']);
            yadjqc(1:r,j)  = S.(pfields{j}).([b_param,'_ADJUSTED_QC']);
            if isfield(S.(pfields{j}), [b_param,'_ADJUSTED_ERROR'])
                yadjerr(1:r,j) = S.(pfields{j}).([b_param,'_ADJUSTED_ERROR']);
            else
                yadjerr(1:r,j) = S.(pfields{j}).(b_param) * NaN;
            end
        end
        y       = y(:,asort);
        yqc     = yqc(:,asort);
        yadj    = yadj(:,asort);
        yadjqc  = yadjqc(:,asort);
        yadjerr = yadjerr(:,asort);




        % extract data
        %         y=S.(ubgcparams{i}).value(:,asort);
        %         yqc=S.([ubgcparams{i} '_QC']).value(:,asort);
        %         yadj=S.([ubgcparams{i} '_ADJUSTED']).value(:,asort);
        %         yadjqc=S.([ubgcparams{i} '_ADJUSTED_QC']).value(:,asort);
        %         yadjerr=S.([ubgcparams{i} '_ADJUSTED_ERROR']).value(:,asort);
        % check for overlapping portion of N_PROFs
        if ismember(ubgcparams{i},{'PRES'})
            overlap=true(size(y)); % pressures of different N_PROFs (almost) necessarily overlap
            % get max and min ranges per nprof for adjusted anyway
            xrangeadj=ones(2,1+p_axes_ct)*NaN;
            overlapadj=false(size(yadj)); % and flag portions that don't overlap for adjusted data
            for k=1:p_axes_ct
                % and adjusted
                %indadj=~isnan(yadj(:,k)) & xpresqc(:,k); % use only levels without FillValue and with useful PRES_QC
                indadj=(~isnan(yadj(:,k)) | ismember(yadjqc(:,k),4)) & xpresqc(:,k); % consider that adjusted_QC 4 may have FillValue adjusted data
                if any(indadj) % data in current nprof and with proper pres
                    xrangeadj(:,1+k)=[min(xpres(indadj,k)); max(xpres(indadj,k))];
                    % check if it should be added or not
                    if isnan(xrangeadj(1,1)) % first nprof with data
                        xrangeadj(:,1)=xrangeadj(:,1+k);
                        overlapadj(:,k)=indadj;
                    else % more than one nprof with data: keep only data outside existing range
                        overlapadj(:,k)=indadj & (xpres(:,k)>xrangeadj(2,1) | xpres(:,k)<xrangeadj(1,1));
                        xrangeadj(1,1)=min([xrangeadj(1,1) xrangeadj(1,1+k)]);
                        xrangeadj(2,1)=max([xrangeadj(2,1) xrangeadj(2,1+k)]);
                    end
                end
            end % cycle all nprofs
            clear indadj xrangeadj

        else % not PRES
            % get max and min ranges per nprof
            xrange=ones(2,1+p_axes_ct)*NaN;
            xrangeadj=ones(2,1+p_axes_ct)*NaN;
            overlap=false(size(y)); % and flag portions that don't overlap
            overlapadj=false(size(yadj)); % and flag portions that don't overlap for adjusted data
            for k=1:p_axes_ct
                ind=~isnan(y(:,k)) & xpresqc(:,k);
                if any(ind) % data in current nprof and with proper pres
                    xrange(:,1+k)=[min(xpres(ind,k)); max(xpres(ind,k))];
                    % check if it should be added or not
                    if isnan(xrange(1,1)) % first nprof with data
                        xrange(:,1)=xrange(:,1+k);
                        overlap(:,k)=ind;
                    else % more than one nprof with data: keep only data outside existing range
                        overlap(:,k)=ind & (xpres(:,k)>xrange(2,1) | xpres(:,k)<xrange(1,1));
                        xrange(1,1)=min([xrange(1,1) xrange(1,1+k)]);
                        xrange(2,1)=max([xrange(2,1) xrange(2,1+k)]);
                    end
                end
                % and adjusted
                %indadj=~isnan(yadj(:,k)) & xpresqc(:,k); % use only levels without FillValue and with useful PRES_QC
                indadj=(~isnan(yadj(:,k)) | ismember(yadjqc(:,k),4)) & xpresqc(:,k); % consider that adjusted_QC 4 may have FillValue adjusted data
                if any(indadj) % data in current nprof and with proper pres
                    xrangeadj(:,1+k)=[min(xpres(indadj,k)); max(xpres(indadj,k))];
                    % check if it should be added or not
                    if isnan(xrangeadj(1,1)) % first nprof with data
                        xrangeadj(:,1)=xrangeadj(:,1+k);
                        overlapadj(:,k)=indadj;
                    else % more than one nprof with data: keep only data outside existing range
                        overlapadj(:,k)=indadj & (xpres(:,k)>xrangeadj(2,1) | xpres(:,k)<xrangeadj(1,1));
                        xrangeadj(1,1)=min([xrangeadj(1,1) xrangeadj(1,1+k)]);
                        xrangeadj(2,1)=max([xrangeadj(2,1) xrangeadj(2,1+k)]);
                    end
                end
            end % cycle all nprofs
            clear ind xrange indadj xrangeadj
        end % check for overlapping portion of N_PROFs PRES (int32/FV) or other style (double/NaN)
        % do not use <PARAM>_QC of 8 or FillValue
        yflagqc = ismember(yqc,[0 1 2 3 4 5]);
        yflagnoFV = ~isnan(y);



        yflagadjqc = ismember(yadjqc,[0 1 2 3 4 5]);
        %yflagadjnoFV=~isnan(yadj); % double check FV, use only levels without FillValue
        yflagadjnoFV = ~isnan(yadj) | ismember(yadjqc,4);  % double check FV, consider adjusted QC 4 as exception to FillValue

        % clean up: only data of current parameter without QC 8, only pflag
        x=xpres(xpresqc & yflagqc & yflagnoFV);
        xadj=xpres(xpresqc & yflagadjqc & yflagadjnoFV);
        %xqc=xpresqc(xpresqc & yflagqc);
        y=y(xpresqc & yflagqc & yflagnoFV);
        yqc=yqc(xpresqc & yflagqc & yflagnoFV);
        yadj=yadj(xpresqc & yflagadjqc & yflagadjnoFV);
        yadjqc=yadjqc(xpresqc & yflagadjqc & yflagadjnoFV);
        yadjerr=yadjerr(xpresqc & yflagadjqc & yflagadjnoFV);
        overlapadj=overlapadj(xpresqc & yflagadjqc & yflagadjnoFV);
        overlap=overlap(xpresqc & yflagqc & yflagnoFV);
        yadjerrpresence=~all(isnan(yadjerr)); % error not mandatory for adjusted fields..
        clear yflagqc yflagnoFV yflagadjqc yflagadjnoFV

        % preallocate
        synth.(ubgcparams{i}).value=double(spresaxis)*NaN;

        %         if ismember(ubgcparams{i},{'PRES'})
        %             synth.(ubgcparams{i}).value=int32(ones(size(spresaxis))*FV);
        %         else
        %             synth.(ubgcparams{i}).value=double(spresaxis)*NaN;
        %         end

        synth.([ubgcparams{i} '_QC']).value=double(spresaxis)*NaN;
        synth.([ubgcparams{i} '_ADJUSTED']).value=double(spresaxis)*NaN;
        synth.([ubgcparams{i} '_ADJUSTED_QC']).value=double(spresaxis)*NaN;
        synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value=double(spresaxis)*NaN;
        %synth.([ubgcparams{i} '_dPRES']).value=int32(ones(size(spresaxis))*FV);
        synth.([ubgcparams{i} '_dPRES']).value= ones(size(spresaxis))*NaN;

        if ismember(ubgcparams{i},{'PRES'}) % add synthetic levels (incl. vertical offsets of BGC sensors)
            synth.PRES.value(synthind)=presaxis;
        end

        if ~isempty(y) % not only NaN/FV data
            % use only non-overlapping portion
            x=x(overlap);y=y(overlap);yqc=yqc(overlap);
            clear overlap
            % make monotonic for interpolation (looses nprof priority! But non-overlapping anyway)
            %[~,ind]=sort(x); % monotonic sorting. Commented out 02/12/24 JP

            % make monotonic & unique for interpolation (ss4003, cycle 29
            % breaks the above commented out line -has repeated surface 
            % values). JP 02/12/24
            [~,ind,~] = unique(x); 

            x=x(ind);y=y(ind);yqc=yqc(ind);
            % and make sure that they are column vectors
            x=x(:);y=y(:);yqc=yqc(:);

            % copy data for levels that are part of the synthetic pressure axis:
            % get indices to which data to copy
            % take only 'first' occurence (in nprof-priority sorted record!) and
            % toss away repeated occurences (e.g., MD5904767_004.nc PSAL @ 46.0 dbar)
            [~,fillind,cpind]=intersect(spresaxis,x);
            synth.(ubgcparams{i}).value(fillind)=y(cpind);
            synth.([ubgcparams{i} '_QC']).value(fillind)=yqc(cpind);
            synth.([ubgcparams{i} '_dPRES']).value(fillind)=0;
            %if ismember(ubgcparams{i},{'PRES';'TEMP';'PSAL'}) % keep which BGC PTS were truly measured for later
            %    full.(ubgcparams{i}).fillind=fillind;
            %end
            clear fillind cpind

            % rest of data:
            % interpolate data for other levels of the synthetic pressure axis:
            % toss away repeated occurence of pressures: e.g., MD5904767_004.nc PSAL @ 46.0 dbar
            if ismember(ubgcparams{i},{'PRES'})
                [~,uind]=unique(x);
            else
                uind=1:length(x); % should have been dealt with by overlapping portions already
            end
            if length(x(uind))>1
                % data interpolation: linear, no extrapolation
                synth.(ubgcparams{i}).value(isnan(synth.(ubgcparams{i}).value))=interp1(double(x(uind)),double(y(uind)),double(spresaxis(isnan(synth.(ubgcparams{i}).value))),'linear',NaN);
                % data extrapolation: nearest-neighbour, no limit on extrapolation
                if ~ismember(ubgcparams{i},{'PRES'}) % PRES extrapolation as nearest-neighbour does not work..
                    synth.(ubgcparams{i}).value(isnan(synth.(ubgcparams{i}).value))=interp1(double(x(uind)),double(y(uind)),double(spresaxis(isnan(synth.(ubgcparams{i}).value))),'nearest','extrap');
                end % should not occur for PRES field, either..
                % deal with qc
                % qc interpolation: next and previous, no extrapolation
                qcnext=interp1(double(x(uind)),yqc(uind),double(spresaxis),'next',NaN);
                qcprevious=interp1(double(x(uind)),yqc(uind),double(spresaxis),'previous',NaN);
                % take maximum of QC; order 1 < 2 < 5 < 3 < 4
                qcnext(qcnext==5)=2.5; qcprevious(qcprevious==5)=2.5; % replace QC 5 with 2.5:
                qcfill=max(qcnext,qcprevious); % max for interpolated QC
                synth.([ubgcparams{i} '_QC']).value(isnan(synth.([ubgcparams{i} '_QC']).value))=qcfill(isnan(synth.([ubgcparams{i} '_QC']).value));
                synth.([ubgcparams{i} '_QC']).value(synth.([ubgcparams{i} '_QC']).value==2.5)=5; % and reverse QC 5
                % qc extrapolation: nearest-neighbour, no limit on extrapolation
                synth.([ubgcparams{i} '_QC']).value(isnan(synth.([ubgcparams{i} '_QC']).value))=interp1(double(x(uind)),yqc(uind),double(spresaxis(isnan(synth.([ubgcparams{i} '_QC']).value))),'nearest','extrap');
                clear qcnext qcprevious qcfill
            else % only one value, keep this value as well as its QC, and place it closest to the original pressure
                [~,ifill]=min(abs(spresaxis-x(uind)));
                synth.(ubgcparams{i}).value(ifill)=y(uind);
                synth.([ubgcparams{i} '_QC']).value(ifill)=yqc(uind);
                if addTSeverywhere && ismember(ubgcparams{i},{'TEMP';'PSAL'})
                    % extrapolate T and S context by replication, no limit on extrapolation
                    synth.(ubgcparams{i}).value(:)=y(uind);
                    synth.([ubgcparams{i} '_QC']).value(:)=yqc(uind);
                end
            end
            % and kick out unmatched data
            %inomatch=false(nos,1); inomatch(setdiff(isynth,dsearchn(synth.PRES.value,x)))=1;
            % keep nearest data (can be two points) and remove all other that
            % are further away
            inomatch=true(length(spresaxis),1);
            if length(spresaxis)==1 % just one value on synthetic pressure axis
                %neardp=abs(synth.PRES.value-x(uind)); % on original sampling axis
                neardp=abs(spresaxis-x(uind)); % on original sampling axis
            else % more than one value on synthetic pressure axis
                %neardp=abs(synth.PRES.value(dsearchn(synth.PRES.value/1e3,x(uind)/1e3))-x(uind)); % on original sampling axis
                %neardp=abs(spresaxis(dsearchn(double(spresaxis)/1e3,double(x(uind))/1e3))-x(uind)); % on original sampling axis
                neardp=abs(spresaxis(dsearchn(double(spresaxis),double(x(uind))))-x(uind)); %JP
            end
            for k=1:length(uind)
                %inomatch(abs(synth.PRES.value-x(uind(k)))==neardp(k))=0;
                inomatch(abs(spresaxis-x(uind(k)))==neardp(k))=0;
            end
            % decide which data to really remove
            iremove=inomatch;
            % keep all interpolated core-data within [max(pres)-2 dbar float length; min(pres)+1 dbar antenna length]
            if addTSeverywhere && ismember(ubgcparams{i},{'PRES';'TEMP';'PSAL'})
                if ismember(ubgcparams{i},{'PRES'})
                    iremove(:)=0; % don't remove any interpolated/extrapolated PRES
                else % TEMP, PSAL
                    % if there are valid data before or after, keep those in between, too
                    %iremove(cumsum(~iremove,'forward')>0 & cumsum(~iremove,'reverse')>0)=0;
                    % keep all within pressure range: extrapolate deepest data
                    % point up to 1 float length (2 dbar) deeper and shallowest
                    % data point up to 1 antenna length (1 dbar) shallower
                    %iremove(spresaxis<=max(x)+2*1000 & spresaxis>=min(x)-1*1000)=0;
                    iremove(spresaxis<=max(x)+2 & spresaxis>=min(x)-1)=0; % JP
                end
            end

            % keep an isolated hole between data
            % standard case: 1 isolated hole between two data before and after
            if nos>=5
                iremove(find(inomatch(3:nos-2) & ~inomatch(1:nos-4) & ~inomatch(2:nos-3) & ~inomatch(4:nos-1) & ~inomatch(5:nos))+2)=0;
            end
            if nos>=4
                % special case: 'data point # 2' and 'data point # end-1'
                iremove(find(inomatch(2) & ~inomatch(1) & ~inomatch(3) & ~inomatch(4))+1)=0;
                iremove(find(inomatch(nos-1) & ~inomatch(nos-3) & ~inomatch(nos-2) & ~inomatch(nos))+nos-2)=0;
            end
            if nos>2
                % special case: shallowest data point: up to 1 float length (2
                % dbar) and deepest data point: up to 1 antenna length (1 dbar)
                %iremove(find(inomatch(1) & ~inomatch(2) & ~inomatch(3) & diff(synth.PRES.value([1 2]))<=2*1000))=0;
                %iremove(find(inomatch(nos) & ~inomatch(nos-2) & ~inomatch(nos-1) & diff(synth.PRES.value([nos-1 nos]))<=1*1000)+nos-1)=0;
                %iremove(find(inomatch(1) & ~inomatch(2) & ~inomatch(3) & diff(spresaxis([1 2]))<=2*1000))=0;
                iremove(find(inomatch(1) & ~inomatch(2) & ~inomatch(3) & diff(spresaxis([1 2]))<=2))=0;
                %iremove(find(inomatch(nos) & ~inomatch(nos-2) & ~inomatch(nos-1) & diff(spresaxis([nos-1 nos]))<=1*1000)+nos-1)=0;
                iremove(find(inomatch(nos) & ~inomatch(nos-2) & ~inomatch(nos-1) & diff(spresaxis([nos-1 nos]))<=1)+nos-1)=0; %JP
            end

            if ismember(ubgcparams{i},{'PRES'})
                %synth.(ubgcparams{i}).value(iremove)=FV;
            else
                synth.(ubgcparams{i}).value(iremove)=NaN;
                synth.([ubgcparams{i} '_QC']).value(iremove)=NaN;
                %synth.([ubgcparams{i} '_dPRES']).value(iremove)=FV;
                % assign QC flag to gap in a series points
                %synth.([ubgcparams{i} '_QC']).value(inomatch & ~iremove)=8;
                %fillind=inomatch & ~iremove & ~ismember(synth.([ubgcparams{i} '_QC']).value,[3 4]);
                fillind=inomatch & ~iremove & ~(ismember(synth.([ubgcparams{i} '_QC']).value,[3 4]) | isnan(synth.(ubgcparams{i} ).value));
                synth.([ubgcparams{i} '_QC']).value(fillind)=8;
            end % PRES (int32) or not (double)
            clear y yqc uind inomatch iremove fillind

            %             if strcmp(ubgcparams{i},'NITRATE'), pause,end % TESTING
            %             %if strcmp(ubgcparams{i},'PRES'), pause,end  % TESTING

            % check dPRES assignment
            % and add pressure difference where it is missing (and there are data)
            %fillind=find(~isnan(synth.([ubgcparams{i} '_QC']).value) & synth.([ubgcparams{i} '_dPRES']).value==FV);
            fillind=find(~isnan(synth.([ubgcparams{i} '_QC']).value) & isnan(synth.([ubgcparams{i} '_dPRES']).value));

            % if there are two nearest samples, take the deeper one..
            for k=1:length(fillind)
                dp=x-spresaxis(fillind(k)); % deeper sample of minimum pressure difference
                synth.([ubgcparams{i} '_dPRES']).value(fillind(k))=dp(find(abs(dp)==min(abs(dp)),1,'last'));

                %%synth.([ubgcparams{i} '_dPRES']).value(synth.([ubgcparams{i} '_dPRES']).value==FV)=x(uind(nearind(synth.([ubgcparams{i} '_dPRES']).value==FV)))-spresaxis(synth.([ubgcparams{i} '_dPRES']).value==FV);
                %dpind=abs(x-spresaxis(k))==abs(x(uind(nearind(k)))-spresaxis(k)); % points with same distance and to be used (not in low priority nprofs)
                %synth.([ubgcparams{i} '_dPRES']).value(k)=max(x(dpind))-spresaxis(k); % take deeper value if there are any two samples at +/- the same distance
                %clear dpind
            end
            clear x xqc fillind dp

            if ~isempty(yadj) % not only NaN/FV adjusted data
                % use only non-overlapping portion
                xadj=xadj(overlapadj);yadj=yadj(overlapadj);yadjqc=yadjqc(overlapadj);yadjerr=yadjerr(overlapadj);
                clear overlapadj
                % make monotonic for interpolation (looses nprof priority!)
                %[~,ind]=sort(xadj); % monotonic sorting, Commented out 02/12/24 JP

                % make monotonic & unique for interpolation (ss4003, cycle 29
                % breaks the above commented out line -has repeated surface
                % values). JP 02/12/24
                [~,ind,~] = unique(xadj);

                xadj=xadj(ind);yadj=yadj(ind);yadjqc=yadjqc(ind);yadjerr=yadjerr(ind);
                % and make sure that they are column vectors
                xadj=xadj(:);yadj=yadj(:);yadjqc=yadjqc(:);yadjerr=yadjerr(:);

                % do the same with adjusted fields
                % copy data for levels that are part of the synthetic pressure axis:
                [~,fillindadj,cpindadj]=intersect(spresaxis,xadj);
                synth.([ubgcparams{i} '_ADJUSTED']).value(fillindadj)=yadj(cpindadj);
                synth.([ubgcparams{i} '_ADJUSTED_QC']).value(fillindadj)=yadjqc(cpindadj);
                synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value(fillindadj)=yadjerr(cpindadj);
                %if ismember(ubgcparams{i},{'PRES';'TEMP';'PSAL'}) % keep which BGC PTS were truly measured for later
                %    full.(ubgcparams{i}).fillindadj=fillindadj;
                %end
                clear fillindadj cpindadj

                % rest of data:
                % interpolate data for other levels of the synthetic pressure axis:
                % toss away repeated occurence of pressures: e.g., MD5904767_004.nc PSAL @ 46.0 dbar
                if ismember(ubgcparams{i},{'PRES'})
                    [~,uindadj]=unique(xadj);
                else
                    uindadj=1:length(xadj); % should have been dealt with by overlapping portions already
                end
                if length(xadj(uindadj))>1
                    % data interpolation: linear, no extrapolation
                    synth.([ubgcparams{i} '_ADJUSTED']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED']).value))=interp1(double(xadj(uindadj)),double(yadj(uindadj)),double(spresaxis(isnan(synth.([ubgcparams{i} '_ADJUSTED']).value))),'linear',NaN);
                    % data extrapolation: nearest-neighbour, no limit on extrapolation
                    if ismember(ubgcparams{i},{'PRES'}) % PRES_ADJUSTED extrapolation as nearest-neighbour does not work.. use linear extrapolation
                        synth.([ubgcparams{i} '_ADJUSTED']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED']).value))=interp1(double(xadj(uindadj)),double(yadj(uindadj)),double(spresaxis(isnan(synth.([ubgcparams{i} '_ADJUSTED']).value))),'linear','extrap');
                    else % nearest-neighbour
                        synth.([ubgcparams{i} '_ADJUSTED']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED']).value))=interp1(double(xadj(uindadj)),double(yadj(uindadj)),double(spresaxis(isnan(synth.([ubgcparams{i} '_ADJUSTED']).value))),'nearest','extrap');
                    end
                    if yadjerrpresence % same with errors if any
                        synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value))=interp1(double(xadj(uindadj)),double(yadjerr(uindadj)),double(spresaxis(isnan(synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value))),'linear',NaN);
                        synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value))=interp1(double(xadj(uindadj)),double(yadjerr(uindadj)),double(spresaxis(isnan(synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value))),'nearest','extrap');
                    end
                    % deal with qc
                    % qc interpolation: next and previous, no extrapolation
                    qcnextadj=interp1(double(xadj(uindadj)),yadjqc(uindadj),double(spresaxis),'next',NaN);
                    qcpreviousadj=interp1(double(xadj(uindadj)),yadjqc(uindadj),double(spresaxis),'previous',NaN);
                    % take maximum of QC; order 1 < 2 < 5 < 3 < 4
                    qcnextadj(qcnextadj==5)=2.5; qcpreviousadj(qcpreviousadj==5)=2.5; % replace QC 5 with 2.5:
                    qcfilladj=max(qcnextadj,qcpreviousadj); % max for interpolated QC
                    synth.([ubgcparams{i} '_ADJUSTED_QC']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED_QC']).value))=qcfilladj(isnan(synth.([ubgcparams{i} '_ADJUSTED_QC']).value));
                    synth.([ubgcparams{i} '_ADJUSTED_QC']).value(synth.([ubgcparams{i} '_ADJUSTED_QC']).value==2.5)=5; % and reverse QC 5
                    % qc extrapolation: nearest-neighbour, no limit on extrapolation
                    synth.([ubgcparams{i} '_ADJUSTED_QC']).value(isnan(synth.([ubgcparams{i} '_ADJUSTED_QC']).value))=interp1(double(xadj(uindadj)),yadjqc(uindadj),double(spresaxis(isnan(synth.([ubgcparams{i} '_ADJUSTED_QC']).value))),'nearest','extrap');
                    clear qcnextadj qcpreviousadj qcfilladj
                elseif ~isempty(xadj) % only one value, keep this value, its error, and its QC, and place it closest to the original pressure
                    [~,ifill]=min(abs(spresaxis-xadj));
                    synth.([ubgcparams{i} '_ADJUSTED']).value(ifill)=yadj(uindadj);
                    if yadjerrpresence % same with errors if any
                        synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value(ifill)=yadjerr(uindadj);
                    end
                    synth.([ubgcparams{i} '_ADJUSTED_QC']).value(ifill)=yadjqc(uindadj);
                else % no adjusted value left from overlapping portions
                    % do nothing
                end
                % and kick out unmatched data
                inomatch=true(length(spresaxis),1);
                %neardp=abs(synth.PRES.value(dsearchn(synth.PRES.value/1e3,xadj(uindadj)/1e3))-xadj(uindadj)); % on original sampling axis
                %neardp=abs(spresaxis(dsearchn(double(spresaxis)/1e3,double(xadj(uindadj))/1e3))-xadj(uindadj));
                neardp=abs(spresaxis(dsearchn(double(spresaxis),double(xadj(uindadj))))-xadj(uindadj)); %JP
                % on original sampling axis
                for k=1:length(uindadj)
                    %inomatch(abs(synth.PRES.value-xadj(uindadj(k)))==neardp(k))=0;
                    inomatch(abs(spresaxis-xadj(uindadj(k)))==neardp(k))=0;
                end
                % decide which data to really remove
                iremove=inomatch;
                % keep all interpolated core-data within [max(pres)-2 dbar float length; min(pres)+1 dbar antenna length]
                if addTSeverywhere && ismember(ubgcparams{i},{'PRES';'TEMP';'PSAL'})
                    if ismember(ubgcparams{i},{'PRES'})
                        iremove(:)=0; % don't remove any interpolated/extrapolated PRES_ADJUSTED
                    else % TEMP, PSAL
                        % if there are valid data before or after, keep those in between, too
                        %iremove(cumsum(~iremove,'forward')>0 & cumsum(~iremove,'reverse')>0)=0;
                        % keep all within pressure range: extrapolate deepest data
                        % point up to 1 float length (2 dbar) deeper and shallowest
                        % data point up to 1 antenna length (1 dbar) shallower
                        %iremove(spresaxis<=max(xadj)+2*1000 & spresaxis>=min(xadj)-1*1000)=0;
                        iremove(spresaxis<=max(xadj)+2 & spresaxis>=min(xadj)-1)=0; % JP
                    end
                end

                % but only if it's not an isolated hole between data
                % standard case: 1 isolated hole between two data before and after
                if nos>=5
                    iremove(find(inomatch(3:nos-2) & ~inomatch(1:nos-4) & ~inomatch(2:nos-3) & ~inomatch(4:nos-1) & ~inomatch(5:nos))+2)=0;
                end
                if nos>=4
                    % special case: 'data point # 2' and 'data point # end-1'
                    iremove(find(inomatch(2) & ~inomatch(1) & ~inomatch(3) & ~inomatch(4))+1)=0;
                    iremove(find(inomatch(nos-1) & ~inomatch(nos-3) & ~inomatch(nos-2) & ~inomatch(nos))+nos-2)=0;
                end
                if nos>2
                    % special case: shallowest data point: up to 1 float length (2
                    % dbar) and deepest data point: up to 1 antenna length (1 dbar)
                    %iremove(find(inomatch(1) & ~inomatch(2) & ~inomatch(3) & diff(synth.PRES.value([1 2]))<=2*1000))=0;
                    %iremove(find(inomatch(nos) & ~inomatch(nos-2) & ~inomatch(nos-1) & diff(synth.PRES.value([nos-1 nos]))<=1*1000)+nos-1)=0;
                    % iremove(find(inomatch(1) & ~inomatch(2) & ~inomatch(3) & diff(spresaxis([1 2]))<=2*1000))=0;
                    % iremove(find(inomatch(nos) & ~inomatch(nos-2) & ~inomatch(nos-1) & diff(spresaxis([nos-1 nos]))<=1*1000)+nos-1)=0;
                    iremove(find(inomatch(1) & ~inomatch(2) & ~inomatch(3) & diff(spresaxis([1 2]))<=2))=0;
                    iremove(find(inomatch(nos) & ~inomatch(nos-2) & ~inomatch(nos-1) & diff(spresaxis([nos-1 nos]))<=1)+nos-1)=0;
                end
                %if ismember(ubgcparams{i},{'PRES'})
                %    synth.(ubgcparams{i}).value(iremove)=FV;
                %else
                synth.([ubgcparams{i} '_ADJUSTED']).value(iremove)=NaN;
                %end % PRES (int32) or not (double)
                synth.([ubgcparams{i} '_ADJUSTED_QC']).value(iremove)=NaN;
                % assign QC flag to gap in a series points
                %synth.([ubgcparams{i} '_ADJUSTED_QC']).value(inomatch & ~iremove)=8;
                %fillind=inomatch & ~iremove & ~ismember(synth.([ubgcparams{i} '_ADJUSTED_QC']).value,[3 4]);
                fillind=inomatch & ~iremove & ~(ismember(synth.([ubgcparams{i} '_ADJUSTED_QC']).value,[3 4]) | isnan(synth.([ubgcparams{i} '_ADJUSTED']).value));
                if ~ismember(ubgcparams{i},{'PRES'}) % don't tap with '8' for PRES_ADJUSTED: location, not really data
                    synth.([ubgcparams{i} '_ADJUSTED_QC']).value(fillind)=8;
                end
                synth.([ubgcparams{i} '_ADJUSTED_ERROR']).value(iremove)=NaN;
                clear xadj xadjqc yadj yadjqc yadjerr yadjerrpresence uindadj inomatch fillind
            end % ~isempty(yadj) : some adjusted data available
        end % ~isempty(y) : some data available
        % and convert parameter dPRES back to proper 1/1000 double
        %fillind=synth.([ubgcparams{i} '_dPRES']).value==FV;
        fillind=isnan(synth.([ubgcparams{i} '_dPRES']).value);
        %synth.([ubgcparams{i} '_dPRES']).value=double(synth.([ubgcparams{i} '_dPRES']).value)./1000;
        synth.([ubgcparams{i} '_dPRES']).value=double(synth.([ubgcparams{i} '_dPRES']).value);
        synth.([ubgcparams{i} '_dPRES']).value(fillind)=NaN;
        clear fillind
    end % fill in data in synthetic profile
end



% add synthetic profile to merged data with full dimension, BGC levels interleaved
fnames=fieldnames(synth); % use non-interleaved profile as template; fill gaps with NaN
%for i=2:length(fnames) % start at index 2: first one is PRES=presaxis -> already present
for i=1:length(fnames) % start at index 1: first one is PRES=presmerge -> copy all to be sure
    synthfull.(fnames{i}).value=ones(nosf,1)*NaN;
    if ismember(fnames{i},{'PRES';'TEMP';'PSAL';'PRES_QC';'TEMP_QC';'PSAL_QC';'PRES_ADJUSTED';'TEMP_ADJUSTED';'PSAL_ADJUSTED';'PRES_ADJUSTED_QC';'TEMP_ADJUSTED_QC';'PSAL_ADJUSTED_QC';'PRES_ADJUSTED_ERROR';'TEMP_ADJUSTED_ERROR';'PSAL_ADJUSTED_ERROR';'PRES_dPRES';'TEMP_dPRES';'PSAL_dPRES'}) % core already on full axis
        synthfull.(fnames{i}).value=synth.(fnames{i}).value;
    else % BGC on synthetic pressure axis
        synthfull.(fnames{i}).value(synthind)=synth.(fnames{i}).value;
    end % already on full axis (core)?
end % cycle fields

% and convert PRES back to proper 1/1000 double
% fillind=synthfull.PRES.value==FV;
% synthfull.PRES.value=double(synthfull.PRES.value)./1000;
% synthfull.PRES.value(fillind)=NaN;

fillind=isnan(synthfull.PRES.value); % JP
synthfull.PRES.value=double(synthfull.PRES.value); % JP
synthfull.PRES.value(fillind)=NaN; % JP

% sanity-check ADJUSTED_QC in case of ADJUSTED=FV
for i=1:length(ubgcparams)
    if ~all(isnan(synthfull.([ubgcparams{i} '_ADJUSTED']).value)) && ~all(isnan(synthfull.([ubgcparams{i} '_ADJUSTED_QC']).value))
        fillind=~isnan(synthfull.(ubgcparams{i}).value) & isnan(synthfull.([ubgcparams{i} '_ADJUSTED']).value) & isnan(synthfull.([ubgcparams{i} '_ADJUSTED_QC']).value);
        synthfull.([ubgcparams{i} '_ADJUSTED_QC']).value(fillind)=4;
    end % adjusted data or adjusted_QC are present?
end

% LASTLY ADD INFO STURCTURE FOR MBARI ODV CREATION
INFO.sdn              = S.INFO.sdn;
INFO.cast             = S.INFO.cast;
INFO.gps              = S.INFO.gps;
INFO.INST_ID          = S.INFO.INST_ID;
INFO.WMO_ID           = S.INFO.WMO_ID;
INFO.float_type       = S.INFO.float_type;
INFO.file_name        = S.INFO.file_name;
% INFO.WMO_ID           = S.INFO.WMO_ID;
INFO.BOP_sdn          = S.INFO.BOP_sdn;
INFO.EOP_sdn          = S.INFO.EOP_sdn;
INFO.sensors          = S.INFO.sensors;
INFO.pres_axes_ct     = S.INFO.pres_axes_ct;
INFO.bgc_pres_axes_ct = S.INFO.bgc_pres_axes_ct;

synthfull.INFO        = INFO;

clearvars -except synthfull





return

