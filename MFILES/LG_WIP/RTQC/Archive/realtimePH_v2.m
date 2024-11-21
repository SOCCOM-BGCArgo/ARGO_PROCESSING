warning("off")

%IDENTIFY RECENT FLOATS

%Load table of floats
floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");

%Sort any dead floats and those that haven't been deployed yet (Shouldn't be getting new data from those anyways)
floatList = floatList(floatList.tfDead==0 & ~isnat(floatList.x1stDate),:);
%floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "descend");
%newFloats = floatList(datenum(floatList.maxCycleProcDate)>=datenum(now-hours(24)),:);
floatList = sortrows(floatList(~ismissing(floatList.latestCycleFileDate),:), "latestCycleFileDate", "descend");
newFloats = floatList(datenum(floatList.latestCycleFileDate)<=datenum(now),:);

wmoList = newFloats.WMO;

%Open up a qc log text file to store qc messages
userDir = getenv("USERPROFILE");
qcLog = fopen(userDir+"\Documents\realtime_qc_log.txt","w");

%Set QC parameters
deltaThd = .01;
refAnomThd = .0075;
phMagThd = [7.3; 8.5];

%Write all parameters at the top of the file, %.3g limits numbers to 3 sig figs
fprintf(qcLog, "Cycle delta threshold = %.3g\n" + ...
    "Reference anomaly threshold = %.3g\n" + ...
    "pH magnitude bounds = %.3g-%.3g\n\n", ...
    deltaThd, refAnomThd, phMagThd(1), phMagThd(2));

for c = 1:length(wmoList)
    
    %For keeping track of issues
    disp(wmoList(c))

    %SEARCH FOR QC FLOATVIZ FILE
    if exist("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+wmoList(c)+"QC.txt", 'file')
        d = get_FloatViz_data("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+wmoList(c)+"QC.txt");
    
    %LOAD NON-QC FILE IF IT'S NOT AVAILABLE
    elseif  ~exist("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+wmoList(c)+"QC.txt", 'file') & exist("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\"+wmoList(c)+".txt", 'file')
        disp("No QC floatviz file, loading non-processed file")
        d = get_FloatViz_data("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\"+wmoList(c)+".txt");
    else
        disp("No floatviz files exist for this float, check float arrays to check status")
        continue;
    end
    
    %CHECK IF FLOAT HAS PH (MOVE TO NEXT FLOAT IF NO PH)
    if any(contains(d.hdr,'pHinsitu[Total]'))==0
        disp("Float has no pH")
        continue;
    end
    %CHECK IF FLOAT PH IS MARKED INDEFINITELY ON BSL

    %Load the BSL as a table
    BSL = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT");
    BSL = sortrows(BSL,"WMO_","ascend"); %Sort by WMO's in ascending order to match floatsDir
    BSL= BSL((contains(BSL.SENSOR,"PH")&BSL.FLAG==4),:); %Isolate lines with PH
    
    %Create new column with all bad cycle numbers as doubles
    BSL.BADCYCLENUM = cellfun(@str2double, regexp(BSL.CYCLES,'\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
    
    %IF ALL FUTURE CYCLES ARE ON BSL, MOVE ONTO NEXT FLOAT
    if ~isempty(BSL.BADCYCLENUM(BSL.WMO_==wmoList(c)))
        disp("All future cycles already listed on BSL")
        continue;
    end %CONTINUE WITH QC
        
    %CHECK THAT A FLOATQCLIST FILE EXISTS TO OBTAIN START CYCLE
    if exist("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+wmoList(c)+"_FloatQCList.txt",'file')

        %Identify the last cycle where qc was performed by looking at the floatQClist file
        fid = fopen("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+wmoList(c)+"_FloatQCList.txt", 'r');
        
        %Create empty string array
        lines = {};
        
        %Get the first line of the floatQClist file
        next_line = fgetl(fid);
        
        %Iterate through the file until reaching a blank line, indicating
        %the end of the most recent text block
        while ~isempty(next_line)
            next_line = fgetl(fid); %Read the current line
            if isempty(next_line)|next_line==-1 %If the current line is empty, end the while loop
                break;
            else
                lines{end+1} = split(next_line,',');  %Split the string and store it
            end
        end
        
        %The cycle that the qc process will begin on
        startCycle = str2double(lines{end}{2});

    %IF NO FLOATQCLIST FILE EXISTS, START FROM CYCLE 1
    else
        startCycle = 1;
    end

    %The end cycle is the latestCycleMsgFile on the MBARI float list,
    %being used for diagnostic purposes but relatively useless
    % endCycle = newFloats.latestCycleMsgFile(c);
    % disp(string(endCycle-startCycle)+" profiles since last changepoint")
    
    % GET SOME RAW INDICES
    DATA.iStn  = find(strcmp('Station', d.hdr) == 1);
    DATA.iP    = find(strcmp('Pressure[dbar]', d.hdr)  == 1);
    DATA.iT    = find(strcmp('Temperature[°C]', d.hdr)  == 1);
    DATA.iS    = find(strcmp('Salinity[pss]', d.hdr)  == 1);
    DATA.iZ    = find(strcmp('Depth[m]', d.hdr)  == 1);
    DATA.iO    = find(strcmp('Oxygen[µmol/kg]', d.hdr)  == 1);
    DATA.iN    = find(strcmp('Nitrate[µmol/kg]', d.hdr) == 1);
    DATA.iPH   = find(strcmp('pHinsitu[Total]', d.hdr)  == 1);
    
    endCycle = d.data(end,DATA.iStn);

    if startCycle == endCycle
        disp('No new cycles to process')
        continue;
    end

    %Subsample the deep data so that only the new cycles after previous QC (which we
    %assume are correct) are analyzed
    % dAll = d.data(d.data(:,DATA.iStn)>=startCycle & d.data(:,DATA.iStn)<=endCycle & d.data(:,DATA.iPH)>-1e10,:);
    dAll = d.data(d.data(:,DATA.iStn)>=startCycle & d.data(:,DATA.iPH)>-1e10,:);

    %CHECK IF FLOAT IS ON PUMP OFFSET LIST
    pH_pumpoffset_980 %Call the list of pump offset floats
    
    %IF IT IS, SUBSAMPLE DEEP PH TO 950m
    if ismember(wmoList(c),pH_pumpoffset_980_floats)

        disp("Pump offset float")

        %Subsample the deep data to 950-980m
        dDeep = dAll(dAll(:,DATA.iZ)>=920 & dAll(:,DATA.iZ)<=980 & dAll(:,DATA.iPH)>-1e10,:);

    else %IF NOT THEN SUBSAMPLE DEEP PH FROM 980-1520m

        %Subsample the deep data to 980-1520m
        dDeep = dAll(dAll(:,DATA.iZ)>=980 & dAll(:,DATA.iZ)<=1520 & dAll(:,DATA.iPH)>-1e10,:);
    end
    
    if isempty(dDeep)
        disp('No recent deep data for this float')
        %LOGICAL STRING CODE FROM "JAN" AT: 
        % https://www.mathworks.com/matlabcentral/answers/258287-fprintf-for-logical-statement
        pumpOffset = ismember(wmoList(c),pH_pumpoffset_980_floats);
        LogicalStr = {"", "(pump offset float)"};
        fprintf(qcLog, "\n%d\t%s\tcycles:\t%d-%d\n" + ...
            "No data from 950-1520m, no qc performed\n", ...
            wmoList(c), LogicalStr{pumpOffset + 1}, startCycle, endCycle);
        continue;
    end


    %CHECK FOR GPS ISSUES
    if any(dDeep(abs(dDeep(:,3))>360|abs(dDeep(:,4))>90))
  
        disp('Missing gps locations for some stations')
        
        %List of stns without a gps fix
        noGPSstns = unique(dDeep(abs(dDeep(:,3))>360|abs(dDeep(:,4))>90,DATA.iStn));
        
        %Index of data with a gps fix
        gpsIDX = abs(dDeep(:,3))<360 & abs(dDeep(:,4))<90;
        
        %Subset data that has a gps fix to calculate Esper
        dDeep_noGPS = dDeep(gpsIDX,:);
        
        %Create an empty array of NaNs to fill with Esper data
        newEstimates = NaN(length(dDeep),1);
        
        % CALCULATE ESPER pH & NO3
        DesireVar = [3]; %pH & nitrate (5)
        OutCoords = [dDeep_noGPS(:,3), dDeep_noGPS(:,4), dDeep_noGPS(:,DATA.iZ)]; % lon,lat,Z
        PredictorTypes    = [1 2 6]; % PSAL, TEMP, DOXY_ADJ,
        Measurements = [dDeep_noGPS(:,DATA.iS), dDeep_noGPS(:,DATA.iT), dDeep_noGPS(:,DATA.iO)];
        dvec = datevec(dDeep_noGPS(:,1));
        dec_yr = dvec(:,1) +(dvec(:,2)*30)/365 + dvec(:,3)/365; % very crude decimal year
        Equations    = 7; % S, T, O2
        
        [Estimates,~] = ESPER_Mixed(DesireVar, OutCoords, Measurements,...
                        PredictorTypes, 'Equations', Equations, 'EstDates', dec_yr);
        
        %Fills the NaN array with Esper estimates where GPS fixes exist available
        newEstimates(gpsIDX,1) = Estimates.pH ;
        
        %Replace Estimates.pH with the modified NaN array
        Estimates.pH = newEstimates;
     
    else

        % CALCULATE ESPER pH & NO3
        DesireVar = [3]; %pH & nitrate (5)
        OutCoords = [dDeep(:,3), dDeep(:,4), dDeep(:,DATA.iZ)]; % lon,lat,Z
        PredictorTypes    = [1 2 6]; % PSAL, TEMP, DOXY_ADJ,
        Measurements = [dDeep(:,DATA.iS), dDeep(:,DATA.iT), dDeep(:,DATA.iO)];
        dvec = datevec(dDeep(:,1));
        dec_yr = dvec(:,1) +(dvec(:,2)*30)/365 + dvec(:,3)/365; % very crude decimal year
        Equations    = 7; % S, T, O2
        
        [Estimates,~] = ESPER_Mixed(DesireVar, OutCoords, Measurements,...
                        PredictorTypes, 'Equations', Equations, 'EstDates', dec_yr);
    end
    
    %Create array of unique stations in deep PH dataset
    stns = unique(dDeep(:,DATA.iStn));
    
    %Retrieve indices of the first cycle to start from
    stnInd = dDeep(:,DATA.iStn)==stns(1);
    
    %Mean pH over 980-1520m
    deepPH = NaN(length(stns),1);
    
    %Mean Esper at depth
    EsperPH = NaN(length(stns),1);

    %Mean Float-to-Esper anomaly at depth
    refAnom = NaN(length(stns),1);

    %Test for anomalously high or low pH
    profileMagTest = NaN(length(stns),1); 

    for i = 1:length(stns)
        
        %Assign indices for the current station
        stnInd = dDeep(:,DATA.iStn)==stns(i);
        
        %Find pH at all deep depths
        deepPH_allPoints = dDeep(stnInd,DATA.iPH);
        
        %Find EsperMix pH at depth
        deepEsper = Estimates.pH(stnInd);
        
        %Average. float pH at depth that can be used to calculate the
        %cycle-to-cycle difference in pH
        deepPH(i,1) = mean(deepPH_allPoints);
    
        %Avg. EsperMix pH at depth
        EsperPH(i,1) = mean(deepEsper);
    
        %Avg. anomaly between deep float and EsperMix pH
        refAnom(i,1) = mean(deepPH_allPoints-deepEsper);
        
        %Test whether any part of the profile is anomalously high or low
        profileMagTest(i,1) = any(dAll(dAll(:,DATA.iStn)==stns(i),DATA.iPH) <= phMagThd(1)|dAll(dAll(:,DATA.iStn)==stns(i),DATA.iPH) >= phMagThd(2));
    end
    
    %Append phDelta with a 0 at the beginning to compensate for indices offset
    phDelta = [0; diff(deepPH)];

    %BUILD REALTIME QC LOG FILE

    if exist("noGPSstns", 'var')||any(phDelta>=deltaThd)||any(refAnom>=refAnomThd)||any(profileMagTest==1)
        %LOGICAL STRING CODE FROM "JAN" AT: 
        % https://www.mathworks.com/matlabcentral/answers/258287-fprintf-for-logical-statement
        pumpOffset = ismember(wmoList(c),pH_pumpoffset_980_floats);
        LogicalStr = {"", "(pump offset float)"};
        fprintf(qcLog, "\n%d\t%s\tcycles:\t%d-%d\n", wmoList(c), LogicalStr{pumpOffset + 1}, startCycle, endCycle);
    end

    if exist("noGPSstns", 'var')
        horizontal_str = sprintf('%d\t', noGPSstns);
        fprintf(qcLog, "No GPS fix (and no EsperMix) for cycle(s)\t%s\n", horizontal_str);
        clear noGPSstns
    end

    
    if any(phDelta>=deltaThd)
        horizontal_str = sprintf('%d\t', stns(phDelta>=deltaThd));
        fprintf(qcLog, "Cycle delta threshold exceeded at cycle(s)\t%s\n", horizontal_str);
    end
    
    if any(refAnom>=refAnomThd)
        horizontal_str = sprintf('%d\t', stns(refAnom>=refAnomThd));
        fprintf(qcLog, "Reference anomaly threshold exceeded at cycle(s)\t%s\n", horizontal_str);
    end
    
    if any(profileMagTest==1)
        horizontal_str = sprintf('%d\t', stns(any(profileMagTest==1)));
        fprintf(qcLog, "pH magnitude threshold exceeded at cycle(s)\t%s\n", horizontal_str);
    end

end

fclose(qcLog);