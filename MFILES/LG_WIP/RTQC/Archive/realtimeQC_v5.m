%===========================================================================================
% OBJECTIVE: This code is meant to be a daily autojob that looks over the past 24 hours of floats
% with new profiles and checks if the data is acceptable quality based on a cycle delta test and
% reference anomaly test.
%
% TESTS:
% - CYCLE DELTA TEST: Refers to the difference between consecutive cycles; if the difference
% between consecutive cycles is too high, it can be indicative of an overly noisy or failing
% sensor.
%
% - REFERENCE ANOMALY TEST: Refers to the difference between the adjusted float data and the
% estimated ESPER MIX reference model; if the anomaly between the data and the reference is 
% too high, it can be indicative of a drifting sensor that may need to be re-qc'd.
%
% - CONSECUTIVE OCCURENCES: The previously mentioned tests are often triggered by the properties
% of a dynamic region, biofouling, or a less accurate reference for an underrepresented region, 
% but not a failing or misbehaving sensor. These are regular artifacts that can often cause 1-2
% cycles of noise, but then the sensor reverts back to nominal. To try and mitigate these false
% positives, a consecutive occurence threshold is added to the cycle delta and referenace anomaly
% tests.
%
% 
%===========================================================================================

warning("off")

%===========================================================================================
%===================================== INPUTS ==============================================
%===========================================================================================
deltaThd_ph = .01; %Cycle-to-cycle delta for PH
refAnomThd_ph = .01; %Float-to-Esper anomaly for PH
deltaThd_nitrate = 1; %Cycle-to-cycle delta for NITRATE
refAnomThd_nitrate = 1; %Float-to-Esper anomaly for NITRATE
consecCounter = 3; %Number of consecutive test failures before recording the cycle

email_list = {'lgrady@mbari.org';'tmaurer@mbari.org';"jplant@mbari.org";"sbartoloni@mbari.org";"nicolag@mbari.org" };
subjectline = "TEST: RTQC DAILY UPDATE";
timelag = now-hours(24);
%===========================================================================================


%===========================================================================================
%===================================== LOAD DATA ===========================================
%===========================================================================================

%LOAD FLOATS
floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");

%Sort any dead floats and those that haven't been deployed yet (Shouldn't be getting new data from those anyways)
floatList = floatList(floatList.tfDead==0 & ~isnat(floatList.x1stDate),:);
floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "descend");
floatList = floatList(datenum(floatList.maxCycleProcDate)>=datenum(timelag),:);

%FOR TESTING ON THE WHOLE FLEET AND NOT PAST 24 HOURS
%floatList = floatList(datenum(floatList.maxCycleProcDate)<=datenum(now),:);

%Create list of float WMO's to call on
wmoList = floatList.WMO;
idList = floatList.MBARIID;
% wmoList = 5906533;
% idList = {'ua20728'};

%LOAD BSL AND IDENTIFY CYCLES FOR PH, P, T, AND S:
BSL = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT");
BSL = sortrows(BSL,"WMO_","ascend"); %Sort by WMO's in ascending order to match floatsDir

%Regular expression matches any string that's a number followed by a dash but NOT FOLLOWED by a dash 
%and a number, ensuring that only indefinite failures are accounted for
badCycles.PH_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"PH")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.PH_wmo = BSL(contains(BSL.SENSOR,"PH")&BSL.FLAG==4,"WMO_").WMO_;

%Do the same for all questionable cycles
badCycles.PH_cycles_3 = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"PH")&BSL.FLAG==3), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.PH_wmo_3 = BSL(contains(BSL.SENSOR,"PH")&BSL.FLAG==3,"WMO_").WMO_;

%NITRATE
badCycles.N_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"N")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.N_wmo = BSL(contains(BSL.SENSOR,"N")&BSL.FLAG==4,"WMO_").WMO_;

badCycles.N_cycles_3 = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"N")&BSL.FLAG==3), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.N_wmo_3 = BSL(contains(BSL.SENSOR,"N")&BSL.FLAG==3,"WMO_").WMO_;

%OXYGEN
badCycles.O_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"O")&~contains(BSL.SENSOR,"CDOM")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.O_wmo = BSL(contains(BSL.SENSOR,"O")&~contains(BSL.SENSOR,"CDOM")&BSL.FLAG==4,"WMO_").WMO_;

%LOAD THE LIST OF PUMP OFFSET FLOATS:
pH_pumpoffset_980 

%LOAD THE QC LOG AND POPULATE WITH TEST PARAMETERS:
%Open up a qc log text file to store qc messages
userDir = getenv("USERPROFILE");
qcLog = fopen(userDir+"\Documents\rtqc_triggerLog.txt","w");
msgLog = fopen(userDir+"\Documents\rtqc_msgLog.txt","w");

%Write all parameters at the top of the file, %.3g limits numbers to 3 sig figs
fprintf(qcLog, "Nitrate cycle delta threshold = +/- %.3g\n" + ...
    "Nitrate reference anomaly threshold = +/- %.3g\n" + ...
    "PH cycle delta threshold = +/- %.3g\n" + ...
    "PH reference anomaly threshold = +/- %.3g\n" + ...
    "-------------------------------------------------------\n", ...
    deltaThd_nitrate, refAnomThd_nitrate, deltaThd_ph, refAnomThd_ph);

%LOAD THE LIST OF FLOATS THAT NEED INITIAL QC
%===========================================================================================
if ~isempty(check_1st_QC)
   
    no_QC_floats = str2num(cell2mat(check_1st_QC()));

    horizontal_str = sprintf('%d, ', no_QC_floats);
    fprintf(qcLog, "\nFLOATS THAT REQUIRE INITIAL QC: %s\n" + ...
        "\n-------------------------------------------------------\n", horizontal_str);
else
    no_QC_floats = [];
end
%===========================================================================================

%Print an example line to show the format of the warning messages
fprintf(qcLog, "\nQC MESSAGE FORMAT:\nWMO (MBARIID) (Other descriptors) | N or PH (qc depth): last qc'd cycle - current cycle |\n");

%===========================================================================================
%================================ BEGIN REALTIME QC LOOP ===================================
%===========================================================================================

for c = 1:length(wmoList)
    
    %For keeping track of issues from specific floats
    disp(wmoList(c))
    disp(idList(c))
    fprintf(msgLog, "\n%d\n", wmoList(c));

    %SEARCH FOR QC FLOATVIZ FILE
    %===========================================================================================
    if exist("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+wmoList(c)+"QC.txt", 'file')
        DATA = get_FloatViz_data("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+wmoList(c)+"QC.txt");

        %Convert all -1e10 points (missing) to NaN
        DATA.data(DATA.data==-1e10) = NaN;
   
    %CHECK IF FLOATVIZ FILE EXISTS AND IF ITS UP FOR INITIAL QC
    elseif ismember(wmoList(c),no_QC_floats)
        disp("Float requires initial QC, skipping for now")
        fprintf(msgLog, "Float requires initial QC, skipping for now\n");
        continue;
        
    %IF NO FLOATVIZ FILE EXISTS THEN MOVE ONTO NEXT FLOAT
    else
        disp("No QC floatviz file exists and it's not eligible for initial QC")
        fprintf(msgLog, "No QC floatviz file exists and it's not eligible for initial QC\n");
        continue;
    end
    %===========================================================================================

    %CHECK IF FLOAT HAS PH or NITRATE
    %===========================================================================================

    %Add the various test thresholds to the DATA structure
    DATA.deltaThd_ph = deltaThd_ph; %Cycle-to-cycle delta for PH
    DATA.refAnomThd_ph = refAnomThd_ph; %Float-to-Esper anomaly for PH
    DATA.deltaThd_nitrate = deltaThd_nitrate; %Cycle-to-cycle delta for NITRATE
    DATA.refAnomThd_nitrate = refAnomThd_nitrate; %Float-to-Esper anomaly for NITRATE
    DATA.consecCounter = consecCounter; %Number of consecutive test failures before recording the cycle

    % GET SOME RAW INDICES
    DATA.iStn  = find(strcmp('Station', DATA.hdr) == 1);
    DATA.iP    = find(strcmp('Pressure[dbar]', DATA.hdr)  == 1);
    DATA.iT    = find(strcmp('Temperature[°C]', DATA.hdr)  == 1);
    DATA.iS    = find(strcmp('Salinity[pss]', DATA.hdr)  == 1);
    DATA.iZ    = find(strcmp('Depth[m]', DATA.hdr)  == 1);
    DATA.iO    = find(strcmp('Oxygen[µmol/kg]', DATA.hdr)  == 1);
    DATA.iN    = find(strcmp('Nitrate[µmol/kg]', DATA.hdr) == 1);
    DATA.iPH   = find(strcmp('pHinsitu[Total]', DATA.hdr)  == 1);
    DATA.iLat  = find(strcmp('Lat [°N]', DATA.hdr) == 1);
    DATA.iLon  = find(strcmp('Lon [°E]', DATA.hdr) == 1);

    %Start a data struct for all of the logical tests
    tests.phColTest = any(contains(DATA.hdr,'pHinsitu[Total]'));
    tests.nitrateColTest = any(contains(DATA.hdr,'Nitrate[µmol/kg]'));
    
    %If float isn't equipped with either nitrate or ph, skip to next float
    if ~tests.nitrateColTest && ~tests.phColTest
        disp('No NITRATE or PH detected, skipping float')
        fprintf(msgLog, 'No NITRATE or PH detected, skipping float\n');
        continue;
    end
    %===========================================================================================

    %CHECK THAT A FLOATQCLIST FILE EXISTS TO OBTAIN START CYCLE
    %===========================================================================================
    if exist("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+wmoList(c)+"_FloatQCList.txt",'file')

        %Identify the last cycle where qc was performed by looking at the floatQClist file
        fid = fopen("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+wmoList(c)+"_FloatQCList.txt", 'r');
        
        %Create empty string array
        startCycle_ph = {};
        startCycle_nitrate = {};
    
        %Get the first line of the floatQClist file
        next_line = fgetl(fid);
        
        %Iterate through the file until reaching a blank line, indicating
        %the end of the most recent text block
        while ~isempty(next_line)
            next_line = fgetl(fid); %Read the current line
            
            if isempty(next_line)|next_line==-1 %If the current line is empty, end the while loop
                break;
            end
    
            if contains(next_line, 'Nitrate')
                startCycle_nitrate{end+1} = split(next_line,',');
            end
    
            if contains(next_line, 'pH')
                startCycle_ph{end+1} = split(next_line,',');  %Split the string and store it
            end
        end
        
        %The cycle that the qc process will begin on
        if ~isempty(startCycle_nitrate)
            startCycle_nitrate = str2double(startCycle_nitrate{end}{2});
        else
            startCycle_nitrate = 1;
        end
    
        if ~isempty(startCycle_ph)
            startCycle_ph = str2double(startCycle_ph{end}{2});
        else
            startCycle_ph = 1;
        end

    %IF NO FLOATQCLIST FILE EXISTS, START FROM CYCLE 1
    else
        startCycle_nitrate = 1;
        startCycle_ph = 1;
    end

    %CHECK IF THERE ARE ANY NEW CYCLES
    endCycle = DATA.data(end,DATA.iStn); %The final station in the dataset
    
    if (startCycle_nitrate == endCycle) && (startCycle_ph == endCycle)
        disp('No new cycles to process, skipping float')
        fprintf(msgLog,'No new cycles to process, skipping float\n');
        continue;
    end
    %===========================================================================================
  
    %CHECK IF PH OR NITRATE ARE MARKED INDEFINITELY ON BSL
    %===========================================================================================

    %PH
    BSLCycle_ph = cell2mat(badCycles.PH_cycles(badCycles.PH_wmo==wmoList(c)));
    %Check if BSL listing blacklists all future cycles, or only some of the new cycles
    
    %If BSL cycle is earlier than most recent cycle, sensor fails the test
    if BSLCycle_ph <= startCycle_ph + 1
        %Just a placeholder, populate the endCycle variable, but PH will no be processes since it fails this test
        endCycle_ph = BSLCycle_ph;
        startCycle_ph = BSLCycle_ph;
        tests.phBSLtest = 0;

    %If the BSL cycle is placed AFTER our last QC, there may still be cycles to process
    elseif BSLCycle_ph > startCycle_ph + 1
        disp('PH listed on BSL, but there are still unprocessed cycles')
        fprintf(msgLog,'PH listed on BSL, but there are still unprocessed cycles\n');
        endCycle_ph = BSLCycle_ph;
        tests.phBSLtest = 1;
    %If there's any other weirdness, just make the end cycle the last processed cycle from the floatlist
    else
        endCycle_ph = endCycle;
        tests.phBSLtest = 1;
    end
    
    %NITRATE (same steps as PH before)
    BSLCycle_nitrate = cell2mat(badCycles.N_cycles(badCycles.N_wmo==wmoList(c)));

    if BSLCycle_nitrate <= startCycle_nitrate + 1
        endCycle_nitrate = BSLCycle_nitrate;
        startCycle_nitrate = BSLCycle_nitrate;
        tests.nitrateBSLtest = 0;

    elseif BSLCycle_nitrate > startCycle_nitrate + 1
        disp('NITRATE listed on BSL, but there are still unprocessed cycles')
        fprintf(msgLog,'NITRATE listed on BSL, but there are still unprocessed cycles\n');
        endCycle_nitrate = BSLCycle_nitrate;
        tests.nitrateBSLtest = 1;
    else
        endCycle_nitrate = endCycle;
        tests.nitrateBSLtest = 1;
    end
    
    %If ph and nitrate are both indefinitely listed on BSL then skip the float
    if ~tests.phBSLtest && ~tests.nitrateBSLtest
        disp('PH and NITRATE indefinitely listed on BSL, skipping float')
        fprintf(msgLog,'PH and NITRATE indefinitely listed on BSL, skipping float\n');
        continue;
    end

    %Check OXYGEN on the BSL
    BSLCycle_oxygen = cell2mat(badCycles.O_cycles(badCycles.O_wmo==wmoList(c)));
    tests.badO2Test = ismember(wmoList(c),badCycles.O_wmo) & ~isempty(badCycles.O_cycles(BSLCycle_oxygen));
    if tests.badO2Test
        disp('O2 is indefinitely on BSL')
        fprintf(msgLog,'O2 is on BSL\n');
    end
    %===========================================================================================

    %CHECK IF ANY DATA IS INDEFINITELY MARKED AS QUESTIONABLE ON THE BSL
    %===========================================================================================
    %This is not an essential test, but if there are cycles that trigger a warning in the qc log, then 
    %it may be helpful to know that the float has been marked as questionable

    %NITRATE
    questionableCycle_nitrate = cell2mat(badCycles.N_cycles_3(badCycles.N_wmo_3==wmoList(c)));
    if ~isempty(questionableCycle_nitrate)
        tests.questionableNITRATEtest = 0;
    else
        tests.questionableNITRATEtest = 1;
    end

    %PH
    questionableCycle_ph = cell2mat(badCycles.PH_cycles_3(badCycles.PH_wmo_3==wmoList(c)));
    if ~isempty(questionableCycle_ph)
        tests.questionablePHtest = 0;
    else
        tests.questionablePHtest = 1;
    end
    %===========================================================================================
    
    %CHECK FOR ANY GOOD NITRATE AND PH DATA
    %===========================================================================================
    %Now that the data has made it through initial screening, we can begin
    %to prep it for QC analysis by subsetting good data at depth
    
    DATA.PH = DATA.data(DATA.data(:,DATA.iStn)>=startCycle_ph & DATA.data(:,DATA.iStn)<=endCycle_ph ...
        & DATA.data(:,DATA.iPH+1)~=8,:);

    DATA.N = DATA.data(DATA.data(:,DATA.iStn)>=startCycle_nitrate & DATA.data(:,DATA.iStn)<=endCycle_nitrate & ...
        DATA.data(:,DATA.iN+1)~=8,:);

    %Checks that GOOD PH and NITRATE data exists
    tests.phDataTest = ~isempty(DATA.PH);
    tests.nitrateDataTest = ~isempty(DATA.N);
    
    if ~tests.phDataTest && ~tests.nitrateDataTest
        disp('No good PH or NITRATE data, skipping float')
        fprintf(msgLog,'No good PH or NITRATE data, skipping float\n');
        continue;
    end
    %===========================================================================================
    
    %CHECK FOR ANY GOOD DEEP NITRATE AND PH DATA then subsample
    %===========================================================================================

    %SUBSAMPLE DEEP PH
    %If the float is on the pump offset list, subsample to 950m 
    tests.pumpOffset = ismember(wmoList(c),pH_pumpoffset_980_floats);
    if tests.pumpOffset
        disp("Pump offset float")
        fprintf(msgLog,'Pump offset float\n');

        %Subsample the deep data to 950-980m
        depthRange = [920 980];
        DATA.deepPH = DATA.PH(DATA.PH(:,DATA.iZ)>=depthRange(1) & ...
            DATA.PH(:,DATA.iZ)<=depthRange(2),:);
        
        disp("Final PH depth range: 920-980m")
        fprintf(msgLog,'Final PH depth range: 920-980m\n');
        qcDepth_ph = "950";

        tests.shallowFloatTest_ph = 0;
     
    %If not, then subsample pH from 1480-1520m
    else 

        %Start with initial depth range standard for QC
        depthRange = [1480 1520];
        
        DATA.deepPH = DATA.PH(DATA.PH(:,DATA.iZ)>=depthRange(1) & DATA.PH(:,DATA.iZ)<=depthRange(2),:);
        totalStns = unique(DATA.PH(:,DATA.iStn));
        deepStns = unique(DATA.deepPH(:,DATA.iStn));
        tests.shallowFloatTest_ph = 0;

        %If there's any missing profiles, progressively increase the depth range by 100
        while all(length(totalStns) ~= length(deepStns))
            
            depthRange = depthRange - 100;
            DATA.deepPH = DATA.PH(DATA.PH(:,DATA.iZ)>=depthRange(1) & DATA.PH(:,DATA.iZ)<=depthRange(2),:);
            deepStns = unique(DATA.deepPH(:,DATA.iStn));

            %If no consistent data exists below 500m, consider this a shallow, possible Ross Sea float
            if depthRange(1) <= 500
                disp('No PH data below 500m')
                fprintf(msgLog,'No PH data below 500m\n');
                tests.shallowFloatTest_ph = 1;
                break;
            end
        end
        %While loop ends when a depth with no missing data is reached
        
        disp("Final PH depth range: "+depthRange(1)+"-"+depthRange(2))
        fprintf(msgLog,"Final PH depth range: %d-%dm\n", depthRange(1), depthRange(2));
        qcDepth_ph = string(mean(depthRange));
    end
    DATA.qcDepth_ph = depthRange;

    %SUBSAMPLE DEEP NITRATE
    
    %Start with initial depth range standard for QC
    depthRange = [1480 1520];

    DATA.deepN = DATA.N(DATA.N(:,DATA.iZ)>=depthRange(1) & DATA.N(:,DATA.iZ)<=depthRange(2),:);
    totalStns = unique(DATA.N(:,DATA.iStn));
    deepStns = unique(DATA.deepN(:,DATA.iStn));
    tests.shallowFloatTest_nitrate = 0;
    
    %If there's any missing profiles, progressively increase the depth range by 100
    while all(length(totalStns) ~= length(deepStns))
        
        depthRange = depthRange - 100;
        
        DATA.deepN = DATA.N(DATA.N(:,DATA.iZ)>=depthRange(1) & DATA.N(:,DATA.iZ)<=depthRange(2),:);
        deepStns = unique(DATA.deepN(:,DATA.iStn));

        if depthRange(1) <= 500
            disp('No NITRATE data below 500m')
            fprintf(msgLog,'No NITRATE data below 500m\n');
            tests.shallowFloatTest_nitrate = 1;
            break;
        end
    end
    %While loop ends when a depth with no missing data is reached
    
    disp("Final NITRATE depth range: "+depthRange(1)+"-"+depthRange(2))
    fprintf(msgLog,"Final NITRATE depth range: %d-%dm\n", depthRange(1), depthRange(2));
    qcDepth_nitrate = string(mean(depthRange));
    DATA.qcDepth_nitrate = depthRange;

    %FUTURE FEATURE: SHALLOW AND ROSS SEA FLOATS
    %Test if floats are in the Ross Sea
    % tests.RossSeaTest = any(DATA.data(:,DATA.iLat) <= -71 & DATA.data(:,DATA.iLon) >= 160 & DATA.data(:,DATA.iLon) <= 200);
    % if tests.RossSeaTest
    %     disp("Float has profiles in the Ross Sea")
    %     fprintf(msgLog,"Float has profiles in the Ross Sea\n");
    % end

    %Test if both nitrate and pH are shallow to determine if float shoals
    tests.shallowFloatTest = tests.shallowFloatTest_ph & tests.shallowFloatTest_nitrate;

    %Check that there is deep data after subsampling
    tests.phDeepDataTest = ~isempty(DATA.deepPH);
    tests.nitrateDeepDataTest = ~isempty(DATA.deepN);
    
    %If there's no data at depth, skip the float
    if ~tests.phDeepDataTest && ~tests.nitrateDeepDataTest
        disp('No good DEEP PH or DEEP NITRATE data, skipping float')
        fprintf(msgLog,'No good DEEP PH or DEEP NITRATE data, skipping float\n');
        continue;
    end
    
    %NITRATE QC
    %===========================================================================================
    %Only runs through realtimeQC function if data passes ALL tests
    if all([tests.nitrateColTest; tests.nitrateBSLtest; tests.nitrateDataTest; tests.nitrateDeepDataTest])
        disp('Processing NITRATE')
        fprintf(msgLog,'Processing NITRATE\n');
        nitrateDATA = realtimeQC(DATA, 5, tests.badO2Test);
    %If not then display the reason it was terminated
    else
        failureMSG = {"Float has no NITRATE"; "All future NITRATE is on BSL"; "Data has no good NITRATE"; "Data has no DEEP NITRATE"};
        tests.NfailureTests = [tests.nitrateColTest==0 tests.nitrateBSLtest==0 tests.nitrateDataTest==0 tests.nitrateDeepDataTest==0];
        fprintf('%s\n', failureMSG{tests.NfailureTests})
        fprintf(msgLog,'%s\n', failureMSG{tests.NfailureTests});
        clear nitrateDATA %Clear the data structure from previous float that is not replaced by current float
    end
    %===========================================================================================

    %PH QC
    %===========================================================================================
    if all([tests.phColTest; tests.phBSLtest; tests.phDataTest; tests.phDeepDataTest])
        disp('Processing PH')
        fprintf(msgLog,'Processing PH\n');
        phDATA = realtimeQC(DATA, 3, tests.badO2Test);
    else
        failureMSG = {"Float has no PH"; "All future PH is on BSL"; "Data has no good PH"; "Data has no DEEP PH"};
        tests.PHfailureTests = [tests.phColTest==0 tests.phBSLtest==0 tests.phDataTest==0 tests.phDeepDataTest==0];
        fprintf('%s\n', failureMSG{tests.PHfailureTests})
        fprintf(msgLog,'%s\n', failureMSG{tests.PHfailureTests});
        clear phDATA
    end
    %=========================================================================================== 

    %FILL IN REALTIME QC LOG FILE
    %=========================================================================================== 

    %Identify if any tests were triggered for NITRATE data
    if exist('nitrateDATA', 'var')
        tests.noO2test_nitrate = ~isempty(nitrateDATA.noO2stns); %Check if there are any stations without O2 over the whole profile
        tests.noO2Deeptest_nitrate = ~isempty(nitrateDATA.noO2stns_deep);
        tests.noGPStest_nitrate = ~isempty(nitrateDATA.noGPSstns); %Check if there are any stations without GPS fix
        tests.refAnomaly_nitrate = ~isempty(nitrateDATA.refTriggerStns); %Check if there are any stations that breach anomaly/delta thresholds
        tests.cycleDelta_nitrate = ~isempty(nitrateDATA.deltaTriggerStns);

        % tests.refAnomaly_nitrate = nitrateDATA.refAnom>=refAnomThd_nitrate; %Check where anomaly threshold is tripped
        % tests.cycleDelta_nitrate = nitrateDATA.cycleDelta>=deltaThd_nitrate & tests.refAnomaly_nitrate; %Check where cycle delta AND ref anomaly thresholds are tripped
        %If cycle delta is tripped but the sudden change is predicted by EsperMix, then no warning should be raised
        
    else
        tests.noO2test_nitrate = 0;
        tests.noO2Deeptest_nitrate = 0;
        nitrateDATA.noO2stns = []; %In case nitrate data doesn't exist but pH data does, create a dummy variable to union
        nitrateDATA.noO2stns_deep = [];
        tests.noGPStest_nitrate = 0;
        nitrateDATA.noGPSstns = [];
        tests.cycleDelta_nitrate = 0;
        tests.refAnomaly_nitrate = 0;
    end 
    
    %Do the same for PH data
    if exist('phDATA', 'var')
        tests.noO2test_ph = ~isempty(phDATA.noO2stns);
        tests.noO2Deeptest_ph = ~isempty(phDATA.noO2stns_deep);
        tests.noGPStest_ph = ~isempty(phDATA.noGPSstns);
        tests.refAnomaly_ph = ~isempty(phDATA.refTriggerStns);
        tests.cycleDelta_ph = ~isempty(phDATA.deltaTriggerStns);

        % tests.refAnomaly_ph = phDATA.refAnom>=refAnomThd_ph;
        % tests.cycleDelta_ph = phDATA.cycleDelta>=deltaThd_ph & tests.refAnomaly_ph;
        
    else
        tests.noO2test_ph = 0;
        tests.noO2Deeptest_ph = 0;
        phDATA.noO2stns = [];
        phDATA.noO2stns_deep = [];
        tests.noGPStest_ph = 0;
        phDATA.noGPSstns = [];
        tests.cycleDelta_ph = 0;
        tests.refAnomaly_ph = 0;
    end
    
    %List all logical tests for the qc log
    tests.finalTestList = [tests.noGPStest_nitrate;
        tests.noGPStest_ph;
        tests.noO2test_nitrate & ~tests.badO2Test;
        tests.noO2Deeptest_nitrate & ~tests.badO2Test;
        tests.noO2test_ph & ~tests.badO2Test;
        tests.noO2Deeptest_ph & ~tests.badO2Test;
        any(tests.cycleDelta_nitrate);
        any(tests.cycleDelta_ph);
        any(tests.refAnomaly_nitrate);
        any(tests.refAnomaly_ph)];

    %Set up an initial trigger where if ANY of the tests are tripped, create a new section in the log file for the float
    if any(tests.finalTestList)

        %LOGICAL STRING CODE FROM "JAN" AT: 
        % https://www.mathworks.com/matlabcentral/answers/258287-fprintf-for-logical-statement
        %pumpOffset = ismember(wmoList(c),pH_pumpoffset_980_floats);
        LogicalStr = {"(pump offset float)"; "(bad O2 float)"; "(shallow float)"; "(bad NITRATE float)"; "(bad PH float)"; ""}; %True false string that indicates if float is marked with pump offset
        LogicalStrTests = [tests.pumpOffset tests.badO2Test tests.shallowFloatTest tests.nitrateBSLtest==0 tests.phBSLtest==0 ...
            (~tests.pumpOffset & ~tests.badO2Test & ~tests.shallowFloatTest & tests.nitrateBSLtest & tests.phBSLtest)];
        horizontal_str = sprintf('%s ', LogicalStr{LogicalStrTests});
        fprintf(qcLog, "\n%s (%s) %s | N cycles (%sm): %d-%d | PH cycles (%sm): %d-%d\n", ...
            string(wmoList(c)), idList{c}, horizontal_str, qcDepth_nitrate, startCycle_nitrate, endCycle_nitrate, qcDepth_ph, startCycle_ph, endCycle_ph); %Print the string
        
        %List any indefinite questionable listings on the BSL
        if ~tests.questionableNITRATEtest
            fprintf(qcLog, "NITRATE marked indefinitely QUESTIONABLE at cycle: %d\n", ...
                questionableCycle_nitrate);
        end
    
        if ~tests.questionablePHtest
            fprintf(qcLog, "PH marked indefinitely QUESTIONABLE at cycle: %d\n", ...
                questionableCycle_ph);
        end
    end
    
    %Check for cycles with no GPS fix (should be the same for PH and NITRATE)
    if tests.noGPStest_nitrate||tests.noGPStest_ph
        noGPSCycles = union(nitrateDATA.noGPSstns, phDATA.noGPSstns);
        horizontal_str = sprintf('%d, ', noGPSCycles);
        fprintf(qcLog, "No GPS fix for cycle(s): %s\n", horizontal_str);
        clear noGPSCycles %Clear the cycles at the end so they don't carry over to other floats
    end
    
    %Check if only deep oxygen is missing
    if (tests.noO2Deeptest_nitrate||tests.noO2Deeptest_ph) && ~tests.badO2Test
        
        %Check if entire oxygen profiles are missing
        if tests.noO2test_nitrate||tests.noO2test_ph
            noO2Cycles = union(nitrateDATA.noO2stns, phDATA.noO2stns);
            horizontal_str = sprintf('%d, ', noO2Cycles);
            fprintf(qcLog, "Missing OXYGEN for cycle(s): %s\n",horizontal_str);
            clear noO2Cycles %Clear the cycles at the end so they don't carry over to other floats
        else
            noO2Cycles = union(nitrateDATA.noO2stns_deep, phDATA.noO2stns_deep);
            horizontal_str = sprintf('%d, ', noO2Cycles);
            fprintf(qcLog, "Missing OXYGEN at depth for cycle(s): %s\n",horizontal_str);
            clear noO2Cycles %Clear the cycles at the end so they don't carry over to other floats
        end
    end

    %Check if the cycle-to-cycle delta thresholds are tripped for NITRATE
    if any(tests.cycleDelta_nitrate)
        %horizontal_str = sprintf('%d, ', nitrateDATA.stns(tests.cycleDelta_nitrate));
        horizontal_str = sprintf('%d, ', nitrateDATA.nitrateDATA.deltaTriggerStns);
        fprintf(qcLog, "NITRATE cycle delta threshold exceeded at cycle(s): %s\n", horizontal_str);
    end
   
    %Check if the float-to-Esper anomaly thresholds are tripped for NITRATE
    if any(tests.refAnomaly_nitrate)
        %horizontal_str = sprintf('%d, ', nitrateDATA.stns(tests.refAnomaly_nitrate));
        horizontal_str = sprintf('%d, ', nitrateDATA.refTriggerStns);
        fprintf(qcLog, "ESPER NITRATE anomaly threshold exceeded at cycle(s): %s\n", horizontal_str);
    end
    
    %Repeat tests for PH
    if any(tests.cycleDelta_ph)
        %horizontal_str = sprintf('%d, ', phDATA.stns(tests.cycleDelta_ph));
        horizontal_str = sprintf('%d, ', phDATA.deltaTriggerStns);
        fprintf(qcLog, "PH cycle delta threshold exceeded at cycle(s): %s\n", horizontal_str);
    end

    if any(tests.refAnomaly_ph)
        %horizontal_str = sprintf('%d, ', phDATA.stns(tests.refAnomaly_ph));
        horizontal_str = sprintf('%d, ', phDATA.refTriggerStns);
        fprintf(qcLog, "ESPER PH anomaly threshold exceeded at cycle(s): %s\n", horizontal_str);
    end
    
end
%%% END OF MAIN FOR LOOP %%%
%=========================================================================================== 
%===========================================================================================
%================================ CONSTRUCT THE EMAIL ======================================
%===========================================================================================
%This block of code takes the txt file that was written in the previous code block and converts
%it into a string array that can be sent as a nicely formatted email.

%Read the document as a string array 
S = readlines(getenv("USERPROFILE")+"\Documents\rtqc_triggerLog.txt");

%Final check to indicate that no floats triggered a qc warning
if length(S) <= 10
    fprintf(qcLog, "\nNO QC WARNINGS TO REPORT");

    %Update string array with new line
    S = readlines(getenv("USERPROFILE")+"\Documents\rtqc_triggerLog.txt");
end


%Close the qc log file after all floats have been viewed
fclose(qcLog);
fclose(msgLog);

%--------------------------------------------------------------------------------------
% Send the email (This is only) for local testing purposes
%--------------------------------------------------------------------------------------
sender = 'lgrady@mbari.org';
setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
setpref('Internet','E_mail',sender); % define sender
sendmail(email_list,subjectline, S) %Send the email

















%===========================================================================================
%===================================== QC FUNCTION =========================================
%===========================================================================================

function [rtqcDATA] = realtimeQC(DATA, varType, badO2Test)
    
    % INPUTS
    % DATA: A data struct created by the realtime qc code that contains all pH and nitrate data
    % and subsample data at depth for each variable

    % VARTYPES: inputs to the EsperMix function to output a desired variable
    % 3 = PH
    % 5 = NITRATE
    
    %Identify the subsampled deep data from DATA depending on the desired varType
    if varType == 5
        dDeep = DATA.deepN;
        dData = DATA.N;
        varIDX = DATA.iN;
        deltaThd = DATA.deltaThd_nitrate; %Cycle-to-cycle delta for NITRATE
        refAnomThd = DATA.refAnomThd_nitrate; %Float-to-Esper anomaly for NITRATE
    elseif varType == 3
        dDeep = DATA.deepPH;
        dData = DATA.PH;
        varIDX = DATA.iPH;
        deltaThd = DATA.deltaThd_ph; %Cycle-to-cycle delta for PH
        refAnomThd = DATA.refAnomThd_ph; %Float-to-Esper anomaly for PH
    else
        disp("Error location data in DATA structure; check varType")
    end

    %CHECK FOR GPS ISSUES AND CALCULATE ESPERMIX
    %===========================================================================================
    %GPS is a necessary variable to estimate ESPERMIX, so if a float is missing a GPS fix for a cycle
    %the ESPERMIX can't be estimated for that cycle. This block of code creates a NaN array and an
    %array of cycles WITH a gps fix. EsperMix is calculated for the cycles WITH a gps fix and appended
    %to the NaN array. This way, regardless of a float with, without, or with partial gps fixes, Espermix
    %is estimated for everything.

    %List of stns WITHOUT a gps fix
    rtqcDATA.noGPSstns = unique(dDeep(isnan(dDeep(:,DATA.iLon))|isnan(dDeep(:,DATA.iLat)),DATA.iStn));

    %Indices of data WITH a gps fix
    gpsIDX = ~isnan(dDeep(:,DATA.iLon))|~isnan(dDeep(:,DATA.iLat));
    
    %Subset data that has a gps fix to calculate Esper
    dDeep_GPS = dDeep(~isnan(dDeep(:,DATA.iLon))|~isnan(dDeep(:,DATA.iLat)),:);

    %Create an empty array of NaNs to fill with Esper data while NaN'ing data with NO GPS
    RefData = NaN(length(dDeep),1);
    
    %CHECK FOR ANY MISSING OXYGEN
    %===========================================================================================
    %O2 is not required to estimate EsperMix, but does change the equation. This block locates any cycles
    %that are missing O2 throughout the profile and at the qc depth. If any cycles are missing O2 at the 
    %qc depth, then EsperMix is calculated without O2. This is done to keep consistent with how MBARI
    %qc's floats.

    %In case there's missing O2
    stns_deep = unique(dDeep(:,DATA.iStn)); %Find all unique stations at depth
    stns_profile = unique(dData(:,DATA.iStn)); %Find all unique stations at over the entire profile
    
    %Find indices of unique stations
    G = findgroups(dDeep(:,DATA.iStn));

    %Find indices of unique stations for entire profile
    G_profile = findgroups(dData(:,DATA.iStn));

    %Split the array by station number into a cell array
    C = splitapply(@(m){m},dDeep(:,DATA.iO),G);

    %Split the array by station number into a cell array
    C_profile = splitapply(@(m){m},dData(:,DATA.iO),G_profile); 
    
    %Apply a logical test to find which stations contain NO O2 at depth and over the entire profile
    rtqcDATA.noO2stns_deep = stns_deep(cellfun(@(x) all(isnan(x)), C)); 
    rtqcDATA.noO2stns = stns_profile(cellfun(@(x) all(isnan(x)), C_profile));
    %===========================================================================================
    
    %ESTIMATE ESPERMIX
    %=========================================================================================== 
    %This block takes all of the data with GPS fixes and estimates EsperMix either with O2 or without O2.

    %If there's no good gps fixes, then EsperMix = NaN array
    if isempty(dDeep_GPS)
        disp('Missing gps locations for ALL stations')
    elseif any(isnan(dDeep(:,DATA.iLon))|isnan(dDeep(:,DATA.iLat))) & ~isempty(dDeep_GPS)
        disp('Missing gps locations for SOME stations, estimating reference where possible')
    end
    
    %If there's no missing O2 and GPS fixes, calculate ESPERMIX as normal
    if (isempty(rtqcDATA.noO2stns) | isempty(rtqcDATA.noO2stns_deep)) & ~badO2Test & ~isempty(dDeep_GPS)

        % CALCULATE ESPER pH or NO3
        Estimates.DesireVar = [varType];
        Estimates.OutCoords = [dDeep_GPS(:,3), dDeep_GPS(:,4), dDeep_GPS(:,DATA.iZ)]; % lon,lat,Z
        Estimates.PredictorTypes    = [1 2 6]; % PSAL, TEMP, DOXY_ADJ,
        Estimates.Measurements = [dDeep_GPS(:,DATA.iS), dDeep_GPS(:,DATA.iT), dDeep_GPS(:,DATA.iO)];
        Estimates.dvec = datevec(dDeep_GPS(:,1));
        Estimates.dec_yr = Estimates.dvec(:,1) +(Estimates.dvec(:,2)*30)/365 + Estimates.dvec(:,3)/365; % very crude decimal year
        Equations    = 7; % S, T, O2
        
        %Estimate ESPERMIX with no O2 for the desired variable
        [Estimates,~] = ESPER_Mixed(Estimates.DesireVar, Estimates.OutCoords, Estimates.Measurements,...
                    Estimates.PredictorTypes, 'Equations', Equations, 'EstDates', Estimates.dec_yr);

        %Identify the subsampled deep data from DATA depending on the desired varType
        if varType == 5
            %Fills the NaN array with reference estimates where GPS fixes exist
            RefData(gpsIDX,1) = Estimates.nitrate;
        elseif varType == 3
            RefData(gpsIDX,1) = Estimates.pH;
        else
            disp("Error location data in DATA structure; check varType")
        end

    %If there are stations without O2 but still have gps fixes, estimate ESPERMIX without O2 (Equation 8)
    elseif ((~isempty(rtqcDATA.noO2stns) | ~isempty(rtqcDATA.noO2stns_deep)) | badO2Test) & ~isempty(dDeep_GPS)
        
        % CALCULATE ESPER pH or NO3
        Estimates.DesireVar = [varType];
        Estimates.OutCoords = [dDeep_GPS(:,3), dDeep_GPS(:,4), dDeep_GPS(:,DATA.iZ)]; % lon,lat,Z
        Estimates.PredictorTypes    = [1 2 6]; % PSAL, TEMP
        Estimates.Measurements = [dDeep_GPS(:,DATA.iS), dDeep_GPS(:,DATA.iT), dDeep_GPS(:,DATA.iO)];
        Estimates.dvec = datevec(dDeep_GPS(:,1));
        Estimates.dec_yr = Estimates.dvec(:,1) +(Estimates.dvec(:,2)*30)/365 + Estimates.dvec(:,3)/365; % very crude decimal year
        Equations    = 8; % S, T

        %Estimate ESPERMIX with no O2 for the desired variable
        [Estimates,~] = ESPER_Mixed(Estimates.DesireVar, Estimates.OutCoords, Estimates.Measurements,...
                    Estimates.PredictorTypes, 'Equations', Equations, 'EstDates', Estimates.dec_yr);

        %Identify the subsampled deep data from DATA depending on the desired varType
        if varType == 5
            %Fills the NaN array with reference estimates where GPS fixes exist
            RefData(gpsIDX,1) = Estimates.nitrate;
        elseif varType == 3
            RefData(gpsIDX,1) = Estimates.pH;
        else
            disp("Error location data in DATA structure; check varType")
        end
    end
    rtqcDATA.Esper = Estimates;
    %===========================================================================================
    
    %REALTIME QC TESTING
    %===========================================================================================
    %The heart of this function is the realtime qc tests. This block runs the data through the cycle delta
    %and reference anomaly tests. Using logical indices of the cycles that breach the qc thresholds (in the
    %beginning of the code), a third test is run to find any sequences of 3 or more consecutive cycles that 
    %breach the thresholds. This ensures that natural noise during 1 or 2 cycles doesn't trigger a warning.

    %Create array of unique stations in deep PH dataset
    stns = unique(dDeep(:,DATA.iStn))';
    rtqcDATA.stns = stns;

    %Retrieve indices of the first cycle to start from
    stnInd = dDeep(:,DATA.iStn)==stns(1);
    
    %Create blank array for mean deep data at each station
    rtqcDATA.deepData_mean = NaN(length(stns),1);
    
    %Create blank array for mean deep Esper data at each station
    rtqcDATA.RefData_mean = NaN(length(stns),1);

    %Mean Float-to-Esper anomaly at depth
    rtqcDATA.refAnom = NaN(length(stns),1);

    for i = 1:length(stns)
        
        %Assign indices for the current station
        stnInd = dDeep(:,DATA.iStn)==stns(i);
        
        %Find pH over entire deep range
        deepData_all = dDeep(stnInd,varIDX);
        
        %Find EsperMix pH at depth
        deepRefData = RefData(stnInd);
        
        %Average. float pH at depth that can be used to calculate the
        %cycle-to-cycle difference in pH
        rtqcDATA.deepData_mean(i,1) = mean(deepData_all);
    
        %Avg. EsperMix pH at depth
        rtqcDATA.RefData_mean(i,1) = mean(deepRefData);
    
        %Avg. anomaly between deep float and EsperMix pH
        rtqcDATA.refAnom(i,1) = mean(deepData_all-deepRefData);
    end
    
    %Append phDelta with a 0 at the beginning to compensate for indices offset
    rtqcDATA.cycleDelta = [0; diff(rtqcDATA.deepData_mean)];

    %RUN CONSECUTIVE TRIGGER TESTS
    %Code based on example found from:
    %https://www.mathworks.com/matlabcentral/answers/366126-how-many-consecutive-ones

    %Create logical array for delta and ref anom tests (any delta or refAnom beyond the 
    %threshold is flagged)
    deltaConsec = (abs(rtqcDATA.cycleDelta) >= deltaThd) & (abs(rtqcDATA.refAnom) >= refAnomThd);
    refConsec = abs(rtqcDATA.refAnom) >= refAnomThd;
    
    %Create empty zero arrays to populate with indices later
    deltaOut = zeros(size(deltaConsec'));
    refOut = zeros(size(refConsec'));
    
    %Duplicate that zero array but prepend and append a zero for indexing
    deltaOutB = [0,deltaConsec',0];
    refOutB = [0,refConsec',0];
    
    %Find the index of any INITIAL failure, where (i.e. where deltaConsec goes from 0 to 1)
    iDelta = strfind(deltaOutB, [0 1]);
    iRef = strfind(refOutB, [0 1]);

    %Now find the index of any FINAL failure and subtract the INITIAL failure index to
    %calculate how many consecutive points exist between the initial and final failure 
    deltaOut(iDelta) = strfind(deltaOutB, [1 0]) - iDelta; %append the zeros array "delta/ref Out"
    refOut(iRef) = strfind(refOutB, [1 0]) - iRef;

    %i.e. If failures begin at cycle 3 and end at cycle 7, this means cycle 3,4,5,6,7 get flagged. By taking
    %the difference between the initial (3) and final (7), we find that this test was triggered 5 
    %consecutive times from cycle 3 (the initial cycle is counted)

    %Since we're only interested in 3 or more consecutive failures (filters out instances of 1 or 2 failures,
    %that could just be caused by natural noise in the environment)
    iDelta = iDelta(deltaOut(iDelta) >= DATA.consecCounter);
    iRef = iRef(refOut(iRef) >= DATA.consecCounter);
    
    %Create an array of cycle blocks that consists of each INITIAL and FINAL failure
    deltaStnBlocks = [stns(iDelta); stns(iDelta)+deltaOut(iDelta)-1]';
    refStnBlocks = [stns(iRef); stns(iRef)+refOut(iRef)-1]';

    %Create empty arrays to be appended with all failed stations
    rtqcDATA.deltaTriggerStns = [];
    rtqcDATA.refTriggerStns = [];
    
    %For each pair of cycles, populate the previous empty array with ALL cycles between the INITIAL and FINAL failure cycles
    deltaDim = size(deltaStnBlocks);
    refDim = size(refStnBlocks);
    
    if deltaDim(1) == 1 & ~isempty(deltaStnBlocks)
           rtqcDATA.deltaTriggerStns = [rtqcDATA.deltaTriggerStns, deltaStnBlocks(1):1:deltaStnBlocks(2)];
    elseif deltaDim(1) > 1
        for x = 1:deltaDim(1)
            rtqcDATA.deltaTriggerStns = [rtqcDATA.deltaTriggerStns, deltaStnBlocks(x,1):1:deltaStnBlocks(x,2)];
        end
    end
    if refDim(1) == 1 & ~isempty(refStnBlocks)
        rtqcDATA.refTriggerStns = [rtqcDATA.refTriggerStns, refStnBlocks(1):1:refStnBlocks(2)];
    elseif refDim(1) > 1
        for x = 1:refDim(1)
            rtqcDATA.refTriggerStns = [rtqcDATA.refTriggerStns, refStnBlocks(x,1):1:refStnBlocks(x,2)];
        end
    end
end
%=========================================================================================== 



%===========================================================================================
%================================= FIND INITIAL QC FLOATS ==================================
%===========================================================================================
%This function checks MBARI's float list for any floats without a QC file and 6 or more profiles,
%indicating that they require an initial QC by MBARI standards.
function [no_QC_floats] = check_1st_QC()
    % Check_1st_QC.m 
    
    dirs.QCList = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\';
    dirs.cal    = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\';
    
    % ************************************************************************
    % LOAD MBARI MASTER LIST & SUBSET TO pH FLOATS ONLY
    load([dirs.cal, 'MBARI_float_list.mat'])
    fhdr   = d.hdr;
    flist  = d.list;
    
    I.MBA  = strcmp(fhdr,'MBARI ID'); % Define lookup indices for float list
    I.INST = strcmp(fhdr,'INST ID');
    I.WMO  = strcmp(fhdr,'WMO');
    I.TYP  = strcmp(fhdr,'float type');
    I.DIR  = strcmp(fhdr,'msg dir');
    I.PRG  = strcmp(fhdr,'Program');
    I.REG  = strcmp(fhdr,'Region');
    I.PH   = strcmp(fhdr,'tf pH');
    I.CMAX = strcmp(fhdr,'max cycle proc');
    
    % Exclusions from Initial QC Check by MBARI ID
    exclude_str = 'un0948|ua19191|ua19298|ua19314|ua19843';
    tf_exclude = cellfun(@isempty, regexp(flist(:,I.MBA), exclude_str,'once'));
    flist = flist(tf_exclude,:);
    
    rflist = size(flist,1);
    clear d
    
    % ***********************************************************************
    % GET LIST OF QC LIST FILES
    tmp = dir([dirs.QCList,'*FloatQCList.txt']);
    QCfiles = {tmp.name}';
    
    tf_1st_QC = ones(rflist,1)*0;
    for i = 1:rflist
        % Is there a QC List file for float
        cycle_max =flist{i,I.CMAX};
        t1 = ~cellfun(@isempty, regexp(QCfiles, flist{i,I.WMO},'once'));
        if sum(t1) == 1
            tf_1st_QC(i) = 1;
        elseif isnan(cycle_max) || cycle_max < 6 % enough cycles to QC
            tf_1st_QC(i) = 1;
        end
    end
    
    no_QC_floats = flist(~tf_1st_QC,3);
    clearvars -except no_QC_floats
end
%===========================================================================================