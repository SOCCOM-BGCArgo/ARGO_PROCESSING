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
% CHANGELOG
% 
% LG 9/3/2024: Revised the deep data code so that the QC depth is found only using data that is not
% flagged as "fail". This is to try and catch floats such as 5906473, which has bad data at depth
% that we no longer qc to. The actual qc depth is around 1100m, but the code would still see the bad
% data and analyze the deeper values, consistently flagging this float.
%
% LG 9/9/2024: Added 'omitmissing' to the means when calculating refanom and cycle delta. May start to 
% experience an uptick in float catches.
%
% LG 9/11/2024: Added a specific exception for float 6990585, which has Muxer board issues, so that it 
% will always be counted as a "badO2" float. It still receives O2, just at <900m and it doesn't always
% record N, but kept getting flagged as missing O2. Now it doesn't.
%
% LG 9/19/2024: Changed the pump offset list so it loads from CHEM rather than local to prevent the constant
% need to update.
%
% LG 9/23/2024: Floats with cycles between last DMQC and the "failure" point on the BSL no longer
% look at the last cycle. I.E. wn1200 is marked indefinitely as bad pH data from cycle 28 and beyond. The code
% finds the qc depth to analyze deep data based on the depth where ALL stations have data. Since cycle 28 was
% marked bad, there is effectively no data there, since the code only subsets with data marked as "good/questionable".
% In these specific cases, the range of cycles to analyze is now set to the most recent DMQC cycle and the BSL
% cycle - 1, so that only the "good" cycles are analyzed.
%
% LG 10/07/2024: Added a condition to all tests that the last cycle of the breached tests must equal the most recent
% cycle of the float in order to trigger an RTQC warning message. This prevents repeats of warnings from past weeks
% that haven't been looked at yet.
%
% LG 10/24/2024: Starting progress on making a warning "table" to serve as a prototype for global DAC warning table.
%
% LG 10/29/2024: Added a warning catch log table that is created and appended based on new catches. Tidied up the tests
% to be more intuitive when sorting through float descriptions. Added a sequential spike test to the "RTQC_tests"
% function that has yet to be implemented into the warning emails.
%
% LG 11/05/2024: Added a test to check the catch log for any false alarms (flags of 1). These have to be manually changed
% in the catch log after the fact, but once they are added, the exact same catch hopefull doesn't have to be addressed
% again.
%===========================================================================================


%===========================================================================================
%===================================== INPUTS ==============================================
%===========================================================================================

%Turn off warnings
warning("off")

%Set the current user directory to pull local files from
userDir = getenv("USERPROFILE");

%Indicate if an email should be send out to recipients (0 = NO, 1 = YES)
tests.sendEmail = 0;

%RTQC test thresholds to display on the email (NOT LINKED TO ACTUAL RTQC_tests FUNCTION)
tests.rangeChkThd_nitrate = [-10 55]; %Value range check for NITRATE
tests.deltaThd_nitrate = 1; %Cycle-to-cycle delta for NITRATE
tests.refAnomThd_nitrate = 1; %Float-to-Esper anomaly for NITRATE

tests.rangeChkThd_ph = [7.3 8.5]; %Value range check for PH
tests.deltaThd_ph = .01; %Cycle-to-cycle delta for PH
tests.refAnomThd_ph = .012; %Float-to-Esper anomaly for PH

%Inputs for the email
% if tests.sendEmail
%     email_list = {'lgrady@mbari.org';'tmaurer@mbari.org';"jplant@mbari.org";"sbartoloni@mbari.org";"nicolag@mbari.org";"johnson@mbari.org"}; %Recipients
%     subjectline = "TEST: RTQC DAILY UPDATE"; %Subject line
% end

%The range of time to subset floats on the float list (default is from when the script is run to 24 hours ago)
timelag = datetime("now")-hours(24);

%Load the full MBARI float list as a table
floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");

%Sort any dead floats and those that haven't been deployed yet (Shouldn't be getting new data from those anyways)
floatList = floatList(floatList.tfDead==0 & ~isnat(floatList.x1stDate),:);
floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "ascend");

%Subset all floats that have come in from the past 24 hours
floatList = floatList(datenum(floatList.maxCycleProcDate)>=datenum(timelag),:);

%Format the latestCycleDate Column to a simple date without the time extension
floatList.latestCycleFileDate.Format = 'MM/dd/uuuu';

%Create list of float WMOs and MBARI IDs to call on
wmoList = floatList.WMO;
idList = floatList.MBARIID;

%Create a list of dates to ID floats
latestDateList = floatList.latestCycleFileDate;

%Open up the catchlog as a table
catchLog = readtable(userDir+"\Documents\qcCatchLog.csv");
%===========================================================================================

%===========================================================================================
%=============================== INPUTS FOR TESTING ========================================
%===========================================================================================

% %Overwrite any email inputs down here
% if tests.sendEmail
%     email_list = {'lgrady@mbari.org'};
%     subjectline = "TEST";
% end

%FOR TESTING DIFFERENT TIME LAGS
% timelag = datetime("now")-hours(24);

%FOR TESTING ON A SPECIFIC DATE RANGE
dateStart = datenum(datetime(2024,11,12,6,0,0));
dateEnd = datenum(datetime(2024,11,13,6,0,0));
floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");
floatList = floatList(floatList.tfDead==0 & ~isnat(floatList.x1stDate),:);
floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "ascend");
floatList = floatList(datenum(floatList.maxCycleProcDate)>=dateStart & datenum(floatList.maxCycleProcDate)<=dateEnd,:);
floatList.latestCycleFileDate.Format = 'MM/dd/uuuu';
wmoList = floatList.WMO;
idList = floatList.MBARIID;
latestDateList = floatList.latestCycleFileDate;

%FOR TESTING ON SPECIFIC FLFOATS OR INDIVIDUAL FLOATS
% wmoList = [3902557;5907054;5906540];
% wmoList = [5906533;5906487];
% floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");
% floatList = floatList(ismember(floatList.WMO,wmoList),:);
% floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "ascend");
% floatList.latestCycleFileDate.Format = 'MM/dd/uuuu';
% wmoList = floatList.WMO;
% idList = floatList.MBARIID;
% latestDateList = floatList.latestCycleFileDate;

%FOR TESTING ON THE WHOLE FLEET
% floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");
% floatList = floatList(floatList.tfDead==0 & ~isnat(floatList.x1stDate),:);
% floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "ascend");
% floatList.latestCycleFileDate.Format = 'MM/dd/uuuu';
% wmoList = floatList.WMO;
% idList = floatList.MBARIID;
% latestDateList = floatList.latestCycleFileDate;
% sz = [0 16];
% varTypes = ["datetime","double","string","double","double","double","double","double","double","string","string","double","double","double","string","double"];
% varNames = ["DATE","WMO","MBARIID","PUMPOFFSET","N_QBL","N_BAD","PH_QBL","PH_BAD","O2_BAD","SENSOR","CATCH","DEPTH","STARTCYC","ENDCYC","CYCLES_CAUGHT","FLAG"];
% catchLog = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
% catchLog.DATE.Format = 'MM/dd/uuuu';
%===========================================================================================

%LOAD BSL AND IDENTIFY CYCLES FOR PH, P, T, AND S
%===========================================================================================
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
run("\\atlas\Chem\ARGO_PROCESSING\MFILES\FLOATS\pH_pumpoffset_980.m")

%START WRITING THE RTQC EMAIL
%===========================================================================================

%Open up a qc log text file to store qc messages
qcLog = fopen(userDir+"\Documents\rtqc_triggerLog.txt","w");
msgLog = fopen(userDir+"\Documents\rtqc_msgLog.txt","w");

%Write all test parameters from the inputs at the top of the file, %.3g limits numbers to 3 sig figs
fprintf(qcLog, "Nitrate range check thresholds: [%.3g %.3g]\n" + ...
    "Nitrate sensor drift threshold = +/- %.3g\n" + ...
    "Nitrate reference anomaly threshold = +/- %.3g\n" + ...
    "pH range check thresholds: [%.3g %.3g]\n" + ...
    "PH sensor threshold = +/- %.3g\n" + ...
    "PH reference anomaly threshold = +/- %.3g\n" + ...
    "-------------------------------------------------------\n", ...
    tests.rangeChkThd_nitrate(1), tests.rangeChkThd_nitrate(2), tests.deltaThd_nitrate, tests.refAnomThd_nitrate, ...
    tests.rangeChkThd_ph(1), tests.rangeChkThd_ph(2), tests.deltaThd_ph, tests.refAnomThd_ph);

%Load the list of floats that need initial qc using the function "check_1st_QC"
if ~isempty(check_1st_QC_func)
   
    %Convert the output list to a double
    no_QC_floats = str2double(cell2mat(check_1st_QC_func()));
    
    %Horizontally merge the WMO's and print into a line on qcLog
    horizontal_str = sprintf('%d, ', no_QC_floats);
    fprintf(qcLog, "\nFLOATS THAT REQUIRE INITIAL QC: %s\n" + ...
        "\n-------------------------------------------------------\n", horizontal_str);
%If not initial qc floats, make the variable empty
else
    no_QC_floats = [];
end

%Print an example line to show the format of the warning messages
fprintf(qcLog, "\nQC MESSAGE FORMAT:\nWMO (MBARIID) (Other descriptors) | N or PH (qc depth): last qc'd cycle - current cycle |\n");
%===========================================================================================

%===========================================================================================
%================================ BEGIN REALTIME QC LOOP ===================================
%===========================================================================================

%Keep a tally of all floats being analyzed and the floats that trigger RTQC tests
totalFloats = length(wmoList);
triggeredFloats = 0;

%START THE LOOP
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
    
    %Start a data struct for all of the logical tests and add this test that checks if
    %DATA struct has any columns for PH or NITRATE
    tests.hasPH = any(contains(DATA.hdr,'pHinsitu[Total]'));
    tests.hasNITRATE = any(contains(DATA.hdr,'Nitrate[µmol/kg]'));
    
    %If float isn't equipped with either nitrate or ph, skip to next float
    if ~tests.hasNITRATE && ~tests.hasPH
        disp('No NITRATE or PH detected, skipping float')
        fprintf(msgLog, 'No NITRATE or PH detected, skipping float\n');
        continue;
    end
    %===========================================================================================

    %RETRIEVE THE CYCLE NUMBER OF LAST KNOWN DMQC FROM QCLIST FILE
    %===========================================================================================

    %If the qc list file exists for the current float
    if exist("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+wmoList(c)+"_FloatQCList.txt",'file')

        %Identify the last cycle where qc was performed by looking at the floatQClist file
        fid = fopen("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+wmoList(c)+"_FloatQCList.txt", 'r');
        
        %Create empty string array
        startCycle_ph = {};
        startCycle_nitrate = {};
    
        %Get the first line of the floatQClist file
        next_line = fgetl(fid);
        
        %Iterate through the file until reaching a blank line, indicating
        %the end of the most recent qc list file text block
        while ~isempty(next_line)
            next_line = fgetl(fid); %Read the current line
            
            %If the current line is empty, end the while loop
            if isempty(next_line)|next_line==-1 
                break;
            end
            
            %If the line contains the word "Nitrate" or "pH", split it by comma and
            %store the second to last number, which should be the cycle
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
            %If there's no start cycle from the qc list file, then just default to cycle 1
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
    endCycle_ph = DATA.data(end,DATA.iStn); %The final station in the dataset
    endCycle_nitrate = DATA.data(end,DATA.iStn);
    
    %If the end cycle is equal to the start cycle, no new cycles to process and skip the float
    if (startCycle_nitrate == endCycle_nitrate) & (startCycle_ph == endCycle_ph)
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
        tests.phBSL = 1;

    %If the BSL cycle is placed AFTER our last QC, there may still be cycles to process
    elseif BSLCycle_ph > startCycle_ph + 1
        disp('PH listed on BSL, but there are still unprocessed cycles')
        fprintf(msgLog,'PH listed on BSL, but there are still unprocessed cycles\n');
        %end cycle is the cycle JUST BEFORE the cycle listed on the BSL, because the cycle on the BSL should be "bad"
        endCycle_ph = BSLCycle_ph-1;
        tests.phBSL = 0;
    %If there's any other weirdness, just make the end cycle the last processed cycle from the floatlist
    else
        % endCycle_ph = endCycle;
        tests.phBSL = 0;
    end

    %NITRATE (same steps as PH before)
    BSLCycle_nitrate = cell2mat(badCycles.N_cycles(badCycles.N_wmo==wmoList(c)));

    if BSLCycle_nitrate <= startCycle_nitrate + 1
        endCycle_nitrate = BSLCycle_nitrate;
        startCycle_nitrate = BSLCycle_nitrate;
        tests.nitrateBSL = 1;

    elseif BSLCycle_nitrate > startCycle_nitrate + 1
        disp('NITRATE listed on BSL, but there are still unprocessed cycles')
        fprintf(msgLog,'NITRATE listed on BSL, but there are still unprocessed cycles\n');
        endCycle_nitrate = BSLCycle_nitrate - 1;
        tests.nitrateBSL = 0;
    else
        % endCycle_nitrate = endCycle;
        tests.nitrateBSL = 0;
    end

    %If ph and nitrate are both indefinitely listed on BSL then skip the float
    if tests.phBSL && tests.nitrateBSL
        disp('PH and NITRATE indefinitely listed on BSL, skipping float')
        fprintf(msgLog,'PH and NITRATE indefinitely listed on BSL, skipping float\n');
        % continue;
    end

    %Check OXYGEN on the BSL
    BSLCycle_oxygen = cell2mat(badCycles.O_cycles(badCycles.O_wmo==wmoList(c)));
    
    %EXEPTION: this float will always have O2 issues and shouldn't be checked for missing O2
    tests.MUXERfloat = wmoList(c)==6990585; 

    %Create a test that marks a float as having bad O2 if it's on the BSL with bad O2 cycles, AND it isn't the MUXER issue float (6990585)
    tests.badO2 = (ismember(wmoList(c),badCycles.O_wmo) & ~isempty(BSLCycle_oxygen)) | tests.MUXERfloat;
    %IF the badO2Test is TRUE, missing O2 cycles located during the RTQC test process will NOT be marked in the email
    %IF the badO2Test is FALSE, then missing O2 cycles will be marked in the email, because these are presumed to be novel issues

    if tests.badO2
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
        tests.questionableNITRATE = 1;
    else
        tests.questionableNITRATE = 0;
    end

    %PH
    questionableCycle_ph = cell2mat(badCycles.PH_cycles_3(badCycles.PH_wmo_3==wmoList(c)));
    if ~isempty(questionableCycle_ph)
        tests.questionablePH = 1;
    else
        tests.questionablePH = 0;
    end
    %===========================================================================================
    
    %SUBSET THE DATA FROM LAST KNOWN DMQC
    %===========================================================================================
    %Now that the data has made it through initial screening, we can begin
    %to prep it for QC analysis by subsetting GOOD data at depth
    
    %Retrieve RAW data from all depths at the desired cycles
    DATA.PH = DATA.data(DATA.data(:,DATA.iStn)>=startCycle_ph & DATA.data(:,DATA.iStn)<=endCycle_ph,:);
    DATA.N = DATA.data(DATA.data(:,DATA.iStn)>=startCycle_nitrate & DATA.data(:,DATA.iStn)<=endCycle_nitrate,:);
    
    %CHECK THAT NEW NITRATE AND PH DATA EXISTS (bit of a redundancy test for weird floats within our system)
    tests.phDataTest = ~isempty(DATA.PH);
    tests.nitrateDataTest = ~isempty(DATA.N);
    
    %If there is somehow a float that makes it past the new cycles test but is marked on the BSL (for example)
    %and now doesn't have new data past the last DMQC, skip the float
    if ~tests.phDataTest && ~tests.nitrateDataTest
        disp('No good PH or NITRATE data, skipping float')
        fprintf(msgLog,'No good PH or NITRATE data, skipping float\n');
        continue;
    end
    %===========================================================================================
    
    %SUBSAMPLE DEEP PH DATA
    %===========================================================================================
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

        tests.shallowFloat_ph = 0;
     
    %If not, then subsample pH from 1480-1520m
    else 

        %Start with initial depth range standard for QC
        depthRange = [1480 1520];
        
        %Subset data within the depth range AND data that's not flagged as fail (just for this deep subset)
        DATA.deepPH = DATA.PH(DATA.PH(:,DATA.iZ)>=depthRange(1) & DATA.PH(:,DATA.iZ)<=depthRange(2),:);
        totalStnsPH = unique(DATA.PH(:,DATA.iStn));
        deepStns = unique(DATA.deepPH(:,DATA.iStn));
        tests.shallowFloat_ph = 0;

        %If there's any missing profiles, progressively increase the depth range by 100
        while all(length(totalStnsPH) ~= length(deepStns))
            depthRange = depthRange - 100;
            DATA.deepPH = DATA.PH(DATA.PH(:,DATA.iZ)>=depthRange(1) & DATA.PH(:,DATA.iZ)<=depthRange(2),:);
            deepStns = unique(DATA.deepPH(:,DATA.iStn));

            %If no consistent data exists below 500m, consider this a shallow, possible Ross Sea float
            if depthRange(1) <= 500
                %If the depth range is this shallow then widen the range to most of the profile regardless
                depthRange = [100 500];
                disp('No PH data below 500m')
                fprintf(msgLog,'No PH data below 500m\n');
                tests.shallowFloat_ph = 1;
                break;
            end
        end
        %While loop ends when a depth with no missing data is reached
        
        disp("Final PH depth range: "+depthRange(1)+"-"+depthRange(2))
        fprintf(msgLog,"Final PH depth range: %d-%dm\n", depthRange(1), depthRange(2));

        %Make a string of the final qc depth for the email message
        qcDepth_ph = string(mean(depthRange));
    end

    %Populate the DATA structure with the qcDepth range
    DATA.qcDepth_ph = depthRange;

    %Check if there is any GOOD deep data, which will determine if bad data is omitted in the final RTQC
    %tests. If there is no good data at depth, the RTQC code can't evaluate the recent data, so this test
    %tells RTQC_tests to keep bad data so there is something to analyze.
    tests.goodDeepPH = ~isempty(DATA.deepPH(DATA.deepPH(:,DATA.iPH+1)~=8,:));
    %===========================================================================================

    %SUBSAMPLE DEEP NITRATE
    %===========================================================================================
    %Start with initial depth range standard for QC
    depthRange = [1480 1520];

    DATA.deepN = DATA.N(DATA.N(:,DATA.iZ)>=depthRange(1) & DATA.N(:,DATA.iZ)<=depthRange(2),:);
    totalStnsN = unique(DATA.N(:,DATA.iStn));
    deepStns = unique(DATA.deepN(:,DATA.iStn));
    tests.shallowFloat_nitrate = 0;
    
    %If there's any missing profiles, progressively increase the depth range by 100
    while all(length(totalStnsN) ~= length(deepStns))
        
        depthRange = depthRange - 100;
        
        DATA.deepN = DATA.N(DATA.N(:,DATA.iZ)>=depthRange(1) & DATA.N(:,DATA.iZ)<=depthRange(2),:);
        deepStns = unique(DATA.deepN(:,DATA.iStn));
        
        if depthRange(1) <= 500
            %If the depth range is this shallow then widen the range to most of the profile regardless
            depthRange = [100 500];
            disp('No NITRATE data below 500m')
            fprintf(msgLog,'No NITRATE data below 500m\n');
            tests.shallowFloat_nitrate = 1;
            break;
        end
    end
    %While loop ends when a depth with no missing data is reached
    
    disp("Final NITRATE depth range: "+depthRange(1)+"-"+depthRange(2))
    fprintf(msgLog,"Final NITRATE depth range: %d-%dm\n", depthRange(1), depthRange(2));
    qcDepth_nitrate = string(mean(depthRange));
    DATA.qcDepth_nitrate = depthRange;
    tests.goodDeepN = ~isempty(DATA.deepN(DATA.deepN(:,DATA.iN+1)~=8,:));
    %===========================================================================================

    %QC DEPTH TESTS
    %===========================================================================================
    %Test if both nitrate and pH are shallow to determine if float shoals
    tests.shallowFloat = tests.shallowFloat_ph & tests.shallowFloat_nitrate;
    %===========================================================================================

    %QC DEPTH EXCEPTION LIST (TEMPORARY UNTIL UPDATED QC LIST FILE FORMAT IS RELEASED)
    %===========================================================================================
    %Some floats are qc'd to specific depths that cannot be categorized by the qc depth finding code above
    %RTQC_tests function allows for manual input of qcDepth, so this exception list can be used to properly
    %analyze floats with odd qc depths.

    if wmoList(c)==5906473
        %For this specific float, change the depth ranges
        DATA.qcDepth_nitrate = [1480 1520];
        qcDepth_nitrate = string(mean(DATA.qcDepth_nitrate));

        DATA.qcDepth_ph = [1100 1200];
        qcDepth_ph = string(mean(DATA.qcDepth_ph));
        
        %Make sure the shallow float test is reset to NOT TRUE
        tests.shallowFloat_ph = 0;
    end
    
    %Deprecating this specific exception after qc'ing to 1500m; float should be treated as normal
    % if wmoList(c)==2903854
    %     DATA.qcDepth_nitrate = [600 800];
    %     qcDepth_nitrate = string(mean(DATA.qcDepth_nitrate));
    %     DATA.qcDepth_ph = [600 800];
    %     qcDepth_ph = string(mean(DATA.qcDepth_ph));
    %     tests.shallowFloat_ph = 0;
    % end
    %===========================================================================================
    
    %RUN NITRATE RTQC TESTS
    %===========================================================================================
    %Only runs through realtimeQC function if data passes ALL tests
    if all([tests.hasNITRATE; ~tests.nitrateBSL])
        disp('Processing NITRATE')
        fprintf(msgLog,'Processing NITRATE\n');

        DATA.cycleRange_nitrate = [startCycle_nitrate endCycle_nitrate];

        %RTQC_TESTS(full dataset, varType for ESPER, qc depth, cycle range, badO2test, omit bad data, use reference anomaly test)
        nitrateDATA = RTQC_tests(DATA, 5, DATA.qcDepth_nitrate, DATA.cycleRange_nitrate, tests.badO2, tests.goodDeepN, 1);

    %If not all tests are passed then display the reason it was terminated
    else
        failureMSG = {"Float has no NITRATE"; "All future NITRATE is on BSL"; "No NITRATE data below 500m"};
        tests.NfailureTests = logical([~tests.hasNITRATE tests.nitrateBSL tests.shallowFloat_nitrate]);
        fprintf('%s\n', failureMSG{tests.NfailureTests})
        fprintf(msgLog,'%s\n', failureMSG{tests.NfailureTests});
        clear nitrateDATA %Clear the data structure from previous float that is not replaced by current float
    end
    %===========================================================================================

    %RUN PH RTQC TESTS
    %===========================================================================================
    if all([tests.hasPH; ~tests.phBSL])
        disp('Processing PH')
        fprintf(msgLog,'Processing PH\n');
        DATA.cycleRange_ph = [startCycle_ph endCycle_ph];
        phDATA = RTQC_tests(DATA, 3, DATA.qcDepth_ph, DATA.cycleRange_ph, tests.badO2, tests.goodDeepPH, 1);
    else
        failureMSG = {"Float has no PH"; "All future PH is on BSL"; "No PH data below 500m"};
        tests.PHfailureTests = logical([~tests.hasPH tests.phBSL tests.shallowFloat_ph]);
        fprintf('%s\n', failureMSG{tests.PHfailureTests})
        fprintf(msgLog,'%s\n', failureMSG{tests.PHfailureTests});
        clear phDATA
    end
    %=========================================================================================== 

    %FILL IN REALTIME QC LOG FILE
    %=========================================================================================== 
    %Identify if any tests were triggered for NITRATE data
    if exist('nitrateDATA', 'var')
        %Check if there are any stations without O2 AND that O2 isn't already on BSL (&& used to prevent issues where noO2stns is empty)
        tests.noO2_nitrate = ~tests.badO2 & ~isempty(nitrateDATA.noO2stns);% && (nitrateDATA.noO2stns(end)==totalStnsN(end));

        %Check if there are any stations without GPS fix (&& used to prevent issues where noGPSstns is empty)
        tests.noGPStest_nitrate = ~isempty(nitrateDATA.noGPSstns);% && (nitrateDATA.noGPSstns(end)==totalStnsN(end));

        %Check if there are any stations that breach anomaly/delta thresholds
        tests.refAnomaly_nitrate = ~isempty(nitrateDATA.refAnomTriggerStns);% && (nitrateDATA.refAnomTriggerStns(end)==totalStnsN(end));

        %Check if there are any stations that breach anomaly/delta thresholds
        tests.cycleDelta_nitrate = ~isempty(nitrateDATA.cycleDeltaTriggerStns);% && (nitrateDATA.cycleDeltaTriggerStns(end)==totalStnsN(end));
        
        %Check for any stations that breach the sequential spike test threshold
        tests.seqSpike_nitrate = ~isempty(nitrateDATA.seqSpikeTriggerStns);% && (nitrateDATA.seqSpikeTriggerStns(end)==totalStnsN(end));

        %Check if there are stations that breach the range check AND if they are on the recent cycle (&& used to prevent issues where rangeTriggerstns is empty)
        % if ~isempty(nitrateDATA.rangeTriggerStns)
        tests.rangeChk_nitrate = ~isempty(nitrateDATA.rangeChkTriggerStns);% && (nitrateDATA.rangeChkTriggerStns(end)==totalStnsN(end));
    
    %In case nitrate data doesn't exist but pH data does, create dummy variables for if statements below to continue working
    else
        tests.noO2_nitrate = 0;
        nitrateDATA.noO2stns = []; 
        tests.noGPStest_nitrate = 0;
        nitrateDATA.noGPSstns = [];
        tests.cycleDelta_nitrate = 0;
        tests.seqSpike_nitrate = 0;
        tests.refAnomaly_nitrate = 0;
        tests.rangeChk_nitrate = 0;
    end 
    
    %Do the same for PH data
    if exist('phDATA', 'var')
        tests.noO2_ph = ~tests.badO2 & ~isempty(phDATA.noO2stns);% && (phDATA.noO2stns(end)==totalStnsPH(end));
        tests.noGPStest_ph = ~isempty(phDATA.noGPSstns);% && (phDATA.noGPSstns(end)==totalStnsPH(end));
        tests.refAnomaly_ph = ~isempty(phDATA.refAnomTriggerStns);% && (phDATA.refAnomTriggerStns(end)==totalStnsPH(end));
        tests.cycleDelta_ph = ~isempty(phDATA.cycleDeltaTriggerStns);% && (phDATA.cycleDeltaTriggerStns(end)==totalStnsPH(end));
        tests.seqSpike_ph = ~isempty(phDATA.seqSpikeTriggerStns);% && (phDATA.seqSpikeTriggerStns(end)==totalStnsN(end));
        tests.rangeChk_ph = ~isempty(phDATA.rangeChkTriggerStns);% && (phDATA.rangeChkTriggerStns(end)==totalStnsPH(end));
        
    else
        tests.noO2_ph = 0;
        phDATA.noO2stns = [];
        tests.noGPStest_ph = 0;
        phDATA.noGPSstns = [];
        tests.cycleDelta_ph = 0;
        tests.seqSpike_ph = 0;
        tests.refAnomaly_ph = 0;
        tests.rangeChk_ph = 0;
    end
    
    %CHECK CATCHLOG FOR ANY PREVIOUS DRIFT LISTINGS THAT HAVE BEEN MARKED AS FALSE ALARMS
    tests.floatOnCatchLog = ismember(catchLog.WMO, wmoList(c));

    if any(tests.floatOnCatchLog) & tests.cycleDelta_ph
        repeatFloats = catchLog(tests.floatOnCatchLog,:);
        repeatFloats = repeatFloats(contains(repeatFloats.CATCH,'DRIFT'),:);
        cyclesCaught = cell2mat(cellfun(@str2double,regexp(repeatFloats.CYCLES_CAUGHT,'\d+','match'),'UniformOutput',false));
        tests.falseDrift_ph = any(ismember(phDATA.cycleDeltaTriggerStns,cyclesCaught)) & ...
            contains(repeatFloats.CATCH,'DRIFT') & ...
            ((repeatFloats.FLAG==1)|(repeatFloats.FLAG==6));
    else
        tests.falseDrift_ph = 0;
    end

    if any(tests.floatOnCatchLog) & tests.cycleDelta_nitrate
        repeatFloats = catchLog(tests.floatOnCatchLog,:);
        repeatFloats = repeatFloats(contains(repeatFloats.CATCH,'DRIFT'),:);
        cyclesCaught = cell2mat(cellfun(@str2double,regexp(repeatFloats.CYCLES_CAUGHT,'\d+','match'),'UniformOutput',false));
        tests.falseDrift_nitrate = any(ismember(nitrateDATA.cycleDeltaTriggerStns,cyclesCaught)) & ...
            contains(repeatFloats.CATCH,'DRIFT') & ...
            ((repeatFloats.FLAG==1)|(repeatFloats.FLAG==6));
    else
        tests.falseDrift_nitrate = 0;
    end

    %List all logical tests for the qc log
    tests.finalTestList = [tests.noGPStest_nitrate;
        tests.noGPStest_ph;
        tests.noO2_nitrate & ~tests.badO2;
        tests.noO2_ph & ~tests.badO2;
        tests.cycleDelta_nitrate & ~tests.falseDrift_nitrate;
        tests.cycleDelta_ph & ~tests.falseDrift_ph;
        tests.refAnomaly_nitrate;
        tests.refAnomaly_ph;
        tests.seqSpike_nitrate;
        tests.seqSpike_ph;
        tests.rangeChk_nitrate;
        tests.rangeChk_ph];
    
    %Set up an initial trigger where if ANY of the tests are tripped, create a new section in the log file for the float
    if any(tests.finalTestList)

        %For tracking number of floats that get triggered
        triggeredFloats = triggeredFloats + 1;

        %LOGICAL STRING CODE FROM "JAN" AT: 
        % https://www.mathworks.com/matlabcentral/answers/258287-fprintf-for-logical-statement
        %pumpOffset = ismember(wmoList(c),pH_pumpoffset_980_floats);
        LogicalStr = {"(pump offset float)"; "(bad O2 float)"; "(shallow float)"; "(bad NITRATE float)"; "(bad PH float)"; ""}; %True false string that indicates if float is marked with pump offset
        
        %Now create a logical array that will indicate which descriptors from LogicalStr to add to the listing for this float on the email
        LogicalStrTests = logical([tests.pumpOffset tests.badO2 tests.shallowFloat tests.nitrateBSL tests.phBSL ...
            (~tests.pumpOffset & ~tests.badO2 & ~tests.shallowFloat & ~tests.nitrateBSL & ~tests.phBSL)]); %This line is a series of tests that indicate whether the current float has
                                                                                                                  % NO extra descriptors

        %Combine all of the descriptors (if there are any)
        horizontal_str = sprintf('%s ', LogicalStr{LogicalStrTests});

        %Print everything into one nice line in qcLog
        fprintf(qcLog, "\n%s (%s) %s | N cycles (%sm): %d-%d | PH cycles (%sm): %d-%d\n", ...
            string(wmoList(c)), idList{c}, horizontal_str, qcDepth_nitrate, startCycle_nitrate, endCycle_nitrate, qcDepth_ph, startCycle_ph, endCycle_ph); %Print the string
        
        %List any indefinite questionable listings on the BSL
        if tests.questionableNITRATE
            fprintf(qcLog, "NITRATE marked indefinitely QUESTIONABLE at cycle: %d\n", ...
                questionableCycle_nitrate);
        end
    
        if tests.questionablePH
            fprintf(qcLog, "PH marked indefinitely QUESTIONABLE at cycle: %d\n", ...
                questionableCycle_ph);
        end
    end
    
    %Check for cycles with no GPS fix (should be the same for PH and NITRATE)
    if tests.noGPStest_nitrate||tests.noGPStest_ph
        noGPSCycles = union(nitrateDATA.noGPSstns, phDATA.noGPSstns);
        horizontal_str = sprintf('%d,', noGPSCycles);
        fprintf(qcLog, "No GPS fix for cycle(s): %s\n", horizontal_str);

        %Table = {datetime, wmo, 'ID', pumpOffset, N_QBL, N_BAD, PH_QBL, PH_BAD, O2_BAD, 'SENSOR', 'CATCH', DEPTH, STARTCYC, ENDCYC, 'CYCLES_CAUGHT', FLAG}
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            '' 'NO GPS' NaN NaN NaN horizontal_str 5 ''};

        %Append the catchLog with the new information
        catchLog = [catchLog; newRow];
    
        %Clear the cycles at the end so they don't carry over to other floats
        clear noGPSCycles 
    end
    
    %Check if O2 is missing from certain stations
    if tests.noO2_nitrate||tests.noO2_ph
        noO2Cycles = union(nitrateDATA.noO2stns, phDATA.noO2stns);
        horizontal_str = sprintf('%d,', noO2Cycles);
        fprintf(qcLog, "Missing OXYGEN for cycle(s): %s\n",horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            '' 'NO O2' NaN NaN NaN horizontal_str 5 ''};
        catchLog = [catchLog; newRow];
        
        %Clear the cycles at the end so they don't carry over to other floats
        clear noO2Cycles
    end
    
    if tests.rangeChk_nitrate
       horizontal_str = sprintf('%d,', nitrateDATA.rangeChkTriggerStns); 
       fprintf(qcLog, "NITRATE range check exceeded at cycle(s): %s\n", horizontal_str);
       newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'N' 'RNGCHK' qcDepth_nitrate startCycle_nitrate endCycle_nitrate horizontal_str 4 ''};
        catchLog = [catchLog; newRow];
    end

    %Check if the cycle-to-cycle delta thresholds are tripped for NITRATE
    if tests.cycleDelta_nitrate & ~tests.falseDrift_nitrate
        horizontal_str = sprintf('%d,', nitrateDATA.cycleDeltaTriggerStns);
        fprintf(qcLog, "NITRATE sensor drift threshold exceeded at cycle(s): %s\n", horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'N' 'DRIFT' qcDepth_nitrate startCycle_nitrate endCycle_nitrate horizontal_str 2 ''};
        catchLog = [catchLog; newRow];
    end

    if tests.seqSpike_nitrate
        horizontal_str = sprintf('%d,', nitrateDATA.seqSpikeTriggerStns);
        fprintf(qcLog, "NITRATE jump threshold exceeded at cycle(s): %s\n", horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'N' 'JUMP' qcDepth_nitrate startCycle_nitrate endCycle_nitrate horizontal_str 4 ''};
        catchLog = [catchLog; newRow];
    end
   
    %Check if the float-to-Esper anomaly thresholds are tripped for NITRATE
    if tests.refAnomaly_nitrate
        horizontal_str = sprintf('%d,', nitrateDATA.refAnomTriggerStns);
        fprintf(qcLog, "NITRATE reference anomaly threshold exceeded at cycle(s): %s\n", horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'N' 'REFANOM' qcDepth_nitrate startCycle_nitrate endCycle_nitrate horizontal_str 2 ''};
        catchLog = [catchLog; newRow];
    end
    
    %Repeat tests for PH
    if tests.rangeChk_ph
       horizontal_str = sprintf('%d,', phDATA.rangeChkTriggerStns); 
       fprintf(qcLog, "PH range check exceeded at cycle(s): %s\n", horizontal_str);
       newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'PH' 'RNGCHK' qcDepth_ph startCycle_ph endCycle_ph horizontal_str 4 ''};
       catchLog = [catchLog; newRow];
    end

    if tests.cycleDelta_ph & ~tests.falseDrift_ph
        horizontal_str = sprintf('%d,', phDATA.cycleDeltaTriggerStns);
        fprintf(qcLog, "PH sensor drift threshold exceeded at cycle(s): %s\n", horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'PH' 'DRIFT' qcDepth_ph startCycle_ph endCycle_ph horizontal_str 2 ''};
        catchLog = [catchLog; newRow];
    end

    if tests.seqSpike_ph
        horizontal_str = sprintf('%d,', phDATA.seqSpikeTriggerStns);
        fprintf(qcLog, "PH jump threshold exceeded at cycle(s): %s\n", horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'PH' 'JUMP' qcDepth_ph startCycle_ph endCycle_ph horizontal_str 4 ''};
        catchLog = [catchLog; newRow];
    end

    if tests.refAnomaly_ph
        horizontal_str = sprintf('%d,', phDATA.refAnomTriggerStns);
        fprintf(qcLog, "PH reference anomaly threshold exceeded at cycle(s): %s\n", horizontal_str);
        newRow = {latestDateList(c) wmoList(c) idList{c} ...
            tests.pumpOffset tests.questionableNITRATE tests.nitrateBSL tests.questionablePH tests.phBSL tests.badO2 ...
            'PH' 'REFANOM' qcDepth_ph startCycle_ph endCycle_ph horizontal_str 2 ''};
        catchLog = [catchLog; newRow];
    end   
end
%===========================================================================================
%================================ END OF REALTIME QC LOOP ==================================
%===========================================================================================

%CONSTRUCT THE EMAIL
%===========================================================================================
%This block of code takes the txt file that was written in the previous code block and converts
%it into a string array that can be sent as a nicely formatted email.

%Read the document as a string array 

fprintf(qcLog, "\n%d / %d floats breached RTQC tests", triggeredFloats, totalFloats);

S = readlines(getenv("USERPROFILE")+"\Documents\rtqc_triggerLog.txt");


%Final check to indicate that no floats triggered a qc warning
if triggeredFloats == 0
    fprintf(qcLog, "\nNO QC WARNINGS TO REPORT");

    %Update string array with new line
    S = readlines(getenv("USERPROFILE")+"\Documents\rtqc_triggerLog.txt");
end

%Write the catch log as a csv
writetable(catchLog,userDir+"\Documents\qcCatchLog.csv", 'Delimiter',',')

%Close the qc log file after all floats have been viewed
fclose(qcLog);
fclose(msgLog);
%===========================================================================================

% Send the email
%===========================================================================================
if tests.sendEmail
    sender = 'lgrady@mbari.org';
    setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
    setpref('Internet','E_mail',sender); % define sender
    sendmail(email_list,subjectline, S) %Send the email
end
%===========================================================================================