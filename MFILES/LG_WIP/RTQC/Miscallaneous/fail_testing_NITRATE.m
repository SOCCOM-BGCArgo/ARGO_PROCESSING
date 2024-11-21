%===========================================================================================
% OBJECTIVE: This is a specific version of "realtimeQC" that looks at floats that have failed
% Nitrate and PH sensors and checks where the realtimeQC tests would mark the failure point
% compared to what is on the bad sensor list.

%===========================================================================================

warning("off")

%===========================================================================================
%===================================== INPUTS ==============================================
%===========================================================================================
% deltaThd_ph = .01; %Cycle-to-cycle delta for PH
% refAnomThd_ph = .01; %Float-to-Esper anomaly for PH
% rangeChk_ph = [7.3 8.5]; %Value range check for PH
% deltaThd_nitrate = 1; %Cycle-to-cycle delta for NITRATE
% refAnomThd_nitrate = 1; %Float-to-Esper anomaly for NITRATE
% rangeChk_nitrate = [-1 55]; %Value range check for NITRATE
% consecCounter = 3; %Number of consecutive test failures before recording the cycle


% email_list = {'lgrady@mbari.org'};
% subjectline = "SENSOR FAIL TESTING";
% timelag = now-hours(24);
%===========================================================================================

%===========================================================================================
%===================================== LOAD DATA ===========================================
%===========================================================================================

%LOAD FLOATS
% floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");

% %Sort any dead floats and those that haven't been deployed yet (Shouldn't be getting new data from those anyways)
% floatList = floatList(floatList.tfDead==0 & ~isnat(floatList.x1stDate),:);
% floatList = sortrows(floatList(~ismissing(floatList.maxCycleProcDate),:), "maxCycleProcDate", "descend");
% floatList = floatList(datenum(floatList.maxCycleProcDate)>=datenum(timelag),:);
% 
% %FOR TESTING ON THE WHOLE FLEET AND NOT PAST 24 HOURS
% %floatList = floatList(datenum(floatList.maxCycleProcDate)<=datenum(now),:);
% 
% %Create list of float WMO's to call on
% wmoList = floatList.WMO;
% idList = floatList.MBARIID;
% wmoList = [3902557;5907054;5906540];
% idList = {'ua23596';'ua21977';'un1452'};

% wmoList = 5906524;
% idList = {'ua20644'};

% wmoList = 5906218;
% idList = {'ua18081'};

% wmoList = 5904984;
% idList = {'un0569'};

%LOAD BSL AND IDENTIFY CYCLES FOR PH, P, T, AND S:
BSL = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT");
BSL = sortrows(BSL,"WMO_","ascend"); %Sort by WMO's in ascending order to match floatsDir

%Regular expression matches any string that's a number followed by a dash but NOT FOLLOWED by a dash 
%and a number, ensuring that only indefinite failures are accounted for
badCycles.PH_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"PH")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.PH_wmo = BSL(contains(BSL.SENSOR,"PH")&BSL.FLAG==4,"WMO_").WMO_;
badCycles.PH_mid = BSL(contains(BSL.SENSOR,"PH")&BSL.FLAG==4,"MBARIIDSTR").MBARIIDSTR;

%NITRATE
badCycles.N_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"N")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.N_wmo = BSL(contains(BSL.SENSOR,"N")&BSL.FLAG==4,"WMO_").WMO_;
badCycles.N_mid = BSL(contains(BSL.SENSOR,"N")&BSL.FLAG==4,"MBARIIDSTR").MBARIIDSTR;

%OXYGEN
badCycles.O_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"O")&~contains(BSL.SENSOR,"CDOM")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.O_wmo = BSL(contains(BSL.SENSOR,"O")&~contains(BSL.SENSOR,"CDOM")&BSL.FLAG==4,"WMO_").WMO_;

%LOAD THE LIST OF PUMP OFFSET FLOATS:
pH_pumpoffset_980 

wmoList = badCycles.N_wmo;
% wmoList = 5906439;
% wmoList = 5904676;
wmoList = 5904660;

driftTestDiff = NaN(length(wmoList),1);
rangeTestDiff = NaN(length(wmoList),1);

%===========================================================================================
%================================ BEGIN REALTIME QC LOOP ===================================
%===========================================================================================

for c = 1:length(wmoList)
    
    badCycleIDX = badCycles.N_wmo==wmoList(c);
    if isempty(badCycles.N_cycles{badCycleIDX})
        continue;
    end
    
    %For keeping track of issues from specific floats
    disp(wmoList(c))
    
    if isnan(wmoList(c))
        dips('WMO is NaN');
        continue;
    else
        DATA = get_FloatViz_data("\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+wmoList(c)+"QC.txt");
    end    


    %Convert all -1e10 points (missing) to NaN
    DATA.data(DATA.data==-1e10) = NaN;
    %===========================================================================================

    % GET SOME VARIABLE INDICES
    %===========================================================================================
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
    %===========================================================================================
    
    %CHECK FOR ANY GOOD DEEP NITRATE AND PH DATA then subsample
    %===========================================================================================

    depthRange = [1480 1520];

    %Check OXYGEN on the BSL
    BSLCycle_oxygen = cell2mat(badCycles.O_cycles(badCycles.O_wmo==wmoList(c)));
    
    %EXEPTION: this float will always have O2 issues and shouldn't be checked for missing O2
    tests.MUXERfloat = wmoList(c)==6990585; 

    %Create a test that marks a float as having bad O2 if it's on the BSL with bad O2 cycles, AND it isn't the MUXER issue float (6990585)
    badO2Test = (ismember(wmoList(c),badCycles.O_wmo) & ~isempty(BSLCycle_oxygen)) | tests.MUXERfloat;
    
    startCycle = DATA.data(1,DATA.iStn);
    endCycle = DATA.data(end,DATA.iStn);
    cycleRange = [startCycle endCycle];

    %NITRATE QC
    %===========================================================================================
    %Only runs through realtimeQC function if data passes ALL tests
    disp('Processing NITRATE')
    nitrateDATA = RTQC_tests(DATA,5,depthRange,cycleRange,badO2Test,0,0);
    %===========================================================================================
    
    bslCyc = badCycles.N_cycles{badCycleIDX};
    
    if ~isempty(nitrateDATA.rangeChkTriggerStns)
        rngChkDiff = nitrateDATA.rangeChkTriggerStns - bslCyc;
        [rangeMin,rangeIDX] = min(abs(rngChkDiff));
        rangeTestDiff(c) = rngChkDiff(rangeIDX);
    end
    if ~isempty(nitrateDATA.cycleDeltaTriggerStns)
        driftDiff = nitrateDATA.cycleDeltaTriggerStns - bslCyc;
        [driftMin,driftIDX] = min(abs(driftDiff));
        driftTestDiff(c) = driftDiff(driftIDX);
    end

    if isempty(nitrateDATA.stnsDEEP)
        continue;
    end

    %PH QC
    %===========================================================================================
    % disp('Processing PH')
    % %RTQC_tests((DATA, varType, qcDepth, badO2Test, omitBadDeepData,refAnomTest))
    % phDATA = RTQC_tests(DATA, 3, depthRange, cycleRange, badO2Test, 0, 0); %Run the tests without omitting bad data (for fail testing)
    %                                                            %and without the refAnomTest
    %=========================================================================================== 
    
%Check if a directory to store plots already exists
    userDir = getenv("USERPROFILE");
    fileDir = userDir+"\Documents\MATLAB\ARGO_PROCESSING\MFILES\LG_WIP\RTQC\NITRATE_RTQC_PLOTS\";
    % fileDir = userDir+"\Documents\data_pics\"+wmoList(c);
    % 
    % if exist(fileDir,"dir")~=7
    %     disp('Making new directory')
    %     mkdir(fileDir)
    % end

    %PLOTTING
    %===========================================================================================
    
    %Basic parameters to set for each plot
    ms = 8; %Markersize control for all plots
    lw = 2; %Linewidth control
    fs = 20; %Fontsize control
    
    %Subset the first failed profile using the BSL
    BSL_prof = DATA.data(DATA.data(:,DATA.iStn)==badCycles.N_cycles{badCycleIDX},:);

    %Subset the first failed profile using the RTQC tests
    if ~isempty(nitrateDATA.cycleDeltaTriggerStns)
        RTQC_prof = DATA.data(DATA.data(:,DATA.iStn)==nitrateDATA.cycleDeltaTriggerStns(1),:);
    end

    %PLOT 1: pH at 1500m with only the range check cycle marked
    %------------------------------------------------------------------------------------------------
    % f1 = figure;

    % plot(nitrateDATA.stnsDEEP,nitrateDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    %             'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette
    % 
    % %Check that there are any range check trigger stations at all, then plot as an xline
    % if ~isempty(nitrateDATA.rangeChkTriggerStns)
    %     xline(nitrateDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    % end 
    % 
    % set(gca,'TickDir','out');
    % xlim([badCycles.N_cycles{badCycleIDX}-2 inf])
    % if ~isempty(nitrateDATA.rangeChkTriggerStns) & (nitrateDATA.rangeChkTriggerStns(1)+2>badCycles.N_cycles{badCycleIDX}-2)
    %     xlim([badCycles.N_cycles{badCycleIDX}-2 nitrateDATA.rangeChkTriggerStns(1)+2])
    % end 
    % title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
    % ylabel("Nitrate adjusted");
    % xlabel("Station number");
    % set(gca, 'FontSize', fs);
    % savefig(f1, fileDir+"\float_RNGCHK\"+wmoList(c)+"_rtqcPlot_w_RNGCHK.fig");
    % close;
    %------------------------------------------------------------------------------------------------   

    %PLOT 2: pH at 1500m with range check cycle and BSL cycle marked
    %------------------------------------------------------------------------------------------------
    % f2 = figure;
    % plot(nitrateDATA.stnsDEEP,nitrateDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    %             'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette
    % 
    % %Make a single line where the data has been marked as "bad" during DMQC on the BSL
    % xline(badCycles.N_cycles{c}, '--k','LineWidth',lw)
    % 
    % %Check that there are any range check trigger stations at all, then plot as an xline
    % if ~isempty(nitrateDATA.rangeChkTriggerStns)
    %     xline(nitrateDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    % end 
    % 
    % set(gca,'TickDir','out');
    % xlim([badCycles.N_cycles{c}-2 inf])
    % if ~isempty(nitrateDATA.rangeChkTriggerStns) & (nitrateDATA.rangeChkTriggerStns(1)+2>badCycles.N_cycles{c}-2)
    %     xlim([badCycles.N_cycles{c}-2 nitrateDATA.rangeChkTriggerStns(1)+2])
    % end 
    % title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
    % ylabel("Nitrate adjusted");
    % xlabel("Station number");
    % set(gca, 'FontSize', fs);
    % savefig(f2, fileDir+"\float_RNGCHK_BSL\"+wmoList(c)+"_rtqcPlot_w_RNGCHKandBSL.fig");
    % close;
    % %------------------------------------------------------------------------------------------------
    % 
    % %PLOT 3: pH at 1500m with range check cycle, BSL cycle, and sensor drift cycle marked
    % %------------------------------------------------------------------------------------------------
    f3 = figure;
    subplot(2,1,1)
    plot(nitrateDATA.stnsDEEP,nitrateDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette

    %Make a single line where the data has been marked as "bad" during DMQC on the BSL
    xline(badCycles.N_cycles{badCycleIDX}, '--k','LineWidth',lw)

    %Mark the lines where the sensor drift test flags the data
    if ~isempty(nitrateDATA.cycleDeltaTriggerStns)
        xline(nitrateDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
    end

    %Check that there are any range check trigger stations at all, then plot as an xline
    if ~isempty(nitrateDATA.rangeChkTriggerStns)
        xline(nitrateDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    end 

    set(gca,'TickDir','out');
    xlim([badCycles.N_cycles{badCycleIDX}-2 inf])
    if ~isempty(nitrateDATA.rangeChkTriggerStns) & (nitrateDATA.rangeChkTriggerStns(1)+2>badCycles.N_cycles{badCycleIDX}-2)
        xlim([badCycles.N_cycles{badCycleIDX}-2 nitrateDATA.rangeChkTriggerStns(1)+2])
    end 
    title("Float "+wmoList(c)+": Adjusted nitrate at "+mean(depthRange)+" m");
    ylabel("Nitrate adjusted [\mumol/kg]");
    xlabel("Station number");
    set(gca, 'FontSize', fs);
    % savefig(f3, fileDir+"\float_RNGCHK_BSL_DRIFT\"+wmoList(c)+"_rtqcPlot_w_RNGCHKandBSLandDRIFT.fig");
    % close;
    %------------------------------------------------------------------------------------------------

    %PLOT 4: change in sequential pH data with all flags marked
    %------------------------------------------------------------------------------------------------
    if ~isempty(nitrateDATA.cycleDeltaTriggerStns)
        % f4 = figure;
        subplot(2,1,2)
        plot(nitrateDATA.stnsDEEP,nitrateDATA.cycleDelta,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
        
        %Make a single line where the data has been marked as "bad" during DMQC on the BSL
        xline(badCycles.N_cycles{badCycleIDX}, '--k','LineWidth',lw)
        
        %Mark the lines where the sensor drift test flags the data
        if ~isempty(nitrateDATA.cycleDeltaTriggerStns)
            xline(nitrateDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
        end
        
        %Add dashed lines to denote the thresholds of the test
        yline([-1 1], '--k','LineWidth', lw)
    
        %Check that there are any range check trigger stations at all, then plot as an xline
        if ~isempty(nitrateDATA.rangeChkTriggerStns)
            xline(nitrateDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
        end 
        
        set(gca,'TickDir','out');
        xlim([badCycles.N_cycles{badCycleIDX}-2 inf])
        if ~isempty(nitrateDATA.rangeChkTriggerStns) & (nitrateDATA.rangeChkTriggerStns(1)+2>badCycles.N_cycles{badCycleIDX}-2)
            xlim([badCycles.N_cycles{badCycleIDX}-2 nitrateDATA.rangeChkTriggerStns(1)+2])
        end 
        title("Float "+wmoList(c)+": difference between cycles at "+mean(depthRange)+" m");
        ylabel("\Delta Nitrate adjusted");	
        xlabel("Station number");
        set(gca, 'FontSize', fs);
        % savefig(f4, fileDir+"\cycleDelta\"+wmoList(c)+"_rtqcPlot_cycleDelta.fig");
        % close;
    end
end



