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

% wmoList = 5906439;
% idList = {'ua19531'};

%LOAD BSL AND IDENTIFY CYCLES FOR PH, P, T, AND S:
% BSL = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT");
BSL = readtable("C:\Users\lgrady\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.txt");
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

%OXYGEN
badCycles.O_cycles = cellfun(@str2double, regexp(BSL.CYCLES(contains(BSL.SENSOR,"O")&~contains(BSL.SENSOR,"CDOM")&BSL.FLAG==4), ...
    '\d+(?!\-\d+)(?=-)','match'), "UniformOutput",false);
badCycles.O_wmo = BSL(contains(BSL.SENSOR,"O")&~contains(BSL.SENSOR,"CDOM")&BSL.FLAG==4,"WMO_").WMO_;

%LOAD THE LIST OF PUMP OFFSET FLOATS:
pH_pumpoffset_980 

wmoList = badCycles.PH_wmo;
% wmoList = 1902380;
% wmoList = 1902304;
% wmoList = 5906439;
% wmoList = 1902303;
% wmoList = 4903365;
% wmoList = 5905108;
% wmoList = 5906533;
wmoList = 5906561;

%Do you want to plot the refAnom test?
refAnomTest = 1;

%Do you to only want to look at floats with BSL entries?
lookAtBSL = 0;

driftTestDiff = NaN(length(wmoList),1);
rangeTestDiff = NaN(length(wmoList),1);

%===========================================================================================
%================================ BEGIN REALTIME QC LOOP ===================================
%===========================================================================================

for c = 1:length(wmoList)

    onBSL = ismember(wmoList(c),badCycles.PH_wmo);
    if ~onBSL & (lookAtBSL==1)
        continue;
    end
    badCycleIDX = badCycles.PH_wmo==wmoList(c);

    %For keeping track of issues from specific floats
    disp(wmoList(c))
    
    if isnan(wmoList(c))
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

    %SUBSAMPLE DEEP PH
    %If the float is on the pump offset list, subsample to 950m 
    tests.pumpOffset = ismember(wmoList(c),pH_pumpoffset_980_floats);
    if tests.pumpOffset
        disp("Pump offset float")

        %Subsample the deep data to 950-980m
        depthRange = [920 980];
    else
        depthRange = [1480 1520];
    end

    %Check OXYGEN on the BSL
    BSLCycle_oxygen = cell2mat(badCycles.O_cycles(badCycles.O_wmo==wmoList(c)));
    
    %EXEPTION: this float will always have O2 issues and shouldn't be checked for missing O2
    tests.MUXERfloat = wmoList(c)==6990585; 

    %Create a test that marks a float as having bad O2 if it's on the BSL with bad O2 cycles, AND it isn't the MUXER issue float (6990585)
    badO2Test = (ismember(wmoList(c),badCycles.O_wmo) & ~isempty(BSLCycle_oxygen)) | tests.MUXERfloat;

    %NITRATE QC
    %===========================================================================================
    %Only runs through realtimeQC function if data passes ALL tests
    % disp('Processing NITRATE')
    % nitrateDATA = RTQC_tests(DATA, 5, depthRange, badO2Test, 0);
    %===========================================================================================
    startCycle = DATA.data(1,DATA.iStn);
    endCycle = DATA.data(end,DATA.iStn);
    cycleRange = [startCycle endCycle];

    %PH QC
    %===========================================================================================
    disp('Processing PH')
    %RTQC_tests((DATA, varType, qcDepth, badO2Test, omitBadDeepData,refAnomTest))
    phDATA = RTQC_tests(DATA, 3, depthRange, cycleRange, badO2Test, 0, refAnomTest); %Run the tests without omitting bad data (for fail testing)
                                                               %and without the refAnomTest
    %===========================================================================================

    % %Calculating the different between test cycles and DMQC
    % bslCyc = badCycles.PH_cycles{badCycleIDX};
    % 
    % if ~isempty(phDATA.rangeChkTriggerStns)
    %     rngChkDiff = phDATA.rangeChkTriggerStns - bslCyc;
    %     [rangeMin,rangeIDX] = min(abs(rngChkDiff));
    %     rangeTestDiff(c) = rngChkDiff(rangeIDX);
    % end
    % if ~isempty(phDATA.cycleDeltaTriggerStns)
    %     driftDiff = phDATA.cycleDeltaTriggerStns - bslCyc;
    %     [driftMin,driftIDX] = min(abs(driftDiff));
    %     driftTestDiff(c) = driftDiff(driftIDX);
    % end
    
    %Check if a directory to store plots already exists
    % userDir = getenv("USERPROFILE");
    % fileDir = userDir+"\Documents\MATLAB\ARGO_PROCESSING\MFILES\LG_WIP\RTQC\PH_RTQC_PLOTS\";
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
    % 
    % %Subset the first failed profile using the BSL
    % BSL_prof = DATA.data(DATA.data(:,DATA.iStn)==badCycles.PH_cycles{c},:);
    % 
    % %Subset the first failed profile using the RTQC tests
    % if ~isempty(phDATA.cycleDeltaTriggerStns)
    %     RTQC_prof = DATA.data(DATA.data(:,DATA.iStn)==phDATA.cycleDeltaTriggerStns(1),:);
    % end

    % %PLOT 1: pH at 1500m with only the range check cycle marked
    % %------------------------------------------------------------------------------------------------
    % % f1 = figure;
    % % % % subplot(2,1,1)
    % % plot(phDATA.stnsDEEP,phDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    % %             'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette
    % % 
    % % %Check that there are any range check trigger stations at all, then plot as an xline
    % % % if ~isempty(phDATA.rangeChkTriggerStns)
    % % %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    % % % end 
    % % 
    % % set(gca,'TickDir','out');
    % % % xlim([badCycles.PH_cycles{badCycleIDX}-2 inf])
    % % % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{badCycleIDX}-2)
    % % %     xlim([badCycles.PH_cycles{badCycleIDX}-2 phDATA.rangeChkTriggerStns(1)+2])
    % % % end
    % % xlim([0 18]);
    % % % ylim([7.4 8]);
    % % title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
    % % ylabel("pH adjusted");
    % % xlabel("Station number");
    % % set(gca, 'FontSize', fs);
    % % savefig(f1, fileDir+"\float_RNGCHK\"+wmoList(c)+"_rtqcPlot_w_RNGCHK.fig");
    % % close;
    % % % %------------------------------------------------------------------------------------------------   
    % % % 
    % % % %PLOT 2: pH at 1500m with range check cycle and BSL cycle marked
    % % % %------------------------------------------------------------------------------------------------
    % % f2 = figure;
    % % plot(phDATA.stnsDEEP,phDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    % %             'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette
    % % 
    % % %Make a single line where the data has been marked as "bad" during DMQC on the BSL
    % % xline(badCycles.PH_cycles{badCycleIDX}, '--k','LineWidth',lw)
    % % 
    % % %Check that there are any range check trigger stations at all, then plot as an xline
    % % if ~isempty(phDATA.rangeChkTriggerStns)
    % %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    % % end 
    % % 
    % % set(gca,'TickDir','out');
    % % xlim([badCycles.PH_cycles{badCycleIDX}-2 inf])
    % % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{badCycleIDX}-2)
    % %     xlim([badCycles.PH_cycles{badCycleIDX}-2 phDATA.rangeChkTriggerStns(1)+2])
    % % end 
    % % title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
    % % ylabel("pH adjusted");
    % % xlabel("Station number");
    % % set(gca, 'FontSize', fs);
    % % savefig(f2, fileDir+"\float_RNGCHK_BSL\"+wmoList(c)+"_rtqcPlot_w_RNGCHKandBSL.fig");
    % % close;
    % % % %------------------------------------------------------------------------------------------------
    % % % 
    % % % %PLOT 3: pH at 1500m with range check cycle, BSL cycle, and sensor drift cycle marked
    % % % %------------------------------------------------------------------------------------------------
    % % f3 = figure;
    % % plot(phDATA.stnsDEEP,phDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    % %             'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette
    % % 
    % % %Make a single line where the data has been marked as "bad" during DMQC on the BSL
    % % xline(badCycles.PH_cycles{badCycleIDX}, '--k','LineWidth',lw)
    % % 
    % % %Mark the lines where the sensor drift test flags the data
    % % if ~isempty(phDATA.cycleDeltaTriggerStns)
    % %     xline(phDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
    % % end
    % % 
    % % %Check that there are any range check trigger stations at all, then plot as an xline
    % % if ~isempty(phDATA.rangeChkTriggerStns)
    % %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    % % end 
    % % 
    % % set(gca,'TickDir','out');
    % % xlim([badCycles.PH_cycles{badCycleIDX}-2 inf])
    % % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{badCycleIDX}-2)
    % %     xlim([badCycles.PH_cycles{badCycleIDX}-2 phDATA.rangeChkTriggerStns(1)+2])
    % % end 
    % % title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
    % % ylabel("pH adjusted");
    % % xlabel("Station number");
    % % set(gca, 'FontSize', fs);
    % % savefig(f3, fileDir+"\float_RNGCHK_BSL_DRIFT\"+wmoList(c)+"_rtqcPlot_w_RNGCHKandBSLandDRIFT.fig");
    % % close;
    % % %------------------------------------------------------------------------------------------------
    % % 
    % % %PLOT 4: change in sequential pH data with all flags marked
    % % %------------------------------------------------------------------------------------------------
    % if ~isempty(phDATA.cycleDeltaTriggerStns) %& ~isempty(badCycles.PH_cycles{badCycleIDX}) & ~isempty(phDATA.rangeChkTriggerStns)
    %     f4 = figure;
    %     subplot(2,1,1);
    %     plot(phDATA.stnsDEEP,phDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    %         'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette
    % 
    %     %Make a single line where the data has been marked as "bad" during DMQC on the BSL
    %     % xline(badCycles.PH_cycles{badCycleIDX}, '--k','LineWidth',lw)
    % 
    %     %Mark the lines where the sensor drift test flags the data
    %     if ~isempty(phDATA.cycleDeltaTriggerStns)
    %         xline(phDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
    %     end
    % 
    %     %Check that there are any range check trigger stations at all, then plot as an xline
    %     % if ~isempty(phDATA.rangeChkTriggerStns)
    %     %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    %     % end 
    % 
    %     set(gca,'TickDir','out');
    %     % xlim([badCycles.PH_cycles{badCycleIDX}-2 inf])
    %     % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{badCycleIDX}-2)
    %     %     xlim([badCycles.PH_cycles{badCycleIDX}-2 phDATA.rangeChkTriggerStns(1)+2])
    %     % end 
    %     xlim([50 inf]);
    %     title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
    %     ylabel("pH adjusted");
    %     xlabel("Station number");
    %     set(gca, 'FontSize', fs);
    % 
    %     subplot(2,1,2);
    %     plot(phDATA.stnsDEEP,phDATA.cycleDelta,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    %                 'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880])
    % 
    %     %Make a single line where the data has been marked as "bad" during DMQC on the BSL
    %     % xline(badCycles.PH_cycles{c}, '--k','LineWidth',lw)
    % 
    %     %Mark the lines where the sensor drift test flags the data
    %     if ~isempty(phDATA.cycleDeltaTriggerStns)
    %         xline(phDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
    %     end
    % 
    %     %Add dashed lines to denote the thresholds of the test
    %     yline([.01 -.01], '--k','LineWidth', lw)
    % 
    %     %Check that there are any range check trigger stations at all, then plot as an xline
    %     % if ~isempty(phDATA.rangeChkTriggerStns)
    %     %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
    %     % end 
    % 
    %     set(gca,'TickDir','out');
    %     % xlim([badCycles.PH_cycles{c}-2 inf])
    %     % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{c}-2)
    %     %     xlim([badCycles.PH_cycles{c}-2 phDATA.rangeChkTriggerStns(1)+2])
    %     % end 
    %     xlim([50 inf]);
    %     % ylim([-.15 .1]);
    %     title("Float "+wmoList(c)+": difference between cycles at "+mean(depthRange)+" m");
    %     ylabel("\Delta pH adjusted");	
    %     xlabel("Station number");
    %     set(gca, 'FontSize', fs);
    %     % savefig(f4, fileDir+"\cycleDelta\"+wmoList(c)+"_rtqcPlot_cycleDelta.fig");
    %     % close;
    % end
    % % ------------------------------------------------------------------------------------------------
    % 
    % % %PLOT 5: Reference anomaly
    % % %------------------------------------------------------------------------------------------------
    if refAnomTest==1
        disp("we're getting here")
        f5 = figure;
        subplot(2,1,1);
        plot(phDATA.stnsDEEP,phDATA.deepDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
            'MarkerEdgeColor','k','MarkerFaceColor', [0 130 245]/255) %Plots light blue of SAGE color pallette

        hold on;

        plot(phDATA.stnsDEEP,phDATA.RefDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
            'MarkerEdgeColor','k','MarkerFaceColor', [1 .2 .2]) %Plots light blue of SAGE color pallette

        %Make a single line where the data has been marked as "bad" during DMQC on the BSL
        % xline(badCycles.PH_cycles{badCycleIDX}, '--k','LineWidth',lw)

        %Mark the lines where the sensor drift test flags the data
        if ~isempty(phDATA.cycleDeltaTriggerStns)
            xline(phDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
        end

        %Check that there are any range check trigger stations at all, then plot as an xline
        % if ~isempty(phDATA.rangeChkTriggerStns)
        %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
        % end 

        set(gca,'TickDir','out');
        % xlim([badCycles.PH_cycles{badCycleIDX}-2 inf])
        % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{badCycleIDX}-2)
        %     xlim([badCycles.PH_cycles{badCycleIDX}-2 phDATA.rangeChkTriggerStns(1)+2])
        % end
        xlim([50 inf]);
        title("Float "+wmoList(c)+": Adjusted pH at "+mean(depthRange)+" m");
        ylabel("pH adjusted");
        xlabel("Station number");
        set(gca, 'FontSize', fs);

        subplot(2,1,2);
        plot(phDATA.stnsDEEP,phDATA.refAnom,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[0.4940 0.1840 0.5560])

        %Make a single line where the data has been marked as "bad" during DMQC on the BSL
        % xline(badCycles.PH_cycles{c}, '--k','LineWidth',lw)

        %Mark the lines where the sensor drift test flags the data
        if ~isempty(phDATA.cycleDeltaTriggerStns)
            xline(phDATA.cycleDeltaTriggerStns(1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',lw)
        end

        %Add dashed lines to denote the thresholds of the test
        yline([-.012 .012], '--k','LineWidth', lw)

        % %Check that there are any range check trigger stations at all, then plot as an xline
        % if ~isempty(phDATA.rangeChkTriggerStns)
        %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw)
        % end 

        set(gca,'TickDir','out');
        % xlim([badCycles.PH_cycles{c}-2 inf])
        % if ~isempty(phDATA.rangeChkTriggerStns) & (phDATA.rangeChkTriggerStns(1)+2>badCycles.PH_cycles{c}-2)
        %     xlim([badCycles.PH_cycles{c}-2 phDATA.rangeChkTriggerStns(1)+2])
        % end 
        % xlim([90 106]);
        ylim([-.014 .014]);
        xlim([50 inf]);
        title("Float "+wmoList(c)+": reference anomaly at "+mean(depthRange)+" m");
        ylabel("Float - ESPER pH");	
        xlabel("Station number");
        set(gca, 'FontSize', fs);
    %     % savefig(f4, fileDir+"\cycleDelta\"+wmoList(c)+"_rtqcPlot_cycleDelta.fig");
    %     % close;
    end
    %------------------------------------------------------------------------------------------------

    % Plot pH profiles with the bad cycle marked in red
    % [0.5 1 0.2]
    % subplot(2,2,[1,3]);
    % scatter(DATA.data(:,DATA.iPH), DATA.data(:,DATA.iP), 16, DATA.data(:,DATA.iStn), 'filled')
    % hold on
    % % scatter(BSL_prof(:,DATA.iPH),BSL_prof(:,DATA.iP), ms,'k','filled')
    % % hold on
    % % if ~isempty(phDATA.cycleDeltaTriggerStns)
    % %     scatter(RTQC_prof(:,DATA.iPH),RTQC_prof(:,DATA.iP), 16,'g','filled')
    % % end
    % 
    % if ~isempty(phDATA.rangeChkTriggerStns)
    %     scatter(RTQC_prof(:,DATA.iPH),RTQC_prof(:,DATA.iP), 16,'r','filled')
    % end
    % set(gca, 'YDir','reverse');
    % %set(gca, 'xlim',[4 8.5]);
    % set(gca, 'ylim',[0 2000]);
    % set(gca,'TickDir','out');
    % title("Profiles for float: "+wmoList(c)+" ("+badCycles.PH_mid{c}+")");
    % ylabel("Pressure [dbar]");
    % xlabel("pH");
    % col = colorbar;
    % col.Label.String = 'Profile #';
    % col.Label.FontSize = 10;
    % col.Label.Rotation = -90;
    % set(gca, 'FontSize', fs);
        % xline(badCycles.PH_cycles{c}, '--k','LineWidth', lw,'DisplayName', 'DMQC failure point')
    % if ~isempty(phDATA.cycleDeltaTriggerStns)
    %     xline(phDATA.cycleDeltaTriggerStns(1), '.-g','LineWidth', lw,'DisplayName', '1st Sensor drift flag')
    % end
    %Plot 
    % subplot(2,2,2);


    %Plot the difference in pH per cycle
    % subplot(2,2,4);
    % 
    % %set(gca, 'ylim',[-.05 .05]);
    % 
    % 
    % % if ~isempty(phDATA.rangeChkTriggerStns)
    % %     xline(phDATA.rangeChkTriggerStns(1),'--r','LineWidth',lw,'DisplayName','1st Range check flag')
    % % end
    % set(gca,'TickDir','out');
    % title("Mean drift between cycles from "+depthRange(1)+"-"+depthRange(2)+"m");
    % ylabel("Difference between cycles");
    % xlabel("Station number");
    % % legend({'Float' 'DMQC failure point' '1st Sensor drift flag' '1st Range check flag'}, 'Location','bestoutside');
    % set(gca, 'FontSize', fs);
    % 
    % 
    % savefig(f, userDir+"\Documents\MATLAB\ARGO_PROCESSING\MFILES\LG_WIP\RTQC\PH_RTQC_PLOTS\"+string(wmoList(c))+"_rtqcPlot.fig");
    % % saveas(f,userDir+"\Documents\MATLAB\ARGO_PROCESSING\MFILES\LG_WIP\RTQC\PH_RTQC_PLOTS\"+string(wmoList(c))+"_rtqcPlot.png")
    % close;
    % 
    % 
    % f2 = figure;
    % plot(phDATA.stnsDEEP,phDATA.RefDATA_mean,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
    %         'MarkerEdgeColor','k','MarkerFaceColor',[1 .2 0.2])
    % % xline(badCycles.PH_cycles{c}, '--k','LineWidth', lw,'DisplayName', 'DMQC failure point')
    % if ~isempty(phDATA.cycleDeltaTriggerStns)
    %     xline(phDATA.cycleDeltaTriggerStns(1), '.-g','LineWidth', lw,'DisplayName', '1st Sensor drift flag')
    % end
    % % if ~isempty(phDATA.rangeChkTriggerStns)
    % %     xline(phDATA.rangeChkTriggerStns(1), '--r','LineWidth',lw,'DisplayName','1st Range check flag')
    % % end
    % %set(gca, 'ylim',[7.0 8.5]);
    % set(gca,'TickDir','out');
    % title("Mean ESPER pH from "+depthRange(1)+"-"+depthRange(2)+"m");
    % ylabel("pH");
    % xlabel("Station number");
    % % legend({'Float' 'DMQC failure point' '1st Sensor drift flag' '1st Range check flag'}, 'Location','bestoutside');
    % set(gca, 'FontSize', fs);
end



