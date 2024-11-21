%Get user profile
userDir = getenv("USERPROFILE");

%Get the FLOATS directory for .mat files for calculating delta

%Load the BSL as a table
BSL = readtable(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT");
BSL_srt = sortrows(BSL,"WMO_","ascend"); %Sort by WMO's in ascending order to match floatsDir
BSL_ph= BSL_srt(contains(BSL_srt.SENSOR,"PH"),:); %Isolate lines with PH
BSL_ph_fail = BSL_ph(BSL_ph.FLAG==4,:); %Isolate PH flagged as fail (4)

%Create new column with all bad cycle numbers as doubles
BSL_ph_fail.BADCYCLENUM = cellfun(@str2double, regexp(BSL_ph_fail.CYCLES,'^\d+|(?<=,)\d+','match'), "UniformOutput",false);
%BSL_ph_fail.BADCYCLENUM = cellfun(@str2double, regexp(BSL_ph_fail.CYCLES,'\d+(?=-)','match'), "UniformOutput",false);

deepThd = .01;

% BEGIN LOOP
for f = 1:length(BSL_ph_fail.WMO_)
    try
        [a, MSGID] = lastwarn();
        warning('off', MSGID)
        
        %Call float WMO from FLOATS directory
        wmoID = BSL_ph_fail.WMO_(f);
        wmoDir = dir(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\"+wmoID); %Load WMO directory
        wmoTable = struct2table(wmoDir); %Convert to table for easy indexing
        firstFile = 3; %First .mat file always starts on 3
        lastFile = length(find(~contains(wmoTable.name,'cal'))); %Last file is the profile before calxxx file
        
        %Subsample the directory so only the profile files are included
        wmoDir = wmoDir(firstFile:lastFile,:);
        
        %Initialize empty matrix to store points that go beyond delta threshold
        anomMat = NaN(length(wmoDir),0); 
        anomCnt = 0;
        
        %Load in the floatViz data for the WMO
        floatVizDir = userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\QC\"+string(wmoID)+"QC.TXT";
        floatVizData = readtable(floatVizDir);
        % 
        %Subsample all good and failing ph, pressure, and station number data
        %at 1500m for future plots
        pH_qc = floatVizData.pHinsitu_Total_((floatVizData.QF_12==0|floatVizData.QF_12==8)& ...
            (floatVizData.pHinsitu_Total_<99999&floatVizData.pHinsitu_Total_>-99999));
        pres_qc = floatVizData.Pressure_dbar_((floatVizData.QF_12==0|floatVizData.QF_12==8)& ...
            (floatVizData.pHinsitu_Total_<99999&floatVizData.pHinsitu_Total_>-99999));
        stn_qc = floatVizData.Station((floatVizData.QF_12==0|floatVizData.QF_12==8)& ...
            (floatVizData.pHinsitu_Total_<99999&floatVizData.pHinsitu_Total_>-99999));
        ph1500 = pH_qc(pres_qc<=1500 & pres_qc>=950);
        stn1500 = stn_qc(pres_qc<=1500 & pres_qc>=950);

        phMat = NaN(length(wmoDir));

        %Cycle through all .mat files in the directory up to the
        %second-to-last profile
        for i=2:length(wmoDir)
            
            %Load the "current" profile and subsample data point at 1500m
            currentProf = load ([wmoDir(1).folder,'\',wmoDir(i).name]);
            currentPH_depth = currentProf.LR.PH_IN_SITU_TOTAL_ADJUSTED(currentProf.LR.PRES_ADJUSTED>=950 & currentProf.LR.PRES_ADJUSTED<=1500);
            currentPH_depth = currentPH_depth(currentPH_depth<99999);
            
            %Do the same for the profile before
            prevProf = load ([wmoDir(1).folder,'\',wmoDir(i-1).name]);
            prevPH_depth = prevProf.LR.PH_IN_SITU_TOTAL_ADJUSTED(prevProf.LR.PRES_ADJUSTED>=950 & prevProf.LR.PRES_ADJUSTED<=1500);
            prevPH_depth = prevPH_depth(prevPH_depth<99999);

            %Store the previous mean of deep PH into the matrix
            phMat(i-1) = mean(prevPH_depth);

            %Find the difference between the two
            delta_depth = abs(mean(prevPH_depth)-mean(currentPH_depth));
            
            delta_threshold_test = delta_depth>=deepThd;
            high_ph_test = mean(currentPH_depth)>=8.5;
            low_ph_test = mean(currentPH_depth)<=7.5;

            test_sum = delta_threshold_test+high_ph_test+low_ph_test;

            %If any of the tests fail increase anomaly counter by 1
            if test_sum >= 1 
                anomCnt = anomCnt + 1;
                anomMat(i,1) = anomCnt;

            else %Otherwise reset anomaly counter back to 0
                anomCnt = 0;
                anomMat(i,1) = anomCnt;
            end
        end
        
        %Fill the last slot in the ph matrix 
        phMat(i) = mean(currentPH_depth);
        
        %Check whether the anomaly matrix has consecutive triggers up to 3
        if ~isempty(find(anomMat == 3))
    
            % Estimate the bad cycle based on the anomaly counter reaching 3,
            % then choose the cycle just before that data point
            % estBadCycle = find(anomMat == 3)-1;
            BSL_ph_fail.EST_CYCLE_NUM(f,:) = {find(anomMat == 3)-2};
            %Use the floatViz data to isolate the data at the estimated bad cycle
            %EstBadProf = floatVizData(floatVizData.Station==estBadCycle(1),:);

            %Do the same using the bad cycles marked in the BSL
            BadProf = pH_qc(stn_qc==(BSL_ph_fail.BADCYCLENUM{f}(end)));
            BadPres = pres_qc(stn_qc==(BSL_ph_fail.BADCYCLENUM{f}(end)));

            %Attempt to make the plots
            fig = figure;
            set(fig,'Visible','off')

            % Plot pH profiles with the bad cycle marked in red
            subplot(2,2,[1,3]);
            scatter(pH_qc, pres_qc, 8, stn_qc, 'filled');
            hold on;
            scatter(BadProf,BadPres, 20,'k','filled');
            % hold on;
            % scatter(EstBadProf.pHinsitu_Total_,EstBadProf.Pressure_dbar_, 10,'r','filled');
            set(gca, 'YDir','reverse');
            %set(gca, 'xlim',[4 8.5]);
            set(gca, 'ylim',[0 2000]);
            set(gca,'TickDir','out');
            title("pH Profiles");
            ylabel("Pressure [dbar]");
            xlabel("pH");
            c = colorbar;
            c.Label.String = 'Profile #';
            c.Label.FontSize = 10;
            c.Label.Rotation = -90;

            %Plot pH at 1500m
            subplot(2,2,2);
            plot(phMat,"b.-",'LineWidth',1,'MarkerSize',10);
            xline(BSL_ph_fail.BADCYCLENUM{f}(:), 'k-',LineWidth=2);
            xline(BSL_ph_fail.EST_CYCLE_NUM{f}(:), '--r',LineWidth=2)
            set(gca,'TickDir','out');
            title("pH");
            ylabel("pH");
            xlabel("Station number");

            %Plot the difference in pH per cycle
            subplot(2,2,4);
            plot(diff(phMat),"b.-",'LineWidth',1,'MarkerSize',10);
            set(gca, 'ylim',[-.05 .05]);
            xline(BSL_ph_fail.BADCYCLENUM{f}(:), 'k-',LineWidth=2);
            xline(BSL_ph_fail.EST_CYCLE_NUM{f}(:), '--r',LineWidth=2)
            yline([deepThd -deepThd], '--k');
            set(gca,'TickDir','out');
            title("pH delta");
            ylabel("pH");
            xlabel("Station number");
            
        else

            BSL_ph_fail.EST_CYCLE_NUM(f,:) = NaN;
            %Do the same using the bad cycles marked in the BSL
            BadProf = pH_qc(stn_qc==(BSL_ph_fail.BADCYCLENUM{f}(end)));
            BadPres = pres_qc(stn_qc==(BSL_ph_fail.BADCYCLENUM{f}(end)));

            %Attempt to make the plots as well
            fig = figure;
            set(fig,'Visible','off')

            % Plot pH profiles with the bad cycle marked in red
            subplot(2,2,[1,3]);
            scatter(pH_qc, pres_qc, 8, stn_qc, 'filled');
            hold on;
            scatter(BadProf,BadPres, 20,'k','filled');
            set(gca, 'YDir','reverse');
            %set(gca, 'xlim',[4 8.5]);
            set(gca, 'ylim',[0 2000]);
            set(gca,'TickDir','out');
            title("pH Profiles");
            ylabel("Pressure [dbar]");
            xlabel("pH");
            c = colorbar;
            c.Label.String = 'Profile #';
            c.Label.FontSize = 10;
            c.Label.Rotation = -90;

            %Plot pH at 1500m
            subplot(2,2,2);
            plot(phMat,"b.-",'LineWidth',1,'MarkerSize',10);
            xline(BSL_ph_fail.BADCYCLENUM{f}(:), 'k-',LineWidth=2);
            set(gca,'TickDir','out');
            title("pH");
            ylabel("pH");
            xlabel("Station number");

            %Plot the difference in pH per cycle
            subplot(2,2,4);
            plot(diff(phMat),"b.-",'LineWidth',1,'MarkerSize',10);
            xline(BSL_ph_fail.BADCYCLENUM{f}(:), 'k-',LineWidth=2);
            yline([.01 -.01], '--k');
            set(gca,'TickDir','out');
            title("pH delta");
            ylabel("pH");
            xlabel("Station number");
        end
    
        saveas(fig, userDir+"\Documents\realtime_ph_qc\Plots\"+string(wmoID)+".png");

    catch ME
        disp(wmoID)
        disp(ME.message)
    end

end

writetable(BSL_ph_fail,userDir+"\Documents\realtime_ph_qc\realtimePH.TXT")
% fclose(ph_log);
% 
% cycleDat = readtable ("C:\Users\lgrady\Documents\realtime_ph_qc\realtimePH.TXT");
% 
% deltaDat = cycleDat.bad_cycles - cycleDat.est_bad_cycles;

% fig2 = figure;
% histogram(deltaDat, 20);
% title("Fail cycle delta");
% xlabel("BSL-estimated cycle");
% saveas(fig2, userDir+"\Documents\realtime_ph_qc\FailCycleHist.png");