function[] = plot_GDAC_profile_timeseries(PLOT_type)
    %=====================================================================================================================
    % OBJECTIVE: This code uses the SPROF index text file, available on the GDAC, and counts the number of profiles
    %            from floats off of the MBARI_float_list from every year. There are options to plot only the total
    %            profiles, and also the profiles by GO-BGC, SOCCOM, and OTHER.
    %
    % INPUTS:
    %            SPROF_path: the filepath for "argo_synthetic-profile_index.txt" on the users machine.
    %
    %            PLOT_type: A binary string that's either "ALL" to plot the cumulative yearly profiles for each program,
    %                       or "TOTAL" to only plot the cumulative yearly sum all programs together.
    %
    % OUTPUTS:   A plot with a yearly cumulative sum of all programs together or each program individually.
    %=====================================================================================================================
    
    %TO USE
    %=====================================================================================================================
    %Input the desired plot type; in this case, all programs will be plotted
    %   IF PLOT_type = TOTAL: all programs will be plotted
    %   IF PLOT_type = BY_PROGRAM: the yearly cumulative sum for each program will be plotted.
    %=====================================================================================================================
    %Retrieve directory that the file will go into
    userDir = getenv('USERPROFILE');
    userDir = userDir+"\Documents\";
    
    %Identify the directory of the desired file on the GDAC
    webtarget = "https://data-argo.ifremer.fr";
    web_dir  = "/";
    fname = "argo_synthetic-profile_index.txt";
    
    %Full name of the file url
    dataUrl = webtarget+web_dir+fname;
    
    %Name of the file to be put into your local directory
    sProfFile = userDir+"argo_synthetic-profile_index.txt";
    
    %Save the file into userDir as sProfFile
    SPROF_path = websave(sProfFile,dataUrl);
    
    %LOAD FLOATLIST
    floatList = readtable("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.txt");
    
    %PARSE SPROF INDEX
    %--------------------------------------------------------------------------------------------------------
    %LOAD SPROF INDEX (find a way to automate this)
    argoDat = readtable(SPROF_path);
    
    %Split the SPROF INDEX filenames by '/' and isolate the 7-digit WMO_ID
    fileSplit = split(argoDat.file,'/');
    
    %Create a seperate column in the data table for WMOs that uses the second column of the split string
    argoDat.WMO=str2double(fileSplit(:,2));
    
    %Now that SPROF index dataset has WMO column, we can subset all profiles that share a WMO with the MBARI floatlist
    isMBARIfloat = ismember(argoDat.WMO,floatList.WMO); %Tests whether each row in argoDat contains the WMO from MBARI floatlist
    MBARIfloats = argoDat(isMBARIfloat,:); %Subset argoDat using the above logical array as the index
    
    %Now that the dataset is much smaller and more manageable, create a column for datetimes so the profile counts can be split by year
    MBARIfloats.date = num2str(MBARIfloats.date); %Converts the original date column to str
    MBARIfloats.datetime = datetime(MBARIfloats.date, 'InputFormat', 'yyyyMMddHHmmss'); %Creates a new column that converts date str into datetime
    %---------------------------------------------------------------------------------------------------------
    
    %If the desired plot is "BY_PROGRAM", then go through the following processes of splitting the floatlist and SPROF index by program and finding
    %the yearly cumulative sum for each. If PLOT_type is "TOTAL", then these steps will be skipped, and only the total sum will be calculated.
    if PLOT_type == "BY_PROGRAM"
    
        %INDEX WMOs BASED ON FLOAT PROGRAM
        %---------------------------------------------------------------------------------------------------------
        
        %Create indices of WMOs from GOBGC, SOCCOM, and OTHER programs from the MBARI floatlist
        gobgcIDX = contains(floatList.Program,'GO-BGC');
        soccomIDX = contains(floatList.Program,'SOCCOM');
        otherIDX = ~contains(floatList.Program,'GO-BGC')&~contains(floatList.Program,'SOCCOM');
        
        %Call the WMOs from the program indices
        GOBGCfloats = floatList.WMO(gobgcIDX);
        SOCCOMfloats = floatList.WMO(soccomIDX);
        OTHERfloats = floatList.WMO(otherIDX);
        %---------------------------------------------------------------------------------------------------------
        
        %RETREIVE ALL PROGRAM SPECIFIC PROFILES FROM SPROF INDEX USING ABOVE INDICES
        %---------------------------------------------------------------------------------------------------------
        %GO-BGC
        isGOBGCfloat = ismember(MBARIfloats.WMO,GOBGCfloats); %Tests whether each row in argoDat contains the WMO from GOBGC floats
        GOBGCprofs = MBARIfloats(isGOBGCfloat,:); %Subset argoDat using the above logical array as the index
        
        %SOCCOM: does the same as above but using SOCCOM floats
        isSOCCOMfloat = ismember(MBARIfloats.WMO,SOCCOMfloats);
        SOCCOMprofs = MBARIfloats(isSOCCOMfloat,:);
        
        %OTHER: does the same as above but using floats from non-GOBGC and non-SOCCOM floats
        isOTHERfloat = ismember(MBARIfloats.WMO,OTHERfloats);
        OTHERprofs = MBARIfloats(isOTHERfloat,:);
        %---------------------------------------------------------------------------------------------------------
        
        %CALCULATE YEARLY CUMULATIVE SUM OF PROFILES FOR EACH PROGRAM AND THE CUMULATIVE TOTAL
        %---------------------------------------------------------------------------------------------------------
        %GOBGC:
        %Group the profiles by year
        [GOBGCprofs_byYr, GOBGC_YR] = groupcounts(year(GOBGCprofs.datetime));
        
        %Find the cumulative sum of floats for every year
        GOBGCprofs_byYr = cumsum(GOBGCprofs_byYr);
        
        %SOCCOM:
        [SOCCOMprofs_byYr, SOCCOM_YR] = groupcounts(year(SOCCOMprofs.datetime));
        SOCCOMprofs_byYr = cumsum(SOCCOMprofs_byYr);
        
        %OTHER:
        [OTHERprofs_byYr, OTHER_YR] = groupcounts(year(OTHERprofs.datetime));
        OTHERprofs_byYr = cumsum(OTHERprofs_byYr);
    end
    
    %CALCULATE THE YEARLY CUMULATIVE SUM FOR ALL PROGRAMS TOGETHER
    %TOTAL:
    [TOTALprofs_byYr, TOTAL_YR] = groupcounts(year(MBARIfloats.datetime));
    TOTALprofs_byYr = cumsum(TOTALprofs_byYr);
    
    %PLOTTING
    %---------------------------------------------------------------------------------------------------------
    fig = figure;
    
    hold on;
    
    ms = 10; %Markersize control for all plots
    lw = 2; %Linewidth control
    fs = 16; %Fontsize control
    if PLOT_type == "BY_PROGRAM"
        %TOTAL
        bar(TOTAL_YR,TOTALprofs_byYr,'FaceColor',[.75 .75 .75])
        
        %GOBGC
        plot(GOBGC_YR,GOBGCprofs_byYr,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[0.5 1 0.2]) %rgb
        
        %SOCCOM
        plot(SOCCOM_YR,SOCCOMprofs_byYr,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[0.2 0.5 1])
        %OTHER
        plot(OTHER_YR,OTHERprofs_byYr,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.5])
        
        title("Profile counts by program (from GDAC)");
        legend({'TOTAL' 'GO-BGC', 'SOCCOM', 'OTHER'},'Location','bestoutside');
        ax.XLim = [2006 GOBGC_YR(end)+1];
    
    elseif PLOT_type == "TOTAL"
        plot(TOTAL_YR,TOTALprofs_byYr,'-o','MarkerSize',ms,'LineWidth',lw,'Color','k', ...
                    'MarkerEdgeColor','k','MarkerFaceColor','k') %rgb
        title("Total profile counts (from GDAC)");
        ax.XLim = [2006 TOTAL_YR(end-1)+1];
    end
    
    
    ax = gca;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.FontSize = fs;
    ax.YAxis.Exponent = 3;
    xline(max(xlim));
    yline(max(ylim));
    
    ylabel("Profiles (in thousands)",'FontSize', fs);
    
    clearvars
    %---------------------------------------------------------------------------------------------------------
end
