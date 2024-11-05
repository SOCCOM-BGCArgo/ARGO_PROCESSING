function [allInfo] = MBARI_floatCount_summary()
    % MBARI_floatCount_summary.m
    %
    % OBJECTIVES OF THIS SCRIPT:
    %   Using the MBARI-float-list, generate a summary including the following:
    %   -Total deployed floats (for each SOCCOM, GOBGC, ALL)
    %   -Total active floats (for each SOCCOM, GOBGC, ALL)
    %   -Total floats across all platforms (APEX, NAVIS, SOLO)
    %   -Total floats across all regions (NAtl, NPac, Eq, etc.)
    %
    %--------------------------------------------------------------------------
    % Temporary script started by T. Maurer 4/24/24 for L. Grady
    % 
    % LG 5/16/2024 Program counts all floats (total and active only) by program
    % and region, combines the values, and assembles a formatted string array
    % called "allInfo". This string array is sent out as
    % 
    % LG 8/1/2024 Added a column to the email string that counts the total 
    % number of profiles for SOCCOM, GO-BGC, and all other programs
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------------------
    % FUNCTIONS
    %--------------------------------------------------------------------------------------
    
    %For counting ALL floats by program
    %INPUT(S): A program string such as 'SOCCOM' or 'GO-BGC'
    %OUTPUT(S): A count of apex, navis, solo, and combined total floats
    function [apex,navis,solo,all,apexAct,navisAct,soloAct,allAct,totalProfs] = floatCount_byProg(var)
        
        fltData = load('\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.mat'); %load 'MBARI_float_list' (contains full listing of all MBARI-processed floats)
        idx = find(strcmp('Program',fltData.d.hdr) == 1);
        iFLT  = find(strcmp('float type',fltData.d.hdr) == 1);
        iDATE = find(strcmp('1st date',fltData.d.hdr) == 1);
        iSTATUS  = find(strcmp('tf Dead',fltData.d.hdr) == 1);
        iCYC = find(strcmp('max cycle proc',fltData.d.hdr) == 1);
        depData = fltData.d.list(~isnan(cell2mat(fltData.d.list(:,iDATE))),:);
       
        apex = length(depData(matches(depData(:,idx),var) & matches(depData(:,iFLT),"APEX")));
        navis = length(depData(matches(depData(:,idx),var) & matches(depData(:,iFLT),"NAVIS")));
        solo = length(depData(matches(depData(:,idx),var) & matches(depData(:,iFLT),"SOLO")));
        all = length(depData(matches(depData(:,idx),var)));
        totalProfs = nansum(cell2mat(depData(matches(depData(:,idx),var),iCYC)));
             
        dActive = fltData.d.list(cell2mat(fltData.d.list(:,iSTATUS))==0,:); %subset list to all active floats across all programs.
        dActive = dActive(~isnan(cell2mat(dActive(:,iDATE))),:); %Exclude floats that haven't been deployed (without a 1st date)
    
        apexAct = length(dActive(matches(dActive(:,idx),var) & matches(dActive(:,iFLT),"APEX")));
        navisAct = length(dActive(matches(dActive(:,idx),var) & matches(dActive(:,iFLT),"NAVIS")));
        soloAct = length(dActive(matches(dActive(:,idx),var) & matches(dActive(:,iFLT),"SOLO")));
        allAct = length(dActive(matches(dActive(:,idx),var)));
    end
    
    %For counting ALL floats by region
    %INPUT(S): A program string such as 'NAtl' or 'NPac'
    %OUTPUT(S): A count of apex, navis, solo, and combined total floats
    function [apex,navis,solo,all,apexAct,navisAct,soloAct,allAct] = floatCount_byRegion(var)
    
        fltData = load('\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.mat'); %load 'MBARI_float_list' (contains full listing of all MBARI-processed floats)
        idx = find(strcmp('Region',fltData.d.hdr) == 1);
        iFLT  = find(strcmp('float type',fltData.d.hdr) == 1);
        iDATE = find(strcmp('1st date',fltData.d.hdr) == 1);
        iSTATUS  = find(strcmp('tf Dead',fltData.d.hdr) == 1);

        depData = fltData.d.list(~isnan(cell2mat(fltData.d.list(:,iDATE))),:);
        dActive = fltData.d.list(cell2mat(fltData.d.list(:,iSTATUS))==0,:); %subset list to all active floats across all programs.
        dActive = dActive(~isnan(cell2mat(dActive(:,iDATE))),:); %Exclude floats that haven't been deployed (without a 1st date)
    
        apex = length(depData(matches(depData(:,idx),var) & matches(depData(:,iFLT),"APEX")));
        navis = length(depData(matches(depData(:,idx),var) & matches(depData(:,iFLT),"NAVIS")));
        solo = length(depData(matches(depData(:,idx),var) & matches(depData(:,iFLT),"SOLO")));
        all = length(depData(matches(depData(:,idx),var)));

        apexAct = length(dActive(matches(dActive(:,idx),var) & matches(dActive(:,iFLT),"APEX")));
        navisAct = length(dActive(matches(dActive(:,idx),var) & matches(dActive(:,iFLT),"NAVIS")));
        soloAct = length(dActive(matches(dActive(:,idx),var) & matches(dActive(:,iFLT),"SOLO")));
        allAct = length(dActive(matches(dActive(:,idx),var)));
    end
    
    %--------------------------------------------------------------------------------------
    % LOADING DATA AND COUNTING FLOATS
    %--------------------------------------------------------------------------------------
    
    %Load the MBARI float list from chem for a directory of all floats
    data = load('\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.mat');
    %Identify the columns with program and region names
    iPROG = find(strcmp('Program',data.d.hdr) == 1);
    iREG = find(strcmp('Region',data.d.hdr) == 1);
    
    %Create a unique list of all program names
    varsProg = unique(data.d.list(:,iPROG)).';
    
    %Run the program float count function for each program in 'varsProg'
    [aProg,nProg,sProg,tProg,aProgActive,nProgActive,sProgActive,tProgActive,totalProfs] = cellfun(@floatCount_byProg,varsProg);
    %[aProgActive,nProgActive,sProgActive,tProgActive] = cellfun(@floatCount_byProgActive,varsProg);
    
    %Calculate the sum of all MBARI floats (dead and active)
    allFloats = sum(tProg);
    allFloatsActive = sum(tProgActive);
    
    %Create indices for identifying "other" and "SOCCOM and GOBGC" floats
    otherIDX = [1 3 4 5 7 8];
    SOCandGBC = [2 6];

    %Non SOCCOM and GO-BGC programs are lumped into 'OTHER' category
    a1= sum(aProg(otherIDX)); %All apex
    a2 = sum(aProgActive(otherIDX)); %Active apex
    a3 = sum(nProg(otherIDX)); %Navis
    a4 = sum(nProgActive(otherIDX));
    a5 = sum(sProg(otherIDX)); %Solo
    a6 = sum(sProgActive(otherIDX));
    a7 = sum(tProg(otherIDX)); %All floats
    a8 = sum(tProgActive(otherIDX));
    a9 = sum(totalProfs(otherIDX));

    %Repeat with a list of regions
    varsReg = unique(data.d.list(:,iREG)).';
    [aReg,nReg,sReg,tReg,aRegActive,nRegActive,sRegActive,tRegActive] = cellfun(@floatCount_byRegion,varsReg);
    %[aRegActive,nRegActive,sRegActive,tRegActive] = cellfun(@floatCount_byRegionActive,varsReg);
    
    %--------------------------------------------------------------------------------------
    % LIST CONFIGURATION
    %--------------------------------------------------------------------------------------
    
    %Create empty list to be populated with float counts
    FltCntList = {};

    %Create the initial text to be displayed at the top of the message
    FltCntList = [FltCntList; ["TOTAL Floats Managed by MBARI:   "+ allFloats+"("+allFloatsActive+")", "", "", "", "",""]];
    FltCntList = [FltCntList; ["TOTAL Profiles Recorded from Floats:   "+ sum(totalProfs), "", "", "", "",""]]; 
    FltCntList = [FltCntList; ["*** PARENTHESES INDICATE ACTIVE FLOATS ***", "", "", "", "",""]];
    FltCntList = [FltCntList; ["--------------------------------------------------", "", "", "", "",""]];
    
    %Generate headers
    FltCntList= [FltCntList; ["PROGRAM", "APEX", "NAVIS", "SOLO", "TOTAL", "PROFILES"]]; %First line is headers
    FltCntList= [FltCntList; [" ", " ", " ", " ", " "," "]]; %Adds a line of blank space (clunky, revise in the future)
    
    %Loop for writing program float counts
    for i = [2 6] %2 nad 6 pertain to SOCCOM and GO-BGC, respectively
    
        %Writes the line of float counts (Very long and clunky, should be a
        %more efficient way to write this)
        FltCntList = [FltCntList; [varsProg(i)+":" aProg(i)+"("+aProgActive(i)+")",nProg(i)+"("+nProgActive(i)+")", ...
            sProg(i)+"("+sProgActive(i)+")", tProg(i)+"("+tProgActive(i)+")",totalProfs(i)]];
        FltCntList= [FltCntList; [" ", " ", " ", " ", " "," "]];
    end
    
    %Manually add 'OTHER' category of floats
    FltCntList = [FltCntList; ["OTHER:" a1+"("+a2+")", a3+"("+a4+")", a5+"("+a6+")", a7+"("+a8+")",a9]];
    FltCntList= [FltCntList; [" ", " ", " ", " ", " "," "]];
    FltCntList = [FltCntList; ["--------------------------------------------------", "", "", "", ""," "]];
    
    %Create headers for Regions section
    FltCntList= [FltCntList; ["REGION", "APEX", "NAVIS", "SOLO", "TOTAL", ""]];
    FltCntList= [FltCntList; [" ", " ", " ", " ", " "," "]];
    
    %Loop for writing region float counts (again, a little clunky)
    for i = 1:length(varsReg)
        FltCntList = [FltCntList; [varsReg(i)+":" aReg(i)+"("+aRegActive(i)+")", nReg(i)+"("+nRegActive(i)+")", sReg(i)+"("+sRegActive(i)+")", tReg(i)+"("+tRegActive(i)+")"," "]];
        FltCntList= [FltCntList; [" ", " ", " ", " ", " "," "]];
    end
    
    %For each of the five columns, calculate the number of spaces needed to form a string of equal size, 
    %then create a string with required number of blanks
    C1_spacing = arrayfun(@blanks, 10 - cellfun('length', FltCntList(:,1)), UniformOutput=false); %Number of spaces after the category variables
    C2_spacing = arrayfun(@blanks, 10 - cellfun('length', FltCntList(:,2)), UniformOutput=false); %Apex counts
    C3_spacing = arrayfun(@blanks, 10 - cellfun('length', FltCntList(:,3)), UniformOutput=false); %Navis counts
    C4_spacing = arrayfun(@blanks, 10 - cellfun('length', FltCntList(:,4)), UniformOutput=false); %Solos counts
    C5_spacing = arrayfun(@blanks, 10 - cellfun('length', FltCntList(:,5)), UniformOutput=false); %Total counts
    C6_spacing = arrayfun(@blanks, 10 - cellfun('length', FltCntList(:,6)), UniformOutput=false); %Profile counts
    
    %Add the blanks to each column of FltCntList, creating uniformly spaced columns
    categories = FltCntList(:,1)+C1_spacing;
    apex = FltCntList(:,2)+C2_spacing;
    navis = FltCntList(:,3)+C3_spacing;
    solo = FltCntList(:,4)+C4_spacing;
    total = FltCntList(:,5)+C5_spacing;
    profiles = FltCntList(:,6)+C6_spacing;
    
    %Combine all of the spaced out columns into a single cell array
    allInfo = categories+apex+navis+solo+total+profiles;
    
    % %--------------------------------------------------------------------------------------
    % % Send the email (This is only) for local testing purposes
    % %--------------------------------------------------------------------------------------
    % 
    % sender = 'lgrady@mbari.org';
    % %email_list ={'tmaurer@mbari.org';'jplant@mbari.org';'johnson@mbari.org';'nicolag@mbari.org';'lgrady@mbari.org'};
    % email_list = {'lgrady@mbari.org'};%;'tmaurer@mbari.org'}
    % setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
    % setpref('Internet','E_mail',sender); % define sender
    % sendmail(email_list,"Logan's float count test", allInfo) %Send the email
end    