% Float qc command center
% DESCRIPTION: running this code will create a seperate directory for this
% round of DMQC that synthesizes all of the relevant data from the DMQC
% operator assignment list
%
%
% INPUTS:
% - Location of the qc directory where the information will be stored
% after the user profile (e.g. \Documents\dmqc_may_2024) 
% ***RECOMMEND CREATING A SEPERATE FOLDER FOR THE OUTPUTS OF THIS CODE***
%
% - The location of the assigned floats table, which should be the correctly formatted
% .csv file with all of your assigned floats 
% (e.g. column names should be "WMO UWID PROGRAM REGION CRUISE", and
% then filled with the relevant information)
%
%
% OUTPUTS:
% - List of cruise directories within the specified qc directory
%
% - Within each cruise directory:
% -- List of WMO directories for that cruise
% -- A .lst file made of the WMOs from that cruise
%
% - Within each WMO directory:
% -- Copy of the floatviz file
% -- A qc history .txt file with all previous dmqc notes for the float
% -- An array of plots relevant for qc (still working on)
%---------------------------------------------------------------------------

userDir = getenv("USERPROFILE"); %Retreives the user directory

%--------------------------------- INPUTS-----------------------------------
tableDir = userDir+"\Documents\assignedFloatsforDMQC.csv"; %LOCATION OF THE TABLE


qcDir = userDir+"\Documents\dmqc_Floats"; %LOCATION OF ALL OUTPUTS IN "dmqc_Floats' directory
%---------------------------------------------------------------------------

%----------------------- START OF CODE -----------------------------------
mkdir(qcDir) %Make the qc directory where all of the information will be stored

floatVizDir = userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ"; %Location of floatviz files

%Load the dmqc float list
assignedFloats = readtable(tableDir);

%Load unique cruise ID's and the relevant directory paths for each
cruiseID = string(unique(assignedFloats.Cruise)); %Creates string array of each unique cruise
cruiseDirs = qcDir+"\"+cruiseID; %Creates filepath for each cruise within the qc directory

%Begin loop
disp("Collecting floats by cruise")
for i = 1:length(cruiseID) %For each unique Cruise ID:
    
    %If updating your qc directories, a warning will pop up that the
    %directory already exists, this code turns that off
    [a, MSGID] = lastwarn();
    warning('off', MSGID)

    %Create a folder for the cruise within the same directory as the table
    mkdir(cruiseDirs(i));
    
    %Create a list of WMO's from that specific cruise
    WMOID = assignedFloats.WMOID(matches(assignedFloats.Cruise,cruiseID(i))); 

    %Create a .lst file to append the floatviz filepaths to for ODV
    listFileID = fopen(cruiseDirs(i)+"\"+cruiseID(i)+"_list.lst",'w');

    %Loop through the WMO's
    for j = 1:length(WMOID) %For each WMO:
        %keyboard
        %Create a WMO directory for each float in the cruise
        wmoDir = cruiseDirs(i)+"\"+string(WMOID(j));
        mkdir(wmoDir);

        if isfile(floatVizDir+"\QC\"+WMOID(j)+"QC.TXT") %Check if a QC file exists

            %---------------------Cruise list file----------------------------
            %Call the location of the specific QC floatviz file
            listFileStr = floatVizDir+"\QC\"+WMOID(j)+"QC.TXT";

            %Copy it from chem to the table directory
            copyfile(listFileStr,wmoDir);
            
            %Write it as a new line to the listFile opened before
            fprintf(listFileID,'%s\n',listFileStr);
            %-----------------------------------------------------------------

          
            %---------------------QC HISTORY FILE----------------------------
            % Open the text file for reading
            floatQCListLog = fopen(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\FloatQCList_log.TXT", 'r');
            badSensorList = fopen(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT", 'r');
            badSampleList = fopen(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\bad_sample_list.TXT", 'r');

            %Create a new text file to append with qc notes for the current float
            qc_log = fopen(wmoDir+"\"+string(WMOID(j))+'_qcLog.TXT','w');
            fprintf(qc_log, '%s\n', 'QC log entries:');

            % Read lines until the end of the file
            while ~feof(floatQCListLog)
                % Read the next line from the file
                line = fgetl(floatQCListLog);
                
                % Check if the line contains specific characters
                if contains(line, string(WMOID(j)))
                    % If the line contains the specific characters, store it in the cell array
                    fprintf(qc_log, '%s\n', line);
                    %lines{end+1} = line;
                end
            end

            %Add new lines for bad sensors to the qcLog
            fprintf(qc_log, '\n%s\n', 'Bad sensors:');
            while ~feof(badSensorList)
                % Read the next line from the file
                line = fgetl(badSensorList);
                
                % Check if the line contains specific characters
                if contains(line, string(WMOID(j)))
                    % If the line contains the specific characters, store it in the cell array
                    fprintf(qc_log, '%s\n', line);
                    %lines{end+1} = line;
                end
            end

            %Add new lines for bad samples to the qcLog
            fprintf(qc_log, '\n%s\n', 'Bad samples:');
            while ~feof(badSampleList)
                % Read the next line from the file
                line = fgetl(badSampleList);
                
                % Check if the line contains specific characters
                if contains(line, string(WMOID(j)))
                    % If the line contains the specific characters, store it in the cell array
                    fprintf(qc_log, '%s\n', line);
                    %lines{end+1} = line;
                end
            end

            % Close the files
            fclose(floatQCListLog);
            fclose(badSensorList);
            fclose(badSampleList);
            fclose(qc_log);

            %Copy your local qc list file for each WMO and copy it to the proper cruise directory
            copyfile(userDir +"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\QC_LISTS\"+WMOID(j)+"_FloatQCList.TXT", cruiseDirs(i)+"\"+WMOID(j));
            %-----------------------------------------------------------------

        %Check if float doesn't exist in QC folder and repeat the process
        %with the non QC floatviz file
        elseif ~isfile(floatVizDir+"\"+WMOID(j)+"QC.TXT") && isfile(floatVizDir+"\"+WMOID(j)+".TXT") %If it doesn't exist in QC, check the regular folder
            disp("No QC file exists for float "+WMOID(j)+", using non-QC data")
            listFileStr = floatVizDir+"\"+WMOID(j)+".TXT";
            copyfile(floatVizDir+"\"+WMOID(j)+".TXT",wmoDir);
            fprintf(listFileID,'%s\n',listFileStr);
            
            %---------------------QC HISTORY FILE----------------------------
            % Open the text file for reading
            floatQCListLog = fopen(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\FloatQCList_log.TXT", 'r');
            badSensorList = fopen(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT", 'r');
            badSampleList = fopen(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\bad_sample_list.TXT", 'r');
            
            %Create a new text file to append with qc notes for the current float
            qc_log = fopen(wmoDir+"\"+string(WMOID(j))+'_qcLog.TXT','w');
            fprintf(qc_log, '%s\n', 'QC log entries:');
            % Read lines until the end of the file
            while ~feof(floatQCListLog)
                % Read the next line from the file
                line = fgetl(floatQCListLog);

                % Check if the line contains specific characters
                if contains(line, string(WMOID(j)))
                    % If the line contains the specific characters, store it in the cell array
                    fprintf(qc_log, '%s\n', line);
                    %lines{end+1} = line;
                end
            end
            
            %Add new lines for bad sensors to the qcLog PROBLEM AREA
            fprintf(qc_log, '\n%s\n', 'Bad sensors:');
            
            while ~feof(badSensorList)
                % Read the next line from the file
                line = fgetl(badSensorList);
               
                % Check if the line contains specific characters
                if contains(line, string(WMOID(j)))
                    % If the line contains the specific characters, store it in the cell array
                    fprintf(qc_log, '%s\n', line);
                    %lines{end+1} = line;
                end
            end
            
            %Add new lines for bad samples to the qcLog
            fprintf(qc_log, '\n%s\n', 'Bad samples:');
            while ~feof(badSampleList)
                % Read the next line from the file
                line = fgetl(badSampleList);

                % Check if the line contains specific characters
                if contains(line, string(WMOID(j)))
                    % If the line contains the specific characters, store it in the cell array
                    fprintf(qc_log, '%s\n', line);
                    %lines{end+1} = line;
                end
            end
            

            % Close the files
            fclose(floatQCListLog);
            fclose(badSensorList);
            fclose(badSampleList);
            fclose(qc_log);

            %-----------------------------------------------------------------

        else ~isfile(floatVizDir+"\"+WMOID(j)+"QC.TXT") && ~isfile(floatVizDir+"\"+WMOID(j)+".TXT"); %If it still doesn't exist in either, move on
            disp(WMOID(j)+" does not exist in regular or QC directories")
        end
    end
    %fclose(floatQCListLog);
    fclose(listFileID);
end
disp("All information stored in: " + qcDir)
clearvars
