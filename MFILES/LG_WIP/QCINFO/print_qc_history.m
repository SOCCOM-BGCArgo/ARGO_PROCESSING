function print_qc_history(WMO)
% OBJECTIVE: Print out a list of all qc list log comments, bad sensor listings, and bad sample
%            listings for a specified float WMO. All results are pulled from CHEM directory.
%
% ====================== INPUTS ======================
%
% WMO:A numeric or string WMO of a float.

% ====================== OUTPUTS ======================
%
% A list of all qc list log comments, bad sensor listings, and bad sample
% listings for the input WMO.
%
% EXAMPLE:
% QC log entries:
% 4903273	wn1200	04/29/22 14:12:33	jplant	BIC = -10.279	Navis & no bottle data. Initial O2 gain vs WOA (but already 35 cycles. O2 died 28-
% 4903273	wn1200	04/29/22 14:58:36	jplant	Sage profile view doesn’t work to full depth for NO3 because of NaN’s. NO3 vs LINR @ 1500m cycles 1-27, NO3 vs LINR(no O2) 28-. pH vs LIPHR @ 1500m cycles 1-27, PH vs LIPHR(no O2) 28-.
% 4903273	wn1200	05/09/22 13:38:53	bgcargovm	Update pH, NO3 to LIR no O2 (for full record as cycle 1 is missing O2 data as well and uncertainty analysis == small impact to this float).
% 4903273	wn1200	12/18/22 19:23:08	jplant	NO3 sensor failed – dark counts shoot through the roof & Intensities erratic starting on cycle 32. Put  NO3 on BSL 32- . NO3 vs LINR-NOo2 @ 1500m. pH vs LIHPR-NOo2 @ 1500m. Looks to me like pH is fading 42-. I would put on BSL 1- QF =3 because of LIHPR-NOo2
% 4903273	wn1200	03/27/23 14:49:26	bgcargovm	Update pH adj to LIRnoO2.  Sensor dies at cycle 64 but starts drifting significantly at cycle 42.
% 4903273	wn1200	04/19/23 21:44:04	bgcargovm	testing new pH adj; some shallow profiles...
% 4903273	wn1200	05/29/24 16:29:46	bgcargovm	Super dynamic region, updated pH and nitrate to ESPERmix, not trusting references much here in the gulf stream, especially LIR noO2.  Data looks reasonable, prior to sensor failures.
% 
% Bad sensors:
% 4903273	wn1200	O	28-	4
% 4903273	wn1200	N	28-	4
% 4903273	wn1200	PH	42-	4
% 4903273	wn1200	BBP	77,	4
% 
% Bad samples:
% 4903273	wn1200	PH	1	1400-1600	4
% 4903273	un1200	PH	1	1400-1500	4
%
% ====================== USAGE ======================
%
% WMO = 4903273;
% print_qc_history(WMO)
%
% ===================================================================================================================================

    % Open the text file for reading
    floatQCListLog = fopen("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\FloatQCList_log.TXT", 'r');
    badSensorList = fopen("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\bad_sensor_list.TXT", 'r');
    badSampleList = fopen("\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\bad_sample_list.TXT", 'r');
    
    %Create a new text file to append with qc notes for the current float
    % userDir = getenv("USERPROFILE");
    % wmoDir = userDir+"\Documents";
    % qc_log = fopen(wmoDir+"\qc_history_log.TXT",'w');
    
    %Write first line of QC log
    %fprintf(qc_log, '%s\n', 'QC log entries:');
    fprintf('%s\n', 'QC log entries:');
    
    % Read lines until the end of the file
    while ~feof(floatQCListLog)
    
        % Read the next line from the file
        line = fgetl(floatQCListLog);
        
        % Check if the line contains the WMO
        if contains(line, string(WMO))
    
            % If the line contains the specific characters, store it in the cell array
            %fprintf(qc_log, '%s\n', line);
            fprintf('%s\n', line);
        end
    end
    
    %Add new lines for bad sensors to the qcLog PROBLEM AREA
    %fprintf(qc_log, '\n%s\n', 'Bad sensors:');
    fprintf('\n%s\n', 'Bad sensors:');
    
    while ~feof(badSensorList)
        % Read the next line from the file
        line = fgetl(badSensorList);
       
        % Check if the line contains specific characters
        if contains(line, string(WMO))
            % If the line contains the specific characters, store it in the cell array
            %fprintf(qc_log, '%s\n', line);
            fprintf('%s\n', line);
        end
    end
    
    %Add new lines for bad samples to the qcLog
    %fprintf(qc_log, '\n%s\n', 'Bad samples:');
    fprintf('\n%s\n', 'Bad samples:');
    while ~feof(badSampleList)
        % Read the next line from the file
        line = fgetl(badSampleList);
    
        % Check if the line contains specific characters
        if contains(line, string(WMO))
            % If the line contains the specific characters, store it in the cell array
            %fprintf(qc_log, '%s\n', line);
            fprintf('%s\n', line);
        end
    end
    
    
    % Close the files
    fclose(floatQCListLog);
    fclose(badSensorList);
    fclose(badSampleList);
    %fclose(qc_log);
end