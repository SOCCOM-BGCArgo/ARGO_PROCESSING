userDir = getenv("USERPROFILE");
floatsDir = dir(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS");
ph_log = fopen(userDir+"\Documents\realtimePH.TXT",'w');

for f = 4:length(floatsDir)
    wmoID = floatsDir(f).name;
    wmoDir = dir(userDir+"\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\"+wmoID);
    wmoTable = struct2table(wmoDir);
    firstFile = 3;
    lastFile = length(find(~contains(wmoTable.name,'cal')));

    wmoDir = wmoDir(firstFile:lastFile,:);

    anomMat = double.empty(length(wmoDir),0); %Initialize empty matrix
    anomCnt = 0;
    
    try
        disp(wmoID)

        for i=1:length(wmoDir)-1
            currentProf = load ([wmoDir(1).folder,'\',wmoDir(i+1).name]);
            currentPH = currentProf.LR.PH_IN_SITU_TOTAL(currentProf.LR.PRES>=1480 & ...
                currentProf.LR.PRES<=1520);
        
            prevProf = load ([wmoDir(1).folder,'\',wmoDir(i).name]);
            prevPH = prevProf.LR.PH_IN_SITU_TOTAL(prevProf.LR.PRES>=1480 & ...
                prevProf.LR.PRES<=1520);
        
            delta = abs(currentPH-prevPH);
          
            if delta >= .0075
                anomCnt = anomCnt + 1;
                anomMat(i,1) = anomCnt;
            else
                
                anomCnt = 0;
                anomMat(i,1) = anomCnt;
            end
        end
 
        estBadCycle = find(anomMat == 3)-2;

        try
            %estBadCycle(1)
            fprintf(ph_log, '%s\t%o\n', wmoID, estBadCycle(1));

        catch
            %msg = "No bad cycles detected";
            %fprintf(ph_log, '%s\t%s\n', wmoID, msg);
        end
    catch ME %e is an MException struct
  
        %fprintf(ph_log,'There was an error for %s:\t%s\n',wmoID, ME.message);
    end
    
end
fclose(ph_log);