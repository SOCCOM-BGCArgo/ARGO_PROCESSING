%Load the corrected CDOM float list
CDOMdata = readtable('C:\Users\lgrady\Documents\MATLAB\ARGO_PROCESSING\Testing\mcoms_CDOM_corr_20240508');
CDOMdata_sorted = sortrows(CDOMdata,["WMOID"],"descend");
%keyboard
outPath = "C:\Users\lgrady\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATS\";

errFloats = {}

%35 starts from WMO 5906304
for i = 6:9%length(CDOMdata_sorted.WMOID)
    pullPath = ["\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATS\"+num2str(CDOMdata_sorted.WMOID(i))];
    pullDir = dir(pullPath);
    disp(['Currently on row: ',num2str(i), ' out of ',num2str(length(CDOMdata_sorted.WMOID))])
    %keyboard
    if isnan(CDOMdata_sorted.ScaleFactorCorrected(i))
        disp([num2str(CDOMdata_sorted.WMOID(i)),": no scale factor correction, data not corrected"])
    
    elseif (i~=1) && (CDOMdata_sorted.WMOID(i)==CDOMdata_sorted.WMOID(i-1))
        disp([num2str(CDOMdata_sorted.WMOID(i)),": duplicate WMO, data not corrected (again)"])

    elseif (i==1) || (CDOMdata_sorted.WMOID(i)~=CDOMdata_sorted.WMOID(i-1))
    
        disp(["Correcting CDOM for WMO:",num2str(CDOMdata_sorted.WMOID(i))])

        %Call the scale factor and corrected scale factor
        scaleFactor = CDOMdata_sorted.ScaleFactor(i);
        scaleFactorCorrected = CDOMdata_sorted.ScaleFactorCorrected(i);
    
        %TEMPORARY: make a directory to store data
        mkdir(outPath+num2str(CDOMdata_sorted.WMOID(i)))
       
        %Loop through the profiles in the directory (first 2 files are blank, so first profile starts at 3)
        for n = 3:length(pullDir)-1
            
            %Load the profile
            load ([pullDir(n).folder,'\',pullDir(n).name]);
            
            %NaN out all bad values
            LR.CDOM(LR.CDOM==99999) = NaN;
            HR.CDOM(HR.CDOM==99999) = NaN;
    
            %Apply the correction factor to the data
            LR.CDOM = (LR.CDOM/scaleFactor)*scaleFactorCorrected;
            HR.CDOM = (HR.CDOM/scaleFactor)*scaleFactorCorrected;
    
            %Reapply 99999 as the bad value
            LR.CDOM(isnan(LR.CDOM)) = 99999;
            HR.CDOM(isnan(HR.CDOM)) = 99999;
    
            %Save a test version of the corrected CDOM data
            save(outPath+num2str(CDOMdata_sorted.WMOID(i))+"\"+pullDir(n).name, "INFO","LR","HR")
        end
        try
            argo2odv_LIAR(num2str(CDOMdata_sorted.WMOID(i)),[],'all',0);
        catch ME %e is an MException struct
            
            fprintf(1,'The identifier was:\n%s',ME.identifier);
            fprintf(1,'There was an error! The message was:\n%s',ME.message);
            errFloats = [errFloats; num2str(CDOMdata_sorted.WMOID(i))];
            % more error handling...
        end
    end
end