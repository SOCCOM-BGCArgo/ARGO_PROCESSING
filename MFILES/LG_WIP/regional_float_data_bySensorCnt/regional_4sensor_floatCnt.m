function [floatCount] = regional_4sensor_floatCnt(filepath)
    
    %-------------------------------------------------------------------------------------------------------------------
    % INPUTS
    % SProf_filepath: a STRING of the filepath for the SProf index .txt file from either GDAC
    % EXAMPLE: filepath = "Documents\argo_synthetic-profile_index.txt";
    %
    % OUTPUTS
    % floatCount: a data structure with all of the regional and total float counts from ACTIVE floats off of the GDAC
    %
    % USAGE:
    % filepath = "Documents\argo_synthetic-profile_index.txt";
    % regional_4sensor_floatCnt(filepath)
    %-------------------------------------------------------------------------------------------------------------------

    %LOAD SPROF INDEX
    argoDat = readtable(getenv('USERPROFILE')+"\"+string(filepath));
    
    %Split the SPROF INDEX filenames by '/' and isolate the 7-digit WMO_ID
    fileSplit = split(argoDat.file,'/');
    
    %Create a seperate column in the data table for WMOs that uses the second column of the split string
    argoDat.WMO=str2double(fileSplit(:,2));
    
    %SUBSET FLOATS WITH 4+ SENSORS ONBOARD
    hasCTD = contains(argoDat.parameters,"PRES TEMP PSAL");
    hasDOXY = contains(argoDat.parameters,"DOXY");
    hasNIT = contains(argoDat.parameters,"NITRATE");
    hasPH = contains(argoDat.parameters,"PH_IN_SITU_TOTAL");
    hasCHLA = contains(argoDat.parameters,"CHLA")|contains(argoDat.parameters,"BBP700");
    hasIRR = contains(argoDat.parameters,"IRRADIANCE");
    argoDat.totalSensors = hasCTD+hasDOXY+hasNIT+hasPH+hasCHLA+hasIRR;
    
    %Subset all data where only 4+ sensors are present
    argoDat_4Sensor = argoDat(argoDat.totalSensors>=4,:);
    
    %Subset ACTIVE floats
    
    %Create a column for the date that the float was last updated
    argoDat_4Sensor.date_update = num2str(argoDat_4Sensor.date_update);
    
    %Format column to datetime that's readable
    argoDat_4Sensor.datetime_update = datetime(argoDat_4Sensor.date_update, 'InputFormat', 'yyyyMMddHHmmss');
    
    %Create a time threshold for any datetime between now and 6 months ago for NON POLAR floats
    active_time_threshold_nonPolar = datetime("now") - calmonths(6);
    activeFloats_notPolar = unique(argoDat_4Sensor.WMO(argoDat_4Sensor.datetime_update>=active_time_threshold_nonPolar ...
        & argoDat_4Sensor.latitude>=-60 & argoDat_4Sensor.latitude<=70));
    
    %Make the threshold a year for polar floats
    active_time_threshold_Polar = datetime("now") - calmonths(12);
    activeFloats_Polar = unique(argoDat_4Sensor.WMO(argoDat_4Sensor.datetime_update>=active_time_threshold_Polar ...
        & (argoDat_4Sensor.latitude<-60 | argoDat_4Sensor.latitude>70)));
    
    activeFloatsList = vertcat(activeFloats_Polar,activeFloats_notPolar);
    
    isActive = ismember(argoDat_4Sensor.WMO,activeFloatsList);
    argoDat_4Sensor = argoDat_4Sensor(isActive,:);
    
    %SUBSET DATA TO FIRST KNOWN PROFILE
    
    %Some extra steps have to be taken since a lot of floats don't have a "001.nc" file, but start
    %on "000.nc" or even "004.nc".
    
    %Create an empty array to store logical values and populate the first value as "true",
    %indicating that this is a "1st profile" for the float
    firstProf = NaN(length(argoDat_4Sensor.WMO),1);
    firstProf(1) = true;
    
    %Loop through each row
    for i = 2:length(argoDat_4Sensor.WMO)
    
        %Retrieve the current float and the previous float
        currentWMO = argoDat_4Sensor.WMO(i);
        previousWMO = argoDat_4Sensor.WMO(i-1);
    
        %Run a test to check if the current row is the same as the previous. If not, then 
        %the current row is the first profile of the current WMO
        firstProf(i) = currentWMO~=previousWMO;
    end
    
    %Subset all first profiles
    argoDat_1stProf = argoDat_4Sensor(logical(firstProf),:);
    
    %Check that length of first profile dataset matches the unique WMOs of entire dataset
    floatCount.allFloatsPresent = length(argoDat_1stProf.WMO)==length(unique(argoDat_4Sensor.WMO));
    
    %Take care of specific cases where no first gps is available (verified using later cycles)
    
    %1902457: 68.38E, 6.92N
    argoDat_1stProf(argoDat_1stProf.WMO==1902457,3).latitude = 6.92;
    argoDat_1stProf(argoDat_1stProf.WMO==1902457,4).longitude = 68.38;
    argoDat_1stProf(argoDat_1stProf.WMO==1902457,5).ocean = {'I'};
    %5905634: 170W, 71S
    argoDat_1stProf(argoDat_1stProf.WMO==5905634,3).latitude = -71;
    argoDat_1stProf(argoDat_1stProf.WMO==5905634,4).longitude = -170;
    argoDat_1stProf(argoDat_1stProf.WMO==5905634,5).ocean = {'P'};
    %2902882: 23N, 149E
    argoDat_1stProf(argoDat_1stProf.WMO==2902882,3).latitude = 23;
    argoDat_1stProf(argoDat_1stProf.WMO==2902882,4).longitude = 149;
    argoDat_1stProf(argoDat_1stProf.WMO==2902882,5).ocean = {'P'};
    
    %Check that there are no longer any empty coordinates
    floatCount.noMissingCoords = isempty(argoDat_1stProf(isnan(argoDat_1stProf.longitude),:));
    
    %CREATE REGIONAL INDICES
    %NPac: >10N, 133W, 90W
    isNPac = argoDat_1stProf.latitude>10 & argoDat_1stProf.latitude<=70 & contains(argoDat_1stProf.ocean,"P");
    %EqPac: >=10S, <=10N, 133W, 80W 
    isEqPac = argoDat_1stProf.latitude>=-10 & argoDat_1stProf.latitude<=10 & contains(argoDat_1stProf.ocean,"P");
    %SPac 70W, 
    isSPac = argoDat_1stProf.latitude<=-10 & argoDat_1stProf.latitude>=-30 & contains(argoDat_1stProf.ocean,"P");
    %NAtl
    isNAtl = argoDat_1stProf.latitude>0 & argoDat_1stProf.latitude<=70 & contains(argoDat_1stProf.ocean,"A");
    %SAtl
    isSAtl = argoDat_1stProf.latitude<=0 & argoDat_1stProf.latitude>=-30 & contains(argoDat_1stProf.ocean,"A");
    %IO
    isIO = argoDat_1stProf.latitude>=-30 & argoDat_1stProf.latitude<= 70 & contains(argoDat_1stProf.ocean,"I");
    %SO <30S
    isSO = argoDat_1stProf.latitude<-30;
    %ARC >70N
    isARC = argoDat_1stProf.latitude>70;
    
    floatCount.NPac = length(argoDat_1stProf.WMO(isNPac));
    floatCount.EqPac = length(argoDat_1stProf.WMO(isEqPac));
    floatCount.SPac = length(argoDat_1stProf.WMO(isSPac));
    floatCount.NAtl = length(argoDat_1stProf.WMO(isNAtl));
    floatCount.SAtl = length(argoDat_1stProf.WMO(isSAtl));
    floatCount.IO = length(argoDat_1stProf.WMO(isIO));
    floatCount.SO = length(argoDat_1stProf.WMO(isSO));
    floatCount.ARC = length(argoDat_1stProf.WMO(isARC));
    floatCount.total = floatCount.NPac+floatCount.EqPac+floatCount.SPac+ ...
        floatCount.NAtl+floatCount.SAtl+ ...
        floatCount.IO+ ...
        floatCount.SO+ ...
        floatCount.ARC;
end
