function [A] = sprof2tbl(WMO)
    
    %READ IN THE CORRECT SPROF FILE FROM THE GDAC
    %=====================================================================================
    %Retrieve directory that the file will go into
    userDir = getenv('USERPROFILE');
    userDir = userDir+"\Documents\";
    
    %Identify the directory of the desired file on the GDAC
    webtarget = "https://data-argo.ifremer.fr";

    %Full file url of the synthetic profile index from the GDAC
    indexURL = webtarget+"/argo_synthetic-profile_index.txt";
    
    %Name of the file to be put into your local directory
    indexFile = userDir+"argo_synthetic-profile_index.txt";

    %Save the file path of the index file from the url
    indexPath = websave(indexFile, indexURL);
    
    %Read the index in as a table
    sprofIndex = readtable(indexPath);
    
    %Split the SPROF INDEX filenames column by '/' and isolate the 7-digit WMO_ID
    %This will allow parsing of the desired WMO
    fileSplit = split(sprofIndex.file,'/');
    
    %Create an array of WMO's
    indexWMO = str2double(fileSplit(:,2));
    
    %Create an array of the DAC's from the filename column of the index file
    indexDAC = fileSplit(:,1);
    
    %Retrieve an index of rows where the WMO matches the desired WMO and retrieve the DAC name
    DAC = indexDAC(ismember(indexWMO,WMO));
    DAC = string(DAC{1}); %Convert it to a string
    
    %Using the DAC and the WMO, identify the url to pull the correct SPROF.nc from the GDAC
    SProfURL  = webtarget+"/dac/"+DAC+"/"+WMO+"/"+WMO+"_Sprof.nc";
    
    %Name of the file to be put into your local directory
    sProfFile = userDir+"SProf.nc";
    
    %Save the file into userDir as sProfFile
    fullPath = websave(sProfFile,SProfURL);
    %=====================================================================================

    %LOAD IN THE NETCDF FILE AND RETRIEVE VARIABLE INFORMATION
    %=====================================================================================
    %Read in the nc file info
    ncDATA = ncinfo(fullPath);
    
    %Get a list of variable names
    ncNames = {ncDATA.Variables.Name}';
    ncSizes = {ncDATA.Variables.Size}';
    
    %Define the dimensions of the matrices
    dims = [ncDATA.Dimensions.Length];

    %Get indices for the number of profiles and samples per profile
    iProfs = contains({ncDATA.Dimensions.Name},"N_PROF");
    iLevels = contains({ncDATA.Dimensions.Name},"N_LEVELS");

    %Get the number of profiles and samples per profile from the file
    N_Profs = dims(iProfs);
    N_Levels = dims(iLevels);

    %Get the size of the data arrays for all float parameters (not the metadata)
    dataArraySize = [N_Levels N_Profs];

    %Get the variable types ('char', 'double', 'int32'...etc)
    dType = {ncDATA.Variables.Datatype};
    
    %Create a blank table with the names and the data types
    A = table('Size',[N_Levels*N_Profs,length(ncNames)], ...
        'VariableNames',ncNames', ...
        'VariableTypes',dType);
    %=====================================================================================

    %POPULATE THE TABLE
    %=====================================================================================
    for i = 1:length(ncNames)

        %For each of the variable names, load in some key meta information
        varName = ncNames{i};                 %Variable name
        Var = ncread(fullPath,ncNames{i});    %Variable data
        varSize = ncSizes{i};                 %Variable size (i.e. [3 x 2] array)
        varDims = length(ncSizes{i});         %Number of dimensions
        varType = dType{i};                   %Data type ('char', 'double', 'int32'...etc)
        
        %Start filling in the table columns

        %For any scientific comments or equations that are not needed
        if varDims>3
            continue; %Continue so that no other if statements are run
        end
    
        %For the METADATA that is a 2D character matrix
        if varDims==2 && varSize(1) < N_Levels
            %Variable is repeated over all cycles and transposed
            A.(varName) = repmat(Var,1,dataArraySize(1))';
            continue;
        end
        
        %For any 3D character matrices
        if varDims==3 && varSize(1) < N_Levels
            %Isolate the first matrix
            Var = string(Var(:,:,1)');

            %Pull out any individual words (no spaces)
            Var = regexp(Var,'\w+','match');

            %Morph all of the words into a "horizontal string"
            horizontalStr = sprintf('%s,', Var{:});

            %Replicate the horizontal string to fill out the entire column
            A.(varName) = string(repmat(horizontalStr,N_Levels*N_Profs,1));
            continue;
        end 
    
        %For Julian Date (JULD), but not the JULD_QC
        if contains(varName,'JULD') & ~contains(varName,'QC')

            %Create the reference date that JULD relies on
            dateRef = datenum(datetime(1950,1,1,0,0,0));

            %Convert it to a modern date
            Var = dateRef+Var;
            
            %Replicate the array for each cycle
            Var = repmat(Var,dataArraySize(1),1);

            %Squeeze it all to 1D to fill out the column
            A.(varName) = Var(:);
            continue;
        end
        
        %For the JULD and POSITION QC columns
        if contains(varName,'JULD_QC') | contains(varName,'POSITION_QC')
            %Scan the character arrays into integers, turning blanks into '7'
            QC = sscanf(regexprep(Var(:)',' ','7'),'%1f')';

            %Make all 7's NaN so that the column can maintain integer type
            QC(QC==7) = NaN;

            %Replicate array for each cycle and squeeze to 1D
            QC = repmat(QC,dataArraySize(1),1);
            A.(varName) = QC(:);
            continue;
        end
        
        %For all QC flag columns
        if contains(varName,'_QC') & ~contains(varName,'PROFILE_')
            QC = sscanf(regexprep(Var(:)',' ','7'),'%1f'); 
            QC(QC==7) = NaN;
            A.(varName) = QC(:);
            continue;
        end
    
        %For initial METADATA that is just a single character string
        if varDims==1 && varSize <= N_Profs & contains(varType,'char')

            %Single dimension character array gets replicated, transposed, and converted to a string
            A.(varName) = string(repmat(Var',dataArraySize(1)*dataArraySize(2),1));
            continue;
        end
        
        %For the CYCLE_NUMBER, LONGITUDE, and LATITUDE
        if  contains(varName,'CYCLE_NUMBER') | contains(varName,'LONGITUDE') | contains(varName,'LATITUDE')
            Var = repmat(Var',dataArraySize(1),1);
            A.(varName) = Var(:);
        end
    
        %For all regular variables with the correct array size
        if isequal(varSize, dataArraySize) & contains(varType,'single')
            %These ones are easy, just squeeze to 1D
            A.(varName) = Var(:);
        end
    end

    clearvars -except A
end


