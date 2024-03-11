function [Estimates,Uncertainties]= ...
    ESPER_NN(DesiredVariables,OutputCoordinates,PredictorMeasurements,PredictorTypes, ... % Required inputs
            varargin)                                     % Optional inputs
%  Version 1
%
%  Empirical Seawater Property Estimation Routines: Estimates seawater
%  properties and estimate uncertainty from combinations of other parameter
%  measurements.  ESPER_NN refers specifically to code that uses neural
%  networks as opposed to collections of interpolated linear regressions
%  (LIRs).  
%
%  Reserved for version update notes: (no updates, first version)
%
%
%
%
%
%
%  Documentation and citations:
%  LIARv1: Carter et al., 2016, doi: 10.1002/lom3.10087
%  LIARv2, LIPHR, LINR citation: Carter et al., 2018, https://doi.org/10.1002/lom3.10232
%  LIPR, LISIR, LIOR, first described/used: Carter et al., 2021, https://doi.org/10.1029/2020GB006623
%  LIRv3 and ESPER_NN, submitted.
%
%  ESPER_NN is inspired by CANYON-B, which also uses neural networks: 
%  Bittig et al. 2018: https://doi.org/10.3389/fmars.2018.00328
%
%  This function needs the CSIRO seawater package to run if measurements
%  are povided in molar units or if potential temperature or AOU are
%  needed but not provided by the user.  Scale differences from TEOS-10 are
%  a negligible component of alkalinity estimate error. 
% 
% *************************************************************************
% Input/Output dimensions:
% .........................................................................
% p: Integer number of desired property estimate types (e.g., TA, pH, NO3-)
% n: Integer number of desired estimate locations
% e: Integer number of equations used at each location
% y: Integer number of parameter measurement types provided by the user.
% n*e: Total number of estimates returned as an n by e array
% 
% *************************************************************************
% Required Inputs:
%
% DesiredVariables: (required 1 by p vector)
    % Specifies which variables will be returned, excepting unitless pH,
    % all outputs are in umol/kg.
    % 1. Total Titration Seawater Alkalinity (TA)
    % 2. Total Dissolved Inorganic Carbon (DIC)
    % 3. in situ pH on the total scale
    % 4. Phosphate
    % 5. Nitrate
    % 6. Silicate
    % 7. Dissolved Oxygen (O2)
%
% OutputCoordinates (required n by 3 array): 
    % Coordinates at which estimates are desired.  The first column should
    % be longitude (degrees E), the second column should be latitude
    % (degrees N), and the third column should be depth (m).
% 
% PredictorMeasurements (required n by y array): 
    % Parameter measurements that will be used to estimate alkalinity.  The
    % column order (y columns) is arbitrary, but specified by
    % PredictorTypes. Concentrations should be expressed as micromol per kg
    % seawater unless PerKgSwTF is set to false in which case they should
    % be expressed as micromol per L, temperature should be expressed as
    % degrees C, and salinity should be specified with the unitless
    % convention.  NaN inputs are acceptable, but will lead to NaN
    % estimates for any equations that depend on that parameter.
    % 
% PredictorTypes (required 1 by y vector): 
    % Vector indicating which parameter is placed in each column of the
    % 'PredictorMeasurements' input.  Note that salinity is required for
    % all equations.  If O2 is provided then temperature or potential
    % temperature must also be provided to convert O2 to AOU. For example,
    % if the first three columns of 'PredictorMeasurements' contain
    % salinity, silicate, and temperature, then PredictorTypes should equal
    % [1 5 7].
    % 
    % Input Parameter Key: 
    % 1. Salinity
    % 2. Temperature
    % 3. Phosphate
    % 4. Nitrate
    % 5. Silicate
    % 6. O2
%
% *************************************************************************
% Optional inputs:  All remaining inputs must be specified as sequential
% input argument pairs (e.g. "..., 'Equations',[1:16], 'OAAdjustTF', false,
% etc...")
%
% EstDates (Semi-optional n by 1 array or 1 by 1 value, default 2002.0): 
    % A vector of decimal dates for the estimates (e.g. July 1 2020 would
    % be "2020.5").  If only a single date is supplied that value is used
    % for all estimates.  It is highly recommended that date(s) be provided
    % for estimates of DIC and pH.  This version of the code will accept 1
    % by n inputs as well.
    %
% Equations (optional 1 by e vector, default []): 
    % Vector indicating which equations will be used to estimate
    % alkalinity. This input also determines the order of the columns in
    % the 'Estimates' output. If [] is input or the input is not
    % specified then all 16 equations will be used and only the outputs
    % from the equation with the lowest uncertainty estimate will be
    % returned.
    % 
    % (S=salinity, Theta=potential temperature, N=nitrate, Si=silicate,
    % T=temperature, O2=dissolved oxygen molecule... see
    % 'PredictorMeasurements' for units).
    % ...........................................................
    % Output Equation Key (See below for explanation of A, B, and C):
    % 1.  S, T, A, B, C
    % 2.  S, T, A, C
    % 3.  S, T, B, C
    % 4.  S, T, C
    % 5.  S, T, A, B
    % 6.  S, T, A
    % 7.  S, T, B
    % 8.  S, T
    % 9.  S, A, B, C
    % 10. S, A, C
    % 11. S, B, C
    % 12. S, C
    % 13. S, A, B
    % 14. S, A
    % 15. S, B
    % 16. S
    %
    % DesiredVar   | A             B             C
    % _____________|_____________________________________
    % TA           | Nitrate       Oxygen        Silicate
    % DIC          | Nitrate       Oxygen        Silicate
    % pH           | Nitrate       Oxygen        Silicate
    % phosphate    | Nitrate       Oxygen        Silicate
    % nitrate      | Phosphate     Oxygen        Silicate
    % silicate     | Phosphate     Oxygen        Nitrate
    % O2           | Phosphate     Nitrate       Silicate
    %
% MeasUncerts (Optional n by y array or 1 by y vector, default: [0.003 S,
    % 0.003 degrees C (T or theta), 1% P, 1% AOU or O2, 1% Si]: Array of
    % measurement uncertainties (see 'PredictorMeasurements' for units).
    % Uncertainties should be presented in order indicated by
    % PredictorTypes. Providing these estimates will improve ESPER_NN
    % estimate uncertainties in 'UncertaintyEstimates'. Measurement
    % uncertainties are a small part of ESPER_NN estimate uncertainties for
    % WOCE-quality measurements. However, estimate uncertainty scales with
    % measurement uncertainty, so it is recommended that measurement
    % uncertainties be specified for sensor measurements.  If this optional
    % input argument is not provided, the default WOCE-quality uncertainty
    % is assumed.  If a 1 by y array is provided then the uncertainty
    % estimates are assumed to apply uniformly to all input parameter
    % measurements.
    %
% pHCalcTF (Optional boolean, default false): 
    % If set to true, LIPHR will recalculate the pH to be a better estimate
    % of what the seawater pH value would be if calculated from TA and DIC
    % instead of measured with purified m-cresol dye. This is arguably also
    % a better estimate of the pH that would be obtained from pre-2011
    % measurements with impure dyes.  See the LIPHR paper for details
    %
% PerKgSwTF (Optional boolean, default true): 
    % Many sensors provide measurements in micromol per L (molarity)
    % instead of micromol per kg seawater. Indicate false if provided
    % measurements are expressed in molar units (concentrations must be
    % micromol per L if so).  Outputs will remain in molal units
    % regardless.
    %
% VerboseTF (Optional boolean, default true): 
    % Setting this to false will make ESPER_NN stop printing updates to the
    % command line.  This behavior can also be permanently disabled below.
    % Warnings and errors, if any, will be given regardless.
    %
% *************************************************************************
% Outputs:
% 
% Estimates: 
    % A n by e array of estimates specific to the coordinates and parameter
    % measurements provided as inputs.  Units are micromoles per kg
    % (equivalent to the deprecated microeq per kg seawater).
	%
% UncertaintyEstimates: 
    % A n by e array of uncertainty estimates specific to the coordinates,
    % parameter measurements, and parameter uncertainties provided.
    % 
% *************************************************************************
% Missing data: should be indicated with a NaN.  A NaN coordinate will
% yield NaN estimates for all equations at that coordinate.  A NaN
% parameter value will yield NaN estimates for all equations that require
% that parameter.
% 
% *************************************************************************
% Please send questions or related requests to brendan.carter@gmail.com.
% *************************************************************************

% Determining whether the user requested command-line update text.
a=strcmpi(varargin,'VerboseTF');
if any(a)
    VerboseTF=varargin{1,logical([0 a(1:end-1)])};
else
    VerboseTF=true;
end
% Uncomment following line beginning with "VerboseTF" and save the function
% if you want less command line spam and you don't want to have to keep
% telling the code to be quiet.

% VerboseTF=false;

% *************************************************************************
% Parsing inputs, setting defaults, and sanity checking inputs.
%
% Starting timer
tic

% Verifying required inputs are provided
if nargin<4; error('ESPER_NN called with too few input arguments.'); end

% Checking whether specific equations are specified.
a=strcmpi(varargin,'Equations');
if any(a)
    Equations=varargin{1,logical([0 a(1:end-1)])};
else
    Equations=1:16;
end
% Making [] argument for Equations equivalent to no argument.
if isempty(Equations);Equations=1:16; end
% Making 0 argument for Equations equivalent to no argument.
if Equations==0; Equations=1:16; end

% Checking for PerKgSwTF input and setting default if not given
a=strcmpi(varargin,'PerKgSwTF');
if any(a)
    PerKgSwTF=varargin{1,logical([0 a(1:end-1)])};
else
    PerKgSwTF=true;
end

% Checking for MeasUncerts input and setting default if not given
a=strcmpi(varargin,'MeasUncerts');
if any(a)
    MeasUncerts=varargin{1,logical([0 a(1:end-1)])};
    UseDefaultUncertainties=false;
    % Sanity checking the MeasUncerts argument.  This also deals with the
    % possibility that the user has provided a single set of uncertainties for
    % all estimates.
    if not(max(size(MeasUncerts)==size(PredictorMeasurements))) && not(min(size(MeasUncerts)==size(PredictorMeasurements))) && not(max(size(MeasUncerts))==0)
        error('MeasUncerts must be undefined, a vector with the same number of elements as PredictorMeasurements has columns, [] (for default values), or an array of identical dimension to PredictorMeasurements.')
    elseif not(min(size(MeasUncerts)==size(PredictorMeasurements))) && not(max(size(MeasUncerts))==0)
        if ~(size(MeasUncerts,2)==size(PredictorMeasurements,2))
            error('There are different numbers of columns of input uncertainties and input measurements.')
        end
        MeasUncerts=ones(size(PredictorMeasurements(:,1)))*MeasUncerts;              % Copying uncertainty estimates for all estimates if only singular values are provided
    end
    if ~(size(PredictorTypes,2)==size(PredictorMeasurements,2))            % Making sure all provided predictors are identified.
        error('The PredictorTypes input does not have the same number of columns as the PredictorMeasurements input.  This means it is unclear which measurement is in which column.');
    end
else
    UseDefaultUncertainties=true;
end

% Checking for pHCalcTF input and setting default if not given
a=strcmpi(varargin,'pHCalcTF');
if any(a)
    pHCalcTF=varargin{1,logical([0 a(1:end-1)])};
else
    pHCalcTF=false;
end

% Checking for EstDates input and setting default if not given
a=strcmpi(varargin,'EstDates');
if any(a)
    EstDates=varargin{1,logical([0 a(1:end-1)])};
    UserProvidedDates=true;
else
    EstDates=2002.0;                                                       % Setting default date estimate.  This is unneccessary due to the way this is coded below, but might help with debugging.
    UserProvidedDates=false;
end
if size(EstDates,1)==1 
    EstDates=EstDates*ones(size(OutputCoordinates,1),1);                   % Reshaping EstDates to have the same length as Coordinates if only a single date is provided.
end

% ESPER_NN requires non-NaN coordinates to provide an estimate.  This step
% eliminates NaN coordinate combinations prior to estimation.  NaN
% estimates will be returned for these coordinates.
NaNGridCoords=max(isnan(OutputCoordinates),[],2);

% Doing a size check for the coordinates.
if ~(size(OutputCoordinates,2)==3)
    error('OutputCoordinates has either too many or two few columns.  This version only allows 3 columns with the first being longitude (deg E), the second being latitude (deg N), and the third being depth (m).');
end

% Figuring out how many estimates are required
n=sum(~NaNGridCoords,1);
e=size(Equations,2);
p=size(DesiredVariables,2);

% Checking for common missing data indicator flags and warning if any are
% found.  Consider adding your commonly used flags here.
if max(ismember([-999 -9 -1*10^20],PredictorMeasurements))==1
    warning('ESPER_NN: A common non-NaN missing data indicator (e.g. -999, -9, -1e20) was detected in the input measurements provided.  Missing data should be replaced with NaNs.  Otherwise, ESPER_NN will interpret your inputs at face value and give terrible estimates.')
end

% Book-keeping with coordinate inputs and adjusting negative longitudes.
C=OutputCoordinates(~NaNGridCoords,:);
C(C(:,1)>360,1)=rem(C(C(:,1)>360,1),360);
C(C(:,1)<0,1)=rem(C(C(:,1)<0,1,1),360)+360;

% Throwing an error if latitudes are out of the expected range.
if max(abs(C(:,2)))>90
    error('ESPER_NN: A latitude >90 degrees (N or S) has been detected.  Verify latitude is in the 2nd column of the coordinate input.');
end

% Preparing full PredictorMeasurement uncertainty grid
DefaultUncertainties=diag([1 1 0.02 0.02 0.02 0.01]);
DefaultUAll=zeros(size(PredictorMeasurements,1),6);
DefaultUAll(:,PredictorTypes)=PredictorMeasurements* ...
    DefaultUncertainties(PredictorTypes,PredictorTypes);                   % Setting multiplicative default uncertainties for P, N, O2, and Si.
DefaultUAll(:,ismember(PredictorTypes,[1 2]))=0.003;                       % Then setting additive default uncertainties for T and S.
DefaultUAll=DefaultUAll(~NaNGridCoords,:);
if UseDefaultUncertainties==false
    InputUAll=zeros(size(PredictorMeasurements,1),6);                      % Initializing measurement uncertainties
    InputUAll(:,PredictorTypes)=MeasUncerts;                               % Rearranging user inputs to be in expected order per input parameter key (above)
    InputUAll=max(cat(3,InputUAll, DefaultUAll),[],3);                     % Overriding user provided uncertainties that are smaller than the (minimum) default uncertainties (see paper for rationale)
else
    InputUAll=DefaultUAll;
end    

% Making sure all provided predictors are identified.
if ~(size(PredictorTypes,2)==size(PredictorMeasurements,2))
    error('The PredictorTypes input does not have the same number of columns as the Measurements input.  This means it is unclear which measurement is in which column.');
end

if size(EstDates,2)>size(EstDates,1); EstDates=EstDates'; end              % Sometimes people give EstDates as a 1 by n vector.

% Putting all provided measurement inputs in standard order
MAll=NaN(n,6);                                                             % PredictorMeasurements
UAll=NaN(n,6);                                                             % Uncertainties
MAll(:,PredictorTypes)=PredictorMeasurements(~NaNGridCoords,:);            % Reordering and limiting to viable coordinates
UAll(:,PredictorTypes)=InputUAll(:,PredictorTypes);                        % This was already limited to viable coordinates for later use.
EstDates=EstDates(~NaNGridCoords,1);                                       % Limiting Dates to viable coordinates as well
YouHaveBeenWarnedCanth=false;                                              % Calculating Canth is slow, so this flag is used to make sure it only happens once.

% Light input range checking
if (any(MAll(:,1)<5) | any(MAll(:,1)>50)) && VerboseTF==true; disp('Warning: Salinities greater than 50 or less than 5 have been found. ESPER is not intended for seawater with these properties.'); end
if (any(MAll(:,2)<-5) | any(MAll(:,2)>50)) && VerboseTF==true; disp('Warning: Temperatures greater than 50C or less than -5C have been found. ESPER is not intended for seawater with these properties. Note that ESPER expects temperatures in Centigrade.'); end

% Beginning the iterations through the requested properties
for PIter=1:p                                                              
    Property=DesiredVariables(1,PIter);
    % Specifying which variables are required for each property estimate.
    % Each row is for a different property. 
    NeededForProperty=[1 2 4 6 5
        1 2 4 6 5
        1 2 4 6 5
        1 2 4 6 5
        1 2 3 6 5
        1 2 3 6 4
        1 2 3 4 5];
    % Which of the 5 predictors are required for the equations specified?
    % Each row is for a different equation.
    VarVec=logical([1 1 1 1 1
        1 1 1 0 1    
        1 1 0 1 1    
        1 1 0 0 1    
        1 1 1 1 0    
        1 1 1 0 0   
        1 1 0 1 0   
        1 1 0 0 0   
        1 0 1 1 1    
        1 0 1 0 1  
        1 0 0 1 1
        1 0 0 0 1
        1 0 1 1 0
        1 0 1 0 0
        1 0 0 1 0
        1 0 0 0 0]);
    NeedVars=any(VarVec(Equations,:),1);
    HaveVars=false(1,6);  HaveVars(1,PredictorTypes)=1;
    
    % Temperature is required if O2 is used as a predictor (to convert to
    % AOU)... equation 15 is therefore only valid under some circumstances.
    if ismember(6,PredictorTypes)==1; NeedVars(1,2)=true; end
    
    % Ignoring variables that aren't used for the current property
    HaveVarsSub=HaveVars(1,NeededForProperty(Property,:));

    % Temperature is required if measurements are provided in molar units
    % or if CO2sys will be used to adjust pH values
    if PerKgSwTF==false; NeedVars(1,2)=1; end 
    if (UserProvidedDates==true) && (Property==3) && (HaveVars(1,2)==0)
        warning('Carbonate system calculations will be used to adjust the pH, but no temperature is specified.  10 C will be assumed.  If this is a poor estimate for your region, consider supplying an approximate temperature value in the PredictorMeasurements input.');
        Temperature=10*ones(n,1);
    elseif (UserProvidedDates==true) && (Property==3) && (HaveVars(1,2)==1)
        Temperature=MAll(:,2);        
    end

    % Making sure all needed variables are present
    if ~all(HaveVarsSub(1,NeedVars)) && VerboseTF==true  % We are assuming only power users turn VerboseTF off, and hopefully they understand the function call
        disp('Warning: One or more regression equations for the current property require one or more input parameters that are either not provided or not labeled correctly. These equations will return NaN for all estimates.  All 16 equations are used by default unless specific equations are specified.  Temperature is also required when density or carbonate system calculations are called.'); % They should already be NaN
    end
    
    % Limiting measurements to the subset needed for this property estimate
    M=MAll(:,NeededForProperty(Property,:));

    % Checking to see whether temperature is needed, and subbing in
    % potential temp if it is.  As of the version 3 update the  routines
    % now request temperature, but convert it to potential temperature
    % before using it for predictions.
    if NeedVars(1,2)
        M(:,2)=sw_ptmp(M(:,1),M(:,2),sw_pres(C(:,3),C(:,2)),0);
    end
    % Checking to see whether O2 is needed. Defining AOU and subbing in for
    % O2 if yes (see above).
    if any(ismember(NeededForProperty(Property,:),6))
        M(:,4)=sw_satO2(M(:,1),M(:,2))*44.6596-M(:,4);
    end
    % Converting units to molality if they are provided as molarity.
    if PerKgSwTF==false
        densities=sw_dens(M(:,1),M(:,2),sw_pres(C(:,3),C(:,2)))/1000;
        M(:,3)=M(:,3)./densities;
        M(:,4)=M(:,4)./densities;
        M(:,5)=M(:,5)./densities;
    end
    
    % *********************************************************************
    % Beginning treatment of inputs and calculations
    if     Property==1; VName='TA';
    elseif Property==2; VName='DIC';
    elseif Property==3; VName='pH';
    elseif Property==4; VName='phosphate';
    elseif Property==5; VName='nitrate';
    elseif Property==6; VName='silicate';
    elseif Property==7; VName='oxygen';
    else; error('A property identifier >7 or <1 was supplied, but this routine only has 7 possible property estimates.  The property identifier is the first input.')
    end

    % Loading the data, with an error message if not found
    FN=horzcat('NN_files_',VName,'_Unc_Poly.mat');
    % Making sure you downloaded the needed file and put it somewhere it
    % can be found
    if exist(FN,'file')<2; error('ESPER_NN could not find the file(s) needed to run.  These mandatory file(s) should be distributed from the same website as ESPER_NN.  Contact the corresponding author if you cannot find it there.  If you do have it then make sure all of the contents of the ESPER.zip extract are on the MATLAB path or in the active directory.  This will require adding several subfolders for ESPER.'); end
    L=load(FN,'UncGrid','Polys');
    
    %Preallocating for speed
    OutputEstimates=NaN(size(OutputCoordinates,1),e);                      % Using size instead of n since we want to preallocate for NaN coordinate locations as well                 
    UncertEst=NaN(size(OutputCoordinates,1),e); 
    EMLR=NaN(n,e);
    EstAtl=NaN(sum(~NaNGridCoords),e,4);
    EstOther=NaN(sum(~NaNGridCoords),e,4);

    % This loop performs the estimation
    % Disambiguation:
    % Eq... a counter for which equation LIRv3 is on
    % e... the number of equations that will be used
    % Equation... the specific equation number currently being used
    % Equations... the user provided list of equations that will be used.
    for Eq=1:e                                                             % Iterating over the (up to 16) equations
        Equation=Equations(1,Eq);
        NeedVarsNow=any(VarVec(Equation,:),1);
        if all(HaveVarsSub(1,NeedVarsNow))                                 % Skipping if we don't have what we need
            P=horzcat(cosd(C(:,1)-20),sind(C(:,1)-20),C(:,2:3),M(:,VarVec(Equations(Eq),:)));
            for Net=1:4 % A committee of 4 neural networks is used.
                % Separate neural networks are used for the Arctic/Atlantic and
                % the rest of the ocean (more info below).  This functional
                % form of ESPER_NN is designed to avoid calling the neural
                % network toolbox.  It will require the (many) ESPER_NN
                % functions to be in the 'private' subfolder at the ESPER_NN.m location.
                FunctionOther=str2func(horzcat('ESPER_',VName,'_',num2str(Equation),'_Other_',num2str(Net))); % Basically an eval statement.  Sorry Matlab power users.
                FunctionAtl=str2func(horzcat('ESPER_',VName,'_',num2str(Equation),'_Atl_',num2str(Net)));
                EstAtl(:,Eq,Net)=FunctionAtl(P')';
                EstOther(:,Eq,Net)=FunctionOther(P')';
            end
        end
    end
    % Averaging across neural network committee members
    EstAtl=nanmean(EstAtl,3);
    EstOther=nanmean(EstOther,3);

    % We do not want information to propagate across the Panama Canal (for
    % instance), so data is carved into two segments... the Atlantic/Arctic
    % (AA) and everything else.
    AAInds=or(inpolygon(C(:,1),C(:,2),L.Polys.LNAPoly(:,1),L.Polys.LNAPoly(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LSAPoly(:,1),L.Polys.LSAPoly(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LNAPolyExtra(:,1),L.Polys.LNAPolyExtra(:,2)), ...
        or(inpolygon(C(:,1),C(:,2),L.Polys.LSAPolyExtra(:,1),L.Polys.LSAPolyExtra(:,2)), ...
        inpolygon(C(:,1),C(:,2),L.Polys.LNOPoly(:,1),L.Polys.LNOPoly(:,2))))));
    % We'd like a smooth transition in the Bering Strait and in the South
    % Atlantic.  This linearly interpolates between the networks at their
    % boundaries.
    BeringInds=inpolygon(C(:,1),C(:,2),L.Polys.Bering(:,1),L.Polys.Bering(:,2));
    SoAtlInds=(C(:,1)>290 | C(:,1)<20) & (C(:,2)>-44) & (C(:,2)<-34);
    SoAfrInds=(C(:,1)>19 & C(:,1)<27) & (C(:,2)>-44) & (C(:,2)<-34);
    Est=EstOther;                                                          % Pulling out Indo-Pacific estimates
    Est(AAInds,:)=EstAtl(AAInds,:,:);                                      % Pulling out all other estimates
    Est(BeringInds,:)=EstAtl(BeringInds,:).* ...                           % Blending estimates in the Bering Strait vicinity
        repmat((C(BeringInds,2)-62.5)/(7.5),[1,size(Equations,2)]) ...
        +EstOther(BeringInds,:).* ...
        repmat((70-C(BeringInds,2))/(7.5),[1,size(Equations,2)]);
    Est(SoAtlInds,:)=EstAtl(SoAtlInds,:).* ...                             % Blending estimates in the Southern Atlantic to Southern Ocean transition
        repmat((C(SoAtlInds,2)+44)/(10),[1,size(Equations,2)])+ ...
        EstOther(SoAtlInds,:).* ...
        repmat((-34-C(SoAtlInds,2))/(10),[1,size(Equations,2)]);
    Est(SoAfrInds,:)=Est(SoAfrInds,:).* ...                                % Blending estimates South of South Africa toward Indo-Pac Estimates
        repmat((27-C(SoAfrInds,1))/(8),[1,size(Equations,2)])+ ...
        EstOther(SoAfrInds,:).* ...
        repmat((C(SoAfrInds,1)-19)/(8),[1,size(Equations,2)]);
    
    if nargout>1                                                           % Skipping this step if the user does not request uncertainty estimates as outputs
        % Estimating methodological error from depth and salinity
        UGridArray=cat(2,L.UncGrid.DepthGrid(:),L.UncGrid.SalGrid(:),L.UncGrid.EqnGrid(:),L.UncGrid.RMSEs(:));
        UGridArray(any(isnan(UGridArray),2),:)=[];                         % Getting rid of some NaN Values
        UInterpolant=scatteredInterpolant(UGridArray(:,1),UGridArray(:,2),UGridArray(:,3),UGridArray(:,4));
        % LIRv3 can simplistically extract error sensitivities from the
        % coefficients, but ESPER_NN has to call itself in perturbation
        % tests to derive these uncertainties.  These sensitivity tests are
        % conducted without the Canth adjustment for speed.
        for Eq=1:e
            EMLR(:,Eq)=UInterpolant(C(:,3),M(:,1),Equations(1,Eq)*ones(size(M(:,1))));
        end
        %  Scaling estimate uncertainty from LCs
        PertDiff=NaN(n,e,size(PredictorTypes,2));                          % Preallocating
        DefaultPertDiff=NaN(n,e,size(PredictorTypes,2));                   % Preallocating
        for Pred=1:size(PredictorTypes,2)
            Pert=zeros(size(PredictorMeasurements(~NaNGridCoords,:)));
            DefaultPert=zeros(size(PredictorMeasurements(~NaNGridCoords,:)));
            Pert(:,Pred)=InputUAll(:,PredictorTypes(1,Pred));
            DefaultPert(:,Pred)=DefaultUAll(:,PredictorTypes(1,Pred));
            [PertEst]=ESPER_NN(Property,OutputCoordinates(~NaNGridCoords,:),PredictorMeasurements(~NaNGridCoords,:)+Pert,PredictorTypes,'Equations',Equations,'VerboseTF',false);  % It is okay to omit EstDates because, by this point in the code, the benchmark values are specific to 2002.
            [DefaultPertEst]=ESPER_NN(Property,OutputCoordinates(~NaNGridCoords,:),PredictorMeasurements(~NaNGridCoords,:)+DefaultPert,PredictorTypes,'Equations',Equations,'VerboseTF',false);
            PertDiff(:,:,Pred)=Est-PertEst.(VName);
            DefaultPertDiff(:,:,Pred)=Est-DefaultPertEst.(VName);
        end
        Unc=real(...
            (sum(PertDiff.^2,3)-sum(DefaultPertDiff.^2,3)+...              % Adding the user-provided uncertainties and subtracing the smaller default uncertainties
            EMLR.^2).^(1/2));                                              % Adding baseline methodological uncertainties
        UncertEst(~NaNGridCoords,:)=Unc;                                   % Inserting new estimates into the full array
        Uncertainties.(VName)=UncertEst;                                   % Inserting estimates into output structure
    end
    
    % DIC and pH estimates are appropriate for 2002 by default.  Here we
    % scale these estimates to the user provided dates.
    if UserProvidedDates==true && (Property==2 || Property==3)
        if YouHaveBeenWarnedCanth==false
            if VerboseTF==true
                disp('Estimating anthropogenic carbon for ESPER_NN.')
            end
            [Cant,Cant2002]=SimpleCantEstimateLR(C(:,1),C(:,2),C(:,3),EstDates); % Technically this function estimates anthropogenic carbon, but it does a sloppy job of it.  It is okay here only because we need a fast, repeatable, and approximate estimate. Please contact Carter about alternatives before using this function for any real science.
            YouHaveBeenWarnedCanth=true;
        end
        if Property==2
            % Adjusting DIC
            Est=Est+Cant*ones(1,e)-Cant2002*ones(1,e);                     
        elseif Property==3
            % Adjusting pH
            % Getting approximate inputs for CO2SYS
            EstAlk=ESPER_NN(1,C,M(:,1),1,'Equations',16,'VerboseTF',false);                  % Approximating alkalinity for the calculation of impacts of anthropogenic CO2 on pH.  Using equation 16 because it is the only one that is viable regardless of inputs.
            EstSi=zeros(size(EstAlk)); EstP=zeros(size(EstAlk));           % Not very sensitive to these so using a poor estimate
            Pressure=sw_pres(C(:,3),C(:,2));
            for Eq=1:e
                Out=CO2SYS_ESPER(EstAlk.TA,Est(:,Eq),1,3,M(:,1),Temperature,Temperature,Pressure,Pressure,EstSi,EstP,1,10,3);                    % Calculating DIC for estimated TA and pH
                OutAdj=CO2SYS_ESPER(EstAlk.TA,Out(:,2)+Cant-Cant2002,1,2,M(:,1),Temperature,Temperature,Pressure,Pressure,EstSi,EstP,1,10,3);   % Recalculating pH after increasing DIC for additional Canth
                Est(:,Eq)=real(OutAdj(:,18));
                if any(isnan(OutAdj(:,18)))
                    disp('Warning: CO2SYS took >20 iterations to converge. The corresponding estimate(s) will be NaN. This typically happens when ESPER_NN is very poorly suited for estimating water with the given properties (e.g., very high or low salinity or estimates in marginal seas).');
                end
            end
        end
    elseif (Property==2 || Property==3) && UserProvidedDates==false && YouHaveBeenWarnedCanth==0 && VerboseTF==true
        disp('Warning: DIC or pH is a requested output but the user did not provide dates for the desired estimates.  The estimates will be specific to 2002.0 unless the optional EstDates input is provided (recommended).')
        YouHaveBeenWarnedCanth=true;
    end
    
    % Adjusting the pH to be appropriate for values calculated from TA and
    % DIC, if desired.  Otherwise estimates will be appropriate for
    % spectrophotometric pH.
    if pHCalcTF==true && Property==3
        if VerboseTF==true
            disp('Recalculating the pH to be appropriate for pH values calculated from TA and DIC.')
        end
        % Solving equation 1 in the LIRv2 paper to counter 3->4 adjustment.
        Est=(Est+0.3168)/1.0404;
    end

    % Filling the outputs with the estimates at viable locations.
    OutputEstimates(~NaNGridCoords,:)=Est;
    Estimates.(VName)=OutputEstimates;
    % Updating, off by default but can be turned on here.
%     if VerboseTF==true 
%         disp(horzcat('ESPER_NN finished ',VName)); 
%     end
end