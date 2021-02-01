function [pHEstimates,UncertaintyEstimates,MinUncertaintyEquation, ...
    EstimateCsDate]= ...
    LIPHR(Coordinates,Measurements,MeasIDVec, ...         % Required inputs
            varargin)                                     % Optional inputs
% Version 2.0.2 
% *************************************************************************
% Version 2.0 version at manuscript submission
% Version 2.0.1 (2017.10.17) version at manuscript acceptance
%       -"Molality" changed to "PerKgSw," 
%       -retrained with updated cruise flags
%       -retrained with linear adjustments, and changed pHCalcTF to match
%       -enabled user-defined OARates
% Version 2.0.2 (2018.10.10) Fixed a bug in the training routines.
%        New errors should be closer to estimated errors in the 2nd
%        citation. Older errors were averaged 0.002 higher than estimated.
%        See Readme for details.
% *************************************************************************
%  Locally Interpolated pH Regression (LIPHR): Estimates pH
%  and pH estimate uncertainty from combinations of other parameter
%  measurements.
%
%  Citations:
%  LIARv1: Carter et al., 2016, doi: 10.1002/lom3.10087
%  LIARv2, LIPHR, LINR citation: Carter et al. 2018, 
%         https://doi.org/10.1002/lom3.10232
%
%  This function needs the CSIRO seawater package to run if measurements
%  are povided in molar units or if potential temperature or AOU are
%  needed but not provided by the user.  Scale differences from TEOS-10 are
%  a negligible component of pH estimate error.  The seawater
%  package is distributed freely by CSIRO (website has citation info):
% 
%  http://www.cmar.csiro.au/datacentre/ext_docs/seawater.htm
% 
% *************************************************************************
% Input/Output dimensions:
% n: Integer number of desired estimate locations
% e: Integer number of equations used at each location
% y: Integer number of parameter measurement types provided by the user.
% n*e: Total number of estimates returned as an n by e array
% 
% *************************************************************************
% Required Inputs:
%
% Coordinates (required n by 3 array): 
    % Coordinates at which estimates are desired.  The first column should
    % be longitude (degrees E), the second column should be latitude
    % (degrees N), and the third column should be depth (m).
% 
% Measurements (required n by y array): 
    % Parameter measurements that will be used to estimate pH.  The column
    % order (y columns) is arbitrary, but specified by MeasIDVec.
    % Concentrations (including AOU) should be expressed as micromol per kg
    % seawater unless PerKgSwTF is set to false in which case they should
    % be expressed as micromol per L, temperature should be expressed as
    % degrees C, and salinity should be specified with the unitless
    % convention.  NaN inputs are acceptable, but will lead to NaN
    % estimates for any equations that depend on that parameter.
    % 
% MeasIDVec (required 1 by y vector): 
    % Vector indicating which parameter is placed in each column of the
    % 'Measurements' input.  Note that salinity is required for all
    % equations.  If O2 is provided instead of AOU, then temperature or
    % potential temperature must also be provided to convert O2 to AOU.
    % For example, if the first three columns of 'Measurements' contain
    % salinity, silicate, and temperature, then MeasIDVec should equal [1 5
    % 7].
    % 
    % Input Parameter Key: 
    % 1. Salinity
    % 2. Potential temperature 
    % 3. Nitrate
    % 4. AOU
    % 5. Silicate
    % 6. O2 (recommended over AOU if temperature is included so the user
    %        can be confident AOU will be calculated the same way it was
    %        calculated for the training dataset)
    % 7. Temperature (recommended over theta for same reason as O2 vs. AOU)
%
% LIPHR only! Two "optional" inputs are actually required for LIPHR.  You
% must either provide estimate date information with 'EstDates' (and then
% provide the dates), or provide ...,'OAAdjustTF',false,... to tell the
% routine to skip this adjustment.
% *************************************************************************
% Optional inputs:  All remaining inputs must be specified as sequential
% input argument pairs with the argument specifying the input type in
% single quotes (e.g. "...,'Equations',[1:16],'OAAdjustTF',false, etc...")
%
% Equations (optional 1 by e vector, default []): 
    % Vector indicating which equations will be used to estimate pH. This
    % input also determines the order of the columns in the 'pHEstimates'
    % output. If [] is input or the input is not specified then all 16
    % equations will be used and only the outputs from the equation with
    % the lowest uncertainty estimate will be returned.
    % 
    % (S=salinity, Theta=potential temperature, N=nitrate, Si=silicate,
    % AOU=apparent oxygen utilization... see 'Measurements' for units).
    % 
    % Output Equation Key: 
    % 1.  S, Theta or T, N, AOU or O2, Si
    % 2.  S, Theta or T, N, Si
    % 3.  S, Theta or T, AOU, Si
    % 4.  S, Theta or T, Si
    % 5.  S, Theta or T, N, AOU or O2
    % 6.  S, Theta or T, N
    % 7.  S, Theta or T, AOU or O2
    % 8.  S, Theta
    % 9.  S, N, AOU or O2, Si
    % 10. S, N, Si
    % 11. S, AOU or O2, Si
    % 12. S, Si
    % 13. S, N, AOU or O2
    % 14. S, N
    % 15. S, AOU or O2
    % 16. S
    % 
% OAAdjustTF (Semi-optional boolean, default true): 
    % Ocean acidification continues to decrease seawater pH at varying
    % rates, making LIPHR estimates increasingly inappropriate as time
    % passes.  A true value for this boolean applies an adjustment to the
    % pH estimate to make it more appropriate for the user provided date of
    % the pH estimate.  The default adjustment is density dependent only
    % and will therefore be regionally too large or too small. See the
    % paper for details.  Alternately, the user may provide her/his own OA
    % rates for this adjustments as the 'OARate' optional input.  If the
    % OAAdjustTF boolean is set to true, then an estimate date or dates
    % must also be provided as EstDates.  This input or OAAdjustTF must be
    % supplied.  If you wish to use your own OA rate estimates, see line
    % 585 below for instructions for how to output the mean date of the
    % training data.
    %   
% OARate (Optional n by 1 vector, default determined by interpolation
% against potential density):
    % Users who prefer to provide their own OA rate estimates may do so
    % with this optional input.  Units are pH unit change per year, and
    % should in most cases be negative.
    %  
% EstDates (Semi-optional n by 1 array or 1 by 1 value, no default, but
% required if OAAdjustTF is true): 
    % A vector of decimal dates for the estimates (e.g. July 1 2020 would
    % be "2020.5").  If only a single date is supplied that value is used
    % for all estimates.  This input or OAAdjustTF must be supplied.
    % 
% MeasUncerts (Optional n by y array or 1 by y vector, default: [0.003 S,
% 0.003 degrees C (T or theta), 1% N, 1% AOU or O2, 1% Si]: 
    % Array of measurement uncertainties (see 'Measurements' for units).
    % Uncertainties should be presented in order indicated by MeasIDVec.
    % Providing these estimates will improve LIPHR estimate uncertainties
    % in 'UncertaintyEstimates'. Measurement uncertainties are a small part
    % of LIPHR estimate uncertainties for WOCE-quality measurements.
    % However, estimate uncertainty scales with measurement uncertainty, so
    % it is recommended that measurement uncertainties be specified for
    % sensor measurements.  If this optional input argument is not
    % provided, the default WOCE-quality uncertainty is assumed.  If a 1 by
    % y array is provided then the uncertainty estimates are assumed to
    % apply uniformly to all input parameter measurements.
    % 
% PerKgSwTF (Optional boolean, default true): 
    % Many sensors provide measurements in micromol per L (molarity)
    % instead of micromol per kg seawater. Indicate false if provided
    % measurements are expressed in molar units (concentrations must be
    % micromol per L if so).  Outputs will remain in molal units
    % regardless.
% 
% pHCalcTF (Optional boolean, default false): 
    % If set to true, LIPHR will recalculate the pH to be a better estimate
    % of what the seawater pH value would be if calculated from TA and DIC
    % instead of measured with purified m-cresol dye. This is arguably also
    % a better estimate of the pH that would be obtained from pre-2011
    % measurements with impure dyes.  See the LIPHR paper for details.
%    
% VerboseTF (Optional boolean, default true): 
    % Setting this to false will make LIPHR stop printing updates to the
    % command line.  This behavior can also be permanently disabled below.
% *************************************************************************
% Outputs:
%   
% pHEstimates: 
    % A n by e array of LIPHR estimates specific to the coordinates and
    % parameter measurements provided.  Units are pH_total units.
	%
% UncertaintyEstimates: 
    % A n by e array of LIPHR uncertainty estimates specific to the
    % coordinates, parameter measurements, and parameter uncertainties
    % provided.
    % 
% MinUncertaintyEquation: 
    % A n by 1 array of the index of the equation that returns the smallest
    % uncertainty for each estimate location.
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

% Determining how chatty you want the function to be.
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
if nargin<3; error('LIPHR called with too few input arguments.'); end

% Checking whether specific equations are specified.
a=strcmpi(varargin,'Equations');
if any(a)
    Equations=varargin{1,logical([0 a(1:end-1)])};
    SpecifiedEqn=true;
else
    Equations=[1:16];
    SpecifiedEqn=false;
end
% Making [] argument for Equations equivalent to no argument.
if isempty(Equations);Equations=1:16; end
% Making 0 argument for Equations equivalent to no argument.
if Equations==0; Equations=1:16; end

% Checking whether the minimum uncertainty estimate is desired vs. outputs
% from all used equations.
a=strcmpi(varargin,'MinUncEstTF');
if any(a)
    MinUncEst=varargin{1,logical([0 a(1:end-1)])};
else
    MinUncEst=true;
    if VerboseTF==true & SpecifiedEqn==false; disp('By Default, all equations with enough input variables will be used, but only the estimate with the lowest uncertainty will be returned.  If outputs from multiple equations are desired, specify the desired equations with the Equations input and include the input argument pair setting MinUncEstTF to false.');end;
end

% Checking for PerKgSwTF input and setting default if not given
a=strcmpi(varargin,'PerKgSwTF');
if any(a)
    PerKgSwTF=varargin{1,logical([0 a(1:end-1)])};
else
    PerKgSwTF=true;
end


% Checking for OAAdjustTF input and setting default if not given
a=strcmpi(varargin,'OAAdjustTF');
if any(a)
    OAAdjustTF=varargin{1,logical([0 a(1:end-1)])};
else
    OAAdjustTF=1;
end

% Checking whether OARate is provided, if not, default will be determined
% below
a=strcmpi(varargin,'OARate');
if any(a)
    Gamma=varargin{1,logical([0 a(1:end-1)])};
    CustomOA=true;
else
    CustomOA=false;
end

% Checking for MeasUncerts input and setting default if not given
a=strcmpi(varargin,'MeasUncerts');
if any(a)
    MeasUncerts=varargin{1,logical([0 a(1:end-1)])};
else
    MeasUncerts=[]; % This will be overriden immediately below
end
if nargin<5 || max(size(MeasUncerts))==0;
    % Setting multiplicative uncertainties for N, AOU, O2, and Si.
    DefaultUncertainties=diag([1 1 0.02 0.01 0.02 0.01 1]);
    MeasUncerts=Measurements*DefaultUncertainties(MeasIDVec,MeasIDVec);
    % Then setting additive uncertainties for T, S, and theta.
    MeasUncerts(:,ismember(MeasIDVec,[1 2 7]))=0.003;
end
% Sanity checking the MeasUncerts argument.  This also deals with the
% possibility that the user has provided a single set of uncertainties for
% all estimates.
if not(max(size(MeasUncerts)==size(Measurements))) && not(min(size(MeasUncerts)==size(Measurements))) && not(max(size(MeasUncerts))==0);
    error('MeasUncerts must be undefined, a vector with the same number of elements as Measurements has columns, [] (for default values), or an array of identical dimension to Measurements.')
elseif not(min(size(MeasUncerts)==size(Measurements))) && not(max(size(MeasUncerts))==0);
    if ~(size(MeasUncerts,2)==size(Measurements,2));
        error('There are different numbers of columns of input uncertainties and input measurements.')
    end
    MeasUncerts=ones(size(Measurements(:,1)))*MeasUncerts;
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
else
    EstDates=-9999;
end

% Making sure all provided predictors are identified.
if ~(size(MeasIDVec,2)==size(Measurements,2));
    error('The MeasIDVec input does not have the same number of columns as the Measurements input.  It is unclear to the routine which measurement is in which column.');
end

% Requriring date vector if OAAdjust is used.
if OAAdjustTF==true & EstDates==-99999;
    error('OAAdjustTF is set to true (by default) enabling an adjustment for OA, but the measurement dates needed for the adjustment are currently supplied by the user.  Either specify measurement dates as EstDates as a single year date (e.g. 2010.55 for July 1 2010) or as a vector of year dates for each measurement. Alternately, set OAAdjustTF to false to diable the adjustment.')
end

% Reshaping EstDates to have the same length as Coordinates if only a
% single date is provided.
if ~EstDates==-99999 & size(EstDates,1)==1; 
    EstDates=EstDates*ones(size(Coordinates,1),1);
end

% Making sure the user downloaded the needed file and put it somewhere it
% can be found
if exist('LIPHR_files.mat','file')<2; error('LIPHR could not find the LIPHR_files file needed to run.  This mandatory file should be distributed from the same website as LIPHR.  Contact the corresponding author if you cannot find it there.  If you do have it then make sure it is on the MATLAB path or in the active directory.'); end
 
% LIPHR requires non-NaN coordinates to provide an estimate.  This step
% eliminates NaN coordinate combinations prior to estimation.  NaN
% estimates will be returned for these coordinates.
NaNCoords=max(isnan(Coordinates),[],2);
if size(NaNCoords)==size(EstDates)
    EstDates(NaNCoords,:)=[];
end
[NumCoords,~]=size(Coordinates);
n=sum(~NaNCoords);
e=size(Equations,2);


% Checking for common missing data indicator flags and warning if any are
% found.  Consider adding other common ones here separated by spaces.
if max(ismember([-999 -9 -1*10^20],Measurements))==1
    warning('LIPHR: A common non-NaN missing data indicator (e.g. -999, -9, -1e20) was detected in the input measurements provided.  Missing data should be replaced with NaNs.  Otherwise, LIPHR will interpret your inputs at face value and give terrible estimates.')
end

% Checking to make sure all required parameters are provided
VarVec=[1 1 1 1 1
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
    1 0 0 0 0];
NeedVars=any(VarVec(Equations,:),1);
HaveVars=false(1,6);  HaveVars(1,MeasIDVec)=1; 

%Temperature can sub for potential temperature
if ismember(7,MeasIDVec)==1; HaveVars(1,2)=true; end

% Temperature or potential temperature is required if O2 is provided and
% AOU is not (unless AOU is not used).
if ismember(6,MeasIDVec)==1 && ismember(4,MeasIDVec)==0; NeedVars(1,2)=true; end

% If oxygen and temperature or potential temperature are present, then AOU
% is also calculate-able and effectively present
if ismember(6,MeasIDVec)==1 && or(ismember(7,MeasIDVec)==1,ismember(2,MeasIDVec)==1);
    HaveVars(1,4)=true; 
end

% Temperature or potential temperature is required if
% measurements are provided in molar units
if PerKgSwTF==false; NeedVars(1,2)=1; end 

% Eliminating equations.  LIARv1 used all equations if "Equations" was
% unspecified.  Here we limit this input to the equations for which all
% needed inputs are supplied.
if MinUncEst==true & SpecifiedEqn==false;
    Equations(:,any(VarVec(:,~HaveVars(1,1:5)),2))=[];
    NeedVars=any(VarVec(Equations,:),1);
    e=size(Equations,2);
end

% Making sure all needed variables are present
if ~all(HaveVars(1,NeedVars))
    error('LIPHR error: Check "Equations" and "MeasIDVec" inputs.  One or more regression equations selected require one or more input parameters that are either not provided or not labeled correctly.')
end


% Putting measurement inputs in expected order
M=NaN(n,7);  % Measurements
U=NaN(n,7);  % Uncertainties
M(:,MeasIDVec)=Measurements(~NaNCoords,:);
C=Coordinates(~NaNCoords,:);
U(:,MeasIDVec)=MeasUncerts(~NaNCoords,:);



% Checking to see whether potential temperature is defined and needed.
% Defining it and subbing in temp uncert for theta uncert if necessary.
if ismember(2,MeasIDVec)==0 && NeedVars(1,2);
    M(:,2)=sw_ptmp(M(:,1),M(:,7),sw_pres(C(:,3),C(:,2)),0);
    U(:,2)=U(:,7);
end
% Checking to see whether AOU is defined and needed. Defining it and
% subbing in O2 uncert for AOU uncert if necessary.
if ismember(4,MeasIDVec)==0 && NeedVars(1,4);
    M(:,4)=sw_satO2(M(:,1),M(:,2))*44.64-M(:,6);
    U(:,4)=U(:,6);
end
% Converting units to molality if they are provided as molarity.
if PerKgSwTF==false;
    densities=sw_dens(M(:,1),sw_temp(M(:,1),M(:,2), ...
        sw_pres(C(:,3),C(:,2)),0),sw_pres(C(:,3),C(:,2)))/1000;
    M(:,3)=M(:,3)./densities;
    M(:,4)=M(:,4)./densities;
    M(:,5)=M(:,5)./densities;
end

% *************************************************************************
% Beginning treatment of inputs and calculations

L=load('LIPHR_files.mat');
% Adjusting negative longitudes.
C(C(:,1)>360,1)=rem(C(C(:,1)>360,1),360);
C(C(:,1)<0,1)=rem(C(C(:,1)<0,1,1),360)+360;

% We don't want constants to be interpolated across the Panama Canal (for
% instance), so data is carved into two segments... the Atlantic/Arctic
% (AA) and everything else.  These interpolations are then done separately.
DepthToDegreeConversion=25; % See 2016 LIAR paper for explanation.
C(:,3)=C(:,3)/DepthToDegreeConversion;
L.Coords(:,3)=L.Coords(:,3)/DepthToDegreeConversion; 
AAIndsM=or(inpolygon(C(:,1),C(:,2),L.LNAPoly(:,1),L.LNAPoly(:,2)), ...
    or(inpolygon(C(:,1),C(:,2),L.LSAPoly(:,1),L.LSAPoly(:,2)), ...
    or(inpolygon(C(:,1),C(:,2),L.LNAPolyExtra(:,1),L.LNAPolyExtra(:,2)), ...
    or(inpolygon(C(:,1),C(:,2),L.LSAPolyExtra(:,1),L.LSAPolyExtra(:,2)), ...
    inpolygon(C(:,1),C(:,2),L.LNOPoly(:,1),L.LNOPoly(:,2))))));

% Here we predefine a dummy interpolant.  We fill it with the actual values
% later, but this is the slow step so we don't want it being done in a
% loop.  The appended "Dummy" values are designed to return the
% interpolated value at the nearest location within the domain when the
% desired estimate is outside of the domain.  Note that this is distinct
% from the scatteredInterpolant 'nearest' extrapolation option, which
% simply selectes the nearest neighbor of the input data points, and
% therefore results in discontinuities.  Dummy points accomplish this by
% interpolating between the values within the dataset and the mean value
% within the dataset 'placed' nearly infinitely far away.  Think of the 8
% dummy coordinates as the vertices of a very large cube around the data of
% interest.
Dummy=ones(8,1)*nanmean(nanmean(L.Cs,1),3);
Dummyset=horzcat([10^10; 10^10;10^10; 10^10;-10^10; -10^10;-10^10; ...
    -10^10],[10^10;-10^10;10^10;-10^10;10^10;-10^10;10^10;-10^10], ...
    [10^10;10^10;-10^10;-10^10;10^10;10^10;-10^10;-10^10]);
Dummy=horzcat(Dummyset,Dummy);

% Disambiguation:
% L.Cs (loaded from file)... pre-determined regression constants
% L.Coords (loaded from file)... coordinates of L.Cs
% C... viable subset of user-provided coordinates for LIPHR estimates
% LCs... locally-interpolated L.Cs
InterpolantAA=scatteredInterpolant(vertcat(L.Coords(L.AAIndsCs,1), ...
    Dummyset(:,1)),vertcat(L.Coords(L.AAIndsCs,2),Dummyset(:,2)), ...
    vertcat(L.Coords(L.AAIndsCs,3),Dummyset(:,3)), ...
    (1:size(vertcat(L.Coords(L.AAIndsCs,3),Dummyset(:,3))))');
InterpolantElse=scatteredInterpolant(vertcat(L.Coords(~L.AAIndsCs,1), ...
    Dummyset(:,1)),vertcat(L.Coords(~L.AAIndsCs,2),Dummyset(:,2)), ...
    vertcat(L.Coords(~L.AAIndsCs,3),Dummyset(:,3)), ...
    (1:size(vertcat(L.Coords(~L.AAIndsCs,3),Dummyset(:,3))))');

%Preallocating for speed
if MinUncEst==true & SpecifiedEqn==false;
    pHEstimates=NaN(NumCoords,1);
    UncertaintyEstimates=NaN(NumCoords,1);
else
    pHEstimates=NaN(NumCoords,e);
    UncertaintyEstimates=NaN(NumCoords,e);
end
pHEst=NaN(n,e);
UncertEst=NaN(n,e);

% Disambiguation:
% Eq... a counter for which equation LIPHR is on
% e... the number of equations that will be used
% Equation... the specific equation number currently being used
% Equations... the user provided list of equations that will be used.
for Eq=1:e;
    Equation=Equations(1,Eq);
    clear LCs    
%     selecting variables for this regression
    if Equation==1
        VarList='S Theta N AOU Si';
    elseif Equation==2
        VarList='S Theta N Si';
    elseif Equation==3
        VarList='S Theta AOU Si';
    elseif Equation==4
        VarList='S Theta Si';
    elseif Equation==5
        VarList='S Theta N AOU';
    elseif Equation==6
        VarList='S Theta N';
    elseif Equation==7
        VarList='S Theta AOU';
    elseif Equation==8
        VarList='S Theta';
    elseif Equation==9
        VarList='S N AOU Si';
    elseif Equation==10
        VarList='S N Si';
    elseif Equation==11
        VarList='S AOU Si';
    elseif Equation==12
        VarList='S Si';
    elseif Equation==13
        VarList='S N AOU';
    elseif Equation==14
        VarList='S N';
    elseif Equation==15
        VarList='S AOU';
    elseif Equation==16
        VarList='S';
    end
    UseVars=[1:5].*VarVec(Equation,:);
    UseVars(UseVars==0)=[];
%     Updating users if VerboseTF is set to true.
    if VerboseTF==true;
        disp(horzcat(num2str(floor(toc)),' seconds elapsed... now using eqn. ', ...
        num2str(Equation),': ',VarList,'.'))
    end
    
%     Determining dimensions
    SizUseVars=size(UseVars);
    LCs=NaN(n,SizUseVars(1,2)+1);
    NumVars=2+SizUseVars(1,2);
    % +2 below is for constant and depth terms
    VarNumVec=horzcat(1,2,UseVars+2);
%     Filling the interpolants with the predetermined regression constant
%     information and then using them to interpret 'local constants' (LCs)
%     at the coordinates of interest
    for Var=1:NumVars;
        InterpolantAA.Values=vertcat(L.Cs(L.AAIndsCs,VarNumVec(1,Var), ...
            Equation),Dummy(:,VarNumVec(Var)+3));
        InterpolantElse.Values=vertcat(L.Cs(~L.AAIndsCs,VarNumVec(1,Var), ...
            Equation),Dummy(:,VarNumVec(Var)+3));
        LCs(AAIndsM,Var)=InterpolantAA(C(AAIndsM,1),C(AAIndsM,2), ...
            C(AAIndsM,3));
        LCs(~AAIndsM,Var)=InterpolantElse(C(~AAIndsM,1), ...
            C(~AAIndsM,2),C(~AAIndsM,3));
    end
    % Estimating methodological error from depth
    EMLR=interp1(L.EMLRrec(:,1),L.EMLRrec(:,1+Equation),C(:,3));
    %  Estimating pH and pH estimate uncertainty from LCs.
    pHEst(:,Eq)=real(sum(LCs.*horzcat(ones(n,1),C(:,3)*DepthToDegreeConversion,M(:,UseVars)),2));
    % Estimating uncertainty
    UncertEst(:,Eq)=real((sum((LCs.*horzcat(zeros(n,2), ...
        U(:,UseVars))).^2,2)+0.004^2+EMLR.^2).^(1/2));
end
% Estimating terms for the OA adjust, if used
if OAAdjustTF==true;
    if VerboseTF==true;
        disp('Applying the ocean acidification adjustment.')
    end
    if min(HaveVars(1,1:2))==1;
        % Calculating potential density anomaly to determine rate of OA.
        EstimateDensity(:,1)=sw_dens0(M(:,1),M(:,2))-1000;
    else
        % Estimating potential density anomaly if temperature is not
        % provided.
        InterpolantAA.Values=vertcat(L.Coords(L.AAIndsCs,4), ...
            27.5*ones(size(Dummy(:,VarNumVec(Var)+3))));
        InterpolantElse.Values=vertcat(L.Coords(~L.AAIndsCs,4), ...
            27.5*ones(size(Dummy(:,VarNumVec(Var)+3))));
        EstimateDensity(AAIndsM,1)=InterpolantAA(C(AAIndsM,1),C(AAIndsM,2), ...
            C(AAIndsM,3));
        EstimateDensity(~AAIndsM,1)=InterpolantElse(C(~AAIndsM,1), ...
            C(~AAIndsM,2),C(~AAIndsM,3));
    end
    % Determining the mean date of the data used to estimate the Cs
    InterpolantAA.Values=vertcat(L.Dates(L.AAIndsCs,1), ...
        1995*ones(size(Dummy(:,VarNumVec(Var)+3))));
    InterpolantElse.Values=vertcat(L.Dates(~L.AAIndsCs,1), ...
        1995*ones(size(Dummy(:,VarNumVec(Var)+3))));
    EstimateCsDate(AAIndsM,1)=InterpolantAA(C(AAIndsM,1),C(AAIndsM,2), ...
        C(AAIndsM,3));
    EstimateCsDate(~AAIndsM,1)=InterpolantElse(C(~AAIndsM,1), ...
        C(~AAIndsM,2),C(~AAIndsM,3));
    if CustomOA==false;
        % Estimating default rate of OA by density if none is provided.
        Gamma=interp1(L.OAEstMat(:,1),L.OAEstMat(:,2:end),EstimateDensity);
    else
        Gamma=Gamma(~NaNCoords,1);
    end
    % Multiplying rate by time between measurement of data used for
    % estimating Cs and the estimate.
    OAAdjust=(EstDates-EstimateCsDate)*ones(1,(size(pHEst,2))).*repmat(Gamma(:,1),1,size(pHEst,2));
    pHEst=pHEst+OAAdjust;
end
% Determines which regression has the lowest uncertainty estimate.  If
% no regression was specified for outputs, this lowest-uncertainty
% pH and uncertainty estimate is returned.
[Temp,MinUncertaintyEquation]=min(UncertEst,[],2);
if MinUncEst==true & SpecifiedEqn==false;
    if VerboseTF==true;
        disp('Picking the outputs with the lowest uncertainties.')
    end
    pHEst=pHEst(sub2ind(size(pHEst),[1:1:size(MinUncertaintyEquation,1)]',MinUncertaintyEquation));
    UncertEst=UncertEst(sub2ind(size(UncertEst),[1:1:size(MinUncertaintyEquation,1)]',MinUncertaintyEquation));
end
% Recalculating pH values to be appropriate for pH calculated from TA and
% DIC if requested by the user.  See paper for coefficient sets.
if pHCalcTF==true
    if VerboseTF==true;
        disp('Recalculating the pH to be appropriate for pH values calculated from TA and DIC.')
    end
    % Solving equation 1 in the paper to counter adjustment.
    pHEst=(pHEst+0.3168)/1.0404;
end
% Filling the outputs with the estimates at viable locations.
pHEstimates(~NaNCoords,:)=pHEst;
UncertaintyEstimates(~NaNCoords,:)=UncertEst;
if VerboseTF==true; 
    disp(horzcat('LIPHR finished after ',num2str(round(toc,0)),' seconds.')); 
end
end
