function [NitrateEstimates,UncertaintyEstimates,MinUncertaintyEquation]= ...
    LINR(Coordinates,Measurements,MeasIDVec, ...          % Required inputs
            varargin)                                     % Optional inputs
%  Version 2.0.1 (first LINR version, released alongside LIARv2.0)
%  Updated 2017.10.12: 
%       -"Molality" changed to "PerKgSw," 
%  Updated LINR_files 2018.10.10: Fixed a bug in the training routines.
%        New errors should be closer to estimated errors in the 2nd
%        citation. Older errors were averaged 0.01 umol/kg higher than
%        estimated. See Readme for details.
%
%  Locally Interpolated nitrate Regression (LINR): Estimates nitrate
%  and nitrate estimate uncertainty from combinations of other parameter
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
%  a negligible component of nitrate estimate error.  The seawater
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
    % Parameter measurements that will be used to estimate nitrate.  The
    % column order (y columns) is arbitrary, but specified by MeasIDVec.
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
    % 3. Phosphate
    % 4. AOU
    % 5. Silicate
    % 6. O2 (recommended over AOU if temperature is included so the user
    %        can be confident AOU will be calculated the same way it was
    %        calculated for the training dataset)
    % 7. Temperature (recommended over theta for same reason as O2 vs. AOU)
%
% *************************************************************************
% Optional inputs:  All remaining inputs must be specified as sequential
% input argument pairs (e.g. "...,'Equations',[1:16],'OAAdjustTF',false,
% etc...")
%
% Equations (optional 1 by e vector, default []): 
    % Vector indicating which equations will be used to estimate nitrate.
    % This input also determines the order of the columns in the
    % 'NitrateEstimates' output. If [] is input or the input is not
    % specified then all 16 equations will be used and only the outputs
    % from the equation with the lowest uncertainty estimate will be
    % returned.
    % 
    % (S=salinity, Theta=potential temperature, P=phosphate, Si=silicate,
    % AOU=apparent oxygen utilization, T=temperature, O2=dissolved oxygen
    % molecule... see 'Measurements' for units).
    % 
    % Output Equation Key: 
    % 1.  S, Theta or T, P, AOU or O2, Si
    % 2.  S, Theta or T, P, Si
    % 3.  S, Theta or T, AOU, Si
    % 4.  S, Theta or T, Si
    % 5.  S, Theta or T, P, AOU or O2
    % 6.  S, Theta or T, P
    % 7.  S, Theta or T, AOU or O2
    % 8.  S, Theta or T
    % 9.  S, P, AOU or O2, Si
    % 10. S, P, Si
    % 11. S, AOU or O2, Si
    % 12. S, Si
    % 13. S, P, AOU or O2
    % 14. S, P
    % 15. S, AOU or O2
    % 16. S
    % 
% MeasUncerts (Optional n by y array or 1 by y vector, default: [0.003 S,
% 0.003 degrees C (T or theta), 1% P, 1% AOU or O2, 1% Si]: 
    % Array of measurement uncertainties (see 'Measurements' for units).
    % Uncertainties should be presented in order indicated by MeasIDVec.
    % Providing these estimates will improve LINR estimate uncertainties
    % in 'UncertaintyEstimates'. Measurement uncertainties are a small part
    % of LINR estimate uncertainties for WOCE-quality measurements.
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
% VerboseTF (Optional boolean, default true): 
    % Setting this to false will make LINR stop printing updates to the
    % command line.  This behavior can also be permanently disabled below.
% *************************************************************************
% Outputs:
% 
% NitrateEstimates: 
    % A n by e array of LINR estimates specific to the coordinates and
    % parameter measurements provided.  Units are micromoles per kg
    % seawater.
	%
% UncertaintyEstimates: 
    % A n by e array of LINR uncertainty estimates specific to the
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
if nargin<3; error('LINR called with too few input arguments.'); end

% Checking whether specific equations are specified.
a=strcmpi(varargin,'Equations');
if any(a)
    Equations=varargin{1,logical([0 a(1:end-1)])};
    SpecifiedEqn=true;
else
    Equations=1:16;
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
    if VerboseTF==true && SpecifiedEqn==false; disp('By Default, all equations with enough input variables will be used, but only the estimate with the lowest uncertainty will be returned.  If outputs from multiple equations are desired, specify the desired equations with the Equations input and include the input argument pair setting MinUncEstTF to false.');end;
end

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

% Making sure all provided predictors are identified.
if ~(size(MeasIDVec,2)==size(Measurements,2));
    error('The MeasIDVec input does not have the same number of columns as the Measurements input.  This means it is unclear which measurement is in which column.');
end

% Making sure you downloaded the needed file and put it somewhere it
% can be found
if exist('LINR_files.mat','file')<2; error('LINR could not find the LINR_files file needed to run.  This mandatory file should be distributed from the same website as LINR.  Contact the corresponding author if you cannot find it there.  If you do have it then make sure it is on the MATLAB path or in the active directory.'); end
 
% LINR requires non-NaN coordinates to provide an estimate.  This step
% eliminates NaN coordinate combinations prior to estimation.  NaN
% estimates will be returned for these coordinates.
NaNCoords=max(isnan(Coordinates),[],2);
[NumCoords,~]=size(Coordinates);
n=sum(~NaNCoords);
e=size(Equations,2);

% Checking for common missing data indicator flags and warning if any are
% found.  Consider adding your common ones here
if max(ismember([-999 -9 -1*10^20],Measurements))==1
    warning('LINR: A common non-NaN missing data indicator (e.g. -999, -9, -1e20) was detected in the input measurements provided.  Missing data should be replaced with NaNs.  Otherwise, LINR will interpret your inputs at face value and give terrible estimates.')
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

% Eliminating equations.  Earlier we used all equations if "Equations" was
% unspecified.  Here we limit this input to the equations for which all
% needed inputs are supplied.
if MinUncEst==true && SpecifiedEqn==false;
    Equations(:,any(VarVec(:,~HaveVars(1,1:5)),2))=[];
    NeedVars=any(VarVec(Equations,:),1);
    e=size(Equations,2);
end

% Making sure all needed variables are present
if ~all(HaveVars(1,NeedVars))
    error('LINR error: Check "Equations" and "MeasIDVec" inputs.  One or more regression equations selected require one or more input parameters that are either not provided or not labeled correctly.')
end


% Putting measurement inputs in expected order
M=NaN(n,6);  % Measurements
U=NaN(n,6);  % Uncertainties
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

L=load('LINR_files.mat');
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
% C... viable subset of user-provided coordinates for LINR estimates
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
if MinUncEst==true && SpecifiedEqn==false;
    NitrateEstimates=NaN(NumCoords,1);
    UncertaintyEstimates=NaN(NumCoords,1);
else
    NitrateEstimates=NaN(NumCoords,e);
    UncertaintyEstimates=NaN(NumCoords,e);
end
NitrateEst=NaN(n,e);
UncertEst=NaN(n,e);

% Disambiguation:
% Eq... a counter for which equation LINR is on
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
    LCs=NaN(n,SizUseVars(1,2));
    NumVars=1+SizUseVars(1,2);
%     +1 above and below is for constant term, different in LIPHR
    VarNumVec=horzcat(1,UseVars+1);
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
    %  Estimating nitrate and nitrate estimate uncertainty from LCs.
    NitrateEst(:,Eq)=real(sum(LCs.*horzcat(ones(n,1),M(:,UseVars)),2));
    % Estimating uncertainty
    UncertEst(:,Eq)=real((sum((LCs.*horzcat(zeros(n,1), ...
        U(:,UseVars))).^2,2)+0.3^2+EMLR.^2).^(1/2));
end
% Setting negative LINR values to 0
NitrateEst(NitrateEst<0)=0;
% Determines which regression has the lowest uncertainty estimate.  If
% no regression was specified for outputs, this lowest-uncertainty
% nitrate and uncertainty estimate is returned.
[Temp,MinUncertaintyEquation]=min(UncertEst,[],2);
if MinUncEst==true && SpecifiedEqn==false;
    if VerboseTF==true;
        disp('Picking the outputs with the lowest uncertainties.')
    end
    NitrateEst=NitrateEst(sub2ind(size(NitrateEst),[1:1:size(MinUncertaintyEquation,1)]',MinUncertaintyEquation));
    UncertEst=UncertEst(sub2ind(size(UncertEst),[1:1:size(MinUncertaintyEquation,1)]',MinUncertaintyEquation));
end
% Filling the outputs with the estimates at viable locations.
NitrateEstimates(~NaNCoords,:)=NitrateEst;
UncertaintyEstimates(~NaNCoords,:)=UncertEst;
if VerboseTF==true; 
    disp(horzcat('LINR finished after ',num2str(round(toc,0)),' seconds.')); 
end
end
