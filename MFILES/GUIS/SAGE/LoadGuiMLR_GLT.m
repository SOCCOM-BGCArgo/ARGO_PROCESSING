function MLR = LoadGuiMLR_GLT
% RETURN MLR COEFFICIENTS BASED ON DATE TYPE AND DATA SELECTION

% These algorithms should only be used to adjust float data at 1500m!!
% They will not accurately predict pH elsewhere in the water column.
% They have been trained using data between 1000 and 2100 m and should not
% be applied outside that depth range.

% ************************************************************************
%                     DEFINE NO3 MLR STRUCTURES updataed 11/01/2015
% ************************************************************************
Williams_50Sto80S.cC  =  543.674094604715;  % Constant
Williams_50Sto80S.cO  = -0.108034737; % Oxygen umol/kg
Williams_50Sto80S.cS  = -4.919884289; % Salinity
Williams_50Sto80S.cT  = -2.687795383; % Temperature deg C
Williams_50Sto80S.cST = -11.39462746; % Sigma Theta
Williams_50Sto80S.cP  =  0.000313895; % Pressure
MLR.Williams_50Sto80S.NO3 = Williams_50Sto80S;


Williams_30Sto50S.cC  =  441.543323389259;  % Constant
Williams_30Sto50S.cO  = -0.0906049023210556; % Oxygen umol/kg
Williams_30Sto50S.cS  = -11.1338569605989; % Salinity
Williams_30Sto50S.cT  = -2.00298256515697; % Temperature deg C
Williams_30Sto50S.cST = 0; % Sigma Theta
Williams_30Sto50S.cP  =  -0.00109487287188606; % Pressure
MLR.Williams_30Sto50S.NO3 = Williams_30Sto50S;


% % MLR.NO3.WOA2013 = [];
% % MLR.NO3.CANYON = [];
% % MLR.NO3.LIR_beta = [];

% ************************************************************************
%                     DEFINE PH MLR STRUCTURES
% ************************************************************************
%45S to 80S
Williams_50Sto80S.cC  =  1.380E+00;
Williams_50Sto80S.cO  =  1.802E-03; % Oxygen umol/kg
Williams_50Sto80S.cS  =  1.786E-01; % Salinity
Williams_50Sto80S.cT  =  7.482E-03; % Temperature deg C
Williams_50Sto80S.cST =  0; % Sigma Theta
Williams_50Sto80S.cP  = -3.966E-05; % Pressure
MLR.Williams_50Sto80S.PH = Williams_50Sto80S;

%  updated in email 5/29(19?)/2016 this is for South Pacific 30S to 50S
Williams_30Sto50S.cC  =  6.269E+00;
Williams_30Sto50S.cO  =  -2.787E-04; % Oxygen umol/kg
Williams_30Sto50S.cS  =  4.450E-02; % Salinity
Williams_30Sto50S.cT  =  3.476E-02; % Temperature deg C
Williams_30Sto50S.cST =  0; % Sigma Theta
Williams_30Sto50S.cP  = -2.104E-05; % Pressure
MLR.Williams_30Sto50S.PH = Williams_30Sto50S;


%************************************************************************
%                    DEFINE O2 COMPARISON MODE
%************************************************************************

MLR.Williams_50Sto80S.O2 = [];
MLR.Williams_30Sto50S.O2 = [];


%************************************************************************
%                    DEFINE T AND S COMPARISON MODE
%************************************************************************
MLR.Williams_50Sto80S.S = [];
MLR.Williams_30Sto50S.S = [];

MLR.Williams_50Sto80S.T = [];
MLR.Williams_30Sto50S.T = [];







