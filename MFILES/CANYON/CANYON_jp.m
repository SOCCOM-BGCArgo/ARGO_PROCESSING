function out = CANYON_jp(gtime,lat,lon,pres,temp,psal,doxy,param)
% Multi-layer perceptron to predict oceanographic variables
% Neural network training and original R function by Raphaëlle Sauzede,
% LOV, as Matlab function by Henry Bittig, 24.05.2016, LOV
%
% REFERENCE:
% Sauzède, R., Bittig, H.C., Claustre, H., Pasqueron de Fommervault, O.,
% Gattuso, J.-P., Legendre, L., and Johnson, K.S. (2017). 
% Estimates of water-column nutrient concentrations and carbonate system
% parameters in the global ocean: A novel approach based on neural networks. 
% Front. Mar. Sci. 4:128. doi:10.3389/fmars.2017.00128
% 

% 04/26/17 - Modified by Josh Plant. Combined individual property
%       estimator functions into one function. Only difference is an
%       added input string at end of function call and a look up table to 
%       define files to use for given varaible estimate. IF TRAINING DATA
%       IS UPDATED & FILE NAMES CHANGE, UPDATE TABLE!
%
% INPUTS:
%    gtime - date (UTC) as matlab time (days since 01-Jan-0000)
%    lat   - latitude / °N  [-90 90]
%    lon   - longitude / °E [-180 180] or [0 360]
%    pres  - pressure / dbar
%    temp  - in-situ temperature / °C
%    psal  - salinity
%    doxy  - dissolved oxygen / umol kg-1 (!)
%    param - a char array representing the desired seawater property must
%            be one of the following:
%               NO3  - nitrate                            - µmol / kg
%               PO4  - dissolved inorganic phosphate      - µmol / kg
%               Si   - dissovled inorganic silicate       - µmol / kg
%               PH   - total scale at insitu PTS
%               TA   - total alkalinity                   - µmol / kg
%               DIC  - dissolved inorganic carbon         - µmol / kg
%               pCO2 - partial pressure of carbon dioxide - µatm
%
% OUTPUT:
% out   - prediction of oceanographic variable at lat, lon, time and pres
%
% CHECK VALUES FOR:
%   gtime = datenum(2014,12,9,8,45,00); % 09-Dec-2014 08:45
%   lat   = 17.6;   % 17.6° N
%   lon   = -24.3;  % 24.3° W
%   pres  = 180;    % 180 dbar
%   temp  = 16;     % 16 °C
%   psal  = 36.1;   % 36.1 psu
%   doxy  = 104;    % 104 umol O2 kg-1
%   param = 'pCO2'; % ADJUSTABLE
%
%   NO3  -  18.57849 umol kg-1  [18.57849 JP] 
%   PO4  -  1.062308 umol kg-1  [1.062308  JP] 
%   Si   -  5.603550 umol kg-1  [5.60355  JP]
%   PH   -  7.868775 total scale, insitu PTS [7.868775 JP]
%   TA   -  2358.330 umol kg-1  [2358.330 JP]
%   DIC  -  2197.551 umol kg-1  [2197.551 JP]
%   pCO2 -  647.9297 uatm       [647.9297 JP]
%
% reference:
% Sauzède et al. (2017). Estimates of water-column nutrients and carbonate
% system parameters in the global ocean: A novel approach based on neural
% networks. Frontiers in Marine Science.
%
% Coded for Matlab by Henry Bittig, LOV, 29.03.2017


% No input checks! Assumes informed use, e.g., same dimensions for all
% inputs, ...

% ************************************************************************
% PATHS & LOOKUP TABLE
% ************************************************************************
user_dir = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
data_dir = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\DATA\CANYON\'];

% CHECK TO MAKE SURE DATA DIR EXISTS
dir_chk = ls(data_dir);
if isempty(dir_chk)
    disp('Data directory for canyon_jp can not be found!!')
    disp('Set the path in canyon_jp.m to point to the canyon data files.');
    disp('~line 63, "data_dir" variable');
    out = pres*NaN;
    return
end

% DEFINE LOOKUP TABLE
% CHANGE FILE NAMES FOR NEW VERSIONS
% {VAR, fname wts, fname mean train, fname std train,
% [#inputs, # neurons 1st layer, # neurons 2nd layer, # outputs]}
info(1,:) = {'NO3', 'Fichier_poids_NO3_hidden_ascii20_17.sn', ...
             'moy_NO3.dat', 'std_NO3.dat', [9, 20, 17, 1]};
         
info(2,:) = {'PO4', 'Fichier_poids_PO4_hidden_ascii17_17.sn', ...
             'moy_PO4.dat', 'std_PO4.dat', [9, 17, 17, 1]};
         
info(3,:) = {'Si', 'Fichier_poids_SI_hidden_ascii20_15.sn', ...
             'moy_SI.dat', 'std_SI.dat', [9, 20, 15, 1]};  
         
info(4,:) = {'PH', 'Fichier_poids_PH_hidden_ascii19_8.sn', ...
             'moy_PH.dat', 'std_PH.dat', [10, 19, 8, 1]};  
         
info(5,:) = {'TA', 'Fichier_poids_AT_hidden_ascii19_17.sn', ...
             'moy_AT.dat', 'std_AT.dat', [9, 19, 17, 1]};  
         
info(6,:) = {'DIC', 'Fichier_poids_CT_hidden_ascii18_15.sn', ...
             'moy_CT.dat', 'std_CT.dat', [10, 18, 15, 1]};
         
info(7,:) = {'pCO2', 'Fichier_poids_pCO2_indirect_hidden_ascii18_8.sn', ...
             'moy_pCO2_indirect.dat', 'std_pCO2_indirect.dat', [10, 18, 8, 1]};         
         
% ************************************************************************
% CHECK CHARACTER ARRAY INPUT / LOAD PERTINENT INFO FROM LOOKUP TABLE
%valid_param = {'NO3' 'PO4' 'Si' 'PH' 'TA' 'DIC' 'pCO2'};
% ************************************************************************
t1 = strcmp(param,info(:,1));
if sum(t1) == 1
    wt_tmp = info{t1,2};
    wt_fn  = ls([data_dir,wt_tmp]);
    
    if isempty(wt_fn)
        disp(['Could not find weight file for ',param])
        out = [];
        return
    end
    
    mean_fn = info{t1,3};
    std_fn  = info{t1,4};
    
    ne  = info{t1,5}(1); % Number of inputs
    nc1 = info{t1,5}(2); % Number of neurons of the 1st hidden layer
    nc2 = info{t1,5}(3); % Number of neurons of the 2nd hidden layer
    ns  = info{t1,5}(4); % Number of outputs
    
else
    disp('Param must be one of the following character arrays:')
    disp(info(:,1)')
    out = [];
    return
end

% ************************************************************************
% INPUT PREP
% ************************************************************************
gvec = datevec(gtime);
% only full yearday used; entire year (365 d) mapped to 360°
doy  = floor(datenum(gtime) - datenum(gvec(1),1,0))*360/365; 
year = gvec(:,1); % get year number
lon(lon>180) = lon(lon>180)-360;

% doy sigmoid scaling
presgrid  = dlmread([data_dir, 'CY_doy_pres_limit.csv'],'\t');
[x,y]     = meshgrid(presgrid(2:end,1),presgrid(1,2:end));
prespivot = interp2(x,y,presgrid(2:end,2:end)',lon(:),lat(:)); % Pressure pivot for sigmoid
fsigmoid  = 1./(1+exp((pres(:)-prespivot)./50));

% CHOOSE PROPER INPUT SEQUENCE BASED ON PARAMETER
switch param
    
    case 'NO3' % NO3 independent of year sin(doy) & cos(doy) reversed       
         data = [lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
                 sind(doy(:)).*fsigmoid(:) cosd(doy(:)).*fsigmoid(:), ...
                 temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];
   
    case 'PO4' % PO4 independent of year sin(doy) & cos(doy) reversed
         data=[lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
             sind(doy(:)).*fsigmoid(:) cosd(doy(:)).*fsigmoid(:), ...
             temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];
   
    case 'Si' % Si independent of year sin(doy) & cos(doy) reversed
        data=[lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
            sind(doy(:)).*fsigmoid(:) cosd(doy(:)).*fsigmoid(:), ...
            temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];    
            
    case 'PH'
        data=[lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
            cosd(doy(:)).*fsigmoid(:) sind(doy(:)).*fsigmoid(:), year, ...
            temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];    
            
    case 'TA' % TA independent of year 
       data=[lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
           cosd(doy(:)).*fsigmoid(:) sind(doy(:)).*fsigmoid(:), ...
           temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];
     
    case 'DIC'
        data=[lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
            cosd(doy(:)).*fsigmoid(:) sind(doy(:)).*fsigmoid(:), year, ...
            temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];    
            
    case 'pCO2'
        data=[lat(:)/90 sind(lon(:)) cosd(lon(:)), ...
            cosd(doy(:)).*fsigmoid(:) sind(doy(:)).*fsigmoid(:), year, ...
            temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];        
end

% ************************************************************************
% WEIGHTS, BIAS AND NORMALIZATION
% ************************************************************************
fid   = fopen([data_dir, wt_fn],'r');
poids = fscanf(fid,'%f',[3 inf]);
fclose(fid);
poids = poids(3,:)';

% weight and bias from the input layer to the first hidden layer w1, b1
b1 = poids(1:nc1);    
w1 = reshape(poids(nc1+nc2+ns+(1:ne*nc1)),nc1,ne);
% weight and bias from the first hidden layer to the second hidden layer
b2 = poids(nc1+(1:nc2));
w2 = reshape(poids(nc1+nc2+ns+ne*nc1+(1:nc1*nc2)),nc2,nc1);
% weight and bias from the second hidden layer to the output layer w3, b3    
b3 = poids(nc1+nc2+(1:ns));
w3 = reshape(poids(nc1+nc2+ns+ne*nc1+nc1*nc2+(1:nc2*ns)),ns,nc2);

%%% Mean and standard deviation of the training dataset
%%% These values are used to normalize the inputs parameters
Moy   = load([data_dir, mean_fn]);
Ecart = load([data_dir, std_fn]);

%%% NORMALISATION OF THE INPUT PARAMETERS
[rx,~] = size(data);
data_N = (2./3)*(data-(Moy(1:ne)*ones(1,rx))')./(Ecart(1:ne)*ones(1,rx))';

% TWO HIDDEN LAYERS
% input layer to first hidden layer
a= 1.715905*tanh((2./3)*(data_N*w1'+(b1*ones(1,rx))')); 

% first hidden layer to second hidden layer
b= 1.715905*tanh((2./3)*(     a*w2'+(b2*ones(1,rx))')); 

% second hidden layer to output layer
y=                            b*w3'+(b3*ones(1,rx))';  

% Y is the normalised output value of the neural network, 
% Denormalisation of the output of the NN for getting the true value 
y_rescaled = 1.5*y*Ecart(ne+1)+Moy(ne+1);

if strcmp(param,'pCO2') % recalculate pCO2
    % recalculate pCO2
    outcalc=CO2SYS(2300,y_rescaled,1,2,35,25,NaN,0,NaN,0,0,1,10,1);
    y_rescaled=outcalc(:,4);
end

% and put into same shape as the input variables
out = reshape(y_rescaled,size(pres));