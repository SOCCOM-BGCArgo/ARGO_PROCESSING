function out=CANYON_pCO2(gtime,lat,lon,pres,temp,psal,doxy)
% function out=CANYON_pCO2(gtime,lat,lon,pres,temp,psal,doxy)
% 
% Multi-layer perceptron to predict CO2 partial pressure / uatm 
%
% Neural network training by Raphaëlle Sauzede, LOV; 
% as Matlab function by Henry Bittig, LOV
%
%
% input:
% gtime - date (UTC) as matlab time (days since 01-Jan-0000)
% lat   - latitude / °N  [-90 90]
% lon   - longitude / °E [-180 180] or [0 360]
% pres  - pressure / dbar
% temp  - in-situ temperature / °C
% psal  - salinity
% doxy  - dissolved oxygen / umol kg-1 (!)
%
% output:
% out   - pCO2 / uatm
%
% check value: 647.9297 uatm
% for 09-Dec-2014 08:45, 17.6° N, -24.3° E, 180 dbar, 16 °C, 36.1 psu, 104 umol O2 kg-1
%
% reference:
% Sauzède, R., Bittig, H.C., Claustre, H., Pasqueron de Fommervault, O.,
% Gattuso, J.-P., Legendre, L., and Johnson, K.S. (2017). 
% Estimates of water-column nutrient concentrations and carbonate system
% parameters in the global ocean: A novel approach based on neural networks. 
% Front. Mar. Sci. 4:128. doi:10.3389/fmars.2017.00128 
%
%
% requires CO2SYS-matlab:
% van Heuven et al. (2011). MATLAB Program Developed for CO2 System
% Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis
% Center, Oak Ridge National Laboratory, US Department of Energy,  Oak
% Ridge, Tennessee. http://dx.doi.org/10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1
%
% Henry Bittig, LOV
% 29.03.2017

% No input checks! Assumes informed use, e.g., same dimensions for all
% inputs, ...

basedir='../CANYON_Training/'; % relative or absolute path to CANYON training files

% input preparation
gvec=datevec(gtime);
doy=floor(datenum(gtime)-datenum(gvec(1),1,0))*360/365; % only full yearday used; entire year (365 d) mapped to 360°
year=gvec(:,1); % get year number
lon(lon>180)=lon(lon>180)-360;
% doy sigmoid scaling
presgrid=dlmread([basedir 'CY_doy_pres_limit.csv'],'\t');
[x,y]=meshgrid(presgrid(2:end,1),presgrid(1,2:end));
prespivot=interp2(x,y,presgrid(2:end,2:end)',lon(:),lat(:)); % Pressure pivot for sigmoid
fsigmoid=1./(1+exp((pres(:)-prespivot)./50));
% input sequence
%     lat,      sin(lon),    cos(lon),    cos(day),    sin(day),    year,    temp,   sal,    oxygen, P 
data=[lat(:)/90 sind(lon(:)) cosd(lon(:)) cosd(doy(:)).*fsigmoid(:) sind(doy(:)).*fsigmoid(:) year temp(:) psal(:) doxy(:) pres(:)./2e4+1./((1+exp(-pres(:)./300)).^3)];

temporaire=fopen([basedir 'Fichier_poids_pCO2_indirect_hidden_ascii18_8.sn'],'r');
poids=fscanf(temporaire,'%f',[3 inf]);
fclose(temporaire);
poids=poids(3,:)';

ne=10;  % Number of inputs
nc1=18; % Number of neurons of the first hidden layer
nc2=8; % Number of neurons of the second hidden layer
ns=1;   % Number of outputs

% WEIGHT AND BIAS PARAMETERS   
% weight and bias from the input layer to the first hidden layer w1, b1
b1=poids(1:nc1);    
w1=reshape(poids(nc1+nc2+ns+(1:ne*nc1)),nc1,ne);
% weight and bias from the first hidden layer to the second hidden layer
b2=poids(nc1+(1:nc2));
w2=reshape(poids(nc1+nc2+ns+ne*nc1+(1:nc1*nc2)),nc2,nc1);
% weight and bias from the second hidden layer to the output layer w3, b3    
b3=poids(nc1+nc2+(1:ns));
w3=reshape(poids(nc1+nc2+ns+ne*nc1+nc1*nc2+(1:nc2*ns)),ns,nc2);

%%% Mean and standard deviation of the training dataset
%%% These values are used to normalize the inputs parameters
Moy  =load([basedir 'moy_pCO2_indirect.dat']);
Ecart=load([basedir 'std_pCO2_indirect.dat']);

%%% NORMALISATION OF THE INPUT PARAMETERS
[rx,~]=size(data);
data_N=(2./3)*(data-(Moy(1:ne)*ones(1,rx))')./(Ecart(1:ne)*ones(1,rx))';

% Two hidden layers
a=1.715905*tanh((2./3)*(data_N*w1'+(b1*ones(1,rx))')); % input layer to first hidden layer
b=1.715905*tanh((2./3)*(     a*w2'+(b2*ones(1,rx))')); % first hidden layer to second hidden layer
y=                           b*w3'+(b3*ones(1,rx))';   % second hidden layer to output layer
%%% Y is the normalised output value of the neural network, 

%%% Denormalisation of the output of the NN for getting the true value 
y_rescaled=1.5*y*Ecart(ne+1)+Moy(ne+1);

% recalculate pCO2
outcalc=CO2SYS(2300,y_rescaled,1,2,35,25,NaN,0,NaN,0,0,1,10,1);
y_rescaled=outcalc(:,4);

% and put into same shape as the input variables
out=reshape(y_rescaled,size(pres));
