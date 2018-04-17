function [TR,TA,TW,TO,TU,TAA,TAS,FS,RGD,RGS,RS ] = ...
trans(lat,long,wave,ah2o,ao3,ao2,viz,slvp,relhum,precwat,windspd,windavg,am,zenang,M,MP,MOZ,MPST,OZ)


R = [0.1, 1, 10];
R0 = [0.03, 0.24, 2];

%  CALCULATE SEA-SURFACE REFLECTANCE
%  Initialize the reflectance vectors
%  Set all except scattered specular to zero
RGDSP = zeros(length(windspd),1);	%  direct specular
RGSSP = ones(length(windspd),1);	%  scattered specular
RGF = zeros(length(windspd),1);		%  sea foam
RGD = zeros(length(windspd),1);		%  total direct reflectance
RGS = zeros(length(windspd),1);		%  total scattered reflectance (diffuse)

%  Calculate the direct specular
%  If wind speed is less than 2 m s-1, use fresnel only
if zenang>=40 & windspd>2 
  B=-7.14E-4*windspd+0.0618;
  RGDSP = 0.0253*exp(B.*(zenang-40));
else
  RGDSP = fresnel(zenang);
end

%  Set some constants relating wind stress to foam reflectance
D1=2.2E-5;
D2=4.0E-4;
D3=4.5E-5;
D4=4.0E-5;
RHOA=1.2E3;

%  If wind speed is less than 4 m s-1, specular equal to Burt (1966) for flat sea and no sea foam 
%  If wind speed greater than 4 m s-1, specular equal to Burt (1966) for rough sea and sea foam from
%    drag coefficient
if windspd<=4 
  RGSSP=0.066 .* RGSSP;
  RGF=0.0;
elseif windspd>4 & windspd<=7 
  RGSSP=0.057 .* RGSSP;
  CD=(0.62+1.56/windspd)*1.0E-3;		%  Drag coefficient for sea surface
  RGF=D1*RHOA*CD*(windspd.^2)-D2;
else
  RGSSP=0.057 .* RGSSP;
  CD=(0.49+0.065*windspd)*1.0E-3;
  RGF=(D3*RHOA*CD-D4).*(windspd.^2);
end
RGD=RGDSP+RGF;
RGS=RGSSP+RGF;
clear D1 D2 D3 D4 RHOA B CD RGDSP RGSSP RGF

%  If either of the reflectances are greater than one, set to one
i=find(RGD>1);
RGD(i)=1;
i=find(RGS>1);
RGS(i)=1;

%  CALCULATE ALPHA AND BETA FOR RAYLEIGH TRANSMITTANCE
F=((2-relhum/100)./(6*(1-relhum/100))).^(1/3);
i = find(isinf(F));
F(i) = ones(size(i)) * 5.51;
A=ones(length(windspd),3);
if length(am)==1
  A(:,1)=A(:,1)*2000*(am^2);
else
  A(:,1)=2000*(am.^2);
end
A(:,2)=max(0.5, 5.866*(windavg-2.2));
A(:,3)=max(1.4E-5, 0.01527*(windspd-2.2)*0.05);
for J =1:3
  DNDR=0;
  for I = 1:3
    NR=A(:,I).*exp(-1*(log(R(J)/R0(I).*F)).^2)./F;
    DNDR=DNDR+NR;
  end
  Y(:,J)=log10(DNDR);
end 
xin = log10(R);
for j = 1:size(Y,1)
  yin = Y(j,:);
%   G(j) = fit(xin,yin);
  p = polyfit(xin, yin,1);   % ECR:  See G&C eqn 25, 26
  G(j) = p(1);
end
ALPHA=-1.0*(G+3)';
if length(viz)~=1
  i = find(viz<5);
  viz(i) = ones(size(i)) * 5;
end
BETA = 0.55.^ALPHA .* (3.912./viz - 0.01162) .* (0.02472 * (viz - 5) + 1.132);
clear viz windspd windavg R R0 F DNDR NR A xin yin i j G

%  CALCULATE A WHOLE BUNCH OF STUFF FOR SCATTERING Ts
if ALPHA<0 
  cosT=0.82;
elseif ALPHA>1.2 
  cosT=0.65;
else
  cosT = -0.1417*ALPHA+0.82;
end
ALG = log(1 - cosT);
AFS = ALG .* (1.4590 + ALG .* (0.15950 + ALG .* 0.4129));
BFS = ALG .* (0.0783 + ALG .* (-0.3824 - ALG .* 0.5874));
FS = 1 - 0.5 .* exp((AFS + BFS.*cos(deg2rad(zenang))) .* cos(deg2rad(zenang)));
FSP = 1 - 0.5 .* exp((AFS + BFS./1.8)./1.8);
clear zenang ALG AFS BFS cosT

%  LOOP THROUGH THE SPECTRA AND CALCULATE TRANSMISSIVITIES
TR = zeros(length(precwat),length(wave));
TA = zeros(length(precwat),length(wave));
TW = zeros(length(precwat),length(wave));
TO = zeros(length(precwat),length(wave));
TU = zeros(length(precwat),length(wave));
TAA = zeros(length(precwat),length(wave));
TAS = zeros(length(precwat),length(wave));
RS = zeros(length(precwat),length(wave));
for J = 1:length(wave)
  TR(:,J)=exp(-MP./(wave(J).^4.*(115.6406-(1.3366./wave(J).^2))));
  TRP=exp(-MPST./(wave(J).^4.*(115.6406-(1.3366./wave(J).^2))));
  Tao2 = BETA .* wave(J).^(-1.0.*ALPHA);
  TA(:,J) = exp(-Tao2 .* M);
  TW(:,J)=exp(-0.2385.*ah2o(J).*precwat.*M./(1+20.07.*ah2o(J).*precwat.*M).^0.45);
  TWP=exp(-0.2385.*ah2o(J).*precwat.*1.8./(1+20.07.*ah2o(J).*precwat.*1.8).^0.45);
  TO(:,J) = exp(-ao3(J) .* OZ .* MOZ);
  TOP = exp(-ao3(J) .* OZ .* 1.8);
  TU(:,J) = exp(-1.41 .* ao2(J) .* MP./(1 + 118.3 .*ao2(J) .* MP).^0.45);
  OMEGA = (-0.0032 .* am + 0.977) .* exp(-0.001185 .* relhum .* (log(wave(J)./0.4)).^2);
  TAA(:,J) = exp(-1 .* (1 - OMEGA) .* Tao2 .* M);
  TAAP = exp(-1 .* (1 - OMEGA) .* Tao2 .* 1.8);
  TAS(:,J) = exp(-OMEGA .* Tao2 .* M);
  TASP = exp(-OMEGA .* Tao2 .* 1.8);
  RS(:,J)=TOP.*TWP.*TAAP.*(0.5.*(1.0-TRP)+(1-FSP).*TRP.*(1-TASP));
end
