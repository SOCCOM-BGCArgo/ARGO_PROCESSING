function [M, MP, MOZ, MPST, OZ, precwat] = atmosph(lat,long,dectime,zenang,slvp,dryt,relhum)
%  CALCULATES ATMOSPHERIC VARIABLES

daynum = fix(dectime);

%  CALCULATE AIR MASS, PRESSURE-CORRECTED AIR MASS, STANDARD
%  AIR MASS, & OZONE AIR MASS.
H0 = 22.0;
M = 1.0./(cos(deg2rad(zenang)) + 0.15 .* (93.885 - zenang).^(-1.253));
MP = M .* slvp./1013.25;
MOZ = (1.+(H0/6370.))./(cos(deg2rad(zenang)).^2.+(2.*H0/6370.)).^0.5;
MPST = 1.8 .* slvp./1013.25;

%  SET Y FOR WATER VAPOR CALC DEPENDENT ON MONTH OF YEAR
% if MON >= 4 & MON <= 6
%  Y = -0.0229;
% else
%  Y = 0.02023;
% end

%  CALCULATE WATER VAPOR. UNITS =  CM.
% WP = exp(0.07074 * DP + Y);
%,WV = WP*(slvp/1013.25)^0.75 * (273.0/(DT+273.0))^0.5;

%  Calculate precipitable water from Leckner via Iqbal  UNITS =  CM
precwat = 0.493 * (relhum/100) .* (exp(26.23-5416./(dryt+273))) ./ (dryt + 273);

%  CALCULATE RELATIVE HUMIDITY.  UNITS = %
% A0=6.107799961;
% A1=4.436518521E-1;
% A2=1.428945805E-2;
% A3=2.650648471E-4;
% A4=3.031240396E-6;
% A5=2.034080948E-8;
% A6=6.136820929E-11;
% ES=A0+DT*(A1+DT*(A2+DT*(A3+DT*(A4+DT*(A5+A6*DT)))));
% E=A0+DP*(A1+DP*(A2+DP*(A3+DP*(A4+DP*(A5+A6*DP)))));
% RS=1000.0*0.622*ES/(slvp-ES);
% R=1000.0*0.622*E/(slvp-E);
% RH=(R/RS)*100.0;


%  CALCULATE OZONE. UNITS =  ATM*CM.
%  Northern Hemisphere
OZ=(235+(150+ 40*sin(0.9865*deg2rad(daynum-30))+20*sin(3*deg2rad(-long)).*sin(1.28*deg2rad(lat)).^2))/1000;

%  Southern Hemisphere
% DN=daynum+152.625;
% i=find(DN>365);
% if length(i)~=0, DN(i)=DN(i)-365;, end;
% OZ=(235+(100+ 30*sin(0.9865*deg2rad(DN))+20*sin(2*deg2rad(long-75)).*sin(1.5*deg2rad(lat)).^2))/1000;
