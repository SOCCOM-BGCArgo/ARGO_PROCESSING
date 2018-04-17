function [QT, IT, QD, QS, ID, IS] = ...
speclite(lat,long,dectime,zenang,esd,wave,etr,ah2o,ao3,ao2,viz,slvp,relhum,precwat,windspd,windavg,am, ...
         TR,TA,TW,TO,TU,TAA,TAS,FS,RGD,RGS,RS)

%  LOOP THRU SPECTRA AND CALCULATE IRRADIANCES (J/(m^2 s um))
%  Convert to quantum (mol/(m^2 s um))
IE=zeros(size(TR));
ID=zeros(size(TR));
IS=zeros(size(TR));
IT=zeros(size(TR));
QD=zeros(size(TR));
QS=zeros(size(TR));
QT=zeros(size(TR));

for J=1:length(wave)
  IE(:,J)=etr(J).*esd.*cos(deg2rad(zenang));
  ID(:,J)=IE(:,J).*TR(:,J).*TA(:,J).*TW(:,J).*TO(:,J).*TU(:,J);
  IR=IE(:,J).*TO(:,J).*TU(:,J).*TW(:,J).*TAA(:,J).*(1-(TR(:,J).^0.95)).*0.5;
  IA=IE(:,J).*TO(:,J).*TU(:,J).*TW(:,J).*TAA(:,J).*(TR(:,J).^1.5).*(1.0-TAS(:,J)).*FS;
  IG=ID(:,J).*cos(deg2rad(zenang)).*RS(:,J).*RGD./(1-RS(:,J).*RGD)+(IR+IA).*RS(:,J).*RGS./(1-RS(:,J).*RGS);
  if wave(J)<=0.45
    CS=(wave(J)+0.55).^1.8;
  else
    CS=1;
  end
  IS(:,J)=(IR+IA+IG).*CS;
  
  IS(:,J)=IS(:,J).*(1-RGS);
  ID(:,J)=ID(:,J).*(1-RGD);
  IT(:,J)=ID(:,J)+IS(:,J);
  
  QD(:,J)=ID(:,J).*(wave(J)./(6.6260755E-34*2.99792458E14*6.0221367E23)).*(1-RGD);
  QS(:,J)=IS(:,J).*(wave(J)./(6.6260755E-34*2.99792458E14*6.0221367E23)).*(1-RGS);
  QT(:,J)=QD(:,J)+QS(:,J);
end
