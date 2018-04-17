function [wave dectime zenang IT ID IS QT QD QS Eo] = radmod2
%  radmod
%  This code is a slight variation of the Gregg&Carder adaptation of the Bird&Riordan
%  spectral irradiance model.
%  
%  Output
%    wave   = wavelength (um)
%    dectime = decimal time (days); Jan 1 = 1
%    zenang = solar zenith angle (degrees)
%
%    Downwelling irradiance
%    IT = total spectral irradiance (J/(m^2 s um)) = ID + IS
%    ID = direct spectral irradiance (J/(m^2 s um))
%    IS = scattered (diffuse) spectral irradiance (J/(m^2 s um))
%
%    QT = total spectral irradiance (mole/(m^2 s um)) = QD + QS
%    QD = direct spectral irradiance (mole/(m^2 s um))
%    QS = scattered (diffuse) spectral irradiance (mole/(m^2 s um))
%
%    Scalar (4pi) irradiance
%    Eo = total spectral scalar irradiance (J/(m^2 s um))

%  Written by:     Richard Davis
%  Last Modified:  1-April-98


%  Set the scene
% hold off
% close all
% clear



%  Load the extraterrestrial solar spectrum and atmospheric transmissivities
%  Wavelength is in um
%  ah2o, ao2, ao3 are the spectral absorptionyd coefficients
%    of water, gas, and ozone, repectively (cm-1)
%  etr is the extraterrestrial irradiance (W m-2 um-1)
%    NOTE:  integration of etr does not yield the solar constant (1367 W m-2)
%           this is because etr here stops at 4 um.
load inatmos.mat
wave=awave;
clear awave

%  Load the meteorological data
%  If met data exists, set realmet to 1;
%  otherwise, set realmet to 0 to allow for the creation of fake met vectors
%  input
%    dd = day of year for which you want met data
%  output
%    dectime = decimal time (days); Jan 1 = 1
%    lat = decimal latitude (degrees)
%    long = decimal longitude (degrees)
%    psp = pyranometer readings (W m-2)  (Unused !)
%    slvp = sea-level pressure (mb)
%    dryt = dry air temperature (deg C)
%    relhum = relative humidity (%)
%    windspd = wind speed (m s-1)
%    windavg = 24hr mean wind speed (m s-1)
realmet=1;
if realmet==1
  [dectime,lat,long,psp,slvp,dryt,relhum,windspd,windavg]=metload(0);
else
  [dectime,lat,long,psp,slvp,dryt,relhum,windspd,windavg]=metload('uv:fun305:henlop305:','met.dat');
end

%  Set some constants
% viz = 23;  %  Unlimited horizontal visibility (km)
viz = 15;  %  horizontal visibility (km)
am = 1;    %  The Gregg-Carder factor set to marineam = 1;    %  The Gregg-Carder factor set to marine

%  Calculate zenith angle (deg) and earth-sun distance correction (unitless)
[zenang, esd] = solstuff(lat,dectime);

% ECR - properly calculate zenith angle using dectime in UTC and longitude.
[zalt sorad] = soradna1(dectime-1, 2008, -long, lat);
zenang = 90-zalt;

%  If zenith angle is >= 90 (i.e., sun hasn't risen), set it to 91 for easy indexing
i=find(zenang>=90);
zenang(i)=ones(size(i)) * 91;

%  Calculate various air masses, ozone and precipitable water concentrations
%  input
%    defined above
%  output
%    M = relative optical air mass
%    MP = pressure-corrected relative optical air mass
%    MOZ = relative optical ozone mass
%    MPST = pressure-corrected standard air mass
%    OZ = ozone concentration (atm cm)
%    precwat = precipitable water (cm)
[M, MP, MOZ, MPST, OZ, precwat] = atmosph(lat,long,dectime,zenang,slvp,dryt,relhum);

% NAB08 Day 128 from climatology
% OZ = 369.5/1000;  % 1000 Dobson Units = 1 ATM*CM   

%  Calculate various transmissivities and reflectances
%  input
%    defined above
%  output
%    TR = atmospheric transmittance after Rayleigh scattering
%    TA = aerosol transmittance
%    TW = water vapor transmittance
%    TO = ozone transmittance
%    TU = uniformly mixed gas transmittance
%    TAA = aerosol absorption
%    TAS = aerosol scattering
%    FS = fraction of aerosol scatter that is downward
%    RGD = total direct reflectance
%    RGS = total scatter reflectance
%    RS = sky reflectivity
[TR,TA,TW,TO,TU,TAA,TAS,FS,RGD,RGS,RS ] = ...
trans(lat,long,wave,ah2o,ao3,ao2,viz,slvp,relhum,precwat,windspd,windavg,am,zenang,M,MP,MOZ,MPST,OZ);
% clear M MP MOZ MPST OZ

%  Calculate spectral irradiance
%  input
%    defined above
%  output
%    IT = total spectral irradiance (J/(m^2 s um)) = ID + IS
%    ID = direct spectral irradiance (J/(m^2 s um))
%    IS = scattered (diffuse) spectral irradiance (J/(m^2 s um))
%    QT = total spectral irradiance (mole/(m^2 s um)) = QD + QS
%    QD = direct spectral irradiance (mole/(m^2 s um))
%    QS = scattered (diffuse) spectral irradiance (mole/(m^2 s um))
[QT, IT, QD, QS, ID, IS] = ...
speclite(lat,long,dectime,zenang,esd,wave,etr,ah2o,ao3,ao2,viz,slvp,relhum,precwat,windspd,windavg,am, ...
         TR,TA,TW,TO,TU,TAA,TAS,FS,RGD,RGS,RS);
         
%  If the sun hasn't risen, set the spectra to zero (this removes roundoff problems)         
i=find(zenang==91);
for j=1:length(wave)
  IT(i,j)=zeros(size(i));
  QT(i,j)=zeros(size(i));
end

% Estitmate scalar (4pi) irradiance
[Eo,invmud]=escalarcalc(dectime,ID,IS,IT,zenang);
