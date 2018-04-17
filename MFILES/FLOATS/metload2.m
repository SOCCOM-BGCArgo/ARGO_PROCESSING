function [dectime,lat,long,psp,slvp,dryt,relhum,windspd,windavg]=metload2(dectime);


% North Atlantic
% dectime=95.5

%    dectime = decimal time (days); Jan 1 = 1
%    lat = decimal latitude (degrees)
%    long = decimal longitude (degrees)
%    psp = pyranometer readings (W m-2)   (Unused!)
%    slvp = sea-level pressure (mb)
%    dryt = dry air temperature (deg C)
%    relhum = relative humidity (%)
%    windspd = wind speed (m s-1)
%    windavg = 24hr mean wind speed (m s-1)


%     dectime=[95:0.1:96]
lat=60 * ones(size(dectime));
long=-20  * ones(size(dectime));
% lat=62.5 * ones(size(dectime));
% long=30  * ones(size(dectime));
psp=100 * ones(size(dectime));
slvp=29.9 * 33.863 * ones(size(dectime));  % 1 in = 33.863 mb
dryt=0 * ones(size(dectime));
relhum=80 * ones(size(dectime));
windspd=10 * ones(size(dectime));
windavg=10 * ones(size(dectime));


