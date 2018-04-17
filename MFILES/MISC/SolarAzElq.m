function [Az,El] = SolarAzElq(UTC,Lat,Lon,Alt)
%Sun possition from time and location (matricised)
% [Az,El] = SolarAzEl(UTC,Lat,Lon,Alt)
%UTC: UTC Time (MatLab's datenum or 'yyyy-mm-dd HH:MM:SS' cellstr or char)
%Lat: Latitude [-90 90] (deg)
%Lon: Longitude [-180 180]  or [0 to 360E jp 03/23/17](deg)
%Alt: Altitude above sea level, optional (km)
%Az: Azimuth location of the sun (deg)
%El: Elevation location of the sun (deg)
%
%Example:
% [Az,El] = SolarAzElq('1991-05-19 13:00:00',50,10,0)
%
%Example:
% [UTC,Lat,Lon] = ndgrid(730486:1/24:730487,-90:5:90,-180:5:180);
% tic,[Az,El] = SolarAzElq(UTC,Lat,Lon);toc
%
%References:
% http://stjarnhimlen.se/comp/tutorial.html#5
% http://www.stargazing.net/kepler/altaz.html RA,DEC to Az,Alt

%History:
% Darin C. Koblick 02/17/2009 Authos
% Darin C. Koblick 04/16/2013 Vectorized
% Serge Kharabash 09/02/2016 Metricised

%Test:
% [UTC,Lat,Lon] = ndgrid(730486:10:730852,-90:5:90,-180:5:180);
% tic,[Az1,El1] = SolarAzElq(UTC(:),Lat(:),Lon(:));toc,t1=toc;
% tic,[Az2,El2] = SolarAzEl(UTC(:),Lat(:),Lon(:),0);toc,t2=toc;
% max(abs(Az1-Az2)),max(abs(El1-El2)),t2/t1

%defaults
if nargin<4 || isempty(Alt), Alt = 0; end
d2r = pi/180; %degrees to radiance conversion factor
r2d = 180/pi; %radiance to degrees conversion factor

%checks
if ischar(UTC)
    UTC = cellstr(UTC);
end
if iscell(UTC)
    UTC = reshape(datenum(UTC(:),'yyyy-mm-dd HH:MM:SS'),size(UTC));
end

% IF DEGREES EAST, CONVERT TO +/- 180
t1 = Lon > 180;
Lon(t1) = Lon(t1) - 360;

%julian date
[year,month,day,hour,min,sec] = datevec(UTC);
if ndims(UTC)>2 %#ok<ISMAT>
    year = reshape(year ,size(UTC));
    month = reshape(month,size(UTC));
    day = reshape(day ,size(UTC));
    hour = reshape(hour ,size(UTC));
    min = reshape(min ,size(UTC));
    sec = reshape(sec ,size(UTC));
end
[jd,UTH] = juliandate(year,month,day,hour,min,sec);
day = jd - 2451543.5;

%Keplerian elements for the Sun (geocentric)
w = 282.9404 + 4.70935e-5 * day; %longitude of perihelion degrees
e = 0.016709 - 1.151e-9 * day; %eccentricity
M = mod(356.0470 + 0.9856002585 * day, 360); %mean anomaly degrees
L = w + M; %Sun's mean longitude degrees
oblecl = (23.4393 - 3.563e-7 * day)*d2r; %Sun's obliquity of the ecliptic, rad

%auxiliary angle
E = M + r2d*e.*sin(M*d2r).*(1+e.*cos(M*d2r));

%rectangular coordinates in the plane of the ecliptic (x toward perhilion)
x = cos(E*d2r)-e;
year = sin(E*d2r).*sqrt(1-e.^2);

%distance and true anomaly
r = sqrt(x.^2 + year.^2);
v = atan2(year,x)*r2d;

%longitude of the sun
lon = v + w;

%ecliptic rectangular coordinates
xeclip = r.*cos(lon*d2r);
yeclip = r.*sin(lon*d2r);
zeclip = 0;

%rotate to equitorial rectangular coordinates
xequat = xeclip;
yequat = yeclip.*cos(oblecl) + zeclip*sin(oblecl);
zequat = yeclip.*sin(0.409115648642983) + zeclip*cos(oblecl);

%convert to RA and Dec
r = sqrt(xequat.^2 + yequat.^2 + zequat.^2) - (Alt/149598000); %roll up the altitude correction
RA = atan2(yequat,xequat); %rad
delta = asin(zequat./r); %rad

%local siderial time
GMST0 = mod(L+180,360)/15;
SIDTIME = GMST0 + UTH + Lon/15;

%replace RA with hour angle HA
HA = 15*SIDTIME - RA * r2d;

%convert to rectangular coordinate system
x = cos(HA*d2r).*cos(delta);
year = sin(HA*d2r).*cos(delta);
z = sin(delta);

%rotate along an axis going east-west
xhor = x.*cos((90-Lat)*d2r) - z.*sin((90-Lat)*d2r);
yhor = year;
zhor = x.*sin((90-Lat)*d2r) + z.*cos((90-Lat)*d2r);

%find Az and El
Az = atan2(yhor,xhor) * r2d + 180;
El = asin(zhor) * r2d;

function [jd,UTH] = juliandate(year,month,day,hour,min,sec)
%calculate julian date & J2000 value
UTH = hour + min/60 + sec/3600; %J2000
idx = month <= 2;
year(idx) = year(idx) - 1;
month(idx) = month(idx) + 12;
jd = floor(365.25*(year+4716)) + floor(30.6001*(month+1)) + 2 - ...
    floor(year/100) + floor(floor(year/100)/4) + day - 1524.5 + ...
    UTH/24;