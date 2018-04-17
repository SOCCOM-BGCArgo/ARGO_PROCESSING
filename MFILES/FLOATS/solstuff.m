function [zenang, esd] = solstuff(lat,dectime)
% SOLSTUFF CALCULATES A VARIETY OF EARTH-SUN RELATIONSHIPS.
%  lat is positive N, in decimal degrees
%  dectime is decimal time in days (i.e., 1.5 is noon on Jan 1)

daynum=fix(dectime);
hour=rem(dectime,1).*24;
ly=0;
%  CALCULATE THE DAY ANGLE AS A FUNCTION OF LEAP YEAR
if ly==0
  dayang = 2.0*pi*(daynum-1.0)/365.0;
else
  dayang = 2.0*pi*(daynum-1.0)/366.0;
end

%  CALCULATE THE EARTH-SUN DISTANCE CORRECTION
esd = 1.00011 + 0.034221 * cos(dayang) + 0.00128 * sin(dayang) + 0.000719 * cos(2*dayang) + 0.000077 * sin(2*dayang);
%  CALCULATE sigma FOR SOLAR DEClINATION
sigma = 0.006918 - 0.399912*cos(dayang) + 0.070257*sin(dayang) - 0.006758*cos(2*dayang) + 0.000907*sin(2*dayang) - 0.002697*cos(3*dayang) + 0.00148*sin(3*dayang);
decl = sigma * 180.0/pi;
%  CALCULATE THE CENTERED HOUR ANGLE
hoang = 15.0 * (12.0 - (hour));
%  CALCULATE THE ZENITH ANGLE
cosZ = sin(deg2rad(lat)) .* sin(sigma) + cos(deg2rad(lat)) .* cos(sigma) .* cos(deg2rad(hoang));
zenang = acos(cosZ) .* 180.0/pi;
