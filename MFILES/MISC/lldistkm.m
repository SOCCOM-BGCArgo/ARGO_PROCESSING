function [d1km d2km]=lldistkm(latlon1,latlon2)
% JP modified code created by M Sohrabinia on 30 Oct 2012 - 10/28/2016
% From matlab file exchange File ID: #38812
%
% INPUTS
%     lat1lon1 = n x 2 lat lon array for origin point(s) 
%     lat2lon2 = n x 2 lat lon array for destination point(s) 
%
% OUTPUTS
%     d1km: distance calculated by Haversine formula
%           (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
%
%     d2km: distance calculated based on Pythagoran theorem
%           (see: http://en.wikipedia.org/wiki/Pythagorean_theorem)
%
% After: http://www.movable-type.co.uk/scripts/latlong.html
%
% EXAMPLES:
%     Example 1, short distance:
%         latlon1=[-43 172];
%         latlon2=[-44  171];
%         [d1km d2km]=distance(latlon1,latlon2)
%         d1km = 137.365669065197 (km)
%         d2km = 137.368179013869 (km)
%         d1km is approximately equal to d2km
%
%     Example 2, longer distance:
%         latlon1=[-43 172];
%         latlon2=[20  -108];
%         [d1km d2km]=distance(latlon1,latlon2)
%         d1km = 10734.8931427602 (km)
%         d2km = 31303.4535270825 (km)
%         d1km is significantly different from d2km
%         (d2km is not able to work for longer distances).
%
% First version: 15 Jan 2012
% Updated: 17 June 2012
% Modified by JP 10/28/2016 for vectors
%
% TESTING
% [rr,cc]  = size(Gdata);
% latlon1 = Gdata(:,5:6);
% latlon2 = ones(rr,1)* track(1,2:3);


% ************************************************************************
% ************************************************************************
radius=6371;

lat1 = latlon1(:,1)*pi./180; % Convert to radians
lat2 = latlon2(:,1)*pi./180;
lon1 = latlon1(:,2)*pi./180;
lon2 = latlon2(:,2)*pi./180;

deltaLat = lat2 - lat1;
deltaLon = lon2 - lon1;

a = sin((deltaLat)./2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon./2).^2;
c = 2.*atan2(sqrt(a),sqrt(1-a));
d1km = radius.*c;    %Haversine distance in km

x = deltaLon.*cos((lat1+lat2)./2);
y = deltaLat;
d2km = radius.*sqrt(x.*x + y.*y); %Pythagoran distance in km

%end