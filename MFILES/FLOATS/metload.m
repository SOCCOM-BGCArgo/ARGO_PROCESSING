function [dectime,lat,long,psp,slvp,dryt,relhum,windspd,windavg]=metload(pathname,filename);

if nargin==2; file_load=1; end;
if nargin<2; file_load=0; end;

%    dectime = decimal time (days); Jan 1 = 1
%    lat = decimal latitude (degrees)
%    long = decimal longitude (degrees)
%    psp = pyranometer readings (W m-2)   (Unused!)
%    slvp = sea-level pressure (mb)
%    dryt = dry air temperature (deg C)
%    relhum = relative humidity (%)
%    windspd = wind speed (m s-1)
%    windavg = 24hr mean wind speed (m s-1)

if file_load
    [hdr,data]=hdrload([pathname filename]);
    dectime=data(:,1);
    lat=data(:,2);
    long=data(:,3);
    slvp=data(:,4);
    dryt=data(:,5);
    relhum=data(:,6);
    windspd=data(:,7);
    winddir=data(:,8);
    psp=data(:,9);
    windavg=windspd;
else
    % Manually set the atmopheric parameters
    
    % Richard Davis original data sets
    %   dectime=(249:1/24:250)';
    %   lat=29.5 * ones(size(dectime));
    %   long=87 * ones(size(dectime));
    %   psp=100 * ones(size(dectime));
    %   slvp=1013.25 * ones(size(dectime));
    %   dryt=0 * ones(size(dectime));
    %   relhum=80 * ones(size(dectime));
    %   windspd=0 * ones(size(dectime));
    %   windavg=0 * ones(size(dectime));
    %
    %   dectime=128.710664;
    %   lat=61.0697 * ones(size(dectime));
    %   long=-26.6632 * ones(size(dectime));
    %   psp=100 * ones(size(dectime));
    %   slvp=29.9 * 33.863 * ones(size(dectime));  % 1 in = 33.863 mb
    %   dryt=16.5 * ones(size(dectime));   % With relhum = 80 --> precipwat ~= 2.5
    %   relhum=80 * ones(size(dectime));
    %   windspd=5.5 * ones(size(dectime));
    %   windavg=5.5 * ones(size(dectime));
    %
    %   dectime=141.5927;
    %   lat=61.4932 * ones(size(dectime));
    %   long=-25.0570 * ones(size(dectime));
    %   psp=100 * ones(size(dectime));
    %   slvp=29.9 * 33.863 * ones(size(dectime));  % 1 in = 33.863 mb
    %   dryt=16.5 * ones(size(dectime));   % With relhum = 80 --> precipwat ~= 2.5
    %   relhum=80 * ones(size(dectime));
    %   windspd=12.4980 * ones(size(dectime));
    %   windavg=12.4980 * ones(size(dectime));
    
    % North Atlantic
%     dectime=95.5
% %     dectime=[95:0.1:96]
%     lat=62.5 * ones(size(dectime));
%     long=30  * ones(size(dectime));
%     psp=100 * ones(size(dectime));
%     slvp=29.9 * 33.863 * ones(size(dectime));  % 1 in = 33.863 mb
%     dryt=0 * ones(size(dectime));
%     relhum=80 * ones(size(dectime));
%     windspd=10 * ones(size(dectime));
%     windavg=10 * ones(size(dectime));
    
    
     % North Atlantic
    dectime=95.5
%     dectime=[95:0.1:96]
    lat=62.5 * ones(size(dectime));
    long=30  * ones(size(dectime));
    psp=100 * ones(size(dectime));
    slvp=29.9 * 33.863 * ones(size(dectime));  % 1 in = 33.863 mb
    dryt=0 * ones(size(dectime));
    relhum=80 * ones(size(dectime));
    windspd=10 * ones(size(dectime));
    windavg=10 * ones(size(dectime));
    
end
