function [dd,lt,lg,psp,sp,dt,rh,ws,wa]=bozone_metload(ddd)
%  [dd,lt,lg,psp,sp,dt,rh,ws,wa,cs]=metload(ddd)

if nargin<=1; tstep=5; end;

load uv:bozone:met:metinput

i=find(fix(dectime)==ddd);

%  dstep=tstep/5;
dstep=1;

dd=dectime(i(1):dstep:i(length(i)));
lt=lat(i(1):dstep:i(length(i)));
lg=long(i(1):dstep:i(length(i)));
psp=psp(i(1):dstep:i(length(i)));
sp=slvp(i(1):dstep:i(length(i)));
dt=dryt(i(1):dstep:i(length(i)));
rh=relhum(i(1):dstep:i(length(i)));
ws=windspd(i(1):dstep:i(length(i)));
wa=windavg(i(1):dstep:i(length(i)));
cs=csi(i(1):dstep:i(length(i)));
