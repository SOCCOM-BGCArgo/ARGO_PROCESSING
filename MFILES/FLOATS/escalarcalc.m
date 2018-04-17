function [Es,invmud]=escalarcalc(time,ID,IS,IT,zenang)
%  [Es,invmud]=escalarcalc(time,ID,IS,IT,zenang)

[nr, nc]=size(IT);

zenw=(asin(sin(deg2rad(zenang))./1.33)).*180./pi;

sunfact=cos(deg2rad(zenw*ones(1,nc)));

tempinvmud=(ID./IT)./sunfact + (IS./IT)./0.859;

invmud=zeros(size(IT));

i=find(~isnan(tempinvmud));

invmud(i)=tempinvmud(i);

Es=invmud.*IT;
