function [out] = addTAgpan(in)

out = in;

%% calc geopotential anomaly at 5m referenced to 1000 m (in dynamic meters)
out.gpan=NaN(size(out.Salinity));
out.gpan5to1000=NaN(size(out.Salinity));
for i=1:max(out.Station);
    A=out.Station==i...
        &~isnan(out.Salinity)...
        &~isnan(out.TemperatureC)...
        &~isnan(out.Depthm);
    gpan=flip(sw_gpan(flip(out.Salinity(A)),flip(out.TemperatureC(A)),flip(out.Depthm(A))));
    out.gpan(A)=gpan;
    gpan1000=meanNaN(gpan(out.Depthm(A)>960&out.Depthm(A)<1040));
    gpan50=meanNaN(gpan(out.Depthm(A)>47&out.Depthm(A)<53));
    gpan5to1000=meanNaN(gpan1000-gpan50);
    if ~isnan(gpan5to1000)
        out.gpan5to1000(A)=gpan5to1000;
    else out.gpan5to1000(A)=NaN;
    end
    clear A;
end
out.gpan5to1000dynm=out.gpan5to1000/10;

%% add Williams AlgALKALI (optimized for South of 45S)
% need to update the algorithm coefficients if/when they change

ALLALKALIcoeff=[734.722320773799;-15.4776417841503;-0.111489142191378;...
    59.7516385934224;-2.79046056637350;0.0134602966289807;-37.1209335266590];
    out.AlgALKALI = ALLALKALIcoeff(1)+...
        ALLALKALIcoeff(2)*out.SIGMATHETA+...
        ALLALKALIcoeff(3)*out.OXYGEN+...
        ALLALKALIcoeff(4)*out.Salinity+...
        ALLALKALIcoeff(5)*out.TemperatureC+...
        ALLALKALIcoeff(6)*out.Depthm+...
        ALLALKALIcoeff(7)*out.gpan5to1000dynm;

end