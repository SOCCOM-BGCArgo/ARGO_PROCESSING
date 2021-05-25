%Test_pHvsS.m
% TEST pH sensitivity to salinity
fname = '9313SOOCN';
fnum = regexp(fname,'^\d+', 'match', 'once');
load(['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\', ...
    'cal',fname,'.mat']);
%d = parse_APEXmsg4ARGO('\\atlas\chemwebdata\floats\f9313\9313.001.msg');
target =['\\atlas\chemwebdata\floats\f',fnum,'\',fnum,'.055.msg'];
d = parse_APEXmsg4ARGO(target);
LR = d.lr_d;

Vrs    = LR(:,8);
Press  = LR(:,1);
Temp   = LR(:,2);
Salt   = LR(:,3);
k0     = cal.pH.k0;
k2     = cal.pH.k2;
Pcoefs = cal.pH.pcoefs;

% Salt off sets
S_offsets = [-10 -5 -2 0 2 5 10];
fill1 = ones(size(Press));
data = ones(size(Press,1),size(S_offsets,2))*NaN;
Soff = data;
for i = 1:size(S_offsets,2)
    S = Salt + S_offsets(i);
    [phfree,phtot]= phcalc_jp(Vrs, Press, Temp, S, k0, k2, Pcoefs);
    data(:,i) = phfree;
    Soff(:,i) =fill1 * S_offsets(i);
end


ind0 = S_offsets == 0;
ph_diffs = data - (data(:,ind0) * ones(size(S_offsets)));
Z = Press * ones(size(S_offsets));

figure(1)
set(gcf, 'Position', [488.2000 437.8000 801.6000 420.0000]);

ax(1) = subplot(1,2,1);
scatter(data(:), Z(:), 25, Soff(:), 'filled')
set(gca,'YDir','Reverse')
title(fname)
xlabel('pH in situ')
ylabel('Pressure, mbar')

ax(2) = subplot(1,2,2);
ax(2).Color = [224 224 224]./256;
scatter(ph_diffs(:), Z(:), 25, Soff(:), 'filled')
set(gca,'YDir','Reverse')
title(fname)
xlabel('pH diff from S = S+0')
ylabel('Pressure, mbar')

ax(3) =colorbar;
ax(3).Label.String = 'Salinity offset';
ax(3).Ticks = S_offsets;

ax(1).Position(3) = ax(1).Position(3) -0.05;
%ax(1).Color = [225 229 204]./256;
ax(2).Position(1) = ax(2).Position(1) -0.05;
ax(2).Position(3) = ax(1).Position(3);
%ax(2).Color = [224 224 224]./256;

