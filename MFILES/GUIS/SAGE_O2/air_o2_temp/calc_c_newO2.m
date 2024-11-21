%calc_c_newO2.m

floatID = '9634';
dirs = SetupDirs_sO2;
DATA = getall_APEX_floatdata_sO2(dirs,floatID);
NCEP = getNCEP(DATA.track(:,1),DATA.track(:,3),DATA.track(:,4));
ncepO2 = NCEP.PRES;
pO2air = (ncepO2./100 - mean([DATA.O2air{1}{1}(:,11) DATA.O2air{1}{3}(:,11)],2)).*0.20946; %use mean of pH20 values derived from different air measurements
pO2inf = DATA.O2air{1}{3}(:,10); % air pO2 values (bladder inflated)
pO2def = DATA.O2air{1}{2}(:,10); % surface pO2 values (bladder deflated)
SDN = DATA.O2air{1}{3}(:,1);
M = [SDN pO2inf pO2def pO2air];

defO2sol = DATA.O2air{1}{2}(:,9);
defO2conc = DATA.O2air{1}{2}(:,8);
defO2sat = ((defO2conc).*100)./defO2sol;
% filename = 'f9752_O2calobs.csv';
% 
% fid = fopen(filename, 'w');
% fprintf(fid, 'SDN, pO2inf, pO2def, pO2air\n');
% fclose(fid)
% dlmwrite(filename, M, '-append', 'precision', '%.6f', 'delimiter', ',');

g=pO2air./pO2inf;
avg_g=nanmean(g);

y = pO2inf(1:end-1,:); %there's a nan, exclude last row
x = [pO2def pO2air];
x=x(1:end-1,:); 
B = (inv(transpose(x)*x))*(transpose(x)*y);
c = B(1);
d = B(2);
m = (1-c)/d;