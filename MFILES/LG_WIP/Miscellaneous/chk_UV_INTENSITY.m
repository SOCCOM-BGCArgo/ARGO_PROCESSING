%mbari_id = 'ua19605'
%mbari_id = 'ua20704' % cycle 20 wavy baseline
%mbari_id = 'ua21951' % cycle 11
%mbari_id = 'wn1354' % cycle 11
%mbari_id = 'wn1347' % cycle 11
%mbari_id = 'ua21900' % cycle 11
%mbari_id = 'wn1488'; % cycle 1
%mbari_id = 'wn1351'; % 
%mbari_id = 'wn1476'; %
%mbari_id = 'un1504';
%mbari_id = 'ua21996';%

%mbari_id = 'ua12386';% high baseline slope
%mbari_id = 'ua12398';% onlly 2 samples
%mbari_id = 'ua12549';% large hump @ 230-245
%mbari_id = 'ua12702';% funky profile - fouling
%mbari_id = 'ua12757';% wierd spectral shape
%mbari_id = 'ua18299';% high baseline slope
%mbari_id = 'ua20138';% hump @ 230-245
%mbari_id = 'ua12559';% hump @ 230-245
% %mbari_id = 'ua12757';% 
%mbari_id = 'ua19085';% hump @ 230-245
%mbari_id = 'ua21910'
%mbari_id = 'ua20602'
%mbari_id = 'wn1359'
%mbari_id = 'wn1350'
%mbari_id = 'ua21900';
% mbari_id = 'ua21105';
% mbari_id = 'ua21900';
%mbari_id = 'wn1350';
%mbari_id = 'ua21519';
% mbari_id = 'un1515';
%mbari_id = 'ua21075';
%mbari_id = 'ua17182';
%mbari_id = 'ua19605';
%mbari_id = 'ua21267';
%mbari_id = 'ua20209';
%mbari_id = 'ua21441';
%mbari_id = 'ua19006';
%mbari_id = 'un1509';
%mbari_id = 'ua21818';
%mbari_id = 'ua20579';
%mbari_id = 'ua21792';
%mbari_id = 'ua19085';
% mbari_id = 'ua20103';
%mbari_id = 'ua21105';
% mbari_id = 'un1204';
%mbari_id = 'ua21479';
% mbari_id = 'ua20209';
%mbari_id = 'ua19006';
%mbari_id = 'ua21535';
%mbari_id = 'un1525';
%mbari_id = 'ua19142';
%mbari_id = 'ua21883';
%mbari_id = 'ua21910';
%mbari_id = 'un1501';
%mbari_id = 'un1502';
%mbari_id = 'ua21984';
% mbari_id = 'ua9766';
% mbari_id = 'ua12396';
% mbari_id = 'ua12382';
% mbari_id = 'ua19751';
%mbari_id = 'ua18013';
% mbari_id = 'ua21827';
%mbari_id = 'ua21827';
%mbari_id = 'ua20704';
% mbari_id = 'ua20602';
%mbari_id = 'ua20623';
% mbari_id = 'ua21286';
% mbari_id = 'ua11090';
% mbari_id = 'ua21075';
% mbari_id = 'ua21900';

% mbari_id = 'ua19085';
% mbari_id = 'ua17328';
%mbari_id = 'ua21792';
%mbari_id = 'ua21951';
%mbari_id = 'ua19302';
%mbari_id = 'wn1345';
%mbari_id = 'ua21519';
%mbari_id = 'ua20265';
%mbari_id = 'ua22496';
%mbari_id = 'wn1562';
%mbari_id = 'wn1559';
% mbari_id = 'ua23599';
% mbari_id = 'wn1558';
% mbari_id = 'ua21216';
% mbari_id = 'ua21026';
%mbari_id = 'ua21865';
%mbari_id = 'ua21986';
%mbari_id = 'ua20352';
%mbari_id = 'un1570';
% mbari_id = 'un1572';
% mbari_id = 'un1522';
% mbari_id = 'ua22958';
%mbari_id = 'ua22912';
%mbari_id = 'ua21546';
%mbari_id = 'ua19531';
mbari_id = 'ua22127';
mbari_id = 'ua19389';
mbari_id = 'ua19976';
mbari_id = 'ua19806';
mbari_id = 'ua22725';
mbari_id = 'ua21986';
mbari_id = 'ua21206'; % cycle 2 wonky I values 1 sample
mbari_id = 'ua21811'; % cycle 5 wonky I values 1 sample
mbari_id = 'ua21861'; % cycle 1 wonky I values 1 sample
mbari_id = 'ua21285'; % cycle 9 wonky I values 1 sample

% mbari_id = 'un1119'; % shifting intesity floats
% mbari_id = 'ua20109'; 
% mbari_id = 'ua20492'; 
% mbari_id = 'ua20602';
% mbari_id = 'ua21075';
% mbari_id = 'ua21519';
% mbari_id = 'ua21900';
% mbari_id = 'ua21996';
% mbari_id = 'ua22725';
% mbari_id = 'un1448';

%mbari_id = 'ua23599'; % high fit error deployed since 1/1/2024
% mbari_id = 'ua20075';
% mbari_id = 'ua21206';
% %mbari_id = 'ua21026'; % weird break in abs spectra at 235
% %mbari_id = 'ua21216'; % hump
% mbari_id = 'ua21986';
% mbari_id = 'ua22958'; % weird break in abs spectra at 234
% mbari_id = 'ua22372'; % weird break in abs spectra at 239-240
% mbari_id = 'wn1529'; % weird break in abs spectra at 239-240
% 
% mbari_id = 'wn1556';
% mbari_id = 'ua21206';
% mbari_id = 'ua21811';
% mbari_id = 'wn1356'
mbari_id = 'ua19364'

cycle    = 40;


norm_wl = 250;
% ************************************************************************
% SET UP DIRECTORIES, FILE NAMES & LOAD FILE LIST
user_dir  = getenv('USERPROFILE'); %returns user path,i.e. 'C:\Users\jplant'
user_dir  = [user_dir,'\Documents\MATLAB\ARGO_PROCESSING\'];
% user_dir = '\\atlas\Chem\ARGO_PROCESSING\';

dirs.data = [user_dir,'DATA\FLOATVIZ\'];
dirs.cal  = [user_dir, 'DATA\CAL\'];
dirs.FV   = [user_dir,'DATA\FLOATVIZ\']; % used to get date of failure once cycle known

dirs.save = [user_dir,'DATA\Survival_stats\'];
%dirs.save = ;


% ************************************************************************
% load cal file
calfn = sprintf('cal%s.mat',mbari_id);
load(fullfile(dirs.cal,calfn));
ncal = cal.N;
cts = (1:size(ncal.WL,1))'; % 1 to 255

% **** Play with ref adjustment ****
% ref_wl_adj = 0.075; % nm
% wl_new  = ncal.WL + ref_wl_adj;
% ref_new = interp1(wl_new, ncal.Ref, ncal.WL);
% ncal.Ref = ref_new;

% **** play with NO3 contam on 19142 ****
% sal_tmp = 32.76;
% tmp = (ncal.ESW * sal_tmp - ncal.ENO3 * 2) / sal_tmp; % remove 2 uM NO3 contam
% ncal.ESW =tmp *0.98;

% ncal.ESW = ncal.ESW * 0.95;

% ************************************************************************
% Load master list
list_fn   = 'MBARI_float_list.mat';
load(fullfile(dirs.cal, list_fn));

mlist     = d.list;
mlist_hdr = d.hdr;
rlist    = size(mlist,1);

I.MSG  = find(strcmp(mlist_hdr,'msg dir') == 1);
I.INST = find(strcmp(mlist_hdr,'INST ID') == 1);


tg      = strcmp(mlist(:,1),mbari_id);
msg_dir = mlist{tg,I.MSG};
ID      = mlist{tg,I.INST};
fn      = sprintf('%s.%03.0f.isus',ID,cycle);
fp1      = [msg_dir,fn];
fp2     = regexprep(fp1,'UW\\','UW\\alternate\\');   

if isfile(fp1)
    fp = fp1;
    spec    = parse_NO3msg(fp);
elseif isfile(fp2)
    fp = fp2;
    spec    = parse_NO3msg(fp);
else
    fprintf('% not found', fn)
    return
end

tFIT  = cts >= spec.spectra_pix_range(1) & cts <= spec.spectra_pix_range(2);
wl = ncal.WL(tFIT);
Io = ncal.Ref(tFIT);
UVI = spec.UV_INTEN;
UVI_DC = UVI - spec.DC;
abs = -log10(UVI_DC ./ Io');

tWL = find(wl <= norm_wl, 1,'last');
wl_str = sprintf('abs@%0.1f',wl(tWL));
norm_abs = abs - abs(:,tWL ); % quick and dirty normalization
%norm_abs = abs - abs(:,end); % quick and dirty normalization
%norm_abs = abs - abs(:,end-10); % quick and dirty normalization

% Remove SW abs & normalize
% UPDATED ESW TEMPERATURE CORRECTION ALGORITHM TO dLN[ESW]/dT 03/30/2022 JP
% (03/29/2022 JP)EQUITECH PROBES ONLY with 3 anamolously low float experiments
% removed 05/20/2016, 11/21/2016, 10/09/2017 (P CHAMBER EXPERIMENTS?)
Tcorr_coef  = [1.27353e-07 -7.56395e-06 2.91898e-05 1.67660e-03 1.46380e-02];
% Tcorr       = polyval(Tcorr_coef,(WL - WL_offset)) .* (ctd_temp - cal_temp);
% ESW_in_situ = ESW .* exp(Tcorr);
Tcorr           = (spec.T - ncal.CalTemp) * polyval(Tcorr_coef,(wl - 210)'); % matrix mutliply
ESW             = ones(size(spec.T)) * ncal.ESW(tFIT)';
ESW_abs         = ESW .* exp(Tcorr) .* spec.S;
NO3_BL_abs      = abs - ESW_abs;
norm_NO3_BL_abs = NO3_BL_abs - NO3_BL_abs(:,tWL);
%norm_NO3_BL_abs = NO3_BL_abs - NO3_BL_abs(:,end);
%norm_NO3_BL_abs = NO3_BL_abs - NO3_BL_abs(:,end-10);


% build a cmap matrix yp match sample count for the given cycle
cmap      = colormap; % get colormap & close fig the call generates
cmap_ct   = 1:size(cmap,1); % count array, size of colormao, nx3
f = gcf;
delete(f);


% make count array size of sample count, max value size of cmap
sample_ct = linspace(1, size(cmap,1), size(UVI,1)); 

% equal distant subset of cmp with size of sample count
cmapi     = interp1(cmap_ct, cmap, sample_ct); 


F10 = figure(10);
F10.Units = 'Normalized';
F10.Position = [0.3177 0.2618 0.4109 0.5507];
x_lim = [217 250];

ax(1) = subplot(2,2,1);
p1 = plot(wl,UVI');
xlabel('Wavelength')
ylabel('Intensity')
title(fn)
xlim(x_lim)
hold on
pREF = plot(ncal.WL,ncal.Ref,'k-','LineWidth',2);
hold off
set(p1, {'color'}, num2cell(cmapi,2)); % set color of all plots held in p1
legend(pREF,'Io','Location','SouthEast');
ax(1).FontSize = 14;


ax(2) = subplot(2,2,2);
p2 =plot(wl,norm_abs');
xlabel('Wavelength')
ylabel(['Abs - ',wl_str])
%title(fn)
xlim(x_lim)
set(p2, {'color'}, num2cell(cmapi,2));
colormap(ax(2),cmapi); % make interp colormap the color map for the colorbar
cb            = colorbar;
cb.Location   = 'North';
cb.Direction  = 'Reverse';
cb.Ticks      = [0,1];
cb.TickLabels = fix(spec.P([1,end]));
cb.Position   = cb.Position + [0.04 0 -0.05 0];
cb.Label.String = 'Sample pressure';  
ax(2).FontSize = 14;


ax(3) = subplot(2,2,3);
P = spec.P * ones(1,size(UVI,2));
p3 = plot(UVI',P','o');
xlabel('Intensity')
ylabel('Pressure')
%title(fn)
ax(3).YDir = 'Reverse';
set(p3, {'color'}, num2cell(cmapi,2));
ax(3).FontSize = 14;


ax(4) = subplot(2,2,4);
p4 = plot(wl,norm_NO3_BL_abs');
xlabel('Wavelength')
ylabel(['NO3+BL Abs - NO3+BL ',wl_str])
%title(fn)
xlim(x_lim)
set(p4, {'color'}, num2cell(cmapi,2));
ax(4).FontSize = 14;


% ax(3) = subplot(2,2,3);
% P = spec.P * ones(1,size(UVI,2));
% plot(UVI',P','o')
% xlabel('Intensity')
% ylabel('Pressure')
% title(fn)
% ax(3).YDir = 'Reverse';


