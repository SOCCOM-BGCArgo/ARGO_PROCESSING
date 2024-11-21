

dirList = {'\\atlas\Chem\ARGO_PROCESSING\DATA\CAL';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_GO-BGC_snapshot_23Jun2024\SOCCOM_GO-BGC_LoResQC_CANYONB_23Jun2024\SOCCOM_GO-BGC_LoResQC_CANYONB_23Jun2024';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_GO-BGC_snapshot_20Dec2023\SOCCOM_GO-BGC_LoResQC_CANYONB_20Dec2023\SOCCOM_GO-BGC_LoResQC_CANYONB_20Dec2023';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_GO-BGC_snapshot_28Aug2023\SOCCOM_GO-BGC_LoResQC_CANYONB_28Aug2023\SOCCOM_GO-BGC_LoResQC_CANYONB_28Aug2023';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_GO-BGC_snapshot_26Apr2023\SOCCOM_GO-BGC_LoResQC_CANYON_26Apr2023\SOCCOM_GO-BGC_LoResQC_CANYON_26Apr2023';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_snapshot_21Dec2022\SOCCOM_LoResQC_CANYON_21Dec2022\SOCCOM_LoResQC_CANYON_21Dec2022';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_snapshot_19May2022\SOCCOM_LoResQC_CANYON_19May2022\SOCCOM_LoResQC_CANYON_19May2022';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_snapshot_21Dec2021\SOCCOM_LoResQC_CANYON_21Dec2021\SOCCOM_LoResQC_CANYON_21Dec2021';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_snapshot_01Sep2021\SOCCOM_LoResQC_CANYON_01Sep2021\SOCCOM_LoResQC_CANYON_01Sep2021';
    '\\atlas\Chem\ARGO_PROCESSING\DATA\FloatVIZ_SNAPSHOTS\SOCCOM_snapshot_05May2021\SOCCOM_LoResQC_CANYON_05May2021'};

dateList = [datetime("now");
    datetime(2024,6,23);
    datetime(2023,12,20);
    datetime(2023,8,28);
    datetime(2023,4,26);
    datetime(2022,12,21);
    datetime(2022,5,19);
    datetime(2021,12,21);
    datetime(2021,9,1);
    datetime(2021,5,5)];

%Create an empty array with columns DATE GOBGC SOCCOM OTHER and a row for each snapshot
profCounts = NaN(length(dirList),3);

for i = 1:length(dirList)
    load([dirList{i} '\MBARI_float_list.mat'])
    iCYC = find(strcmp('max cycle proc',d.hdr) == 1);
    iPROG = find(strcmp('Program',d.hdr) == 1);
    
    gobgcIDX = contains(d.list(:,iPROG),'GO-BGC');
    soccomIDX = contains(d.list(:,iPROG),'SOCCOM');
    otherIDX = ~contains(d.list(:,iPROG),'GO-BGC')&~contains(d.list(:,iPROG),'SOCCOM');
    
    gobgc_profs = sum(cell2mat(d.list(gobgcIDX,iCYC)),'omitmissing');
    soccom_profs = sum(cell2mat(d.list(soccomIDX,iCYC)),'omitmissing');
    other_profs = sum(cell2mat(d.list(otherIDX,iCYC)),'omitmissing');

    profCounts(i,1) = gobgc_profs;
    profCounts(i,2) = soccom_profs;
    profCounts(i,3) = other_profs;
end


fig = figure;

% bar(dateList,profCounts)
% 
% hold on;

plot(dateList,sum(profCounts, 2),'-o','MarkerSize',8,'LineWidth',2,'Color','k', ...
            'MarkerEdgeColor','k','MarkerFaceColor','k')

xticks(flip(dateList));
xtickformat("yyyy-MMM")
ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';

%FONT SIZE CONTROL
ax.FontSize = 24;

ax.YAxis.Exponent = 3;
ax.XTickLabelRotation = 45;
xline(max(xlim));
yline(70000);

title("Total profile counts since May 2021");

%Y LABEL FONTSIZE CONTROL
ylabel("Profiles (in thousands)",'FontSize', 24);
% legend({'GO-BGC', 'SOCCOM', 'OTHER', 'TOTAL'},'Location','bestoutside')