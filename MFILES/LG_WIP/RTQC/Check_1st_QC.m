% Check_1st_QC.m 

dirs.QCList = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\QC_LISTS\';
dirs.cal    = '\\atlas\chem\ARGO_PROCESSING\DATA\CAL\';

% ************************************************************************
% LOAD MBARI MASTER LIST & SUBSET TO pH FLOATS ONLY
load([dirs.cal, 'MBARI_float_list.mat'])
fhdr   = d.hdr;
flist  = d.list;

I.MBA  = strcmp(fhdr,'MBARI ID'); % Define lookup indices for float list
I.INST = strcmp(fhdr,'INST ID');
I.WMO  = strcmp(fhdr,'WMO');
I.TYP  = strcmp(fhdr,'float type');
I.DIR  = strcmp(fhdr,'msg dir');
I.PRG  = strcmp(fhdr,'Program');
I.REG  = strcmp(fhdr,'Region');
I.PH   = strcmp(fhdr,'tf pH');
I.CMAX = strcmp(fhdr,'max cycle proc');

% Exclusions from Initial QC Check by MBARI ID
exclude_str = 'un0948|ua19191|ua19298|ua19314|ua19843';
tf_exclude = cellfun(@isempty, regexp(flist(:,I.MBA), exclude_str,'once'));
flist = flist(tf_exclude,:);

rflist = size(flist,1);
clear d



% QC_exception_filter = '';

% ***********************************************************************
% GET LIST OF QC LIST FILES
tmp = dir([dirs.QCList,'*FloatQCList.txt']);
QCfiles = {tmp.name}';

tf_1st_QC = ones(rflist,1)*0;
for i = 1:rflist
    % Is there a QC List file for float
    cycle_max =flist{i,I.CMAX};
    t1 = ~cellfun(@isempty, regexp(QCfiles, flist{i,I.WMO},'once'));
    if sum(t1) == 1
        tf_1st_QC(i) = 1;
    elseif isnan(cycle_max) || cycle_max < 6 % enough cycles to QC
        tf_1st_QC(i) = 1;
    end
end

no_QC_floats = flist(~tf_1st_QC,3);
clearvars -except no_QC_floats