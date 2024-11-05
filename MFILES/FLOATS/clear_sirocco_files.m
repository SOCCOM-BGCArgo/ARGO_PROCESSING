%function clear_sirocco_files(wmo_id)
% CLear FloatViz text files from sirocco to remove QC flag propagation when
% it affects updated float reprocessing

%mbari_id = '7613SOOCN';
%mbari_id = '12559SOOCN';
%mbari_id = '12884SOOCN';
%mbari_id = '18796SOOCN';
%mbari_id = '18013SOOCN';
%wmo_id = '5904677';
%wmo_id = '5906307'; %un1115	1115	
%wmo_id = '4903026'; 
% wmo_id = '5906492'; %wn1359	1359	
% wmo_id = '5906042'; % ua11017 11017
%wmo_id = '5906041'; % ua8482
%wmo_id = '1902459'; %wn1487
%wmo_id = '5906503'; %ua20144
%wmo_id = '1902381'; %ua20144
%wmo_id = '1902381'; %ua20144
%wmo_id = '5904186'; %un1512
%wmo_id = '5906210'; %ua21910
%wmo_id = '5906556'; %ua21910
%wmo_id = '5905991'; 
%wmo_id = '2903472';  %ua21267
%wmo_id = '5906533';  %ua17328 for NO3 reset
% wmo_id = '5906342';  %ua19142 reset before reprocess with NO3 fit window [217 235]
%wmo_id = '5906207';  %ua17534 rremove NaN entries from previous QC
%wmo_id = '5906306'; %wn1483 processed with test msg file
%wmo_id = '1902490'; 
%wmo_id = '7901107';  %JP ua21216 had wrong pH cal making pH bad ( Ihope!)
%wmo_id = '1902371';  %JP ua21026 had wrong pH cal . NO3 cal making pH * NO3 bad ( Ihope!)
% wmo_id = '7902103'; % NG had wrong pH cal files (we hope???)
% wmo_id = '7902102'; % NG had the wrong pH cal files (we hope???) 
%wmo_id = '1902370'; % ss4012
% wmo_id = '1902371'; % ss4013
% wmo_id = '7902103'; % LG 8/27/2024 had wrong pH cal files AGAIN
% wmo_id = '7902102'; % LG 8/27/2024 had wrong pH cal files AGAIN 
% wmo_id = '7902102'; % % JP 09/23//2024 still on BSAL! removed & reprocessed
%wmo_id = '5906448'  %LG 10/14/2024 changed cycles on BSL from 19- to 19-52, pH recovery!
wmo_id = '7902114'

data_dir = '\\sirocco\wwwroot\lobo\Data\FloatVizData\';
FV_dirs  = {'','CANYON_QC\', 'HR\', 'HRQC\','MLR_QC\','QC\'};

% CHECK IF THERE ARE ANY FILES TO REMOVE
file_list = {};
str = ' ';
for ct = 1: size(FV_dirs,2)
    tmp = dir([data_dir, FV_dirs{ct}, wmo_id,'*']);

    if ~isempty(tmp)
        for i = 1:size(tmp,1)
            fp = [tmp(i).folder,'\',tmp(i).name];
            file_list = [file_list; fp];
        end
    end       
end

if ~isempty(file_list)
    disp('The following files will be deleted from SIROCCO:')
    disp(' ')
    disp(file_list)
    disp(' ')
    str = input('Are you sure?','s');
end

if regexpi(str,'^Y', 'once')
    for i = 1: size(file_list,1)
        delete(file_list{i})
    end
end




