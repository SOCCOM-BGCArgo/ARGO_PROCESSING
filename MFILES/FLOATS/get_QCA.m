function QCA = get_QCA(qc_path,float_name)
% RETURN QC ADJUSTMENTS FOR A GIVEN PARAMETER
%
% INPUTS:
%   qc_path    - path to qc adjustment list file
%   float_name - mabri float name
%   data_type  - 'NO3' , 'PH', or 'O2'
%
% OUTPUT:
%   QCA       - A STRUCTURE OF ADJUSTMENTS 
%      .offset - for pH only
%      .data   - a matrix [CYLE GAIN OFFSET DRIFT]


% TESTING
% float_name = '0507SoOcn';
% qc_path = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\TEST\Cal_files\FloatQCList.txt';
% float_name = '12380SOOCN';
% qc_path = 'C:\Users\tmaurer\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\QC_LISTS\\12380SOOCN_FloatQCList.txt';

%CHANGES:
% 08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
% 10/30/2017 - get CHL correction if it exists



% predim ouput
QCA.O2  = [];
QCA.NO3 = [];
QCA.PH_OFFSET = [];
QCA.PH  = [];
QCA.CHL = [];

fid = fopen(qc_path);
if fid < 0 %invalid file identifier
    return
end

% FIND SPECIFIC SECTION
tline = '';
while ischar(tline) % get to Section for float
    if regexpi(tline,float_name,'once')
        break
    end
    tline = fgetl(fid);
end

% FIND SPECIFIC DATA TYPE
while ischar(tline) 
    %disp(tline)
    if regexpi(tline,'PREVIOUS','once') 
        break
    end
    % only gain value & ONLY 1 LINE
    if regexp(tline,'^Oxygen','once') 
        tmp = textscan(tline,'%s%f%f%f%f%s','Delimiter', ',');
        QCA.O2  =[QCA.O2;[ tmp{1,2},tmp{1,3},tmp{1,4}, tmp{1,5}]];
    % cycle gain offset, drift    
    elseif regexp(tline,'^Nitrate','once')
        tmp = textscan(tline,'%s%f%f%f%f%s','Delimiter', ',');
        QCA.NO3  =[QCA.NO3; tmp{1,2},tmp{1,3},tmp{1,4}, tmp{1,5}];
    % for pump ON V shift    
    elseif regexp(tline,'pH,\s+offset','once') 
        tmp = textscan(tline,'%s%s%f%s','Delimiter', ',');
        QCA.PH_OFFSET = tmp{3};
    % cycle, offset, drift    
    elseif regexp(tline,'^pH','once') 
        tmp = textscan(tline,'%s%f%f%f%s','Delimiter', ',');
        QCA.PH =[QCA.PH;[tmp{1,2}, 1, tmp{1,3}, tmp{1,4}]];
    elseif regexp(tline,'^CHL','once') 
        tmp = textscan(tline,'%s%f%f%s','Delimiter', ',');
        QCA.CHL =[QCA.CHL;[1, tmp{1,2}, tmp{1,3},0]];
    end
    tline = fgetl(fid);
end
fclose(fid);

