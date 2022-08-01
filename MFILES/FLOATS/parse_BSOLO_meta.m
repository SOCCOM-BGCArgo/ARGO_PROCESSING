function INFO = parse_BSOLO_meta(meta_file)


% !!!  THIS IS A START TO A META FILE PARSER  !!!
% meta_file = 'C:\temp\09043_000004.meta';

INFO = [];

% CHECK IF FILE EXIST
if ~isfile(meta_file)
    fprints('BGC SOLO meta file  was not found: %s',meta_file);
    return
end

fid = fopen(meta_file);
tline = ' ';

% PREDIM O2 COEFFICIENT ARRAYS
O2.A  = ones(1,3)*NaN;
O2.B  = ones(1,2)*NaN;
O2.C  = ones(1,3)*NaN;
O2.TA = ones(1,4)*NaN;

while ischar(tline)
    
    % *********************************************************************
    %                             OXYGEN
    % *********************************************************************
    if regexp(tline, 'DOXY A\d')
        ind = str2double(regexp(tline,'(?<=DOXY A)\d','once', 'match')) + 1;
        tmp = regexp(tline,'\s+','split');
        O2.A(ind) = str2double(tmp(end));
    end
    
    if regexp(tline, 'DOXY B\d')
        ind = str2double(regexp(tline,'(?<=DOXY B)\d','once', 'match')) + 1;
        tmp = regexp(tline,'\s+','split');
        O2.B(ind) = str2double(tmp(end));
    end
    
    if regexp(tline, 'DOXY C\d')
        ind = str2double(regexp(tline,'(?<=DOXY C)\d','once', 'match')) + 1;
        tmp = regexp(tline,'\s+','split');
        O2.C(ind) = str2double(tmp(end));
    end
    
    if regexp(tline, 'DOXY TEMP TA\d')
        ind = str2double(regexp(tline,'(?<=DOXY TEMP TA)\d','once', 'match')) + 1;
        tmp = regexp(tline,'\s+','split');
        O2.TA(ind) = str2double(tmp(end));
    end
    
    if regexp(tline, '^oxygen sensor type')
        tmp = regexp(tline,'\s+','split');
        O2.Model = regexp(tmp{end},'\w+','match','once');
        O2.Model = regexprep(O2.Model,'_OPTODE','');
    end
    
    if regexp(tline, '^oxygen sensor serial')
        tmp = regexp(tline,'\s+','split');
        O2.SN = regexp(tmp{end},'\w+','match','once');
    end

    tline = fgetl(fid);
end

if sum(isnan(O2.A),2) == 0 % IF O@ cal coefs add to uotput structure
    INFO.O2 = O2;
end
    


