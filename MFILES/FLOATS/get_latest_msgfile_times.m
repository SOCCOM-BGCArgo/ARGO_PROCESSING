% get_latest_msgfile_times.m
%
% Short script to loop through all SOCCOM floats and record the time of
% last msg file receipt.  Sharon will use this information to populate a
% new column in the stats table.  Can later write to txt file if that is
% preferred.  Will likely also expand script to produce a similar file for
% GO-BGC.
%
% Tanya Maurer
% MBARI
% 01/21/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rawPRIMdir = '\\atlas\chemwebdata\floats\';
rawALTdir = '\\atlas\chemwebdata\floats\alternate\';
load \\atlas\Chem\ARGO_MBARI2AOML\float_type_lists.mat
LIST = SOCCOM_list;
for i=1:size(SOCCOM_list,1)
    fltnum = SOCCOM_list{i,2}
    fltID = SOCCOM_list{i,1};
    wmoid = SOCCOM_list{i,3};
    flttype = SOCCOM_list{i,5};
    if strcmp(flttype,'NAVIS') == 1
        dirtag = 'n';
    elseif strcmp(flttype,'APEX') == 1
        dirtag = 'f';
    end
    fltprim = [rawPRIMdir,dirtag,fltnum,'\*.msg'];
    fltalt = [rawALTdir,dirtag,fltnum,'\*.msg'];
    [latestfile_prim,latestfiletime_prim] = getlastfile(fltprim);
    [latestfile_alt,latestfiletime_alt] = getlastfile(fltalt);
    latestFTs = [latestfiletime_prim;latestfiletime_alt];
    latestFs = [latestfile_prim;latestfile_alt];
    [Y,I] = max(latestFTs);
    FLTID{i} = fltID;
    WMOID{i} = wmoid;
    if ~isempty(latestFs)
        latestFILE{i} = latestFs(I,:);
        latestTIMEMAT{i} = Y;
        latestTIME{i} = datestr(Y,0);
    else
        latestFILE{i} = 'NA';
        latestTIMEMAT{i} = NaN;
        latestTIME{i} = 'NA';
    end
end
out(:,1) = WMOID;
out(:,2) = FLTID;
out(:,3) = latestFILE;
out(:,4) = latestTIME;
out(:,5) = latestTIMEMAT;
hdr = {'WMO-ID' 'FLT-ID' 'latest-msg-file' 'time-of-receipt' 'time-of-receipt-MATLAB'};
data = out;
outfilename = 'C:\Users\bgcargovm\Documents\MATLAB\ARGO_PROCESSING\DATA\SOCCOM_last_msgfiles_received.mat';
save(outfilename,'hdr','out')
copyfile(outfilename,'\\atlas\ftp\pub\SOCCOM\RawFloatData\SOCCOM_last_msgfiles_received.mat','f');

%Just Save As mat file for now.  Most useful to Sharon?
% % fid = fopen(['SOCCOM_last_msgfiles_received.txt'], 'w');
% % fmt = '';
% % for i = 1:size(hdr,2)-1
% %     fprintf(fid,'%s\t', hdr{i});
% %     if regexp(hdr{i},'^NC|^tf','once')
% %         fmt = [fmt,'%0.0f\t'];
% %     else
% %         fmt = [fmt,'%s\t'];
% %     end
% % end
% % fmt = [fmt,'%0.0f\r\n'];
% % fprintf(fid,'%s\r\n',hdr{i+1});
% % 
% % D = list';
% % fprintf(fid, fmt, D{:}); %
% % fclose(fid);
