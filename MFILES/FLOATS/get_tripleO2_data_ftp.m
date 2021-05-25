% get_tripleO2_data_ftp.m
%
% Script to grab all files from SBE server for triple-O2 float 1173
% Currently this script just grabs everything (msg, logs et).  Could/should
% revise at some point to get only new files.
%
% Tanya Maurer
% MBARI
% 11/11/20

FTP_target = '209.181.142.198';
FTP_dir = '/phod/incoming/ARGO_FTP/argo/skip_verify';
psword = 'RPl30I12';   % username and password forwarded from Dana Swift on 11/11/20
username = 'mnicholas20130624';
% TRY CONNECTING TO FTP SERVER
count = 1;
err_count = 1;
disp('')
disp(['Attempting to connect to ',FTP_target,'...',char(10),char(10)])
while count == err_count  
    try
%        f = ftp(FTP_target,username,psword); % Connect to FTP server
		f = ftp(FTP_target,username,psword,'LocalDataConnectionMethod','passive'); % Connect to FTP server
        binary(f)
    catch
        if err_count >= 5
            ftperr = [' Connection to: ',FTP_target,' failed after 5 attempts.',char(10),...
            'No Argo Core files were obtained.'];
            disp(ftperr)
            return
        else
            disp(['Could not connect to ftp server at: ',FTP_target])
            disp('No files were pushed')
            disp('Trying again in 5 min')
            pause(3) %pause for 5 min then try again
            err_count = err_count+1;
        end
    end
    count = count+1;
end
disp(['Connection to ',FTP_target,' successful'])
f
cd(f);
%sf = struct(f);
%sf.jobject.enterLocalPassiveMode();
ftp_path = '1173/';
cd(f,ftp_path)
% listing = dir(f);
% thefiles = char(listing.name);
% thesize = [listing(:).bytes];
% X = ~cellfun('isempty',strfind(cellstr(thefiles),'.msg'));
% files = thefiles(X,:);
% sizes = thesize(X);
% mydir = 'C:\Users\tmaurer\Documents\O2\triple_o2\msgfiles_fromftp\';
% D = dir([mydir,'*.msg']);
% myfiles = char(D.name);
% 
% tf = strfind(files,myfiles)
% cd(f, ftp_path);
targetdir = '\\atlas\Chem\tripleO2float\n1173\';
mget(f,'*',targetdir)
close(f)
MATLABfinish

