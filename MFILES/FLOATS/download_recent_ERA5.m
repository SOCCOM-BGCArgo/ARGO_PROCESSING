%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% download_recent_ERA5.m
%
%
% This script will download the most recent year fo ERA5 data for surface
% pressure, U10, and V10 (10 winds for U and V).
%
% If it's Jan/Feb of a given year, it will download the previous year as
% well, to ensure that the data grabbed is QC'd, and not Real-Time.
%
% Yui Takeshita
% MBARI
%
% UPDATES: TM, 2/9/2022 Added error trapping w/ email notification and
%                       modified calls to month() & year() fx.
%          TM, 10/31/22 Modified calls based on recent changes at ECMWF 1)
%                       parameter names have changed slightly...; 2: Avoid errors by
%                       avoiding the 5 day embargo period (otherwise, now generates error)
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


clear
close all


%% ======================= CONFIGURE EMAIL =================================
 
setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
setpref('Internet','E_mail','tmaurer@mbari.org'); % define sender
%email_list ={'tmaurer@mbari.org';'yui@mbari.org'};
email_list = {'tmaurer@mbari.org'};
email_msg = {};

%% ======================= Parameters ======================================

% variable names to grab. First column: name on Copernicus server, Second
% column: name is the subdirectory name in ncroot for the respective
% parameter, and it will be used in the filename too
%vars2get = {'surface_pressure', 'sfcpres'
%    '10m_v_component_of_wind',  'v10'
%    '10m_u_component_of_wind',  'u10'};
% vars2get = {'surface_pressure', 'sp'
%     'northward_wind',  'v'
%     'eastward_wind',  'u'};
vars2get = {'surface_pressure', 'sp'};

% root directory for saving era5 data
ncroot = '//atlas/chem/argo_processing/data/era5/';
% ncroot = '//atlas/chem/maurer/ERAtest/';

% determine which year to grab. If Jan/Feb, array of last year and this
% year. if not, just thi syear
if(month(datetime(datevec(now))) < 3)
    yr = [year(datetime(datevec(now)))-1 year(datetime(datevec(now)))];
else
    yr = year(datetime(datevec(now)));
end
yr = 2023;
% =========================================================================

%% ========================= Main loop starts ==============================

% loop through the variables to download each .nc
for i = 1:size(vars2get,1)
    
    % loop through the years; If Jan/Feb, need to grab previous year too
    for y = 1:length(yr)
        
        % update parameter for specific variable:
        varnow = vars2get{i,1}; % variable name on Copernicus
        yrnow = num2str(yr(y));
        path2save = [ncroot,vars2get{i,2},'/ERA5_',vars2get{i,2},'_',yrnow,'.nc'];
        
        fprintf('Starting download for %s, %s...',varnow,yrnow);
        
        % write a get_era5.py file used for the cdsapi to get the target variable for a certain
        % year. Specify where you want the file to be saved. Embedded function below.
% % %         create_py_for_cds(varnow,yrnow,path2save);
        fprintf(' finished creating .py file...');
        
% % %         pause(3); % give it some time to save .py. probably not necessary
        
        % start downloading ERA5 data by running python through the command prompt
        commandStr = sprintf('python get_era5.py');
        [status, commandOut] = system(commandStr);
        
        if (status ~= 0)
            disp('WARNING: OPERATION ERA5-DOWNLOAD FAILED!')
            disp(commandOut)
            fprintf('\nFinished (with errors), %s,%s!',varnow,yrnow);
            % commandOut is likley way too verbose, so I think this will suffice:
            email_msg = [email_msg;['ERROR downloading ERA ',varnow,' ',yrnow,'!']];
            % remove .py file
            %delete('get_era5.py');
        else
            fprintf('\nFinished %s,%s!',varnow,yrnow);
            email_msg = [email_msg;['SUCCESS downloading ERA ',varnow,' ',yrnow,'!']];
            % remove .py file
            %delete('get_era5.py');
        end
    end
    
end


% Send final notification email:
sendmail(email_list,'ARGOSY: ERA5 DATA DOWNLOAD', email_msg)

return

function [] = create_py_for_cds(variable, yr, path2save)

% 
% Update (10/28/22) TM: exclude the week leading up to current date to prevent ERA5T grab error
Dnow = datevec(now);

if Dnow(1) == str2num(yr)
    maxmo = Dnow(2);
    maxday = Dnow(3)-8;
    if maxday<=0 %early on in the month; go further back
        if maxmo==1 %if in early Jan...just skip this year.
            return
        end
        maxmo = Dnow(2)-1;
        maxday = 20; %keep it easy...
    end
else
    maxmo = 12;
    maxday = 31;
end

% ==================== hard coded inputs =======================
% set input types
product_type = 'reanalysis';
format = 'netcdf';

% input for month; grab all months
mo = '[';
for i = 1:maxmo
    if(i~=maxmo)
        mo = [mo,sprintf('''%i'',',i)];
    else
        mo = [mo,sprintf('''%i''',i)];
    end
end
mo = [mo,']'];

% input for days; grab all days
d = '[';
for i = 1:maxday
    if(i~=maxday)
        d = [d,sprintf('''%i'',',i)];
    else
        d = [d,sprintf('''%i''',i)];
    end
end
d = [d,']'];

% input for time; grab all time
t = '[''00:00'',''06:00'',''12:00'',''18:00'']';

% headers required for .py
headers = sprintf(['#!/usr/bin/env python\r\n',...
    'import cdsapi\r\n',...
    'c = cdsapi.Client()\r\n',...
    'c.retrieve(\r\n']);

% ====================== start writing to .py ==========================

fid = fopen('get_era5.py','w');

fprintf(fid,headers);
fprintf(fid,'    ''reanalysis-era5-single-levels'',\r\n');
fprintf(fid,'    {\r\n');
fprintf(fid,'        ''product_type'': ''%s'',\r\n', product_type);
fprintf(fid,'        ''variable'': ''%s'',\r\n',variable);
fprintf(fid,'        ''year'': ''%s'',\r\n',yr);
fprintf(fid,'        ''month'': %s,\r\n',mo);
fprintf(fid,'        ''day'': %s,\r\n',d);
fprintf(fid,'        ''time'': %s,\r\n',t);
fprintf(fid,'        ''format'': ''%s'',\r\n',format);
fprintf(fid,'    },\r\n');
fprintf(fid,'    ''%s'')',path2save);

%close file
fclose(fid);



end










