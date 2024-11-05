function proc_status = download_ERA5(get_yr, vars)
% ************************************************************************
% download_ERA5.m
% ************************************************************************
% 
% Function to download the specified year(s) of ERA5 data for surface
% pressure and sea level pressure.
% If it's Jan/Feb of a given year, it will download the previous year as
% well, to ensure that the data grabbed is QC'd, and not Real-Time.
%
% GENERAL FUNCTION PROCEDURE:
%   - For each year, write a get_era5.py file.
%   - Run each get_era5.py to download.
%
% USE AS: (1) For download of the current year of ERA5 data, use as:
%             proc_status = download_ERA5('current', {'surface_pressure', 'sp'})
%         (2) For download of multiple years of ERA5 data, 
%             proc_status = download_ERA5('all', {'surface_pressure', 'sp'})
%
% INPUTS:
%    get_yr      - 'current': only download current year
%                   'all'   : download all years since 2014
%                   (ie)[2021, 2022]: single or array of years to download
%    vars        - variable name to grab on Copernicus
%                  First column: name on Copernicus server, 
%                  Second column: name of subdirectory in ncroot for the 
%                  respective parameter, and will be used in the filename 
%                - (ie) {'surface_pressure', 'sp'
%                       'mean_sea_level_pressure', 'msl'}; 
%
% OUTPUTS: 
%    ERA5 data is downloaded to CHEM
%    proc_status = status of job; 1 = success, 0 = failure
%
% SUPPORTING FUNCTIONS:
%    create_py_for_cds (local function)
%   
% AUTHOR:
%   Yui Takeshita
%   MBARI

% UPDATES: TM, 2/9/22   Added error trapping w/ email notification and
%                       modified calls to month() & year() fx.
%          TM, 10/31/22 Modified calls based on recent changes at ECMWF 1)
%                       parameter names have changed slightly...; 2: Avoid errors by
%                       avoiding the 5 day embargo period (otherwise, now generates error)
%          SB, 8/12/24  Turned into a function and automated.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% DEFINE DIRS STRUCTURE____________________________________________________

% root directory for saving era5 data
ncroot = '//atlas/chem/argo_processing/data/era5/';

%% ======================= CONFIGURE EMAIL =================================
 
setpref('Internet','SMTP_Server','mbarimail.mbari.org'); % define server
setpref('Internet','E_mail','tmaurer@mbari.org'); % define sender
% email_list ={'tmaurer@mbari.org';'yui@mbari.org'};
email_list = {'sbartoloni@mbari.org'};
email_msg = {};

%% ======================= Parameters ======================================

% Determine current year. If current month is Jan, grab previous year 
if(month(datetime(datevec(now))) == 1)
    yr = [year(datetime(datevec(now)))-1 year(datetime(datevec(now)))];
else
    yr = year(datetime(datevec(now)));
end

% Determine what years to download.
if strcmp(get_yr, 'all') == 1
    yr = 2014:year(datetime(datevec(now)));
elseif isnumeric(get_yr)
    yr = get_yr;
end

% =========================================================================

%% ========================= Main loop starts ==============================

proc_status = nan(1,length(yr));
    
% loop through the years; If Jan/Feb, need to grab previous year too
for y = 1:length(yr)
        
    % update parameter for specific variable:
    varnow = vars{1,1};
    yrnow = num2str(yr(y));
    path2save = [ncroot,vars{1,2},'/ERA5_',vars{1,2},'_',yrnow,'.nc'];
        
    fprintf('Starting download for %s, %s...',varnow,yrnow);
        
    % write a get_era5.py file used for the cdsapi to get the target variable for a certain
    % year. Specify where you want the file to be saved. Embedded function below.
    create_py_for_cds(varnow,yrnow,path2save);
    fprintf(' finished creating .py file...');
        
    % start downloading ERA5 data by running python through the command prompt
    commandStr = sprintf('python get_era5.py');
    [status, commandOut] = system(commandStr);
        
    if (status ~= 0)
        proc_status(y) = 0;
        disp('WARNING: OPERATION ERA5-DOWNLOAD FAILED!')
        disp(commandOut)
        fprintf('\nFinished (with errors), %s,%s!',varnow,yrnow);
        email_msg = [email_msg;['ERROR downloading ERA ',varnow,' ',yrnow,'!']];
        % remove .py file
        % delete('get_era5.py');
    else
        proc_status(y) = 1;
        fprintf('\nFinished %s,%s!',varnow,yrnow);
        email_msg = [email_msg;['SUCCESS downloading ERA ',varnow,' ',yrnow,'!']];
        % remove .py file
        % delete('get_era5.py');
    end
end

% Send final notification email:
%sendmail(email_list,'ARGOSY: ERA5 DATA DOWNLOAD', email_msg)

return


% Local Function
function [] = create_py_for_cds(variable, yr, path2save)

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

fprintf('Downloading up to %s/%s\n.',num2str(maxmo),num2str(maxday))
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




end
