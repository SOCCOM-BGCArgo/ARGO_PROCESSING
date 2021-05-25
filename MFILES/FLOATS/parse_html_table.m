function d = parse_html_table(url)
% Parse Sharon's SOCCOM sensor stats table to extract Dead floats

% Tmaurer April 14 2020, updated Jan 20, 2021 to clean up code.
% 03/17/2021 JP - updated to return the parsed stats table as a structure
% with hdr & table data fields
% 03/23/2021 JP - minor tweaks for launch of Sharon's new GOBGC formated
% tables

% ************************************************************************ 
% TESTING
%url  = 'http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html'

% url  = ['https://www3.mbari.org/chemsensor/MBARI_float_table/',...
%         'mbari_float_table.html'];
%url  = 'http://soccom.ucsd.edu/SOCCOM_float_performance.html'
% ************************************************************************

% ************************************************************************
% GET SHARON'S DATA REFERENCE TABLE AS A STRING
% ************************************************************************

big_str    = webread(url); % Sharon's HTML page as big string
html_lines = regexp(big_str,'\n','split')'; % cell array of all html lines 

% FIND HDR LINE
thdr     = ~cellfun(@isempty, regexp(html_lines,'^<table class','once')); % table data lines
hdr_line = html_lines(thdr); % only want 1st header definition line
hdr_search_exp = '(?<=<th>)[\w+\<\>-\s*\.]+(?=</th>)'; % may need to expand char set for other tables
tbl_hdr  = regexp(hdr_line{1}, hdr_search_exp,'match');
tbl_hdr  = regexprep(tbl_hdr,'<br>-',''); % end of hdr ID, replace with nothing
tbl_hdr  = regexprep(tbl_hdr,'<br>',' '); % header break code replaced with space


% GET TABLE DATA ROWS
tg         = ~cellfun(@isempty, regexp(html_lines,'^\s*<td','once')); % table data lines
table_rows = html_lines(tg); % rows should equal number of floats

%predim out put cell array & fill in cell array
% tried to avoid loop for table data but too many patterns
html_tbl = cell(size(table_rows,1), size(tbl_hdr,2));
for ct = 1:size(table_rows,1)
    tbl_row  = table_rows{ct};
    td_start = regexp(tbl_row,'<td.','end');
    td_end   = regexp(tbl_row,'/td');
    for i = 1:size(tbl_hdr,2) % step through table data
        str = tbl_row(td_start(i):td_end(i));
        data = regexp(str,'(?<=>)[\w+-\.:\s\(\),;/\~\&\[\]\!]+(?=<)','match');
        if size(data,2)>1
            tmp = cellfun(@num2str,data,'UniformOutput',false);
            html_tbl{ct,i} = sprintf('%s',tmp{:});
        else
        html_tbl{ct,i} = data{1,1};
        end
    end
end

d.hdr  = tbl_hdr;
d.data = html_tbl;
clearvars -except d

