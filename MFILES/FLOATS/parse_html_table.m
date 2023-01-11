function d = parse_html_table(url)
% Parse Sharon's SOCCOM sensor stats table to extract Dead floats

% Tmaurer April 14 2020, updated Jan 20, 2021 to clean up code.
% 03/17/2021 JP - updated to return the parsed stats table as a structure
%     with hdr & table data fields
% 03/23/2021 JP - minor tweaks for launch of Sharon's new GOBGC formated tables
% 10/17/22 JP - updated table header line extraction to account for Sharon's
%    sortable table up date. Minor Modification "hdr_search_exp"
% 01/04/22 JP fix to the 10/17 fix. parsing sortable header data type fix broke 
%    sortable header line parse when no data type. Now test for table type to decide on
%    regexp filter.

% ************************************************************************ 
% TESTING
%url  = 'http://soccom.ucsd.edu/floats/SOCCOM_float_stats.html'
%url  = 'http://soccom.ucsd.edu/SOCCOM_float_performance.html'
%url = 'https://www3.mbari.org/chemsensor/MBARI_float_table/mbari_float_table.html';
% ************************************************************************

% ************************************************************************
% GET SHARON'S DATA REFERENCE TABLE AS A STRING
% ************************************************************************

big_str    = webread(url); % Sharon's HTML page as big string
html_lines = regexp(big_str,'\n','split')'; % cell array of all html lines 

% FIND AND EXTRACT COLUMN ID'sHDR LINE
thdr     = ~cellfun(@isempty, regexp(html_lines,'^<table class','once')); % table data lines
hdr_line = html_lines(thdr); % only want 1st header definition line

% test table type - sortable data type  or just sortable jp 01/04/22
if ~isempty(regexp(hdr_line{1},'<th data', 'once')) % sortable
    hdr_search_exp = '(?<=<th data.*>)[\w+\<\>-\s*\.]+(?=</th>)'; % JP 10/17/22 add ".*" in the look before
else % not sortable
    hdr_search_exp = '(?<=<th>)[\w+\<\>-\s*\.]+(?=</th>)'; % may need to expand char set for other tables
end

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

