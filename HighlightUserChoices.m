function nMarked = HighlightUserChoices(xlsFile, sheetName, params)
%% Highlight User Choices
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 14-11-2026 21:20 IST
%
% Revision history is not kept here. See HighlightUserChoices_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% What This Function Does
%
% Marks, in the Parameters sheet of the output workbook, every setting that is the user's own choice
% and has no basis in the literature. Such a row has a Source entry that begins with "User choice".
% The whole row is filled amber and set in bold, and the header row is set in bold. The Source
% column itself is written by writetable together with the other columns, so the information
% survives even where the highlighting cannot be applied.
%
% Formatting needs Microsoft Excel, which MATLAB drives through its COM server (Windows only).
% When Excel is not installed, the workbook is open in another program, or any COM call fails, the
% function issues a warning and returns without touching the data. The Excel process is always
% closed before the function returns.
% ____________________________________________________________________________________________________________
% Usage
%   nMarked = HighlightUserChoices(xlsFile, sheetName, params)
%
%   xlsFile    full path of the workbook, already written by writetable
%   sheetName  name of the sheet that holds the Parameters table, 'Parameters'
%   params     the Parameters table that was written, with a variable called Source
%
% Output
%   nMarked    number of rows highlighted, 0 when nothing was marked
% ____________________________________________________________________________________________________________

nMarked = 0;
if ~istable(params)||~ismember('Source', params.Properties.VariableNames)
    warning('HighlightUserChoices:NoSource', 'The Parameters table has no Source column. Nothing was highlighted.');
    return
end
isUser = startsWith(string(params.Source), "User choice");
rowsUser = find(isUser)+1;      % Sheet row, since row 1 holds the headings
nCols = width(params);
if isempty(rowsUser), return; end
amber = 255+256*230+65536*153;  % Excel colour is blue*65536 + green*256 + red
xl = [];
try
    xl = actxserver('Excel.Application');
    xl.Visible = false;
    xl.DisplayAlerts = false;
    wb = xl.Workbooks.Open(xlsFile);
    ws = wb.Sheets.Item(sheetName);
    lastCol = ColumnLetter(nCols);
    ws.Range(sprintf('A1:%s1', lastCol)).Font.Bold = true;
    for q = 1:numel(rowsUser)
        rg = ws.Range(sprintf('A%d:%s%d', rowsUser(q), lastCol, rowsUser(q)));
        rg.Interior.Color = amber;
        rg.Font.Bold = true;
    end
    wb.Save;
    wb.Close(false);
    nMarked = numel(rowsUser);
catch ME
    warning('HighlightUserChoices:ExcelFailed', ...
        'Rows were not highlighted (%s). The Source column of the %s sheet still marks the user choices.', ME.message, sheetName);
end
if ~isempty(xl)
    try, xl.Quit; catch, end
    try, delete(xl); catch, end
end
end

%% Column Letter
function s = ColumnLetter(n)
s = '';
while n > 0
    r = mod(n-1, 26);
    s = [char(65+r) s];
    n = floor((n-1)/26);
end
end
