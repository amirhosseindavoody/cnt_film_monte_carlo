function [t, excitonDist] = importExcitonDist(filename, startRow, endRow, startCol, endCol)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   EXCITONDIST1 = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   EXCITONDIST1 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   excitonDist1 = importfile('excitonDist1.csv', 1, 100001);
%
%    See also TEXTSCAN.


%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
formatSpec = char.empty(0,2*(endCol-startCol+1));
for i=1:(endCol - startCol+1)
  formatSpec(i*2-1:i*2) =  '%f'; 
end
formatSpec = [formatSpec '%[^\n\r]'];

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
matrix = [dataArray{1:end-1}];

t = matrix(1:end,1);
excitonDist = matrix(1:end,2:end);
