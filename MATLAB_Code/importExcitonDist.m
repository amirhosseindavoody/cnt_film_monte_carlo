function [numExcitons,Tmax,deltaT,numRegions,xdim,t,dist] = importExcitonDist( resultsFilePath )
% Imports the exciton distribution over time.

file_name = ' '; %#ok<NASGU>
path = ' ';%#ok<NASGU>
if nargin == 1
    [path, file_name, ext] = fileparts(resultsFilePath); %pre-selected file
    path = [path '/'];
    file_name = [file_name ext];
else
    [file_name, path] = uigetfile('.csv','Select exciton distribution file'); %file dialog
    [~,~,ext] = fileparts(file_name);
end

if(~strcmp(ext,'.csv'))
    return; %wrong file type
end

filePath = [path file_name];

% Import header of the file
[numExcitons,Tmax,deltaT,numRegions,xdim] = importExcitonDistHeader(filePath);

%Initialize data matrix
Tsteps = uint64(ceil(Tmax/deltaT));
t = zeros(Tsteps,1);
dist = zeros(Tsteps,numRegions);

% Open file and fast forward to dist data
currFile = fopen(filePath,'r');
fgets(currFile);
fgets(currFile);
fgets(currFile);

%Iterate over all times recording dist information
for i=1:Tsteps
    line = strsplit(strcat(fgets(currFile)),{';' ','});
    t(i) = str2double(line(1));
    dist(i,1:end) = str2double(line(2:end));    
end
    
    

