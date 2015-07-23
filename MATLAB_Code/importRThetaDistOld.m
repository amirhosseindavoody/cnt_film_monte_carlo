function [thetas,rs,hm] = importRThetaDistOld( resultsFilePath )
% Import R & Theta HeatMap files. Plots distribution.

file_name = ' '; %#ok<NASGU>
path = ' ';%#ok<NASGU>
    if nargin == 1
        [path, file_name, ext] = fileparts(resultsFilePath); %pre-selected file
        path = [path '/'];
        file_name = [file_name ext];
    else
        [file_name, path] = ...
            uigetfile('.csv','Select heat map file'); %file dialog
        [~,~,ext] = fileparts(file_name);
    end

    if(~strcmp(ext,'.csv'))
        return; %wrong file type
    end
    
    filePath = [path file_name];
    
    currFile = fopen(filePath,'r');
    % Get r linspace info
    line = strcat(fgets(currFile));
    temp = strsplit(line,{';' ','});
    rlow = str2double(temp(2));
    rhigh = str2double(temp(3));
    rnum = uint64(str2double(temp(4)));
    rs = linspace(rlow,rhigh,rnum);
    
    % Get theta linspace info
    line = strcat(fgets(currFile));
    temp = strsplit(line,{';' ','});
    tlow = str2double(temp(2));
    thigh = str2double(temp(3));
    tnum = uint64(str2double(temp(4)));
    thetas = linspace(tlow,thigh,tnum);
    
    line = fgets(currFile); %#ok<NASGU>
    
    %initialize heatMap matrix
    hm = zeros(length(rs),length(thetas));
    
    for i=1:length(rs)
       line = fgets(currFile);
       temp = str2double(strsplit(line,{';' ','}));
       for j=1:length(thetas)
           hm(i,j) = temp(j);
       end
    end
    
    figure;
    surf(thetas,rs,hm);
    
    
end