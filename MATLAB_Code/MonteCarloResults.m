%% Monte Carlo Results
% Imports the results from file that the Monte Carlo simulation output.

%% Initialize program
clear
close all;

%% Build Workspace from files

folder = uigetdir(pwd);

matFileNameIdx = strfind(folder,'\');
matFileName = folder(matFileNameIdx(end)+1:length(folder));

if(exist([folder '\' matFileName '.mat'],'file') ~= 2) 

    %import supporting information
    [numExcitons,Tmax,deltaT,segLenMin,numRegions,xdim,minBin,rmax,numBins,lowAng,...
        highAng,numAng,numTSteps,regLenMin]= importDetails([folder '/details.csv']);

    %import heat map data
    heatMap = importRThetaDist([folder '/heatMap.csv'],1,inf,1,numAng);

    %Import exciton dist. All rows and columns of T and all num regions
    [t, excitonDist]=importExcitonDist([folder '/excitonDist.csv'],1,inf,1,numRegions+1);

    save([folder '\' matFileName '.mat']);

else
    
    load([folder '\' matFileName '.mat']);
    
end

%% Processing the R and Theta Distribution

rs = linspace(minBin,rmax,numBins); % Possible r values
thetas = linspace(lowAng,highAng,numAng); % Possible theta values

figure;
surf(thetas,rs,heatMap);

%% Processing the Exciton distribution vs time

extra = xdim - segLenMin*numRegions;
regLen = segLenMin + extra/numRegions;
x = linspace(regLen-(xdim/2), xdim/2, numRegions) - regLen/2;

xx = linspace(x(1),x(end),numRegions*4); % Points to be used for spline interpolation

%Figure for all plots. All get captured for video making
capFig = figure('visible','off'); 

numAve = 1000; %number of vectors the results will be averaged over for each plot
numTimeSteps = floor(length(excitonDist)/numAve);
endCount= zeros(numTimeSteps,1); %average exciton count in exit contact
excitonCount = zeros(numTimeSteps,1); %average exciton count in simulation

currDist = zeros(1,length(xx));
exDistAve = zeros(1,length(excitonDist(1,:)));
sparset = zeros(numTimeSteps,1); %[nS] time vector representing time at each average

%Movie structures
writerObj = VideoWriter([folder '/excitonDist']);
open(writerObj);

for i=1:numTimeSteps
    
    %Summation for average
    for j=1:numAve
        exDistAve = exDistAve + excitonDist((i-1)*numAve+j,:);
        %Look at output contact
        endCount(i) = endCount(i) + excitonDist((i-1)*numAve+j,end);
        excitonCount(i) = excitonCount(i) + sum(excitonDist((i-1)*numAve+j,:));
    end
    %Division for average
    exDistAve = exDistAve/numAve;
    endCount(i) = endCount(i)/numAve;
    excitonCount(i) = excitonCount(i)/numAve;
    
    sparset(i) = t(i*numAve)*10^9;
    
    currDist = spline(x,exDistAve,xx);
    plot(xx,currDist,x,exDistAve,'*');
    titleString = sprintf('Time: %f nS',t(i*numAve)*10^9);
    title(titleString);
    axis([x(1) x(end) 0 numExcitons*3]);
    xlabel('Position (nm)');
    ylabel('Exciton Count');
    writeVideo(writerObj,getframe(capFig));
    

end
close(capFig); %Close figure
close(writerObj); %close movie maker

%Plot the num particles at the output contact vs time
plot(sparset,endCount);
title('Exciton count in output contact vs. time');
xlabel('time [nS]');
ylabel('Exciton count');

figure;
plot(sparset,excitonCount);
title('Total exciton count vs. time');
xlabel('time [nS]');
ylabel('Exciton count');

