%% Monte Carlo Results
% Imports the results from file that the Monte Carlo simulation output.

%% Initialize program
clear
close all;

%% Build Workspace from files

folder = uigetdir(pwd);

if(exist([folder '.mat'],'file') ~= 2) 

    %import supporting information
    [numExcitons,Tmax,deltaT,segLenMin,numRegions,xdim,minBin,rmax,numBins,lowAng,...
        highAng,numAng,numTSteps]= importDetails([folder '/details.csv']);

    %import heat map data
    heatMap = importRThetaDist([folder '/heatMap.csv'],1,inf,1,numAng);

    %Import exciton dist. All rows and columns of T and all num regions
    [t, excitonDist]=importExcitonDist([folder '/excitonDist.csv'],1,inf,1,numRegions+1);

    save([folder '.mat']);

else
    
    load([folder '.mat']);
    
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
capFig = figure('Visible','off'); 

numTimeSteps = floor(length(excitonDist)/10);
%F(numTimeSteps) = struct('cdata',[],'colormap',[]); %Array to capture frames
count= zeros(numTimeSteps,1); 

for i=1:floor(numTimeSteps)
    
    currDist = spline(x,excitonDist(i,:),xx);
    plot(xx,currDist,x,excitonDist(i,:),'*');
    titleString = sprintf('Time: %f pS',t(i)*10^12);
    title(titleString);
    axis([x(1) x(end) 0 numExcitons]);
    xlabel('Position (nm)');
    ylabel('Exciton Count');
    F(i) = getframe(capFig);
    
    %Look at output contact
    count(i) = excitonDist(i,end);
%     yaxis([0,100]);
end
close(capFig); %Close figure

%Plot the num particles at the output contact vs time
plot(t(1:numTimeSteps),count);

vidFig = figure; %figure used for displaying movie
movie(vidFig,F,1)

