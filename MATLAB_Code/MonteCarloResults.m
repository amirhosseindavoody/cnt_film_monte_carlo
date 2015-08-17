%% Monte Carlo Results
% Imports the results from file that the Monte Carlo simulation output.

%% Initialize program
clear
close all;
fromTable = true;
%% Build Workspace from files

folder = uigetdir(pwd);

%#ok<*NASGU>
details = '/details.csv'; 
excitonDist = '/excitonDist.csv';
matExt = '.mat';
excitonDistMov = '/excitonDist';
excitonContCountFig = '/ExCountCont.fig';
excitonCountFig = '/ExCountTot.fig';

if(fromTable)
    details = '/details_t.csv';
    excitonDist = '/excitonDist_t.csv';
    matExt = '_t.mat';
    excitonDistMov = '/excitonDist_t';
    excitonContCountFig = '/ExCountCont_t.fig';
    excitonCountFig = '/ExCountTot_t.fig';
end

matFileNameIdx = strfind(folder,'\');
matFileName = folder(matFileNameIdx(end)+1:length(folder));

if(exist([folder '\' matFileName matExt],'file') ~= 2) 

    %import supporting information
    [numExcitons,Tmax,deltaT,segLenMin,numRegions,xdim,minBin,rmax,numBins,lowAng,...
        highAng,numAng,numTSteps,regLenMin]= importDetails([folder details]);
    
    segmentCountPerRegion = ...
        importSegmentCountPerRegion([folder '/segmentCountPerRegion.csv'],1,inf,1,numRegions);

    %import heat map data
    heatMap = importRThetaDist([folder '/heatMap.csv'],1,inf,1,numAng);

    %Import exciton dist. All rows and columns of T and all num regions
    [t, excitonDist]=importExcitonDist([folder excitonDist],1,inf,1,numRegions+1);

    save([folder '\' matFileName matExt]);

else
    
    load([folder '\' matFileName matExt]);
    
end

%% Processing the R and Theta Distribution

rs = linspace(minBin,rmax,numBins); % Possible r values
thetas = linspace(lowAng,highAng,numAng); % Possible theta values

figure;
surf(thetas,rs,heatMap);
xlabel('Theta [rads]','FontSize',20);
ylabel('R [Angstroms]','FontSize',20);
zlabel('Count','FontSize',20);
title('R and Theta Distribution','FontSize',20);
savefig([folder '/heatmap.fig']);

%% Processing the Exciton distribution vs time

extra = xdim - segLenMin*numRegions;
regLen = segLenMin + extra/numRegions;
x = linspace(regLen-(xdim/2), xdim/2, numRegions) - regLen/2;

xx = linspace(x(1),x(end),numRegions*4); % Points to be used for spline interpolation

%Figure for all plots. All get captured for video making
capFig = figure('visible','off'); 

numAve = 500; %number of vectors the results will be averaged over for each plot
numTimeSteps = floor(length(excitonDist)/numAve);

if(length(excitonDist) > numAve)
    endCount= zeros(numTimeSteps,1); %average exciton count in exit contact
    excitonCount = zeros(numTimeSteps,1); %average exciton count in simulation

    currDist = zeros(1,length(xx));
    exDistAve = zeros(1,length(excitonDist(1,:)));
    sparset = zeros(numTimeSteps,1); %[nS] time vector representing time at each average

    %Movie structures
    writerObj = VideoWriter([folder excitonDistMov]);
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
        title(titleString,'FontSize',20);
        axis([x(1) x(end) 0 numExcitons*3]);
        xlabel('Position (nm)','FontSize',20);
        ylabel('Exciton Count','FontSize',20);
        writeVideo(writerObj,getframe(capFig));


    end
    close(capFig); %Close figure
    close(writerObj); %close movie maker

    %Plot the num particles at the output contact vs time
    plot(sparset,endCount);
    title('Exciton count in output contact vs. time','FontSize',20);
    xlabel('time [nS]','FontSize',20);
    ylabel('Exciton count','FontSize',20);
    savefig([folder excitonContCountFig]);
    
    figure;
    plot(sparset,excitonCount);
    title('Total exciton count vs. time','FontSize',20);
    xlabel('time [nS]','FontSize',20);
    ylabel('Exciton count','FontSize',20);
    savefig([folder excitonCountFig]);
else
    disp('Not enough time steps for video making');
end

%% Processing the number segments per region
if(segmentCountPerRegion(1) ~= -1)
    figure;
    bar(x,segmentCountPerRegion,.95,'b');
    title('Number of segments per region in x-direction','FontSize',20);
    ylabel('Number of segments','FontSize',20);
    xlabel('x [Angstroms]','FontSize',20);
    savefig([folder '/SegCountPerRegion.fig']);
end