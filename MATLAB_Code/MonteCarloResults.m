%% Monte Carlo Results
% Imports the results from file that the Monte Carlo simulation output.

folder = uigetdir(pwd);

if(exist([folder '.mat'],'file') ~= 2) 

    [numExcitons,Tmax,deltaT,numRegions,xdim,minBin,rmax,numBins,lowAng,highAng,numAng]...
        = importDetails([folder '/details.csv']);

    heatMap = importRThetaDist([folder '/heatMap.csv']);

    [t, excitonDist] = importExcitonDist([folder '/excitonDist.csv']);

    save([folder '.mat']);

else
    
    load([folder '.mat']);
    
end