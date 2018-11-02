clear
clc
close all

% This 11/01/18 version works to do convergence studies and other kinds of
% plots through various settings.

% This version is considered unstable at this time, possible errors from 
% the 11/01/18 may still be present.

dirList = dir('flowProfileMeshes');
numOfFolders = size(dirList,1);

% Setup what the case requirements are
rResReq = [8 16 32 64];
thetaResReq = 16;
lengthResReq = 16;
pipeWallReq = 0.003;
airGapReq = [0.001 0.002 0.004 0.008];
rPipeReq = 0.0127;
lPipeReq = 0.2;
magORReq = 0.05;
magIRReq = 0.049;
magHReq = 0.03;

numOfrResCases = size(rResReq,2);
numOfthetaResCases = size(thetaResReq,2);
numOflengthResCases = size(lengthResReq,2);
numOfPipeWallCases = size(pipeWallReq,2);
numOfairGapCases = size(airGapReq,2);
numOfrPipeCases = size(rPipeReq,2);
numOflPipeCases = size(lPipeReq,2);
numOfmagORCases = size(magORReq,2);
numOfmagIRCases = size(magIRReq,2);
numOfmagHCases = size(magHReq,2);


for fileIter = 3:numOfFolders
    
    caseParams = importdata(['flowProfileMeshes/',dirList(fileIter).name,...
        '/flowProfile_Parameters.csv']);
    
    if (isempty(rResReq) || sum(rResReq == caseParams(1))) &&...
            (isempty(thetaResReq) || sum(thetaResReq == caseParams(2))) &&...
            (isempty(lengthResReq) || sum(lengthResReq == caseParams(3))) && ...
            (isempty(pipeWallReq) || sum(pipeWallReq == caseParams(4))) &&...
            (isempty(airGapReq) || sum(airGapReq ==  caseParams(5))) &&...
            (isempty(rPipeReq) || sum(rPipeReq == caseParams(6))) &&...
            (isempty(lPipeReq) || sum(lPipeReq == caseParams(7))) &&...
            (isempty(magORReq) || sum(magORReq == caseParams(8))) &&...
            (isempty(magIRReq) || sum(magIRReq == caseParams(9))) &&...
            (isempty(magHReq) || sum(magHReq == caseParams(10)))
    
    rResCases(fileIter-2) = caseParams(1);
    thetaResCases(fileIter-2) = caseParams(2);
    lengthResCases(fileIter-2) = caseParams(3);
    pipeWallCases(fileIter-2) = caseParams(4);
    airGapCases(fileIter-2) = caseParams(5);
    rPipeCases(fileIter-2) = caseParams(6);
    lPipeCases(fileIter-2) = caseParams(7);
    magORCases(fileIter-2) = caseParams(8);
    magIRCases(fileIter-2) = caseParams(9);
    magHCases(fileIter-2) = caseParams(10);
    
    casesOfInterest(fileIter-2) = fileIter-2;
    casesOfInterestParams(:,fileIter-2) = caseParams;
    end
    
end



aCases = [2 4 6 8 10];
numOfaCases = size(aCases,2);

casesOfInterest(casesOfInterest==0) = [];   % Necessary

% rResCases(rResCases==0) = [];
% thetaResCases(thetaResCases==0) = [];
% lengthResCases(lengthResCases==0) = [];
% airGapCases(airGapCases==0) = [];
% pipeWallCases(pipeWallCases==0) = [];
% rPipeCases(pipeWallCases==0) = [];
% lPipeCases(lPipeCases==0) =[];
% magORCases(magORCases==0) = [];
% magIRCases(magIRCases==0) = [];
% magHCases(magHCases==0) = [];

numOfSpacingCases = size(airGapCases,2);
for aIter = 1:numOfaCases
    caseIndex = 0;
    for caseIter = casesOfInterest
        caseIndex = caseIndex + 1;
        % Import the .csv files that contain the mesh data, magnetic field data,
        % and case parameters for the problem, created by pipeMeshForRingMagCalc.m
        Bx = importdata(['flowProfileMeshes/',dirList(caseIter+2).name,...
            '/flowProfile_PipeBFieldX.csv']);
        By = importdata(['flowProfileMeshes/',dirList(caseIter+2).name,...
            '/flowProfile_PipeBFieldY.csv']);
        Bz = importdata(['flowProfileMeshes/',dirList(caseIter+2).name,...
            '/flowProfile_PipeBFieldZ.csv']);
        
        xMesh = importdata(['flowProfileMeshes/',dirList(caseIter+2).name,...
            '/flowProfile_PipeMeshX.csv']);
        yMesh = importdata(['flowProfileMeshes/',dirList(caseIter+2).name,...
            '/flowProfile_PipeMeshY.csv']);
        zMesh = importdata(['flowProfileMeshes/',dirList(caseIter+2).name,...
            '/flowProfile_PipeMeshZ.csv']);
        
        % Define the important values for the problem based on the caseParams file
        rPipe = rPipeCases(caseIter);
        lPipe = lPipeCases(caseIter);
        radialRes = rResCases(caseIter);
        thetaRes = thetaResCases(caseIter);
        lengthRes = lengthResCases(caseIter);
        rOut = magORCases(caseIter);                % ring magnet outer radius
        rIn = magIRCases(caseIter);                 % ring magnet inner radius
        h = magHCases(caseIter);
        rFlowCenter = rOut+airGapCases(caseIter)+pipeWallCases(caseIter)+...
            pipeWallCases(caseIter);   % distance from mag center to pipe center
        
        dr = rPipe/(radialRes-1);   % Radial mesh spacing
        dTheta = 2*pi/thetaRes;     % Theta mesh spacing
        dz = lPipe/(lengthRes-1);   % Length mesh spacing
        
        % Re-shape the imported data to be 3D based on the pipe.
        Bx = reshape(Bx,radialRes,thetaRes,lengthRes);
        By = reshape(By,radialRes,thetaRes,lengthRes);
        Bz = reshape(Bz,radialRes,thetaRes,lengthRes);
        
        xMesh = reshape(xMesh,radialRes,thetaRes,lengthRes);
        yMesh = reshape(yMesh,radialRes,thetaRes,lengthRes);
        zMesh = reshape(zMesh,radialRes,thetaRes,lengthRes);
        
        % Test quiver plot to make sure the magnetic field and mesh came out
        % correctly after manipulation
        quiver3(xMesh,yMesh,zMesh,Bx,By,Bz)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        axis equal
        
        % Pipe flow parameter definition, for now assume u0*(1-r/R)^a
        a = aCases(aIter);
        Q = 5.0671e-4;  % Flow rate such that u0 is 1 m/s for 1" OD pipe
        %u0 = 2; % Use this if you want to fix the velocity in center of pipe
        u0 = Q*(a+2)/(2*pi*a*rPipe^2); % Use this for fixed flowrate
        
        % Omega iteration
        numOfOmegaIters = 20;
        omegaGuess = 0;
        omegaMax = 2*u0/rFlowCenter;
        omegaMin = -2*u0/rFlowCenter;
        omegaIter = 0;
        
        while omegaIter < numOfOmegaIters
            for rIter = 1:radialRes
                u = [u0*(1-((rIter-1)/(radialRes-1))^a),0,0];
                r = rPipe*(rIter-1)/(radialRes-1);
                for thetaIter = 1:thetaRes
                    for zIter = 1:lengthRes
                        B = [Bx(rIter,thetaIter,zIter),By(rIter,thetaIter,zIter),...
                            Bz(rIter,thetaIter,zIter)];
                        momentArm = ...
                            [xMesh(rIter,thetaIter,zIter),yMesh(rIter,thetaIter,zIter),...
                            rFlowCenter+zMesh(rIter,thetaIter,zIter)];
                        uMagnet = cross([0,omegaGuess,0],momentArm);
                        uRel = u-uMagnet;
                        j = cross(uRel,B);
                        jxB = cross(j,B);
                        
                        dV = dr * r*dTheta * dz;
                        
                        Fx(rIter,thetaIter,zIter) = jxB(1)*dV;
                        Fy(rIter,thetaIter,zIter) = jxB(2)*dV;
                        Fz(rIter,thetaIter,zIter) = jxB(3)*dV;
                        FMagnet = -[Fx(rIter,thetaIter,zIter),Fy(rIter,thetaIter,zIter),...
                            Fz(rIter,thetaIter,zIter)];
                        
                        MElement = cross(momentArm,FMagnet);
                        
                        Mx(rIter,thetaIter,zIter) = MElement(1);
                        My(rIter,thetaIter,zIter) = MElement(2);
                        Mz(rIter,thetaIter,zIter) = MElement(3);
                        
                    end
                end
            end
            
            MyTotal = sum(sum(sum(My)));
            FxTotal = sum(sum(sum(Fx)));
            
            if MyTotal > 0
                omegaGuessNew = (omegaGuess+omegaMax)/2;
                omegaMin = omegaGuess;
                % disp('greater')
            elseif MyTotal < 0
                omegaGuessNew = (omegaGuess+omegaMin)/2;
                omegaMax = omegaGuess;
                % disp('less')
            else
                omegaIter = numOfOmegaIters;
                % disp('equal')
            end
            
            omegaIter = omegaIter + 1;
            
%             fprintf('Iteration: %d\n',omegaIter)
%             fprintf('y-Moment: %d\n',MyTotal)
%             fprintf('omega value: %d\n',omegaGuess)
%             fprintf('\n');
            omegaGuess = omegaGuessNew;
        end
        omegaSol(aIter,caseIndex) = omegaGuess;
        FxTotalSol(aIter,caseIndex) = FxTotal;
    end
end
rPipeCases(rPipeCases==0) = [];
airGapCases(airGapCases==0) = [];

% The vector fixing here was added recently--this may be buggy.
rResCases(rResCases==0) = [];
thetaResCases(thetaResCases==0) = [];
lengthResCases(lengthResCases==0) = [];
pipeWallCases(pipeWallCases==0) = [];
lPipeCases(lPipeCases==0) =[];
magORCases(magORCases==0) = [];
magIRCases(magIRCases==0) = [];
magHCases(magHCases==0) = [];
casesOfInterestParams(casesOfInterestParams==0) = [];
numOfCasesOfInterest = size(casesOfInterest,2);
casesOfInterestParams = reshape(casesOfInterestParams,10,numOfCasesOfInterest);

rFlowCenterCases = airGapCases+rOut+pipeWallReq+rPipeReq;

% Choose the x-axis of the plot by setting the desired x-axis value to 1.
% The value of a is handled separately

% NOTE: The choice of the variables to iterate over near the start of the
% code should(?) be the variables that are chosen for x axis and legend.

xAxisParameter = [
    0   % radial resolution
    0   % theta resolution
    0   % length resolution
    0   % pipe wall thickness
    1   % air gap size
    0   % pipe radius
    0   % pipe length
    0   % ring magnet outer radius
    0   % ring magnet inner radius
    0   % ring magnet height/depth
    ];

xAxisLabel = {
    'radial resolution'
    'theta resolution'
    'length resolution'
    'pipe wall thickness [m]'
    'air gap size [m]'
    'pipe radius [m]'
    'pipe length [m]'
    'ring magnet OR [m]'
    'ring magnet IR [m]'
    'ring magnet height [m]'
    };

% Parameter to be iterated on in the legend--a is handled separately.
legendParameter = [
    1   % radial resolution
    0   % theta resolution
    0   % length resolution
    0   % pipe wall thickness
    0   % air gap size
    0   % pipe radius
    0   % pipe length
    0   % ring magnet outer radius
    0   % ring magnet inner radius
    0   % ring magnet height/depth
    ];

legendLabel = {
    'radial resolution'
    'theta resolution'
    'length resolution'
    'pipe wall thickness'
    'air gap size'
    'pipe radius'
    'pipe length'
    'ring magnet OR'
    'ring magnet IR'
    'ring magnet height'
    };

numOfParameterInstances = [
    numOfrResCases
    numOfthetaResCases
    numOflengthResCases
    numOfPipeWallCases
    numOfairGapCases
    numOfrPipeCases
    numOflPipeCases
    numOfmagORCases
    numOfmagIRCases
    numOfmagHCases
    ];
    

copyOfNumOfParameterInstances = numOfParameterInstances;
copyOfNumOfParameterInstances(xAxisParameter==0) = [];
xAxisMaxIter = copyOfNumOfParameterInstances;
copyOfNumOfParameterInstances = numOfParameterInstances;
copyOfNumOfParameterInstances(legendParameter==0) = [];
legendMaxIter = copyOfNumOfParameterInstances;

% Sort through the data

for dataIter = 1:numOfCasesOfInterest
    
    xVarInCaseOrder(dataIter) = sum(xAxisParameter.*casesOfInterestParams(:,dataIter));
    legendVarInCaseOrder(dataIter) = sum(legendParameter.*casesOfInterestParams(:,dataIter));
    
end
% This is where the y axis is decided
yValuesOfInterest = omegaSol(1,:);
yAxisOfInterestLabel = 'omega';

xDataSortIter = zeros(1,legendMaxIter);
xDataSortIter(1,1) = 1;
legendDataSortIter = 1;
previousLegendVal = legendVarInCaseOrder(1);
previousXVal = xVarInCaseOrder(1);
xToPlot(1,1) = previousXVal;
yToPlot(1,1) = yValuesOfInterest(1);

numOfLegendValsEncountered = 1;
orderedLegendVals = legendVarInCaseOrder(1);
for iter = 2:numOfCasesOfInterest
    if sum(orderedLegendVals==legendVarInCaseOrder(iter)) == 0
        numOfLegendValsEncountered = numOfLegendValsEncountered + 1;
        orderedLegendVals(numOfLegendValsEncountered) = legendVarInCaseOrder(iter);
    end
    
end

% sort the x and y data of interest to plot
for caseIter = 2:numOfCasesOfInterest
    if previousLegendVal==legendVarInCaseOrder(caseIter)
        xDataSortIter(legendDataSortIter) = xDataSortIter(legendDataSortIter) + 1;
        xToPlot(xDataSortIter(legendDataSortIter),legendDataSortIter) =...
            xVarInCaseOrder(caseIter);
        yToPlot(xDataSortIter(legendDataSortIter),legendDataSortIter) =...
            yValuesOfInterest(caseIter);

    else

        previousLegendVal = legendVarInCaseOrder(caseIter);
        legendDataSortIter = sum( (1:legendMaxIter).*(previousLegendVal==orderedLegendVals) );
        xDataSortIter(legendDataSortIter) = xDataSortIter(legendDataSortIter) + 1;        
        xToPlot(xDataSortIter(legendDataSortIter),legendDataSortIter) = xVarInCaseOrder(caseIter);
        yToPlot(xDataSortIter(legendDataSortIter),legendDataSortIter) = yValuesOfInterest(caseIter);
    end
end   

figure()
hold on
legendEntries = {};

for legendIter = 1:legendMaxIter
    plot(xToPlot(:,legendIter),yToPlot(:,legendIter),'-o','LineWidth',1.5)
    
    legendEntry = [char(legendLabel(sum((1:10)'.*legendParameter))),...
        sprintf(' = %d',orderedLegendVals(legendIter))];
    legendEntries = [legendEntries;legendEntry];
end
xlabel(char(xAxisLabel(sum((1:10)'.*xAxisParameter))))
ylabel(yAxisOfInterestLabel)
legend(legendEntries)
set(gca,'FontSize',20)

hold off

figure()
hold on
legendEntries = {};
% Used (rFlowCenterCases-rPipeCases).*omegaSol(aIter,:) for Egemen
for aIter = 1:numOfaCases
    plot(rFlowCenterCases,(rFlowCenterCases-rPipeCases).^1.*(omegaSol(aIter,:)),...
        '-o','LineWidth',1.5)
    legendEntry = sprintf('a = %i',aCases(aIter));
    legendEntries = [legendEntries;legendEntry];
end
xlabel('Distance from magnet center to pipe center')
ylabel('Omega * r')
legend(legendEntries)
set(gca,'FontSize',20)

hold off

figure()
hold on
for aIter = 1:numOfaCases
    u = 2*Q*(aCases(aIter)+2)/(2*pi*aCases(aIter)*rPipe^2)*...
        (1-(linspace(0,rPipe,radialRes).^aCases(aIter)/rPipe^aCases(aIter)));
    plot(linspace(0,rPipe,radialRes),u,'-o','LineWidth',1.5)
end
xlabel('radial point in pipe')
ylabel('flow speed [m/s]')
legend(legendEntries)
set(gca,'FontSize',20)

hold off

figure()
hold on
for aIter = 1:numOfaCases
    plot(rFlowCenterCases,FxTotalSol(aIter,:),'-o','LineWidth',1.5)
end
xlabel('Distance from magnet center to pipe center')
ylabel('Force on the magnet ring')
legend(legendEntries)
set(gca,'FontSize',20)

hold off



figure()
hold on
legendEntries = {};
for aIter = 1:numOfaCases
    plot(rFlowCenterCases,(omegaSol(aIter,:))./(omegaSol(1,:)),'-o','LineWidth',1.5)
    legendEntry = sprintf('a = %i',aCases(aIter));
    legendEntries = [legendEntries;legendEntry];
end
xlabel('Distance from magnet center to pipe center')
ylabel('Omega ratio to a = 2 case')
legend(legendEntries)
set(gca,'FontSize',20)

hold off














