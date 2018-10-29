clear
clc
close all

% Adding capabilities as of 10/25/18 12:28:28PM for looping over multiple
% directories created by pipeMeshForRingMagCalc.m

% Made into git repository 10/29/18 ~4PM. Making pointless comment to
% configurate everything..

dirList = dir('flowProfileMeshes');
numOfFolders = size(dirList,1);

for fileIter = 3:numOfFolders
    
    caseParams = importdata(['flowProfileMeshes/',dirList(fileIter).name,...
        '/flowProfile_Parameters.csv']);
    rResCases(fileIter-2) = caseParams(1);
    thetaResCases(fileIter-2) = caseParams(2);
    lengthResCases(fileIter-2) = caseParams(3);
    rFlowCenterCases(fileIter-2) = caseParams(4);
    rPipeCases(fileIter-2) = caseParams(5);
    lPipeCases(fileIter-2) = caseParams(6);
    magORCases(fileIter-2) = caseParams(7);
    magIRCases(fileIter-2) = caseParams(8);
    magHCases(fileIter-2) = caseParams(9);
    
end

numOfSpacingCases = size(rFlowCenterCases,2);

aCases = [2 4 6 8 10];
numOfaCases = size(aCases,2);
for aIter = 1:numOfaCases

    for caseIter = 1:numOfSpacingCases
        
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
        rFlowCenter = rFlowCenterCases(caseIter);   % distance from mag center to pipe center
        
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
        omegaSol(aIter,caseIter) = omegaGuess;
        FxTotalSol(aIter,caseIter) = FxTotal;
    end
end

figure()
hold on
legendEntries = {};
for aIter = 1:numOfaCases
    plot(rFlowCenterCases,rFlowCenterCases.*(omegaSol(aIter,:)),'LineWidth',1.5)
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
    plot(linspace(0,rPipe,radialRes),u,'LineWidth',1.5)
end
xlabel('radial point in pipe')
ylabel('flow speed [m/s]')
legend(legendEntries)
set(gca,'FontSize',20)

hold off

figure()
hold on
for aIter = 1:numOfaCases
    plot(rFlowCenterCases,FxTotalSol(aIter,:),'LineWidth',1.5)
end
xlabel('Distance from magnet center to pipe center')
ylabel('Force on the magnet ring')
legend(legendEntries)
set(gca,'FontSize',20)

hold off

















