clear
clc
close all

rPipe = 0.0127;         % Pipe radius [m], should not exceed 1m
lPipe = 0.2;            % Pipe test length [m], should not exceed 1m
radialRes = 32;    % radial resolution, should not exceed 9999
lengthRes = 16;    % length resolution, should not exceed 9999
thetaRes = 16;       % theta resolution (ideally multiple of 2), should not exceed 9999

rin = 0.049;         % Inner radius of magnet ring
rout = 0.05;        % Outer radius of magnet ring
h = 0.03;            % Depth of magnet ring

pipeWall = 0.003;   % Thickness of pipe walls [m], mm resolution
airGap = 0.002;     % Space between magnet and pipe [m], mm resolution

%for airGap = [0.001 0.002 0.003 0.004 0.005 0.007 0.010 0.015 0.030 0.050]
for airGap = [0.001 0.002]
    
    % Distance from pipe center to magnet center, mm res--THIS SHOULD NOT EXCEED 1m!!!
    rFlowCenter = rout+pipeWall+airGap+rPipe;
    
    %rFlowCenter = 0.060;    % Distance from pipe center to magnet center, mm res--THIS SHOULD NOT EXCEED 1m!!!
    
    pipeRMeshCyl = linspace(0,rPipe,radialRes);                    % Create radial points for pipe mesh
    pipeZMeshCyl = linspace(-lPipe/2,lPipe/2,lengthRes);           % Create lengthwise points for pipe mesh
    pipeThetaMeshCyl = linspace(0,2*pi-2*pi/thetaRes,thetaRes);    % Create theta points for pipe mesh
    
    % Pipe is set up such that in cylindrical coordinates the z-axis is the
    % same as the x-axis in cartesian. Pipe flow travels in the positive
    % x-direction, while the vertical is given as the y-direction, and
    % correspondingly the horizontal-transverse is given as the z-direction.
    % The mesh is laid out corersponding to (r,theta,z).
    
    pipeXMesh = ones(radialRes,thetaRes,lengthRes);
    pipeYMesh = pipeXMesh;
    pipeZMesh = pipeXMesh;
    
    % Create x-coordinates in cartesian
    for zIter = 1:lengthRes
        pipeXMesh(:,:,zIter) = pipeZMeshCyl(zIter)*pipeXMesh(:,:,zIter);
    end
    
    % Create y-coordinates in cartesian
    for thetaIter = 1:thetaRes
        for rIter = 1:radialRes
            pipeYMesh(rIter,thetaIter,:) = pipeRMeshCyl(rIter)*sin(pipeThetaMeshCyl(thetaIter));
        end
    end
    
    % Create z-coordinates in cartesian
    for thetaIter = 1:thetaRes
        for rIter = 1:radialRes
            pipeZMesh(rIter,thetaIter,:) = pipeRMeshCyl(rIter)*cos(pipeThetaMeshCyl(thetaIter));
        end
    end
    
    % Un-comment the following code to plot out the mesh
    
    % figure()
    % hold on
    % for rIter = 1:radialRes
    %     for thetaIter = 1:thetaRes
    %         for zIter = 1:lengthRes
    %             plot3(pipeXMesh(rIter,thetaIter,zIter),pipeYMesh(rIter,thetaIter,zIter),pipeZMesh(rIter,thetaIter,zIter),'bo')
    %         end
    %     end
    % end
    % xlabel('X - flow direction')
    % ylabel('Y - height')
    % zlabel('Z - span')
    % axis([-lPipe/2 lPipe/2 -rPipe rPipe -rPipe rPipe])
    % axis equal
    % hold off
    
    % This is the end of the code that plots out the mesh
    
    % Create r/z values for ring magnet--distance from center of magnet to the
    % pipe center is a necessary parameter, called rFlowEdge
    
    magRMesh = ones(radialRes,thetaRes,lengthRes);
    magZMesh = pipeYMesh;
    
    % Create r-coordinates for ring magnet and angles for each point
    thetaPrime = zeros(radialRes,thetaRes,lengthRes);
    for rIter = 1:radialRes
        for thetaIter = 1:thetaRes
            for zIter = 1:lengthRes
                magRMesh(rIter,thetaIter,zIter) = sqrt(...
                    (rFlowCenter+pipeZMesh(rIter,thetaIter,zIter))^2+pipeXMesh(rIter,thetaIter,zIter)^2);
                thetaPrime(rIter,thetaIter,zIter) = atan(...
                    pipeXMesh(rIter,thetaIter,zIter)/(pipeZMesh(rIter,thetaIter,zIter)+rFlowCenter));
            end
        end
    end
    
    %% Calculate magnetic field values on the magnetic field mesh coordinates
    [Hr,Htheta,Hz] = RadialRingMagnetFieldGen(magRMesh,magZMesh,rin,rout,h);
    
    %% Rotate the magnetic field components to cartesian
    pipeXBField = zeros(radialRes,thetaRes,lengthRes);
    pipeYBField = zeros(radialRes,thetaRes,lengthRes);
    pipeZBField = zeros(radialRes,thetaRes,lengthRes);
    for rIter = 1:radialRes
        for thetaIter = 1:thetaRes
            for zIter = 1:lengthRes
                
                pipeXBField(rIter,thetaIter,zIter) = Hr(rIter,thetaIter,zIter)*sin(thetaPrime(rIter,thetaIter,zIter));
                pipeYBField(rIter,thetaIter,zIter) = Hz(rIter,thetaIter,zIter);
                pipeZBField(rIter,thetaIter,zIter) = Hr(rIter,thetaIter,zIter)*cos(thetaPrime(rIter,thetaIter,zIter));
                
            end
        end
    end
    
    % Plot out resulting magnetic field in mesh
    
    % figure()
    % quiver3(pipeXMesh,pipeYMesh,pipeZMesh,pipeXBField,pipeYBField,pipeZBField)
    % xlabel('X')
    % ylabel('Y')
    % zlabel('Z')
    % axis equal
    
    % End plotting of resultant magnetic field
    
    %% Output mesh file with corresponding magnetic field values
    
    % Parameters used for generating the mesh
    parameterList = [
        radialRes
        thetaRes
        lengthRes
        pipeWall
        airGap
        rPipe
        lPipe
        rout
        rin
        h
        ];
    
    % Output the file with the spacing from the magnet center to pipe center as
    % the end of the file name (0p060 would mean 0.060m--goes down to mm)
    
    directoryNameToPrintTo = sprintf(...
        ['flowProfileMeshes/flowProfile_0p%04dairGap0p%04dpipeWall_',...
        '%04drRes%04dthetaRes%04dzRes0p%04drPipe0p%04dlPipe_0p%03dOR0p%03dIR0p%03dh'],...
        airGap*1e4,pipeWall*1e4,radialRes,thetaRes,lengthRes,rPipe*1e4,lPipe*1e4,rout*1e3,rin*1e3,h*1e3);
    mkdir(directoryNameToPrintTo)
    
    fileName = [directoryNameToPrintTo,sprintf('/flowProfile_Parameters.csv')];
    csvwrite(fileName,parameterList);
    
    % Mesh files of pipe coordinates in cartesian
    csvwrite([directoryNameToPrintTo,'/flowProfile_PipeMeshX.csv'],pipeXMesh);
    csvwrite([directoryNameToPrintTo,'/flowProfile_PipeMeshY.csv'],pipeYMesh);
    csvwrite([directoryNameToPrintTo,'/flowProfile_PipeMeshZ.csv'],pipeZMesh);
    
    % Mesh files of pipe magnetic field in cartesian
    csvwrite([directoryNameToPrintTo,'/flowProfile_PipeBFieldX.csv'],pipeXBField);
    csvwrite([directoryNameToPrintTo,'/flowProfile_PipeBFieldY.csv'],pipeYBField);
    csvwrite([directoryNameToPrintTo,'/flowProfile_PipeBFieldZ.csv'],pipeZBField);
    
end







