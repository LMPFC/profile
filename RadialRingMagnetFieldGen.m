function [HrNew,HthetaNew,HzNew] = RadialRingMagnetFieldGen(r,z,rin,rout,h)

% Based on Babic 2008

%% Constants
mu0 = 4*pi*1e-7;
sigmaStar = 1;      % "surface magnetic pole density"

rSize = size(r,1);
thetaSize = size(r,2);
zSize = size(r,3);


HrNew = zeros(rSize,thetaSize,zSize);
HthetaNew = zeros(rSize,thetaSize,zSize);
HzNew = zeros(rSize,thetaSize,zSize);

% Adjust z-mesh due to magnet going from 0 to h
z = z+h/2;

for rIter = 1:rSize
    for thetaIter = 1:thetaSize
        for zIter = 1:zSize
            
            %% Azimuthal component
            
            % The theta component of magnetic field from a ring magnet is zero.
            
            %% Axial component
            
            HzNew(rIter,thetaIter,zIter) = calcHz(r(rIter,thetaIter,zIter),z(rIter,thetaIter,zIter));  % Using: Babic 2008
            
            %% Radial component
            
            HrNew(rIter,thetaIter,zIter) = calcHr(r(rIter,thetaIter,zIter),z(rIter,thetaIter,zIter));  % Using: Babic 2008
            
        end
    end
end

%% FUNCTIONS

% This is the function according to the paper UNCHANGED
%     function Hz = calcHz(r,z)
%         Hz = sigmaStar/(4*pi*mu0)*(-4*rin*KStarCalc(4*r*rin/(r^2+rin^2-2*r*rin+z^2))/...
%             sqrt(r^2+rin^2-2*r*rin+z^2)+(4*rin*KStarCalc(4*r*rin/(r^2+rin^2-2*r*rin+(z-h)^2)))...
%             /sqrt(r^2+rin^2-2*r*rin+(z-h)^2)); % (21) from Ravaud 2008
%     end
% DETERMINED THAT THERE IS AN EXTRA MINUS SIGN IN THE PAPER IN KSTAR CALCS

% This is the function according to the paper with changes that correct
% what I suspect are mistakes made by the authors... namely, 1 minus
% sign inside the KStarCalc rather than 2 which cancel.
    function KStar = KStarCalc(m)
        KStar = FStarCalc(pi/2,m);  % (22) from Ravaud 2008
    end
    function FStar = FStarCalc(phi,m)
        ellipticIntegral = @(theta) 1./sqrt(1-m*sin(theta).^2);
        FStar = integral(ellipticIntegral,0,phi);   % (23) from Ravud 2008
    end
    function PiStar = PiStarCalc(n,phi,m)
        ellipticIntegral = @(theta) 1./((1-n*sin(theta).^2).*sqrt(1-m*sin(theta).^2));
        PiStar = integral(ellipticIntegral,0,phi);
    end
    function HzNew = calcHz(r,z)     % Relatively confident that this works.
        HzNewPlus = sigmaStar/(2*pi*mu0)*(HzPlusSumCalc(r,z,1)+HzPlusSumCalc(r,z,2));
        HzNewMinus = -sigmaStar/(2*pi*mu0)*(HzMinusSumCalc(r,z,1)+HzMinusSumCalc(r,z,2));
        HzNew = HzNewPlus + HzNewMinus;
    end
    function HzPlusComp = HzPlusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kPlus(1) = sqrt(4*r*rin/((r+rin)^2+t(1)^2));
        kPlus(2) = sqrt(4*r*rin/((r+rin)^2+t(2)^2));
        HzPlusComp = (-1)^(n-1)*kPlus(n)*sqrt(rin/r)*KStarCalc(kPlus(n));
    end
    function HzMinusComp = HzMinusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kMinus(1) = sqrt(4*r*rout/((r+rout)^2+t(1)^2));
        kMinus(2) = sqrt(4*r*rout/((r+rout)^2+t(2)^2));
        HzMinusComp = (-1)^(n-1)*kMinus(n)*sqrt(rout/r)*KStarCalc(kMinus(n));
    end
    function HrNew = calcHr(r,z)
        HrNewPlus = -sigmaStar/(4*pi*mu0)*(HrPlusSumCalc(r,z,1)+HrPlusSumCalc(r,z,2));
        HrNewMinus = sigmaStar/(4*pi*mu0)*(HrMinusSumCalc(r,z,1)+HrMinusSumCalc(r,z,2));
        HrNew = HrNewPlus + HrNewMinus;
    end
    function HrPlusComp = HrPlusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kPlus(1) = sqrt(4*r*rin/((r+rin)^2+t(1)^2));
        kPlus(2) = sqrt(4*r*rin/((r+rin)^2+t(2)^2));
        kMinus(1) = sqrt(4*r*rout/((r+rout)^2+t(1)^2));
        kMinus(2) = sqrt(4*r*rout/((r+rout)^2+t(2)^2));
        hPlus = 4*r*rin/(r+rin)^2;
        hMinus = 4*r*rout/(r+rout)^2;
        epsPlus(1) = asin(sqrt((1-hPlus)/(1-kPlus(1)^2)));
        epsPlus(2) = asin(sqrt((1-hPlus)/(1-kPlus(2)^2)));
        epsMinus(1) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        epsMinus(2) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        
        % "Suitable for regular or singular cases" --doesn't work
%         HrPlusComp = (-1)^(n-1)*sqrt(rin/r)*(2*t(n)*kPlus(n)/(r+rin)*KStarCalc(kPlus(n))+...
%             pi*sqrt(rin/r)*sign(r-rin)*sign(t(n))*(1-Lambda0Calc(epsPlus(n),kPlus(n))));
        
        % Regular Hr
        HrPlusComp = (-1)^(n-1)*t(n)*kPlus(n)/r*sqrt(rin/r)*(KStarCalc(kPlus(n))+...
            (r-rin)/(r+rin)*PiStarCalc(hPlus,pi/2,kPlus(n)));
    end
    function HrMinusComp = HrMinusSumCalc(r,z,n)
        t(1) = z-h;
        t(2) = z;
        kPlus(1) = sqrt(4*r*rin/((r+rin)^2+t(1)^2));
        kPlus(2) = sqrt(4*r*rin/((r+rin)^2+t(2)^2));
        kMinus(1) = sqrt(4*r*rout/((r+rout)^2+t(1)^2));
        kMinus(2) = sqrt(4*r*rout/((r+rout)^2+t(2)^2));
        hPlus = 4*r*rin/(r+rin)^2;
        hMinus = 4*r*rout/(r+rout)^2;
        epsPlus(1) = asin(sqrt((1-hPlus)/(1-kPlus(1)^2)));
        epsPlus(2) = asin(sqrt((1-hPlus)/(1-kPlus(2)^2)));
        epsMinus(1) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        epsMinus(2) = asin(sqrt((1-hMinus)/(1-kMinus(2)^2)));
        
        % "Suitable for regular or singular cases" --doesn't work.
%         HrMinusComp = (-1)^(n-1)*sqrt(rout/r)*(2*t(n)*kMinus(n)/(r+rout)*KStarCalc(kMinus(n))+...
%             pi*sqrt(rout/r)*sign(r-rout)*sign(t(n))*(1-Lambda0Calc(epsMinus(n),kMinus(n))));
        
        % Regular Hr
        HrMinusComp = (-1)^(n-1)*t(n)*kMinus(n)/r*sqrt(rout/r)*(KStarCalc(kMinus(n))+...
            (r-rout)/(r+rout)*PiStarCalc(hMinus,pi/2,kMinus(n)));
        
    end

end




















