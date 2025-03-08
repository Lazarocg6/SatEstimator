function [bistaticRange,bistaticVelocity,R1,R2,baselineLength,baselineECEF] = bistaticParams(LatTX,LonTX,elTX,LatRX,LonRX,elRX, ...
    recefSATS,vecefSATS,f,Re)
%BISTATICPARAMS Generate bistatic parameters
%   Detailed explanation goes here

    wgs72 = referenceEllipsoid('wgs72', 'kilometer');
    
    % Compute geodesic distance in km
    baselineLength = distance(LatTX, LonTX, LatRX, LonRX, wgs72);

    ecefTX = lla2ecef([LatTX LonTX elTX],f,Re);
    ecefRX = lla2ecef([LatRX LonRX elRX],f,Re);

    baselineECEF = sqrt((ecefTX(1)-ecefRX(1))^2+(ecefTX(2)-ecefRX(2))^2+(ecefTX(3)-ecefRX(3))^2);

    bistaticVelocity = zeros(size(recefSATS,1),size(recefSATS,3));
    bistaticRange = zeros(size(recefSATS,1),size(recefSATS,3));
    R1 = zeros(size(recefSATS,1),size(recefSATS,3));
    R2 = zeros(size(recefSATS,1),size(recefSATS,3));

    for j = 1:size(recefSATS,1)
        for i=1:size(recefSATS,3)
            R1(j,i) = sqrt((ecefTX(1)-recefSATS(j,1,i))^2+(ecefTX(2)-recefSATS(j,2,i))^2+(ecefTX(3)-recefSATS(j,3,i))^2);
            R2(j,i) = sqrt((ecefRX(1)-recefSATS(j,1,i))^2+(ecefRX(2)-recefSATS(j,2,i))^2+(ecefRX(3)-recefSATS(j,3,i))^2);
            bistaticRange(j,i)= R1(j,i) + R2(j,i) - baselineECEF;
            bistaticVelocity(j,i) = ((((ecefTX(1)-recefSATS(j,1,i))*vecefSATS(j,1,i))+((ecefTX(2)-recefSATS(j,2,i))* ...
                vecefSATS(j,2,i))+((ecefTX(3)-recefSATS(j,3,i))*vecefSATS(j,3,i)))./R1(j,i))+((((ecefRX(1)-recefSATS(j,1,i))* ...
                vecefSATS(j,1,i))+((ecefRX(2)-recefSATS(j,2,i))*vecefSATS(j,2,i))+((ecefRX(3)-recefSATS(j,3,i))*vecefSATS(j,3,i)))./R2(j,i));
        end
    end
end

