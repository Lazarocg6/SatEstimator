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

    bistaticVelocity = zeros(1,length(recefSATS));
    bistaticRange = zeros(1,length(recefSATS));
    R1 = zeros(1,length(recefSATS));
    R2 = zeros(1,length(recefSATS));

    for i=1:length(recefSATS)
        R1(i) = sqrt((ecefTX(1)-recefSATS(1,i))^2+(ecefTX(2)-recefSATS(2,i))^2+(ecefTX(3)-recefSATS(3,i))^2);
        R2(i) = sqrt((ecefRX(1)-recefSATS(1,i))^2+(ecefRX(2)-recefSATS(2,i))^2+(ecefRX(3)-recefSATS(3,i))^2);
        bistaticRange(i)= R1(i) + R2(i) - baselineECEF;
        bistaticVelocity(i) = -((((ecefTX(1)-recefSATS(1,i))*vecefSATS(1,i))+((ecefTX(2)-recefSATS(2,i))* ...
            vecefSATS(2,i))+((ecefTX(3)-recefSATS(3,i))*vecefSATS(3,i)))./R1(i))-((((ecefRX(1)-recefSATS(1,i))* ...
            vecefSATS(1,i))+((ecefRX(2)-recefSATS(2,i))*vecefSATS(2,i))+((ecefRX(3)-recefSATS(3,i))*vecefSATS(3,i)))./R2(i));
    end
end

