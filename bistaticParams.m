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

    bistaticVelocity = zeros(1,size(recefSATS,2));
    bistaticRange = zeros(1,size(recefSATS,2));
    R1 = zeros(1,size(recefSATS,2));
    R2 = zeros(1,size(recefSATS,2));

    % for i=1:size(recefSATS,2)
    %     R1(i) = sqrt((ecefTX(1)-recefSATS(1,i))^2+(ecefTX(2)-recefSATS(2,i))^2+(ecefTX(3)-recefSATS(3,i))^2);
    %     R2(i) = sqrt((ecefRX(1)-recefSATS(1,i))^2+(ecefRX(2)-recefSATS(2,i))^2+(ecefRX(3)-recefSATS(3,i))^2);
    %     bistaticRange(i)= R1(i) + R2(i) - baselineECEF;
    %     bistaticVelocity(i) = -((((ecefTX(1)-recefSATS(1,i))*vecefSATS(1,i))+((ecefTX(2)-recefSATS(2,i))* ...
    %         vecefSATS(2,i))+((ecefTX(3)-recefSATS(3,i))*vecefSATS(3,i)))./R1(i))-((((ecefRX(1)-recefSATS(1,i))* ...
    %         vecefSATS(1,i))+((ecefRX(2)-recefSATS(2,i))*vecefSATS(2,i))+((ecefRX(3)-recefSATS(3,i))*vecefSATS(3,i)))./R2(i));
    % end

    R1 = sqrt((recefSATS(1,:)-ecefTX(1)).^2 ...
        +(recefSATS(2,:)-ecefTX(2)).^2 ...
        +(recefSATS(3,:)-ecefTX(3)).^2);

    R2 = sqrt((recefSATS(1,:)-ecefRX(1)).^2 ...
        +(recefSATS(2,:)-ecefRX(2)).^2 ...
        +(recefSATS(3,:)-ecefRX(3)).^2);

    bistaticRange= R1 + R2 - baselineECEF;

    bistaticVelocity = ((((recefSATS(1,:)-ecefTX(1)).*vecefSATS(1,:)) ...
    +((recefSATS(2,:)-ecefTX(2)).*vecefSATS(2,:)) ...
    +((recefSATS(3,:)-ecefTX(3)).*vecefSATS(3,:))) ...
    ./sqrt( ...
    (recefSATS(1,:)-ecefTX(1)).^2 ...
    +(recefSATS(2,:)-ecefTX(2)).^2 ...
    +(recefSATS(3,:)-ecefTX(3)).^2) ...
    +...
    (((recefSATS(1,:)-ecefRX(1)).*vecefSATS(1,:)) ...
    +((recefSATS(2,:)-ecefRX(2)).*vecefSATS(2,:)) ...
    +((recefSATS(3,:)-ecefRX(3)).*vecefSATS(3,:))) ...
    ./sqrt( ...
    (recefSATS(1,:)-ecefRX(1)).^2 ...
    +(recefSATS(2,:)-ecefRX(2)).^2 ...
    +(recefSATS(3,:)-ecefRX(3)).^2 ...
    ));

end

