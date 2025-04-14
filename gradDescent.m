function [optimized] = gradDescent(sat, instante, duracion, precision, ...
            filenameEOP, f, Re, latTX, lonTX, elTX, latRX, lonRX, ...
            elRX, lambda, noisyDop)
%GRADDESCENT Summary of this function goes here
%   Detailed explanation goes here

    epochs = 25104.55075651+linspace(-0.05,0.05,50);
    t = epochToUTC(epochs);
    t.TimeZone = 'Local';

    f_doppler = zeros(length(epochs),(2*duracion/precision)+1);

    for i = 1:length(epochs)
    
        [recef, vecef, rlla_out, vlla_out, tsince] = propagar(sat, ...
            instante, duracion, precision, filenameEOP, f, Re, epochs(i));
    
        %Bistatic parameters
        [bistaticRange, bistaticVelocity, R1, R2, llaDIST, ecefDIST] = bistaticParams(latTX, ...
            lonTX, elTX, latRX, lonRX, elRX, recef * 1000, vecef * 1000, f, Re);
            
        %Doppler
        f_doppler(i,:) = (bistaticVelocity*1000)./-lambda;
        

    end

    error = (f_doppler - noisyDop); 
    mse = (vecnorm(error').^2)'./size(error,2);
    figure(10000)
    plot(t,mse)
    optimized = noisyDop;
end

