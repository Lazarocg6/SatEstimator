function [noisy] = noiseGen(f_doppler)

    mean = 0;
    stdDeviation = 0.1;
    l = length(f_doppler);
    n = randn(1,l)*stdDeviation+mean;

    noisy = f_doppler+(n*10^6);
end