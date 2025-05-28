function [noisy,epochd,epoch_og] = noiseGen(time,ID)

    fid = fopen([fullfile('TLEs',int2str(ID)),'.txt'],'rt');
    tline = fgetl(fid);
    tline = fgetl(fid);
    epoch = str2double(tline(19:32));

    epoch_og = epoch;

    epochd = epoch+rand/4000;
    fprintf('-------------------------NoiseGen Data-------------------------\n')
    fprintf('Epoch diff = %s, noisy ->%10.10f, OG ->%10.10f\n', (epochToUTC(epochd)-epochToUTC(epoch)),epochd,epoch)

    f_doppler = paramsToDop(epochd,time);

    cte = 1e3;

    mean = 0;
    stdDeviation = 0.2;
    l = length(f_doppler);
    n = randn(1,l)*stdDeviation+mean;

    mean = 2;
    stdDeviation = 0.5;
    n2 = randn(1)*stdDeviation+mean;

    fprintf('Frequency deviation -> %5.4fHz \n',n2*100)
    fprintf('---------------------------------------------------------------\n')

    noisy = f_doppler+(n*cte)+(n2*100);
end