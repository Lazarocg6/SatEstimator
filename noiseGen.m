function [noisy,epochd,epoch] = noiseGen(time,ID)

    fid = fopen([fullfile('TLEs',int2str(ID)),'.txt'],'rt');
    tline = fgetl(fid);
    tline = fgetl(fid);
    epoch = str2double(tline(19:32));

    epochd = epoch+rand/2000;
    fprintf('Epoch diff = %s, noisy ->%10.10f, OG ->%10.10f\n', (epochToUTC(epochd)-epochToUTC(epoch)),epochd,epoch)

    f_doppler = paramsToDop(epochd,time);

    mean = 0;
    stdDeviation = 0.3;
    l = length(f_doppler);
    n = randn(1,l)*stdDeviation+mean;

    noisy = f_doppler+(n*1e3);
end