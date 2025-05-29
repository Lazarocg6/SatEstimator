function [noisy,epochd,epoch_og,f_dev] = noiseGen(time,ID,RX1,seed,varargin)

    fid = fopen([fullfile('TLEs',int2str(ID)),'.txt'],'rt');
    tline = fgetl(fid);
    tline = fgetl(fid);
    epoch = str2double(tline(19:32));

    epoch_og = epoch; %Epoch usado para generar el ruido

    rng(seed)

    epochd = epoch+rand/4000; %Epoch desplazado

    fprintf('-------------------------NoiseGen Data-------------------------\n')
    fprintf('Epoch diff = %s, noisy ->%10.10f, OG ->%10.10f\n', (epochToUTC(epochd)-epochToUTC(epoch)),epochd,epoch)

    if length(varargin)>=1
        RX2 = varargin{1};
        f_doppler = paramsToDop(epochd,time,RX1,RX2);
    else
        f_doppler = paramsToDop(epochd,time,RX1);
    end

    cte = 1e3;

    mean = 0;
    stdDeviation = 0.2;
    l = length(f_doppler);
    n = randn(1,l)*stdDeviation+mean;

    rng(seed+1)

    mean = 2;
    stdDeviation = 0.5;
    n2 = randn(1)*stdDeviation+mean;

    f_dev = n2*100;

    fprintf('Frequency deviation -> %5.4fHz \n',f_dev)
    fprintf('---------------------------------------------------------------\n')

    % noisy = f_doppler+(n*cte)+(n2*100);
    noisy = f_doppler+(n*cte)+f_dev;
end