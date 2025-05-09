close all
clear 
clc

medidas_dir = "C:\Users\lazar\Desktop\ETSIT\MUIT\TFM\Medidas\";

addpath("C:\Users\lazar\Desktop\ETSIT\MUIT\TFM\Medidas\")

load(medidas_dir+"Procesadas/CSS_070520251050/070520251050IMP.mat")
% load('cropped.mat')

IQ2 = fftshift(flip(fft(IQ'))');
% IQ2 = cropped;

fs = 5e6;
decimation = 80;
colorLims = [62 75];
% dt = datetime(t, 'ConvertFrom', 'posixtime', 'TimeZone', 'Local','Format','dd-MMM HH:mm:ss.SSS');

% figure
% pcolor(linspace(-fs/(decimation*2)+center_freq,fs/(decimation*2)+center_freq,size(IQ2,2)),dt,20.*log10(abs(IQ2)));
% shading flat;
% clim([62 75])
% set(gca, 'YDir', 'reverse');

figure
imagesc(-size(IQ2,2)/2,0,20.*log10(abs(IQ2)));
ax1 = gca;
clim(colorLims)

% step = fs/size(IQ2,2);
% offsetLO = 100e3;
% margin = 200;
% maxSample = floor(offsetLO/step)-margin;
% 
% IQ2 = IQ2(:,(size(IQ2,2)/2)-maxSample:(size(IQ2,2)/2)+maxSample);

%Mapa clutter

orden = 4;

w = 0.02;

cteExp = 0.1;

pesos = [w (1-w).*exp(-cteExp.*(0:orden-2))./(sum(exp(-cteExp.*(0:orden-2))))];
filtered = zeros(size(IQ2));
for i = 1:size(IQ2,1)

    if i<orden
        filtered(i,:) = 0;
    else
        filtered(i,:) = abs(IQ2(i,:)).*w;
        for j = 2:orden
            filtered(i,:) = filtered(i,:)+(abs(IQ2(i-j+1,:)).*pesos(j));
        end
    end

end

% Convolutional CA-CFAR

k = 3.2;
cell = [64,64];
guard = [4,4];

% Create CA-CFAR kernel
kernel_size = cell + guard*2 + 1;
cfar_kernel = ones(kernel_size);
cut_idx = ceil(kernel_size / 2);
cfar_kernel(cut_idx(1)-guard(1):cut_idx(1)+guard(1), ...
            cut_idx(2)-guard(2):cut_idx(2)+guard(2)) = 0;
cfar_kernel = cfar_kernel / sum(cfar_kernel(:));  % Normalize

% Apply 2D convolution
noise_estimate = conv2(abs(filtered), cfar_kernel, 'same');

% Thresholding
detections = abs(filtered) > (k * noise_estimate);


figure
imagesc(-size(IQ2,2)/2,0,20.*log10(abs(filtered)));
ax4 = gca;
title('Clutter map')


figure
imagesc(-size(IQ2,2)/2,0,20.*log10(abs(IQ2))-20.*log10(abs(filtered)));
ax3 = gca;
% clim(colorLims)
title('Decluttered')

figure
imagesc(-size(IQ2,2)/2,0,detections);
ax2 = gca;

figure
imagesc(-size(IQ2,2)/2,0,noise_estimate)
ax5 = gca;

linkaxes([ax1 ax2 ax3 ax4 ax5])