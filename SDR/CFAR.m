close all 
clc

medidas_dir = "C:\Users\lazar\Desktop\ETSIT\MUIT\TFM\Medidas\";

addpath(medidas_dir)

% load(medidas_dir+"Procesadas/CSS_070520251050/070520251050IMP.mat")
load(medidas_dir+"130520250035STARLINK33978CREO.mat")
% load("120520252045NOSE.mat")

% % Comentar si no son medidas antiguas -------------------------------------
% 
% IQ = fftshift(flip(fft(IQ'))');
% BB_sample_rate = 5e6/80;
% NB = size(IQ,1);
% NFB = 1;
% init_t = datetime(t,"ConvertFrom","posixtime");
% t = 0;
% % -------------------------------------------------------------------------


IQ_abs = abs(IQ);

flim = BB_sample_rate/2;

%Mapa clutter

orden = 10;

w = 0.02;

cteExp = 0.1;

pesos = [w (1-w).*exp(-cteExp.*(0:orden-2))./(sum(exp(-cteExp.*(0:orden-2))))];
filtered = zeros(size(IQ_abs));
for i = 1:size(IQ_abs,1)

    if i<orden
        filtered(i,:) = 0;
    else
        filtered(i,:) = abs(IQ_abs(i,:)).*w;
        for j = 2:orden
            filtered(i,:) = filtered(i,:)+(abs(IQ_abs(i-j+1,:)).*pesos(j));
        end
    end

end

decluttered = (IQ_abs./filtered);

% Convolutional CA-CFAR

k = 3;
cell = [16,16];
guard = [1,1];

% Create CA-CFAR kernel ---------------------------------------------------
data = IQ_abs;

kernel_size = cell + guard*2 + 1;
cfar_kernel = ones(kernel_size);
cut_idx = ceil(kernel_size / 2);
cfar_kernel(cut_idx(1)-guard(1):cut_idx(1)+guard(1), ...
            cut_idx(2)-guard(2):cut_idx(2)+guard(2)) = 0;
cfar_kernel = cfar_kernel / sum(cfar_kernel(:));  % Normalize

% Apply 2D convolution
noise_estimate = conv2(abs(data), cfar_kernel, 'same');

% Thresholding
detections = abs(data) > (k * noise_estimate);

% Find clusters -----------------------------------------------------------
threshold = 3;

cell = [16,16];
guard = [0,0];

kernel_size = cell + guard*2 + 1;
cluster_kernel = ones(kernel_size);
cut_idx = ceil(kernel_size / 2);

cluster_kernel(cut_idx(1)-guard(1):cut_idx(1)+guard(1), ...
               cut_idx(2)-guard(2):cut_idx(2)+guard(2)) = 0;

% Apply 2D convolution
convolution = conv2(double(detections), cluster_kernel, 'same');  % Use double

% Thresholding
detections_2 = convolution >= threshold;
detections_out = detections & detections_2;


dynamicRange = [22 27]; % Sin lna dipolo
dynamicRange = [31 40]; % Medidas antiguas
dynamicRange = [31 40]; % Con lna dipolo
dynamicRange = [29 33]; % Con lna dipolo 2

figure(1)
imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),10.*log10(abs(IQ_abs)));
xlim([-6e3,6e3])
clim(dynamicRange)
xlabel('Frequency [Hz]')
ylabel('Time')
title('Spectrogram')
ax1 = gca;


% figure(2)
% imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),10.*log10(abs(filtered)));
% ax4 = gca;
% title('Clutter map')
% clim(dynamicRange)


% figure(3)
% imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),10.*log10(abs(decluttered)));
% ax3 = gca;
% clim(dynamicRange)
% title('Decluttered')

figure(4)
imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),detections);
title('CFAR Detections')
ax2 = gca;

figure(5)
imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),10.*log10(noise_estimate))
title('Noise Estimate')
ax5 = gca;

figure(6)
imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),convolution)
title('Clusters')
ax6 = gca;

figure(7)
imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),detections_2);
title('Cluster Detections')
ax7 = gca;

figure(8)
imagesc(linspace(-flim,flim,size(IQ,1)),init_t+seconds(t(:,1)'),detections_out);
title('Detections')
ax8 = gca;

% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8])
linkaxes([ax1 ax2 ax5 ax6 ax7 ax8])

prompt = input('¿Guardar? (y/n): ','s');

if prompt == 'y'

    save([char(datetime('now','Format','ddMMyyyyHHmm')), '.mat'] ...
        ,'IQ', 'BB_sample_rate','init_t','t');

end
