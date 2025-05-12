
spectrogram_mat = zeros(size(IQ,1),4*size(IQ,2));

fs = 5e6;
DCF = 80; 

% faxis = linspace(-(fs/DCF)/2,(fs/DCF)/2,4*size(IQ,2));
for ii = 1:size(IQ,1)
    spectrogram_mat(ii,:) = 10*log10(abs(fftshift(fft(IQ(ii,:),4*size(IQ,2)))));

    figure(1)
    plot(spectrogram_mat(ii,:))
    xlim([6900 7900])
    ylim([30 38])

end