
spectrogram_mat = zeros(size(IQ,1),4*size(IQ,2));

fs = 5e6;

faxis = linspace(-fs/2,fs/2,4*size(IQ,2));
for ii = 1:size(IQ,1)
    spectrogram_mat(ii,:) = 10*log10(abs(fftshift(fft(IQ(ii,:),4*size(IQ,2)))));

    figure(1)
    plot(faxis,spectrogram_mat(ii,:))
    xlim([-10e3 10e3])
    ylim([40 70])

end