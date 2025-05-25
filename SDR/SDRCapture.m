clc
clear rx

info = findsdru;

MCR = 5e6; % Master clock rate
fft_size = 2^11;
DF = 200; % Decimation factor
BB_sample_rate = MCR/DF;
gain = 31; % dB

NFB = 5000; % Number of frames in burst

f = 143.05e6;
LO_offset = 500e3;

raw_IQ = zeros(NFB,fft_size);
t = raw_IQ;

rx = comm.SDRuReceiver(Platform="B210",SerialNum=info.SerialNum,...
                 SamplesPerFrame=fft_size,MasterClockRate=MCR, ...
                 DecimationFactor=DF,CenterFrequency=f, ...
                 EnableBurstMode=true, NumFramesInBurst=NFB, ...
                 LocalOscillatorOffset=LO_offset,Gain=gain);

for i = 1:NFB
    [raw_IQ(i,:),~,overrun,t(i,:)] = rx();
    if i == 1
        init_t = datetime("now","Format",'dd-MMM HH:mm:ss','TimeZone','local')-seconds((fft_size*NFB)/(BB_sample_rate));
    end
    if overrun == 1
        warning('-----------Warning: Overrun-----------')
    end
end

release(rx);


IQ = flip(fftshift(fft(raw_IQ')',2),2);
