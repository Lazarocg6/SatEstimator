clc
clear

info = findsdru;

MCR = 5e6; % Master clock rate
fft_size = 2^11;
DF = 170; % Decimation factor
BB_sample_rate = MCR/DF;
gain = 31; % dB

NFB = 1000; % Number of frames in burst
NB = 4; %Number of bursts


f = 143.05e6;
LO_offset = 500e3;

raw_IQ = zeros(NFB*NB,fft_size);
t = raw_IQ;

rx = comm.SDRuReceiver(Platform="B210",SerialNum=info.SerialNum,...
                     SamplesPerFrame=fft_size,MasterClockRate=MCR, ...
                     DecimationFactor=DF,CenterFrequency=f, ...
                     EnableBurstMode=true, NumFramesInBurst=NFB, ...
                     LocalOscillatorOffset=LO_offset,Gain=gain);

init_t = datetime("now","Format",'dd-MMM HH:mm:ss','TimeZone','local');

for i = 1:NFB*NB
    [raw_IQ(i,:),~,overrun,t(i,:)] = rx();
    if overrun == 1
        warning('-----------Warning: Overrun-----------')
    end
end

release(rx);

IQ = fftshift(flip(fft(raw_IQ'))');
