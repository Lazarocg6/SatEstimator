clc
clear

rx = comm.SDRuReceiver(...
    'Platform', 'B210', ...
    'SerialNum', '31420EF', ... % Leave blank for auto-detect, or set serial if needed
    'CenterFrequency', 91.3e6, ...
    'Gain', 30, ...
    'MasterClockRate', 20e6, ...
    'DecimationFactor', 200, ...
    'SamplesPerFrame', 4000, ...
    'OutputDataType', 'double');

% Compute sample rate based on decimation
Fs = 20e6 / 200;

% Set up spectrum analyzer
spectrumAnalyzer = spectrumAnalyzer(...
    'SampleRate', Fs, ...
    'Title', 'Live Spectrum - Updated Every 10 Seconds', ...
    'TimeSpanSource', 'Property', ...
    'TimeSpan', 1);  % short span for visual stability

disp('Starting 10-second receiver loop. Press Ctrl+C to stop.');

while true
    numFrames = Fs * 10 / rx.SamplesPerFrame; % Number of frames in 10 sec
    buffer = [];
    
    for i = 1:ceil(numFrames)
        data = rx();
        if ~isempty(data)
            buffer = [buffer; data];
        end
    end

    spectrumAnalyzer(buffer);
end

% Clean up
release(rx);