clc
clear

freq = 143.050e6;
inst = datetime('now', 'TimeZone', 'Local'); %Time origin
inst = datetime('13-May 00:30:22','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 30; % In minutes
precision = 1 / 60; % Precision in minutes

fitter = true;
fitterType = 'real'; % real for real data or sim for generated noise

propagateB4andafter = true; % False propaga hacia delante desde inst
azelFilter = true;

[time, DTtime] = initTimes(inst, duracion, precision,propagateB4andafter);

graphs(time, DTtime, fitter, fitterType, inst, azelFilter, freq)