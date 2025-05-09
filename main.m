clc
clear

freq = 143.050e6;
inst = datetime('now', 'TimeZone', 'Local'); %Time origin
inst = datetime('10-May 09:25:22','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 600; % In minutes
precision = 60 / 60; % Precision in minutes

fitter = false;
propagateB4andafter = true; % False propaga hacia delante desde inst
azelFilter = true;

[time, DTtime] = initTimes(inst, duracion, precision,propagateB4andafter);

graphs(time, DTtime, fitter,inst,azelFilter,freq)