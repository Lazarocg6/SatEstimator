clc
clear

freq = 143.050e6;
inst = datetime('now', 'TimeZone', 'Local'); %Time origin
inst = datetime('28-May 13:28:26','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
% inst = datetime('24-May 16:05:18','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 60; % In minutes
precision = 10 / 60; % Precision in minutes

RX1 = [40.45206046037957, -3.726407299669201, 670];
RX2 = [37.45206046037957, 14.726407299669201, 670];

fitter = true;
bistat = true; % Para simulaciones (fitter_true & sim)
multistat = true; % Para simulaciones (fitter_true & sim)
fitterType = 'sim'; % real for real data or sim for generated noise

propagateB4andafter = true; % False propaga hacia delante desde inst
azelFilter = true;

[time, DTtime] = initTimes(inst, duracion, precision,propagateB4andafter);

graphs(time, DTtime, fitter, fitterType, inst, azelFilter, multistat, bistat, RX1, RX2);