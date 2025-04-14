clc
clear
close all
%Instante duracion precision
inst = datetime('now', 'TimeZone', 'UTC'); %Time origin
% inst = inst + minutes(134);
% inst = datetime('14-Apr 15:46:16','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 120; % In minutes, 'duracion' before and after the current time
precision = 20 / 60; % Precision in minutes

% Satellites
ISS = struct('Name', 'ISS', 'NORAD', 25544);
Starlink = struct('Name', 'Starlink', 'NORAD', 63512);
OneWeb = struct('Name', 'OneWeb', 'NORAD', 44057);

sats = [ISS];

AzElFilter = false;

% Coordinates

ETSIT = [40.45206046037957, -3.726407299669201, 670];
Vill = [40.6223011985758, -4.010124224894723, 930];

graphs(sats,inst,duracion,precision,Vill,AzElFilter)