clc
clear
close all
%Instante duracion precision
inst = datetime('now', 'TimeZone', 'UTC'); %Time origin
% inst = datetime('8-Mar 23:30:17','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 60; % In minutes, 'duracion' before and after the current time
precision = 60 / 60; % Precision in minutes

% Satellites
ISS = struct('Name', 'ISS', 'NORAD', 25544);
Starlink = struct('Name', 'Starlink', 'NORAD', 44737);
OneWeb = struct('Name', 'OneWeb', 'NORAD', 48968);

sats = [ISS OneWeb];

% Coordinates

ETSIT = [40.45206046037957, -3.726407299669201, 670];
Vill = [40.6223011985758, -4.010124224894723, 930];

graphs(sats,inst,duracion,precision,Vill)