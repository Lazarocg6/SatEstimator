%--------------------------------------------------------------------------
%                   SGP4 Orbit Propagator (vectorized)
%
% References:
% Hoots, Felix R., and Ronald L. Roehrich. 1980. Models for Propagation of
% NORAD Element Sets. Spacetrack Report #3. U.S. Air Force: Aerospace Defense
% Command.
% 
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 4th edition (2013).
% 
% Last modified:   2025/03/02   Lázaro Cantos García
%--------------------------------------------------------------------------
clc
clear
close all
format long g

f = 1/298.26; % Params WGS72
Re = 6378135; % Params WGS72
rTierra = 6378.0; %Radio medio de la tierra

addpath('SGP4_Vectorized')

%Coordenadas

latETSIT = 40.45206046037957;
lonETSIT = -3.726407299669201;
elETSIT = 670;%altitud media

latVill = 40.6223011985758;
lonVill = -4.010124224894723;
elVill = 930;%altitud media

latGRAVES = 47.34813145826119;
lonGRAVES = 5.51487314868131;
elGRAVES = 180; %altitud media

% Coordenadas que se van a usar
latTX = latGRAVES;
lonTX = lonGRAVES;
elTX = elGRAVES;
latRX = latVill;
lonRX = lonVill;
elRX = elVill;

duracion = 120;% En minutos, 'duracion' minutos antes y despues del momento actual 
precision = 10/60;% En minutos

% Satellites
ISS = struct('Name','ISS','NORAD',25544);
Starlink = struct('Name','Starlink','NORAD',44737);
OneWeb = struct('Name','OneWeb','NORAD',44057);

sats = [ISS Starlink OneWeb];

n_sats = length(sats);

% filenames
filenameEOP = 'EOP-All.txt';  % Specify the file path

% Actualizar archivo EOP
updateEOP(24,filenameEOP);

% Actualizar TLEs
updateTLE(6,sats);

%SGP4
[recef,vecef,rlla_out,vlla_out,tsince,epoch] = propagar(sats,duracion,precision,filenameEOP,f,Re);

rlla = evitarSaltos(rlla_out); 
vlla = evitarSaltos(vlla_out); 

%bistatic parameters
[bistaticRange,R1,R2,llaDIST,ecefDIST] = bistaticParams(latTX,lonTX,elTX,latRX, ...
    lonRX,elRX,recef*1000,f,Re);

%Figuras
leftX = 200;
rightX = 850;
topY = 550;
botY = 50;


%Figure1-----------------------------------------------------------------------------------
fig1 = figure(1);

plot3(squeeze(recef(:,1,:))',squeeze(recef(:,2,:))',squeeze(recef(:,3,:))')
grid on
pbaspect([1 1 1])
legend({sats.Name}, 'Location', 'best')
set(fig1, 'Position', [leftX, topY, 560, 420]);
title('Órbita [Km]')

%Figure2---------------------------------------------------------------------------------
dist = zeros(n_sats, length(tsince));

% Calculate the distance for each satellite and time point
for j = 1:n_sats
    for i = 1:length(tsince)
        dist(j,i) = norm(recef(j,:,i));
    end
end

fig2 = figure(2);
grid on
hold on
plot(tsince', dist' - rTierra) % Plot vectorized data

legend({sats.Name}, 'Location', 'best') % Use the array of names for the legend
hold off
title('Altura sobre el nivel del mar')
ylabel('Magnitud [Km]')
xlabel('Tiempo [min]')
xlim([-duracion duracion])
set(fig2, 'Position', [rightX, topY, 560, 420]);

%Figure3---------------------------------------------------------------------------------
lim = length(rlla)-(length(rlla)-length(recef));

fig3 = figure(3);
geoplot(squeeze(rlla(1,1,1:lim)),squeeze(rlla(1,2,1:lim)),'DisplayName',sats(1).Name)
geolimits([35 50],[-14 14])
hold on
for i = 2:n_sats
    geoplot(squeeze(rlla(i,1,1:lim)),squeeze(rlla(i,2,1:lim)),'DisplayName',sats(i).Name)
end
geoplot(latTX,lonTX,'xr','DisplayName', 'TX')
geoplot(latRX,lonRX,'xr','DisplayName', 'RX')
geobasemap('darkwater'); 

% Customize the data cursor display
for i = 1:n_sats
    geoplot(squeeze(rlla(i,1,floor(lim/2))),squeeze(rlla(i,2,floor(lim/2))),'ok', ...
        'DisplayName',[sats(i).Name ' position'])
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(obj, event) displayTime(obj, event, ...
    squeeze(rlla(i,1,1:lim)), squeeze(rlla(i,2,1:lim)), tsince(i,:)));
end
hold off
legend('Location','southeast')

% Custom function to display time on hover
function output_txt = displayTime(~, event_obj, latitudes, longitudes, times)
    % Get the index of the point clicked
    idx = event_obj.DataIndex;
    
    % Get the corresponding latitude, longitude, and time
    lat = latitudes(idx);
    lon = longitudes(idx);
    timeStr = num2str(times(idx));

    % Output the time, latitude, and longitude as the label
    output_txt = {['Time: ', timeStr], ...
                  ['Latitude: ', num2str(lat)], ...
                  ['Longitude: ', num2str(lon)]};
end
title('Posición actual y propagada')
set(fig3, 'Position', [leftX, botY, 560, 420]);

%Figure4---------------------------------------------------------------------------------
cte = 1000;
fig4 = figure(4);
plot(tsince',bistaticRange'/cte)
grid on
% hold on
% plot(tsince(1,:),R1(1,:)/cte)
% plot(tsince(1,:),R2(1,:)/cte)
% hold off
legend('Location','best')
title('Bistatic range')
ylabel('Magnitude [Km]')
xlabel('Time [min]')
xlim([-duracion duracion])
legend({sats.Name}, 'Location', 'best')

set(fig4, 'Position', [rightX, botY, 560, 420]);