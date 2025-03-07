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

f = 1/298.26; % WGS72 Parameters
Re = 6378135; % WGS72 Parameters
rTierra = 6378.0; % Average radius of Earth in km

addpath('SGP4_Vectorized')

% Coordinates

latETSIT = 40.45206046037957;
lonETSIT = -3.726407299669201;
elETSIT = 670; % Average altitude

latVill = 40.6223011985758;
lonVill = -4.010124224894723;
elVill = 930; % Average altitude

latGRAVES = 47.34813145826119;
lonGRAVES = 5.51487314868131;
elGRAVES = 180; % Average altitude

% Coordinates to be used
latTX = latGRAVES;
lonTX = lonGRAVES;
elTX = elGRAVES;
latRX = latVill;
lonRX = lonVill;
elRX = elVill;

duracion = 120; % In minutes, 'duracion' before and after the current time
precision = 60 / 60; % Precision in minutes

% Satellites
ISS = struct('Name', 'ISS', 'NORAD', 25544);
Starlink = struct('Name', 'Starlink', 'NORAD', 44737);
OneWeb = struct('Name', 'OneWeb', 'NORAD', 44057);

sats = [ISS Starlink OneWeb];

n_sats = length(sats);

%EOP.txt filename
filenameEOP = 'EOP-All.txt';  

updateEOP(24, filenameEOP); % Update EOP every 24 hours

updateTLE(6, sats); % Update TLE every 6 hours


% SGP4 propagation
[recef, vecef, rlla_out, vlla_out, tsince, epoch] = propagar(sats, ...
    duracion, precision, filenameEOP, f, Re);

rlla = evitarSaltos(rlla_out); 
vlla = evitarSaltos(vlla_out); 

%Bistatic parameters
[bistaticRange, R1, R2, llaDIST, ecefDIST] = bistaticParams(latTX, ...
    lonTX, elTX, latRX, lonRX, elRX, recef * 1000, f, Re);

% Figures
leftX = 200;
rightX = 850;
topY = 550;
botY = 50;

% Figure 1---------------------------------------------------------------------------------
fig1 = figure(1);

plot3(squeeze(recef(:,1,:))',squeeze(recef(:,2,:))',squeeze(recef(:,3,:))')
grid on
pbaspect([1 1 1])
legend({sats.Name}, 'Location', 'best')
set(fig1, 'Position', [leftX, topY, 560, 420]);
title('Orbit [Km]')

% Figure 2---------------------------------------------------------------------------------
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

legend({sats.Name}, 'Location', 'best')
hold off
title('Altitude above sea level')
ylabel('Magnitude [Km]')
xlabel('Time [min]')
xlim([-duracion duracion])
set(fig2, 'Position', [rightX, topY, 560, 420]);

% Figure 3---------------------------------------------------------------------------------
lim = length(rlla) - (length(rlla) - length(recef));

fig3 = figure(3);
ax = geoaxes; % Create Geographic Axes
hold(ax, 'on') % Enable hold on for the geographic axes
geolimits([35 50], [-14 14])

% Plot satellite trajectories
for i = 1:n_sats
    geoplot(ax, squeeze(rlla(i, 1, 1:lim)), squeeze(rlla(i, 2, 1:lim)), ...
        'DisplayName', sats(i).Name, 'Tag', 'trajectory');
end

% Plot TX and RX positions
geoplot(ax, latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX');
geoplot(ax, latRX, lonRX, 'xr', 'DisplayName', 'RX', 'Tag', 'RX');

% Plot current satellite positions
for i = 1:n_sats
    index = find(tsince(i, :) == 0);
    geoplot(ax, squeeze(rlla(i, 1, index)), squeeze(rlla(i, 2, index)), ...
        'ok','DisplayName', [sats(i).Name ' position'], 'Tag', ...
        'current', 'HandleVisibility', 'off');
end

geobasemap(ax, 'darkwater');
hold off
legend('Location', 'southeast')

% Set up the data cursor mode with a custom function
dcm = datacursormode(fig3);
set(dcm, 'UpdateFcn', @(obj, event) displayTime(event, rlla, tsince, ...
    sats, latTX, lonTX, elTX, latRX, lonRX, elRX));

% Custom function to display time on hover (excluding static markers)
function output_txt = displayTime(event_obj, rlla, tsince, sats, latTX, ...
    lonTX, elTX, latRX, lonRX, elRX)

    % Get the clicked point's tag
    targetObj = get(event_obj.Target, 'Tag');
    targetID = find(strcmp({sats.Name}, erase(get(event_obj.Target ...
        , 'DisplayName'), ' position')));

    % Case 1: Trajectory points (show corresponding time)
    if strcmp(targetObj, 'trajectory')
        index = event_obj.DataIndex;
        % Return time, latitude, and longitude
        output_txt = {['Time: ', [num2str(tsince(targetID, index)) ' min']], ...
                      ['Latitude: ', [num2str(rlla(targetID, 1, index))] 'º'], ...
                      ['Longitude: ', [num2str(rlla(targetID, 2, index))] 'º'], ...
                      ['Elevation: ', [num2str(rlla(targetID, 3, index) / 10^3, '%.2f') ' Km']]};

    % Case 2: Current satellite position (show time = 0)
    elseif strcmp(targetObj, 'current')
        index = find(tsince(targetID, :) == 0);
        output_txt = {['Time: ', [num2str(tsince(targetID, index)) ' min']], ...
                      ['Latitude: ', [num2str(rlla(targetID, 1, index))] 'º'], ...
                      ['Longitude: ', [num2str(rlla(targetID, 2, index))] 'º'], ...
                      ['Elevation: ', [num2str(rlla(targetID, 3, index) / 10^3, '%.2f') ' Km']]};

    % Case 3: Static points (TX, RX) - Show lat/lon but hide time
    elseif strcmp(targetObj, 'TX')
        output_txt = {['Latitude: ', [num2str(latTX)] 'º'], ...
                      ['Longitude: ', [num2str(lonTX) 'º']], ...
                      ['Elevation: ', [num2str(elTX)] ' m']};

    elseif strcmp(targetObj, 'RX')
        output_txt = {['Latitude: ', [num2str(latRX)] 'º'], ...
                      ['Longitude: ', [num2str(lonRX)] 'º'], ...
                      ['Elevation: ', [num2str(elRX)] ' m']};

    % Default: Ignore other objects
    else
        output_txt = {};
    end
end

title('Current and propagated position')
set(fig3, 'Position', [leftX, botY, 560, 420]);

% Figure 4---------------------------------------------------------------------------------
cte = 1000;
fig4 = figure(4);
plot(tsince', bistaticRange' / cte)
grid on
% hold on
% plot(tsince(1,:), R1(1,:) / cte)
% plot(tsince(1,:), R2(1,:) / cte)
% hold off
legend('Location', 'best')
title('Bistatic range')
ylabel('Magnitude [Km]')
xlabel('Time [min]')
xlim([-duracion duracion])
legend({sats.Name}, 'Location', 'best')

set(fig4, 'Position', [rightX, botY, 560, 420]);
