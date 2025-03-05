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

addpath('SGP4_Vectorized')

duracion = 60;% En minutos, 'duracion' minutos antes y despues del momento actual 
precision = 10/60;% En minutos

ISS = 25544;
Starlink = 44737;
OneWeb = 44057;

sats = [ISS Starlink OneWeb];

n_sats = length(sats);

% filenames
filenameEOP = 'EOP-All.txt';  % Specify the file path
filenameTLEs = 'TLEs.txt';  % Specify the file path

% Actualizar archivo EOP
updateEOP(24,filenameEOP);

% Actualizar TLEs
updateTLE(6,sats,filenameTLEs);

%SGP4

[recef,vecef,rlla,vlla,tsince] = propagar(sats,duracion,precision,filenameTLEs,filenameEOP);

% figure(1)
% 
% plot3(squeeze(recef(1,1,:)),squeeze(recef(1,2,:)),squeeze(recef(1,3,:)))
% grid on
% hold on 
% for i = 2:n_sats
%     plot3(squeeze(recef(i,1,:)),squeeze(recef(i,2,:)),squeeze(recef(i,3,:)))
% end
% hold off
% pbaspect([1 1 1])
% legend('ISS','Starlink','OneWeb','Location','best')
% 
% dist = zeros(n_sats,length(tsince));
% 
% for j = 1:n_sats
%     for i = 1:length(tsince)
%         dist(j,i) = norm(recef(j,:,i));
%     end
% end
% 
% rTierra = 6378.0; %Radio medio de la tierra
% 
% figure(2)
% plot(tsince,dist(1,:)-rTierra)
% grid on
% hold on
% for i = 2:n_sats
%     plot(tsince,dist(i,:)-rTierra)
% end
% hold off
% legend('ISS','Starlink','OneWeb','Location','best')
% 
lim = length(rlla)-(length(rlla)-length(recef));

figure(3)
geoplot(squeeze(rlla(1,1,1:lim)),squeeze(rlla(1,2,1:lim)))
hold on
for i = 2:n_sats
    geoplot(squeeze(rlla(i,1,1:lim)),squeeze(rlla(i,2,1:lim)))
end
for i = 1:n_sats
    geoplot(squeeze(rlla(i,1,floor(lim/2))),squeeze(rlla(i,2,floor(lim/2))),'ok')
end
hold off
legend('ISS','Starlink','OneWeb','Location','northeast')