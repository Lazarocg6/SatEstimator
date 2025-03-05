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

global const
SAT_Const

% duracion = 90;% En minutos, 'duracion' minutos antes y despues del momento actual 

ge = 398600.8; % Earth gravitational constant [km3/s2]
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

ISS = 25544;
Starlink = 44737;
OneWeb = 44057;

sats = [ISS Starlink OneWeb];

n_sats = length(sats);

% satData = cell(2,9,n_sats);

lines1 = cell(1,length(n_sats));
lines2 = cell(1,length(n_sats));

% Actualizar archivo EOP

filename = 'EOP-All.txt';  % Specify the file path
filenameTLEs = 'TLEs.txt';  % Specify the file path

if isfile(filename)  % Check if the file exists
    fid = fopen(filename, 'r');  % Open for reading only
    if fid == -1
        error('Error opening the file.');
    else
        version = fgetl(fid);
        fecha = fgetl(fid);
    end

    fechaSplitted = split(fecha);
    horaSplitted = split(fechaSplitted{5},':');
    
    months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
          
    % Find the index of the input month in the list
    monthNum = find(strcmpi(fechaSplitted{3}, months));

    fechaUpdate = datetime(str2double(fechaSplitted{2}),monthNum ...
        ,str2double(fechaSplitted{4}),str2double(horaSplitted{1}) ...
        ,str2double(horaSplitted{2}),str2double(horaSplitted{3}));
    fechaActual = datetime("now");

    if hours(fechaActual-fechaUpdate) > 24
        websave(filename,'https://celestrak.org/SpaceData/EOP-All.txt')
        disp('EOP llamado \n')
    end

    fclose(fid);
else
    websave(filename,'https://celestrak.org/SpaceData/EOP-All.txt')
    disp('EOP llamado \n')
end 

% Actualizar TLEs

if isfile(filenameTLEs)  % Check if the file exists

    fid = fopen(filenameTLEs,'rt');
    fechaUpdate = datetime(fgetl(fid));
    fechaActual = datetime("now");
    diff = hours(fechaActual-fechaUpdate);
    fclose(fid);

    if diff > 6
        fid = fopen(filenameTLEs,'wt');
        fprintf(fid,'%s\n',datetime("now"));
        for n_sat = 1:n_sats
            url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(sats(n_sat)),'&FORMAT=2le'];
            data = webread(url);  % Read the content from the URL
            lines = splitlines(data);  % Split the content into lines
            fprintf(fid,'%s\n',lines{1});
            fprintf(fid,'%s\n',lines{2});
        end
        disp('TLEs llamados')
        fclose(fid);
    end
    
else
    fid = fopen(filenameTLEs,'wt');
    fprintf(fid,'%s\n',datetime("now"));
    for n_sat = 1:n_sats
        url = ['https://celestrak.org/NORAD/elements/gp.php?CATNR=',num2str(sats(n_sat)),'&FORMAT=2le'];
        data = webread(url);  % Read the content from the URL
        lines = splitlines(data);  % Split the content into lines
        fprintf(fid,'%s\n',lines{1});
        fprintf(fid,'%s\n',lines{2});
    end
    disp('TLEs llamados')
    fclose(fid);
end 

tsince = 0:1:180; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes] 

num = length(tsince);
rteme = zeros(n_sats,3,num);
vteme = zeros(n_sats,3,num);
reci = zeros(n_sats,3,num);
veci = zeros(n_sats,3,num);
recef = zeros(n_sats,3,num);
vecef = zeros(n_sats,3,num);
rtod = zeros(n_sats,3,num);
vtod = zeros(n_sats,3,num);
% rlla = zeros(n_sats,3,num);
% vlla = zeros(n_sats,3,num);

%SGP4

for n_sat = 1:n_sats

    fid = fopen(filenameTLEs,'rt');
    for i = 1:(2*n_sat)-1
        tline = fgetl(fid);
    end

    % read first line
    tline = fgetl(fid);
    Cnum = tline(3:7);      			        % Catalog Number (NORAD)
    SC   = tline(8);					        % Security Classification
    ID   = tline(10:17);			            % Identification Number
    year = str2num(tline(19:20));               % Year
    doy  = str2num(tline(21:32));               % Day of year
    epoch = str2num(tline(19:32));              % Epoch
    TD1   = str2num(tline(34:43));              % first time derivative
    TD2   = str2num(tline(45:50));              % 2nd Time Derivative
    ExTD2 = str2num(tline(51:52));              % Exponent of 2nd Time Derivative
    BStar = str2num(tline(54:59));              % Bstar/drag Term
    ExBStar = str2num(tline(60:61));            % Exponent of Bstar/drag Term
    BStar = BStar*1e-5*10^ExBStar;
    Etype = tline(63);                          % Ephemeris Type
    Enum  = str2num(tline(65:end));             % Element Number
    
    % read second line
    tline = fgetl(fid);
    i = str2num(tline(9:16));                   % Orbit Inclination (degrees)
    raan = str2num(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e = str2num(strcat('0.',tline(27:33)));     % Eccentricity
    omega = str2num(tline(35:42));              % Argument of Perigee (degrees)
    M = str2num(tline(44:51));                  % Mean Anomaly (degrees)
    no = str2num(tline(53:63));                 % Mean Motion
    a = ( ge/(no*2*pi/86400)^2 )^(1/3);         % semi major axis (km)
    rNo = str2num(tline(65:end));               % Revolution Number at Epoch
    
    fclose(fid);
    
    satdata.epoch = epoch;
    satdata.norad_number = Cnum;
    satdata.bulletin_number = ID;
    satdata.classification = SC; % almost always 'U'
    satdata.revolution_number = rNo;
    satdata.ephemeris_type = Etype;
    satdata.xmo = M * (pi/180);
    satdata.xnodeo = raan * (pi/180);
    satdata.omegao = omega * (pi/180);
    satdata.xincl = i * (pi/180);
    satdata.eo = e;
    satdata.xno = no * TWOPI / MINUTES_PER_DAY;
    satdata.xndt2o = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
    satdata.xndd6o = TD2 * 10^ExTD2 * TWOPI / MINUTES_PER_DAY_CUBED;
    satdata.bstar = BStar;

    % tsince = 0:1:375; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes] 
    % 
    % 
    % fraction = epoch-floor(epoch);
    % total_hours = fraction * 24;  % Convert to total hours
    % epochH = floor(total_hours);   % Get whole hours
    % epochM = floor((total_hours - epochH) * 60);  % Get whole minutes
    % epochS = floor((total_hours - epochH) * 60 - epochM) * 60;  % Get seconds
    % 
    % dateVal = datetime(year, 1, 1) + days(doy - 1); % Convert DOY to date
    % monthNum = month(dateVal);  % Extract the month
    % dayNum = day(dateVal, 'dayofmonth');
    % 
    % epochDaytime = datetime(year+2000,monthNum,dayNum,epochH,epochM,epochS);
    % 
    % tSinceEpoch = floor(minutes(datetime("now")-epochDaytime));

    % tsince = tSinceEpoch-duracion:1:tSinceEpoch+duracion; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes] 
    % 
    % num = length(tsince);
    % rteme = zeros(n_sats,3,num);
    % vteme = zeros(n_sats,3,num);
    % reci = zeros(n_sats,3,num);
    % veci = zeros(n_sats,3,num);
    % recef = zeros(n_sats,3,num);
    % vecef = zeros(n_sats,3,num);
    % rtod = zeros(n_sats,3,num);
    % vtod = zeros(n_sats,3,num);
    % % rlla = zeros(n_sats,3,num);
    % % vlla = zeros(n_sats,3,num);

    [rteme(n_sat,:,:), vteme(n_sat,:,:)] = sgp4(tsince, satdata);
    % read Earth orientation parameters
    fid = fopen('EOP-All.txt','r');
    %  ----------------------------------------------------------------------------------------------------
    % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    % |(0h UTC)           "         "          s          s          "        "          "         "     s 
    %  ----------------------------------------------------------------------------------------------------
    while ~feof(fid)
        tline = fgetl(fid);
        k = strfind(tline,'NUM_OBSERVED_POINTS');
        if (k == 1)
            numrecsobs = str2double(tline(21:end));
            tline = fgetl(fid);
            for i=1:numrecsobs
                eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            for i=1:4
                tline = fgetl(fid);
            end
            numrecspred = str2double(tline(22:end));
            tline = fgetl(fid);
            for i=numrecsobs+1:numrecsobs+numrecspred
                eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            break
        end
    end
    fclose(fid);
    
    
    if (year < 57)
        year = year + 2000;
    else
        year = year + 1900;
    end
    [mon,day2,hr,minute,sec] = days2mdh(year,doy);
    MJD_Epoch = Mjday(year,mon,day2,hr,minute,sec);

    cont = 1;
    for i = 1:num
        MJD_UTC = MJD_Epoch+tsince(i)/1440;
        % Earth Orientation Parameters
        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
        [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
        MJD_UT1 = MJD_UTC + UT1_UTC/86400;
        MJD_TT  = MJD_UTC + TT_UTC/86400;
        T = (MJD_TT-const.MJD_J2000)/36525;
        [reci(n_sat,:,i),veci(n_sat,:,i)] = teme2eci(rteme(n_sat,:,i)',vteme(n_sat,:,i)',T,dpsi,deps);
        [recef(n_sat,:,i),vecef(n_sat,:,i)] = teme2ecef(rteme(n_sat,:,i)',vteme(n_sat,:,i)',T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2);
        [rtod(n_sat,:,i), vtod(n_sat,:,i)] = ecef2tod(recef(n_sat,:,i)',vecef(n_sat,:,i)',T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2,dpsi,deps);
        % LatLonAlt
        rtemp = ecef2lla(recef(n_sat,:,i).*10e3,'WGS84');
        vtemp = ecef2lla(vecef(n_sat,:,i).*10e3,'WGS84');
        if i == 1 
            rlla(n_sat,:,cont) = rtemp;
            vlla(n_sat,:,cont) = vtemp;
        else
            if (abs(rtemp(2)-rlla(n_sat,2,cont-1))>90)
                rlla(n_sat,:,cont) = [NaN NaN NaN];
                vlla(n_sat,:,cont) = [NaN NaN NaN];
                cont = cont+1;
                rlla(n_sat,:,cont) = rtemp;
                vlla(n_sat,:,cont) = vtemp;
            else
                rlla(n_sat,:,cont) = rtemp;
                vlla(n_sat,:,cont) = vtemp;
            end
            rlla(n_sat,:,cont) = rtemp;
            vlla(n_sat,:,cont) = vtemp;
        end
        cont = cont+1;
    end
end

figure(1)

plot3(squeeze(recef(1,1,:)),squeeze(recef(1,2,:)),squeeze(recef(1,3,:)))
grid on
hold on 
for i = 2:n_sats
    plot3(squeeze(recef(i,1,:)),squeeze(recef(i,2,:)),squeeze(recef(i,3,:)))
end
hold off
pbaspect([1 1 1])
legend('ISS','Starlink','OneWeb','Location','best')

dist = zeros(n_sats,length(tsince));

for j = 1:n_sats
    for i = 1:length(tsince)
        dist(j,i) = norm(recef(j,:,i));
    end
end

rTierra = 6378.0; %Radio medio de la tierra

figure(2)
plot(tsince,dist(1,:)-rTierra)
grid on
hold on
for i = 2:n_sats
    plot(tsince,dist(i,:)-rTierra)
end
hold off
legend('ISS','Starlink','OneWeb','Location','best')

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


% NUEVO COMENTARIO