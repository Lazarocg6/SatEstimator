function [f_doppler, recef, vecef, rlla, bistaticRange, bistaticVelocity, ...
    R1, R2, snr, NAME, ID, latTX, latRX, lonTX, lonRX, elTX, elRX] = paramsToDop(epoch_in,time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    addpath('SGP4_Vectorized')
    format long g
    global const
    SAT_Const
%                      _              _            
%                     | |            | |           
%   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___  ___ 
%  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ _ \/ __|
% | (_| (_) | | | \__ \ || (_| | | | | ||  __/\__ \
%  \___\___/|_| |_|___/\__\__,_|_| |_|\__\___||___/

    % NAME = 'ISS';
    % ID = 25544; % NORAD ID de la ISS
    % NAME = 'CSS';
    % ID = 48274;
    NAME = 'STARLINK33984';
    ID = 63830;

    % RX = [40.45206046037957, -3.726407299669201, 670]; % Coords ETSIT
    RX = [40.6223011985758, -4.010124224894723, 930]; % Coords Villalba

    freq = 143.050e6;% frequency in hertz
    f = 1/298.26; % WGS72 Parameters
    Re = 6378135; % WGS72 Parameters
    lambda = 3e8/freq;

    % RADAR eq params

    eirp_tx = 1e6; % EIRP power TX [Watts]
    RCSb = 5; % bistatic RCS [meters^2]

    Grx = 14; % Gain RX [dB]
    Lsys = 3; % System loses RX [dB]
    Fs = 6; % Noise figure RX [dB]
    % Brx = 0.25e6;% BW RX [Hz]
    SR = 5e6; % Sample rate
    fft_size = 2^16;
    Tint = fft_size/SR; % Integration time [s]

    % GRAVES coords

    latGRAVES = 47.34813145826119;
    lonGRAVES = 5.51487314868131;
    elGRAVES = 180; % Average altitude
    
    % Coordinates to be used
    latTX = latGRAVES;
    lonTX = lonGRAVES;
    elTX = elGRAVES;
    latRX = RX(1);
    lonRX = RX(2);
    elRX = RX(3);

    ge = 398600.8; % Earth gravitational constant [km3/s2]
    TWOPI = 2*pi;
    MINUTES_PER_DAY = 1440;
    MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
    MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

%  _____ _      _____  ______  ___ _____ ___  
% |_   _| |    |  ___| |  _  \/ _ \_   _/ _ \ 
%   | | | |    | |__   | | | / /_\ \| |/ /_\ \
%   | | | |    |  __|  | | | |  _  || ||  _  |
%   | | | |____| |___  | |/ /| | | || || | | |
%   \_/ \_____/\____/  |___/ \_| |_/\_/\_| |_/

    updateEOP(36); % Update EOP every 36 hours
    updateTLE(6, ID); % Update TLE every 6 hours
    
    fid = fopen([fullfile('TLEs',int2str(ID)),'.txt'],'rt');
    tline = fgetl(fid);

    % read first line
    tline = fgetl(fid);
    Cnum = tline(3:7);      			        % Catalog Number (NORAD)
    SC   = tline(8);					        % Security Classification
    ID2   = tline(10:17);			            % Identification Number
    % year2 = str2double(tline(19:20));               % Year
    % doy  = str2double(tline(21:32));               % Day of year
    epoch = str2double(tline(19:32));              % Epoch
    TD1   = str2double(tline(34:43));              % first time derivative
    TD2   = str2double(tline(45:50));              % 2nd Time Derivative
    ExTD2 = str2double(tline(51:52));              % Exponent of 2nd Time Derivative
    BStar = str2double(tline(54:59));              % Bstar/drag Term
    ExBStar = str2double(tline(60:61));            % Exponent of Bstar/drag Term
    BStar = BStar*1e-5*10^ExBStar;
    Etype = tline(63);                          % Ephemeris Type
    Enum  = str2double(tline(65:end));             % Element Number
    
    % read second line
    tline = fgetl(fid);
    i = str2double(tline(9:16));                   % Orbit Inclination (degrees)
    raan = str2double(tline(18:25));               % Right Ascension of Ascending Node (degrees)
    e = str2double(strcat('0.',tline(27:33)));     % Eccentricity
    omega = str2double(tline(35:42));              % Argument of Perigee (degrees)
    M = str2double(tline(44:51));                  % Mean Anomaly (degrees)
    no = str2double(tline(53:63));                 % Mean Motion
    a = ( ge/(no*2*pi/86400)^2 )^(1/3);         % semi major axis (km)
    rNo = str2double(tline(65:end));               % Revolution Number at Epoch
    
    fclose(fid);
    if isnan(epoch_in)
        satdata.epoch = epoch;
    else
        satdata.epoch = epoch_in;
    end

    satdata.norad_number = Cnum;
    satdata.bulletin_number = ID2;
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

    epochDaytime = epochToUTC(satdata.epoch);
    
    
    tsince = minutes(datetime(time,'ConvertFrom','datenum','TimeZone','Local')-epochDaytime);

    % tSinceEpoch = floor(minutes(instante-epochDaytime));
    % 
    % tsince = tSinceEpoch-duracion:precision:tSinceEpoch+duracion; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes]

% ______                                  _             
% | ___ \                                | |            
% | |_/ / __ ___  _ __   __ _  __ _  __ _| |_ ___  _ __ 
% |  __/ '__/ _ \| '_ \ / _` |/ _` |/ _` | __/ _ \| '__|
% | |  | | | (_) | |_) | (_| | (_| | (_| | || (_) | |   
% \_|  |_|  \___/| .__/ \__,_|\__, |\__,_|\__\___/|_|   
%                | |           __/ |                    
%                |_|          |___/

    [rteme, vteme] = sgp4(tsince, satdata);

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
    
    num = length(tsince); 
    reci = zeros(3,num);
    veci = zeros(3,num);
    recef = zeros(3,num);
    vecef = zeros(3,num);
    rtod = zeros(3,num);
    vtod = zeros(3,num);
    rlla = zeros(3,num);
    vlla = zeros(3,num);
    
    year2 = year(epochToUTC(satdata.epoch));
    doy = satdata.epoch-((year2-2000)*1000);

    [mon,day2,hr,minute,sec] = days2mdh(year2,doy);
    MJD_Epoch = Mjday(year2,mon,day2,hr,minute,sec);

    for i = 1:num
        MJD_UTC = MJD_Epoch+tsince(i)/1440;
        % Earth Orientation Parameters
        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,~,~,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
        [~,~,~,TT_UTC,~] = timediff(UT1_UTC,TAI_UTC);
        MJD_UT1 = MJD_UTC + UT1_UTC/86400;
        MJD_TT  = MJD_UTC + TT_UTC/86400;
        T = (MJD_TT-const.MJD_J2000)/36525;
        [reci(:,i),veci(:,i)] = teme2eci(rteme(:,i),vteme(:,i),T,dpsi,deps);
        [recef(:,i),vecef(:,i)] = teme2ecef(rteme(:,i),vteme(:,i),T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2);
        [rtod(:,i), vtod(:,i)] = ecef2tod(recef(:,i),vecef(:,i),T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2,dpsi,deps);

        % LatLonAlt
        rlla(:,i) = ecef2lla(recef(:,i)'.*10^3,f,Re);
        vlla(:,i) = ecef2lla(vecef(:,i)'.*10^3,f,Re); 

    end

        rlla = evitarSaltos(rlla);

        % tsince_out = tsince-tSinceEpoch;

% ______ _     _        _   _       ______                             
% | ___ (_)   | |      | | (_)      | ___ \                            
% | |_/ /_ ___| |_ __ _| |_ _  ___  | |_/ /_ _ _ __ __ _ _ __ ___  ___ 
% | ___ \ / __| __/ _` | __| |/ __| |  __/ _` | '__/ _` | '_ ` _ \/ __|
% | |_/ / \__ \ || (_| | |_| | (__  | | | (_| | | | (_| | | | | | \__ \
% \____/|_|___/\__\__,_|\__|_|\___| \_|  \__,_|_|  \__,_|_| |_| |_|___/

        [bistaticRange, bistaticVelocity, R1, R2, llaDIST, ecefDIST] = bistaticParams(latTX, ...
            lonTX, elTX, latRX, lonRX, elRX, recef * 1000, vecef * 1000, f, Re);

% ______                  _           
% |  _  \                | |          
% | | | |___  _ __  _ __ | | ___ _ __ 
% | | | / _ \| '_ \| '_ \| |/ _ \ '__|
% | |/ / (_) | |_) | |_) | |  __/ |   
% |___/ \___/| .__/| .__/|_|\___|_|   
%            | |   | |                
%            |_|   |_|

        f_doppler = (bistaticVelocity)./(-lambda);

%  _____ _____      ______          _            
% |  ___|  _  |     | ___ \        | |           
% | |__ | | | |     | |_/ /__ _  __| | __ _ _ __ 
% |  __|| | | |     |    // _` |/ _` |/ _` | '__|
% | |___\ \/' /     | |\ \ (_| | (_| | (_| | |   
% \____/ \_/\_\     \_| \_\__,_|\__,_|\__,_|_|

        l_sys = 10^(Lsys/10);
        g_rx = 10^(Grx/10);
        fs = 10^(Fs/10);

        snr = snr_in(eirp_tx,g_rx,lambda,RCSb,R1,R2,l_sys,fs,Tint);

end

