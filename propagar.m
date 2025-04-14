function [recef,vecef,rlla,vlla,tsince_out] = propagar(sat,instante,duracion,precision,filenameEOP,f,Re,customEpoch)
%PROPAGAR Funcion que calcula posicion y velocidad tanto en ECEF como en LLA 
%   Precision y duracion en minutos
    format long g
    global const
    SAT_Const

    ge = 398600.8; % Earth gravitational constant [km3/s2]
    TWOPI = 2*pi;
    MINUTES_PER_DAY = 1440;
    MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
    MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

    % n_sats = length(sats);
    
    % for n_sat = 1:n_sats
    
    fid = fopen([fullfile('TLEs',int2str(sat.NORAD)),'.txt'],'rt');
    tline = fgetl(fid);

    % read first line
    tline = fgetl(fid);
    Cnum = tline(3:7);      			        % Catalog Number (NORAD)
    SC   = tline(8);					        % Security Classification
    ID   = tline(10:17);			            % Identification Number
    year = str2double(tline(19:20));               % Year
    doy  = str2double(tline(21:32));               % Day of year
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
    if isnan(customEpoch)
        satdata.epoch = epoch;
    else
        satdata.epoch = customEpoch;
    end
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

    epochDaytime = epochToUTC(satdata.epoch);

    tSinceEpoch = floor(minutes(instante-epochDaytime));

    tsince = tSinceEpoch-duracion:precision:tSinceEpoch+duracion; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes]

    num = length(tsince);
    rteme = zeros(3,num);
    vteme = zeros(3,num);
    reci = zeros(3,num);
    veci = zeros(3,num);
    recef = zeros(3,num);
    vecef = zeros(3,num);
    rtod = zeros(3,num);
    vtod = zeros(3,num);
    rlla = zeros(3,num);
    vlla = zeros(3,num);

    [rteme, vteme] = sgp4(tsince, satdata);
    % read Earth orientation parameters
    fid = fopen(filenameEOP,'r');
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
        [reci(:,i),veci(:,i)] = teme2eci(rteme(:,i),vteme(:,i),T,dpsi,deps);
        [recef(:,i),vecef(:,i)] = teme2ecef(rteme(:,i),vteme(:,i),T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2);
        [rtod(:,i), vtod(:,i)] = ecef2tod(recef(:,i),vecef(:,i),T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2,dpsi,deps);

        % LatLonAlt
        rlla(:,i) = ecef2lla(recef(:,i)'.*10^3,f,Re);
        vlla(:,i) = ecef2lla(vecef(:,i)'.*10^3,f,Re);
        
        % %Out vars
        % recef_out(n_sat,:,:) = recef;
        % vecef_out(n_sat,:,:) = vecef;
        % rlla_out(n_sat,:,:)= rlla;
        % vlla_out(n_sat,:,:)= vlla;            

        end

        tsince_out = tsince-tSinceEpoch;
    % end
end

