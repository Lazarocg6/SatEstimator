clear
close all

addpath('SGP4_Vectorized')
format long g
global const
SAT_Const

% Constants
inst = datetime('30-May 19:37:37','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 95; % In minutes
precision = 10 / 60; % Precision in minutes
propagateB4andafter = true; % False propaga hacia delante desde inst

ge = 398600.8; % Earth gravitational constant in km^3/s^2
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

seed = randi(2^32,1);

fig3 = figure(3);
gx = geoaxes(fig3);                % Create clean geoaxes

fig4 = figure(4);
gx2 = geoaxes(fig4);                % Create clean geoaxes

freq = 143.050e6;% frequency in hertz
f = 1/298.26; % WGS72 Parameters
Re = 6378135; % WGS72 Parameters
lambda = 3e8/freq;

% GRAVES coords

latGRAVES = 47.34813145826119;
lonGRAVES = 5.51487314868131;
elGRAVES = 180; % Average altitude

% RX coords
RX1 = [40.45206046037957, -3.726407299669201, 670];

% Coordinates to be used
latTX = latGRAVES;
lonTX = lonGRAVES;
elTX = elGRAVES;
latRX = RX1(1);
lonRX = RX1(2);
elRX = RX1(3);

[time, DTtime] = initTimes(inst, duracion, precision,propagateB4andafter);

% Open the TLE file
fid = fopen('sat000059513.txt', 'r');

if fid == -1
    error('Failed to open the file.');
end

tle_index = 0;

while ~feof(fid)
    % Read first line of TLE
    tline1 = fgetl(fid);
    if ~ischar(tline1) || length(tline1) < 69
        break;  % Skip incomplete line or end of file
    end

    % Read second line of TLE
    tline2 = fgetl(fid);
    if ~ischar(tline2) || length(tline2) < 69
        break;  % Skip incomplete line or end of file
    end

    % Parse line 1
    Cnum = tline1(3:7);                          % Catalog Number (NORAD)
    SC   = tline1(8);                             % Security Classification
    ID2  = tline1(10:17);                         % Identification Number
    epoch = str2double(tline1(19:32));            % Epoch
    TD1   = str2double(tline1(34:43));            % First time derivative
    TD2   = str2double(tline1(45:50));            % Second time derivative
    ExTD2 = str2double(tline1(51:52));            % Exponent for TD2
    BStar = str2double(tline1(54:59));            % Bstar term
    ExBStar = str2double(tline1(60:61));          % Exponent for Bstar
    BStar = BStar * 1e-5 * 10^ExBStar;
    Etype = tline1(63);                           % Ephemeris Type
    Enum  = str2double(tline1(65:end));           % Element Number

    % Parse line 2
    i     = str2double(tline2(9:16));             % Inclination (deg)
    raan  = str2double(tline2(18:25));            % RA of Ascending Node (deg)
    e     = str2double(['0.', tline2(27:33)]);    % Eccentricity
    omega = str2double(tline2(35:42));            % Argument of Perigee (deg)
    M     = str2double(tline2(44:51));            % Mean Anomaly (deg)
    no    = str2double(tline2(53:63));            % Mean Motion (rev/day)
    a     = (ge / (no * 2 * pi / 86400)^2)^(1/3);  % Semi-major axis (km)
    rNo   = str2double(tline2(65:end));           % Rev number at epoch

    satdata(tle_index+1).norad_number = Cnum;
    satdata(tle_index+1).bulletin_number = ID2;
    satdata(tle_index+1).classification = SC; % almost always 'U'
    satdata(tle_index+1).revolution_number = rNo;
    satdata(tle_index+1).ephemeris_type = Etype;
    satdata(tle_index+1).xmo = M * (pi/180);
    satdata(tle_index+1).xnodeo = raan * (pi/180);
    satdata(tle_index+1).omegao = omega * (pi/180);
    satdata(tle_index+1).xincl = i * (pi/180);
    satdata(tle_index+1).eo = e;
    satdata(tle_index+1).xno = no * TWOPI / MINUTES_PER_DAY;
    satdata(tle_index+1).xndt2o = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
    satdata(tle_index+1).xndd6o = TD2 * 10^ExTD2 * TWOPI / MINUTES_PER_DAY_CUBED;
    satdata(tle_index+1).bstar = BStar;
    satdata(tle_index+1).epoch = epoch;

    % Optional: Do something with the data
    tle_index = tle_index + 1;
    % fprintf('TLE %d: NORAD ID = %s | Epoch = %.8f | SMA = %.2f km\n', ...
    %         tle_index, Cnum, epoch, a);

    

end

fclose(fid);

values = [1 6:6:60];

for tle_index = flip(values)

    epochDaytime = epochToUTC(satdata(tle_index).epoch);
    tsince = minutes(datetime(time,'ConvertFrom','datenum','TimeZone','Local')-epochDaytime);

    % ______                                  _             
    % | ___ \                                | |            
    % | |_/ / __ ___  _ __   __ _  __ _  __ _| |_ ___  _ __ 
    % |  __/ '__/ _ \| '_ \ / _` |/ _` |/ _` | __/ _ \| '__|
    % | |  | | | (_) | |_) | (_| | (_| | (_| | || (_) | |   
    % \_|  |_|  \___/| .__/ \__,_|\__, |\__,_|\__\___/|_|   
    %                | |           __/ |                    
    %                |_|          |___/

    [rteme, vteme] = sgp4(tsince, satdata(tle_index));

    % read Earth orientation parameters
    fid2 = fopen('EOP-All.txt','r');
    %  ----------------------------------------------------------------------------------------------------
    % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    % |(0h UTC)           "         "          s          s          "        "          "         "     s 
    %  ----------------------------------------------------------------------------------------------------
    while ~feof(fid2)
        tline = fgetl(fid2);
        k = strfind(tline,'NUM_OBSERVED_POINTS');
        if (k == 1)
            numrecsobs = str2double(tline(21:end));
            tline = fgetl(fid2);
            for i=1:numrecsobs
                eopdata(:,i) = fscanf(fid2,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            for i=1:4
                tline = fgetl(fid2);
            end
            numrecspred = str2double(tline(22:end));
            tline = fgetl(fid2);
            for i=numrecsobs+1:numrecsobs+numrecspred
                eopdata(:,i) = fscanf(fid2,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            break
        end
    end
    fclose(fid2);
    
    num = length(tsince); 
    reci = zeros(3,num);
    veci = zeros(3,num);
    recef = zeros(3,num);
    vecef = zeros(3,num);
    rtod = zeros(3,num);
    vtod = zeros(3,num);
    rlla = zeros(3,num);
    vlla = zeros(3,num);
    
    year2 = year(epochToUTC(satdata(tle_index).epoch));
    doy = satdata(tle_index).epoch-((year2-2000)*1000);

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

    rlla_t = evitarSaltos(rlla);

    [bistaticRange, bistaticVelocity, R1, R2, llaDIST, ecefDIST] = bistaticParams(latTX, ...
        lonTX, elTX, latRX, lonRX, elRX, recef * 1000, vecef * 1000, f, Re);

    f_doppler = (bistaticVelocity)./(-lambda);

    colors = turbo(max(values));  % or parula(60), jet(60), etc.
    % if ismember(tle_index, [1, 6, 54, 60])
    name = ['TLE ' num2str(tle_index)];
    fitted_values = [1,60];
    if ismember(tle_index, fitted_values)
        figure(3)
            hold(gx, 'on')               % Hold on for plotting
            geoplot(gx, rlla_t(1,:), rlla_t(2,:),'LineWidth', 2, 'DisplayName', name,'LineStyle','--')
            hold(gx, 'off')
        figure(1) 
            hold on
            plot(f_doppler, DTtime, '-', 'LineWidth', 2, 'DisplayName', name,'LineStyle','--');  % normal contrast
            hold off

        %  _                    _   _____                                 
        % | |                  | | /  ___|                                
        % | |     ___  __ _ ___| |_\ `--.  __ _ _   _  __ _ _ __ ___  ___ 
        % | |    / _ \/ _` / __| __|`--. \/ _` | | | |/ _` | '__/ _ \/ __|
        % | |___|  __/ (_| \__ \ |_/\__/ / (_| | |_| | (_| | | |  __/\__ \
        % \_____/\___|\__,_|___/\__\____/ \__, |\__,_|\__,_|_|  \___||___/
        %                                    | |                          
        %
        if tle_index == max(values)
            epoch_og = satdata(tle_index).epoch; %Epoch usado para generar el ruido

            rng(seed)

            epochd = epoch_og+rand/4000; %Epoch desplazado

            fprintf('-------------------------NoiseGen Data-------------------------\n')
            fprintf('Epoch diff = %s, noisy ->%10.10f, OG ->%10.10f\n', (epochToUTC(epochd)-epochToUTC(epoch)),epochd,epoch)

            f_doppler_mod = paramsToDop(epochd,time,RX1,[],satdata(tle_index));

            cte = 1e3;

            mean = 0;
            stdDeviation = 0.2;
            l = length(f_doppler_mod);
            n = randn(1,l)*stdDeviation+mean;

            rng(seed+1)

            mean = 2;
            stdDeviation = 0.5;
            n2 = randn(1)*stdDeviation+mean;

            f_dev = n2*100;

            fprintf('Frequency deviation -> %5.4fHz \n',f_dev)
            fprintf('---------------------------------------------------------------\n')

            % noisy = f_doppler+(n*cte)+(n2*100);
            ydata = f_doppler_mod+(n*cte)+f_dev; 

            og_resnorm = paramsToDop(epoch_og,time,RX1)-ydata;
        end

        x0 = [satdata(tle_index).epoch,0];
        fun = @(optimize,time) paramsToDop(optimize(1),time,RX1,[],satdata(tle_index))+optimize(2);
        [res,resnorm] = lsqcurvefit(fun,x0,time,ydata);
        resnorm_list(tle_index) = resnorm;
        t = (epochToUTC(epochd)-epochToUTC(res(1)));
        t.Format = 'hh:mm:ss.SSS';
        fprintf('-------------------------Bist Data-------------------------tle_index=%i\n',tle_index)
        fprintf('Epoch diff = %s, noisy ->%10.10f, fitted ->%10.10f\n', t,epochd,res(1))
        fprintf('Frequency correction diff -> %5.4fHz, Frequency correction (fitted) -> %5.4fHz\n',res(2)-f_dev,res(2))
        fprintf('ResNorm -> %5.4f\n',resnorm)
        fprintf('-----------------------------------------------------------\n')
        [f_doppler_f, recef_f, vecef_f, rlla_f, bistaticRange_f, bistaticVelocity_f, ...
            R1_f, R2_f, snr_f, name_f, ID_f, latTX_f, latRX_f, lonTX_f, lonRX_f, elTX_f, elRX_f] = paramsToDop(res(1),time,RX1,[],satdata(tle_index));

        name = ['Fitted TLE ' num2str(tle_index)];


        figure(2)
            hold on
            if tle_index == max(values)
                plot(DTtime, ydata,'DisplayName','Noisy')
                plot(DTtime,f_doppler_mod,'DisplayName','Original')
            end
            plot(DTtime,f_doppler_f+res(2),'DisplayName',name)
            hold off
            grid on
            legend
            title('Bistatic Scenario')
        
        figure(4)
            hold(gx2,'on')
            geoplot(rlla_f(1,:),rlla_f(2,:))
            hold(gx2,'off')

    else
        c = colors(tle_index, :);
        pastel = c + 0;  % lighten color
        pastel = pastel / max(pastel);  % normalize to 0–1 range
        figure(1)
        hold on
            plot(f_doppler, DTtime, '-', 'Color', pastel, ...
                 'LineWidth', 0.2, 'DisplayName', name);
            hold off
        figure(3)
            hold(gx,'on')
            geoplot(gx, rlla_t(1,:), rlla_t(2,:),'LineWidth', 0.2, ...
                'DisplayName', name,'LineStyle','-','Color',pastel)
            hold(gx,'off')
    end     

end

figure(1)
    grid on
    ax = gca;
    set(ax, 'YDir','reverse')
    grid on
    hold off
    legend
figure(3)
    hold(gx,"on")
    geolimits(gx, [35 50], [-14 14])
    geobasemap(gx, 'darkwater')
    geoplot(gx,latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX','HandleVisibility','off');
    geoplot(gx,latRX, lonRX, '^r', 'DisplayName', 'RX', 'Tag', 'RX','HandleVisibility','off');
    legend
    hold(gx, 'off')
figure(4)
    hold(gx2,"on")
    geolimits(gx2, [35 50], [-14 14])
    geobasemap(gx2, 'darkwater')
    geoplot(gx2,latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX','HandleVisibility','off');
    geoplot(gx2,latRX, lonRX, '^r', 'DisplayName', 'RX', 'Tag', 'RX','HandleVisibility','off');
    legend
    hold(gx2, 'off')   