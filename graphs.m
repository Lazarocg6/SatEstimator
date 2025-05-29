function graphs(time,DTtime,fitter,fitterType,inst,filter,multistat,bistat,RX1,RX2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all

if multistat
    [f_dopplerMult, recefMult, vecefMult, rllaMult, bistaticRangeMult, bistaticVelocityMult, ...
        R1Mult, R2Mult, snrMult, nameMult, IDMult, latTXMult, latRXMult, lonTXMult, lonRXMult, elTXMult, elRXMult] = paramsToDop(NaN,time,RX1,RX2); 
    f_dopplerMult = [f_dopplerMult(1:end/2) NaN f_dopplerMult((end/2)+1:end)];
end

[f_doppler, recef, vecef, rlla, bistaticRange, bistaticVelocity, ...
    R1, R2, snr, name, ID, latTX, latRX, lonTX, lonRX, elTX, elRX] = paramsToDop(NaN,time,RX1);

%  _____     _  ______      _        
% |  ___|   | | |  _  \    | |       
% | |____  _| |_| | | |__ _| |_ __ _ 
% |  __\ \/ / __| | | / _` | __/ _` |
% | |___>  <| |_| |/ / (_| | || (_| |
% \____/_/\_\\__|___/ \__,_|\__\__,_|

% STARLINK33984
% load('detections130520250037.mat','detections_out','f_axis','t_axis')
% x_range = 900:1100;
% y_range = 2070:2300;

% % STARLINK33984
% load('detections130520252001.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK5978
% load('140520250018STARLINK5979detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 970:1200;
% y_range = 150:600;

% % STARLINK1092
% load('140520250052STARLINK1092detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('140520250212STARLINK3824detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 400:1200;
% y_range = 50:1200;

% % STARLINK3824
% load('140520250245STARLINK1833detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 600:1100;
% y_range = 250:900;

% % STARLINK3824
% load('140520250318STARLINK32303detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('140520251018STARLINK33851detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 800:1200;
% y_range = 1100:1900;

% % STARLINK3824
% load('140520251116STARLINK2030detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 900:1050;
% y_range = 650:1200;

% % STARLINK3824
% load('220520251718STARLINK4671detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1000:1050;
% y_range = 2900:3300;

% % STARLINK3824
% load('220520251240STARLINK1582detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251311STARLINK31963detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251348STARLINK31204detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251416STARLINK30173detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251547STARLINK30232detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251611STARLINK1021detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251658STARLINK2284detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('220520251718STARLINK4671detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('250520251804STARLINK1693detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1:length(f_axis);
% y_range = 1:length(t_axis);

% % STARLINK3824
% load('250520251838STARLINK30191detecciones.mat','detections_out','f_axis','t_axis')
% x_range = 1030:1200;
% y_range = 1200:1740;

% % STARLINK3824
% load('070520251045CSSdetecciones.mat','detections_out','f_axis','t_axis')
% x_range = 7000:8000;
% y_range = 380:450;
% % x_range = 1:length(f_axis);
% % y_range = 1:length(t_axis);

% STARLINK3824
load('280520251304STARLINK4650detecciones.mat','detections_out','f_axis','t_axis')
x_range = 1:length(f_axis);
y_range = 1:length(t_axis);


if strcmp(fitterType,'real')

    [fila,columna] = find(detections_out(y_range,x_range));
    x = f_axis(x_range);
    x = x(columna);
    y = t_axis(y_range);
    y = y(fila);
    y = datenum(y);
    xdata = datetime(y,'ConvertFrom','datenum','TimeZone','local');

end

%  _                    _   _____                                 
% | |                  | | /  ___|                                
% | |     ___  __ _ ___| |_\ `--.  __ _ _   _  __ _ _ __ ___  ___ 
% | |    / _ \/ _` / __| __|`--. \/ _` | | | |/ _` | '__/ _ \/ __|
% | |___|  __/ (_| \__ \ |_/\__/ / (_| | |_| | (_| | | |  __/\__ \
% \_____/\___|\__,_|___/\__\____/ \__, |\__,_|\__,_|_|  \___||___/
%                                    | |                          
%                                    |_|

if fitter
    if strcmp(fitterType,'sim')

        seed = now;

        if multistat
            [ydataMult, epochdMult,epoch_ogMult,f_dev] = noiseGen(time,ID,RX1,seed,RX2);
            x0Mult = [epoch_ogMult,0];
            funMult = @(optimize,time) paramsToDop(optimize(1),time,RX1,RX2)+optimize(2);
            resMult = lsqcurvefit(funMult,x0Mult,time,ydataMult);
            tMult = (epochToUTC(epochdMult)-epochToUTC(resMult(1)));
            tMult.Format = 'hh:mm:ss.SSS';
            fprintf('-------------------------Mult Data-------------------------\n')
            fprintf('Epoch diff = %s, noisy ->%10.10f, fitted ->%10.10f\n', tMult,epochdMult,resMult(1))
            fprintf('Frequency correction diff -> %5.4fHz, Frequency correction (fitted) -> %5.4fHz\n',resMult(2)-f_dev,resMult(2))
            fprintf('-----------------------------------------------------------\n')
            fittedMult = paramsToDop(resMult(1),time,RX1,RX2)+resMult(2);
            fittedMult = [fittedMult(1:end/2) NaN fittedMult((end/2)+1:end)];
            ydataMult_nan = [ydataMult(1:end/2) NaN ydataMult((end/2)+1:end)];
        end
        if bistat

            [ydata, epochd,epoch_og,f_dev] = noiseGen(time,ID,RX1,seed);         
            x0 = [epoch_og,0];
            fun = @(optimize,time) paramsToDop(optimize(1),time,RX1)+optimize(2);
            res = lsqcurvefit(fun,x0,time,ydata);
            t = (epochToUTC(epochd)-epochToUTC(res(1)));
            t.Format = 'hh:mm:ss.SSS';
            fprintf('-------------------------Bist Data-------------------------\n')
            fprintf('Epoch diff = %s, noisy ->%10.10f, fitted ->%10.10f\n', t,epochd,res(1))
            fprintf('Frequency correction diff -> %5.4fHz, Frequency correction (fitted) -> %5.4fHz\n',res(2)-f_dev,res(2))
            fprintf('-----------------------------------------------------------\n')
            fitted = paramsToDop(res(1),time,RX1)+res(2);

        end

    elseif strcmp(fitterType,'real')
        fid = fopen([fullfile('TLEs',int2str(ID)),'.txt'],'rt');
        tline = fgetl(fid);
        tline = fgetl(fid);
        epoch_og = str2double(tline(19:32));
        ydata = x;
        x0 = [epoch_og,0];
        fun = @(optimize,time) paramsToDop(optimize(1),time,RX1)+optimize(2);
        % res = lsqcurvefit(@paramsToDop,x0,y,ydata);
        res = lsqcurvefit(fun,x0,y,ydata);
        t = (epochToUTC(epoch_og)-epochToUTC(res(1)));
        t.Format = 'hh:mm:ss.SSS';
        fprintf('Epoch diff = %s, original ->%10.10f, fitted ->%10.10f\n', t,epoch_og,res(1))
        fprintf('Frequency correction -> %5.4fHz\n',res(2))
        fitted = paramsToDop(res(1),time,RX1)+res(2);
    end

end

% ______                             
% | ___ \                            
% | |_/ /_ _ _ __ __ _ _ __ ___  ___ 
% |  __/ _` | '__/ _` | '_ ` _ \/ __|
% | | | (_| | | | (_| | | | | | \__ \
% \_|  \__,_|_|  \__,_|_| |_| |_|___/

% Para calcular el tiempo de integracion maximo
fs = 5e6; % Sample rate
DF = 200; % Decimation factor

if multistat
    timeMult = [time NaN time];
    salto = datetime(NaT, 'TimeZone', DTtime.TimeZone);
    DTtimeMult = [DTtime salto DTtime];
    recefMult = [recef nan(3,1) recef];
    vecefMult = [vecef nan(3,1) vecef];
    rllaMult = [rlla nan(3,1) rlla];
    % R1 = [R1(1:end/2) NaN R1((end/1)+1:end)];
    % R2 = [R2(1:end/2) NaN R2((end/1)+1:end)];
end

%   ___      _____ _     ______ _ _ _            
%  / _ \    |  ___| |    |  ___(_) | |           
% / /_\ \___| |__ | |    | |_   _| | |_ ___ _ __ 
% |  _  |_  /  __|| |    |  _| | | | __/ _ \ '__|
% | | | |/ /| |___| |    | |   | | | ||  __/ |   
% \_| |_/___\____/|_|    \_|   |_|_|\__\___|_|

%Filter parameters
elevMinTX = 10; % degrees
elevMaxTX = 50;
azMinTX = 90;
azMaxTX = 270;
elevMinRX = 0; % degrees
elevMaxRX = 90;
azMinRX = 0;
azMaxRX = 360;

[azRX,elevRX,slantRangeRX] = ecef2aer(recef(1,:)*10^3,recef(2,:)*10^3, ...
    recef(3,:)*10^3,latRX,lonRX,elRX,referenceEllipsoid('wgs72'));
[azTX,elevTX,slantRangeTX] = ecef2aer(recef(1,:)*10^3,recef(2,:)*10^3, ...
    recef(3,:)*10^3,latTX,lonTX,elTX,referenceEllipsoid('wgs72'));

if filter
    limElTX = (elevTX>elevMinTX) & (elevTX<elevMaxTX);
    limElRX = (elevRX>elevMinRX) & (elevRX<elevMaxRX);

    if azMinTX>azMaxTX
        limAzTX = ~((azTX>azMaxTX) & (azTX<azMinTX));
    else
        limAzTX = (azTX>azMinTX) & (azTX<azMaxTX);
    end

    if azMinRX>azMaxRX
        limAzRX = ~((azRX>azMaxRX) & (azRX<azMinRX));
    else
        limAzRX = (azTX>azMinRX) & (azTX<azMaxRX);
    end

    mask = ~(limElTX & limElRX & limAzTX & limAzRX);
    % mask = ~(limElTX & limElRX);

else
    mask = ones(1,length(elevRX));
end


%  _____           _                     __                  _   _                 
% /  __ \         | |                   / _|                | | (_)                
% | /  \/_   _ ___| |_ ___  _ __ ___   | |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
% | |   | | | / __| __/ _ \| '_ ` _ \  |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
% | \__/\ |_| \__ \ || (_) | | | | | | | | | |_| | | | | (__| |_| | (_) | | | \__ \
%  \____/\__,_|___/\__\___/|_| |_| |_| |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/


function output_txt = displayTime(event_obj, latTX, lonTX, elTX, latRX, lonRX, elRX)

% Get the clicked point's tag
    targetObj = get(event_obj.Target, 'Tag');

    % Case 1: Trajectory points (show corresponding time)
    if strcmp(targetObj, 'trajectory')
        index = event_obj.DataIndex;
        UData = event_obj.Target.UserData;
        output_txt = {['Time: ', char([datetime(UData(2,index),'ConvertFrom','datenum','TimeZone','Local','Format','dd-MMM HH:mm:ss')])], ...
                      ['Latitude: ', [num2str(event_obj.Position(1))] 'º'], ...
                      ['Longitude: ', [num2str(event_obj.Position(2))] 'º'], ...
                      ['AASL: ', [num2str(UData(1,index) / 10^3, '%.2f') ' Km']]};

    % Case 2: Current satellite position (show time = 0)
    elseif strcmp(targetObj, 'current')
        % targetID = find(strcmp({sats.Name}, erase(get(event_obj.Target ...
        % , 'DisplayName'),' position')));
        UData = event_obj.Target.UserData;
        index = UData(3);
        output_txt = {[erase(get(event_obj.Target, 'DisplayName'),' position')], ...
                      ['Time: ', char([datetime(UData(2),'ConvertFrom','datenum','TimeZone','Local','Format','dd-MMM HH:mm:ss')])], ...
                      ['Latitude: ', [num2str(event_obj.Position(1))] 'º'], ...
                      ['Longitude: ', [num2str(event_obj.Position(2))] 'º'], ...
                      ['AASL: ', [num2str(UData(1) / 10^3, '%.2f') ' Km']]};

    % Case 3: Static points (TX, RX) - Show lat/lon but hide time
    elseif strcmp(targetObj, 'TX')
        output_txt = {['Latitude: ', [num2str(latTX)] 'º'], ...
                      ['Longitude: ', [num2str(lonTX) 'º']], ...
                      ['AASL: ', [num2str(elTX)] ' m']};

    elseif strcmp(targetObj, 'RX')
        output_txt = {['Latitude: ', [num2str(latRX)] 'º'], ...
                      ['Longitude: ', [num2str(lonRX)] 'º'], ...
                      ['AASL: ', [num2str(elRX)] ' m']};

    % Default: Ignore other objects
    else
        output_txt = {};
    end
end

%  _____                 _         
% |  __ \               | |        
% | |  \/_ __ __ _ _ __ | |__  ___ 
% | | __| '__/ _` | '_ \| '_ \/ __|
% | |_\ \ | | (_| | |_) | | | \__ \
%  \____/_|  \__,_| .__/|_| |_|___/
%                 | |              
%                 |_|

leftX = 0;
rightX = 1100;
centerX = 550;
topY = 550;
botY = 50;

% Figure 2---------------------------------------------------------------------------------
% Calculate the distance for each satellite and time point
aasl = rlla(3,:)/10^3;
aasl(mask) = NaN;

fig2 = figure(2);
grid on
hold on
if filter
    plot(DTtime, rlla(3,:)/10^3, 'DisplayName',name,'LineStyle',':') % Plot vectorized data
end
plot(DTtime, aasl, 'DisplayName',name) % Plot vectorized data
% legend('Location', 'best')
hold off
title('Altitude above sea level')
ylabel('Magnitude [Km]')
xlabel('Time')
set(fig2, 'Position', [rightX, topY, 560, 420]);

% Figure 3---------------------------------------------------------------------------------
fig3 = figure(3);
rlla_t = rlla;
rlla_t(:,mask) = NaN; 
ax = geoaxes(fig3);
% geolimits([35 50], [-14 14])
geobasemap(ax,'darkwater');

if filter
    % traject = geoplot(ax, rlla(1,:), rlla(2,:),'DisplayName' ...
    %         , name, 'Tag', 'trajectory','LineStyle',':');
    traject = geoplot(ax, rlla(1,:), rlla(2,:), 'Tag', 'trajectory' ...
        ,'LineStyle',':', 'HandleVisibility', 'off');
    traject.UserData = [rlla(3,:);time];
end

hold(ax,'on')

% traject2 = geoplot(ax, rlla_t(1,:), rlla_t(2,:),'DisplayName' ...
%         , name, 'Tag', 'trajectory');
traject2 = geoplot(ax, rlla_t(1,:), rlla_t(2,:) ...
    , 'Tag', 'trajectory', 'HandleVisibility', 'off');
traject2.UserData = [rlla_t(3,:);time];
% Plot TX and RX positions
geoplot(ax,latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX');
if multistat
   geoplot(ax,RX1(1), RX1(2), '^r', 'DisplayName', 'RX', 'Tag', 'RX1');
   geoplot(ax,RX2(1), RX2(2), 'dr', 'DisplayName', 'RX', 'Tag', 'RX2'); 
else
   geoplot(ax,latRX, lonRX, '^r', 'DisplayName', 'RX', 'Tag', 'RX');
end

% Plot current satellite positions
index = find(time == datenum(inst));
curr = geoplot(ax, (rlla(1, index)), (rlla(2, index)), ...
    'ok','DisplayName', [name ' position'], 'Tag', ...
    'current', 'HandleVisibility', 'off');
curr.UserData = [rlla(3,index); time(index);index];

hold(ax,'off')
geolimits([35 50], [-14 14])

legend(ax, 'show', 'Location', 'southeast');

% Set up the data cursor mode with a custom function
dcm = datacursormode(fig3);

set(dcm, 'UpdateFcn', @(obj, event) displayTime(event, latTX, lonTX, ...
    elTX, latRX, lonRX, elRX));

title(ax,'Current and propagated position')
set(fig3, 'Position', [leftX, botY, 560, 420]);

% Figure 5---------------------------------------------------------------------------------
vel = zeros(1, length(time));

% Calculate the distance for each satellite and time point
for i = 1:length(time)
    vel(i) = norm(vecef(:,i));
end

figure(5);
grid on
hold on
vel2 = vel;
vel(mask) = NaN;
if filter
    plot(DTtime, vel2,'DisplayName',name,'LineStyle',':') % Plot vectorized data
end
plot(DTtime, vel,'DisplayName',name) % Plot vectorized data

% legend('Location', 'best')
hold off
title('Satellite speed')
ylabel('Magnitude [Km/s]')
xlabel('Time')

% Figure 8---------------------------------------------------------------------------------
fig8 = figure(8);

azRXp = azRX;
elRXp = elevRX;

azRXp(:,mask) = NaN;
elRXp(:,mask) = NaN;

yyaxis left
ylim([0 360])
hold on
if filter
    plot(DTtime,(azRX(1,:)),'DisplayName',name,'LineStyle',':')
end
plot(DTtime,azRXp,'LineStyle','-')
hold off
ylabel('Azimuth [º]')
yyaxis right
ylim([-40 90])
hold on
if filter
    plot(DTtime,(elevRX(1,:)),'DisplayName',name,'LineStyle',':')
end
plot(DTtime,elRXp,'LineStyle','-')
hold off
ylabel('Elevation [º]')

title('Position relative to RX')
xlabel('Time')
grid on
set(fig8, 'Position', [rightX, botY, 560, 420]);

% Figure 9---------------------------------------------------------------------------------
fig9 = figure(9);

azTXp = azTX;
elTXp = elevTX;

azTXp(:,mask) = NaN;
elTXp(:,mask) = NaN;

yyaxis left
ylim([0 360])
hold on
if filter
    plot(DTtime,(azTX(1,:)),'DisplayName',name,'LineStyle',':')
end
plot(DTtime',azTXp,'LineStyle','-')
hold off
ylabel('Azimuth [º]')
yyaxis right
ylim([-40 90])
hold on
if filter
    plot(DTtime,(elevTX(1,:)),'DisplayName',name,'LineStyle',':')
end
plot(DTtime',elTXp,'LineStyle','-')
hold off
ylabel('Elevation [º]')



title('Position relative to TX')
xlabel('Time [min]')

grid on
set(fig9, 'Position', [leftX, topY, 560, 420]);


% Figure 7---------------------------------------------------------------------------------
fig7 = figure(7);
subplot(1,2,1)
% yyaxis left
divDop = f_doppler;
divDop(mask) = NaN;

hold on
if filter
    plot(DTtime, f_doppler,'-','Tag','Left','DisplayName',name,'LineStyle',':');
end
plot(DTtime, divDop,'-');
hold off
ylabel('Doppler frequency [Hz]')
% yyaxis right
% ylim([min(min(f_doppler)) max(max(f_doppler))])
% ytickformat('%.4f') 
% ax2 = gca;
% axCol = 'black';
% ax2.YAxis(1).Color = axCol;
% ax2.YAxis(2).Color = axCol;
% ylim((ylim + freq)/10^6);
ylabel('Doppler frequency [Hz]')
grid on
sgtitle('Doppler')
xlabel('Time')
% legend('Location', 'best')

subplot(1,2,2)
divDop2 = f_doppler;
divDop2(mask) = NaN;
hold on
if filter
    plot(f_doppler,DTtime,'-','Tag','Right','DisplayName',name,'LineStyle',':');
end
plot(divDop2,DTtime,'-','Tag','Right','DisplayName',name);

if strcmp(fitterType,'real')
    scatter(x,xdata,'.')
    if fitter == true
        plot(fitted,DTtime,'-','Tag','Right','DisplayName',name,'LineStyle',':');
        legend('Original','Filtered','Received','Fitted')
    else
        legend('Original','Filtered','Received')
    end
end

hold off
grid on;
xlabel('Doppler Frequency [Hz]')
ylabel('Time')
ylim([inst-minutes(5) inst+minutes(5)])
ax3 = gca;
set(ax3, 'YDir','reverse') 
    
set(fig7, 'Position', [centerX-280, topY, 2*560, 420]);

% Figure 11---------------------------------------------------------------------------------
fig11 = figure(11);
snr2 = snr;
snr(mask) = NaN;
hold on
if filter 
plot(DTtime,snr2, 'LineStyle',':')
end
plot(DTtime,snr, 'LineStyle','-')
hold off
ylabel('SNR [dB]')
title('Input signal to noise ratio')
xlabel('Time')
% legend('Location','best')
grid on
set(fig11, 'Position', [centerX, botY, 560, 420]);

% ------------------------------------------------------------------------------------------
if fitter
    if strcmp(fitterType,'sim')
        if multistat
            MultiFigfit = figure;
            plot(DTtimeMult, ydataMult_nan)
            hold on
            plot(DTtimeMult,fittedMult)
            plot(DTtimeMult,f_dopplerMult)
            hold off
            grid on
            legend('Noisy','Fitted','Original')
            set(MultiFigfit, 'Position', [leftX, topY, 560, 420]);
            title('Multistatic Scenario')
        end
        if bistat

            BistFigfit = figure;
            plot(DTtime, ydata)
            hold on
            plot(DTtime,fitted)
            plot(DTtime,f_doppler)
            hold off
            grid on
            legend('Noisy','Fitted','Original')
            set(BistFigfit, 'Position', [rightX, topY, 560, 420]);
            title('Bistatic Scenario')
        end
    end

end

if ~multistat
    % Figure 10---------------------------------------------------------------------------------
    fig10 = figure(10);
    
    subplot(1,3,1)
    
    hold on
    if filter
        plot(DTtime(2:end), diff(f_doppler),'-','Tag','Left','DisplayName',name,'LineStyle',':');
    end
    plot(DTtime(2:end), diff(divDop),'-');
    hold off
    
    ax3 = gca;
    
    grid on;
    title('Doppler derivative')
    xlabel('Time')
    ylabel('Doppler frequency derivative')
    
    subplot(1,3,2)
    
    ct = 1e3;
    
    tint = @(var) sqrt(1./abs(diff(var)));
    
    hold on
    if filter
        plot(DTtime(2:end), tint(f_doppler).*ct,'-','Tag','Left','DisplayName',name,'LineStyle',':');
    end
    
    plot(DTtime(2:end), tint(divDop)*ct,'-');
    tintmin_i = find(tint(f_doppler) == min(tint(f_doppler)));
    tint_fdop = tint(f_doppler);
    plot(DTtime(tintmin_i), tint_fdop(tintmin_i).*ct,'ok');
    ax1 = gca;
    
    lims = [tint_fdop(tintmin_i)*0.2*ct tint_fdop(tintmin_i)*5*ct];
    ylim(lims);
    
    hold off
    
    grid on;
    title('Max integration time')
    xlabel('Time')
    ylabel('Integration time [ms]')
    
    subplot(1,3,3)
    
    BB_FS = fs/DF;
    
    fft_size = @(tau) floor(log2(tau.*BB_FS));
    
    hold on
    if filter
        plot(DTtime(2:end),fft_size(tint(f_doppler)),'LineStyle',':')
    end
    
    plot(DTtime(2:end),fft_size(tint(divDop)))
    
    fft_size_fdop = fft_size(tint(f_doppler));
    
    plot(DTtime(tintmin_i), fft_size_fdop(tintmin_i),'ok');
    ax2 = gca;
    
    hold off
    try
        ylim(fft_size(lims/ct))
    catch
    end
    grid on;
    title(sprintf('FFT Size (F_s = %1.2e Hz & DF = %i)',fs,DF))
    xlabel('Time')
    ylabel('n (FFT size = 2^n)')
    
    set(fig10, 'Position', [centerX-520, botY, 3*520, 420]);
    
    set(ax2, 'XLim', get(ax1, 'XLim'));
    linkaxes([ax3 ax1 ax2], 'x')
end

% Figure 12---------------------------------------------------------------------------------
% figure(12)
% plot(DTtime,bistaticVelocity)

%Figure 13---------------------------------------------------------------------------------
figure(13)
plot(DTtime,R1./1000)
hold on
plot(DTtime,R2./1000)
hold off
grid on
title('Bistatic Ranges')
xlabel('Time')
ylabel('Ditance [Km]')

end