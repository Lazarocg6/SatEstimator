function graphs(time,DTtime,fitter,inst,filter,freq)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all

[f_doppler, recef, vecef, rlla, bistaticRange, bistaticVelocity, ...
    R1, R2, snr, name, latTX, latRX, lonTX, lonRX, elTX, elRX] = paramsToDop(NaN,time);

% ______                             
% | ___ \                            
% | |_/ /_ _ _ __ __ _ _ __ ___  ___ 
% |  __/ _` | '__/ _` | '_ ` _ \/ __|
% | | | (_| | | | (_| | | | | | \__ \
% \_|  \__,_|_|  \__,_|_| |_| |_|___/

fs = 5e6; % Sample rate
DF = 170; % Decimation factor

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
elevMinRX = 10; % degrees
elevMaxRX = 50;
azMinRX = 270;
azMaxRX = 90;

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
 
else
    mask = ones(1,length(elevRX));
end

if fitter

    [ydata, epochd, epoch_og] = noiseGen(time,25544);
    
    x0 = [epoch_og];
    
    res = lsqcurvefit(@paramsToDop,x0,time,ydata);
    fprintf('Epoch diff = %s, noisy ->%10.10f, fitted ->%10.10f\n', (epochToUTC(epochd)-epochToUTC(res)),epochd,res)
    
    fitted = paramsToDop(res,time);
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
    traject = geoplot(ax, rlla(1,:), rlla(2,:),'DisplayName' ...
            , name, 'Tag', 'trajectory','LineStyle',':');
    traject.UserData = [rlla(3,:);time];
end

hold(ax,'on')

traject2 = geoplot(ax, rlla_t(1,:), rlla_t(2,:),'DisplayName' ...
        , name, 'Tag', 'trajectory');
traject2.UserData = [rlla_t(3,:);time];
% Plot TX and RX positions
geoplot(ax,latTX, lonTX, 'xr', 'DisplayName', 'TX', 'Tag', 'TX');
geoplot(ax,latRX, lonRX, '^r', 'DisplayName', 'RX', 'Tag', 'RX');

% Plot current satellite positions
index = find(time == datenum(inst));
curr = geoplot(ax, (rlla(1, index)), (rlla(2, index)), ...
    'ok','DisplayName', [name ' position'], 'Tag', ...
    'current', 'HandleVisibility', 'off');
curr.UserData = [rlla(3,index); time(index);index];

hold(ax,'off')
geolimits([35 50], [-14 14])

% legend(ax, 'show', 'Location', 'southeast');

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
hold on
if filter
    plot(DTtime,(azRX(1,:)),'DisplayName',name,'LineStyle',':')
end
plot(DTtime,azRXp,'LineStyle','-')
hold off
ylabel('Azimuth [º]')
yyaxis right
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
hold on
if filter
    plot(DTtime,(azTX(1,:)),'DisplayName',name,'LineStyle',':')
end
plot(DTtime',azTXp,'LineStyle','-')
hold off
ylabel('Azimuth [º]')
yyaxis right
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

if true
    load('detections130520250037.mat','detections_out','f_axis','t_axis')
    [fila,columna] = find(detections_out);
    x = f_axis(columna);
    y = t_axis(fila);
    scatter(x,y,'.')
end

hold off
grid on;
xlabel('Doppler Frequency [Hz]')
ylabel('Time')
ylim([inst-minutes(5) inst+minutes(5)])
ax3 = gca;
set(ax3, 'YDir','reverse')
xlabel('Frequency [MHz]') 
    
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

if fitter
    figfit = figure;
    plot(DTtime, ydata)
    hold on
    plot(DTtime,fitted)
    plot(DTtime,f_doppler)
    hold off
    grid on
    legend('Noisy','Fitted','Original')
    set(figfit, 'Position', [leftX, topY, 560, 420]);
end

% Figure 10---------------------------------------------------------------------------------
fig10 = figure(10);

subplot(1,3,1)

hold on
if filter
    plot(DTtime(2:end), diff(f_doppler),'-','Tag','Left','DisplayName',name,'LineStyle',':');
end
plot(DTtime(2:end), diff(divDop),'-');
hold off

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
linkaxes([ax1 ax2], 'x')

% % Figure 12---------------------------------------------------------------------------------
% figure(12)
% plot(DTtime,bistaticVelocity)

end