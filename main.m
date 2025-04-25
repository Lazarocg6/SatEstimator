clc
clear

inst = datetime('now', 'TimeZone', 'Local'); %Time origin
% inst = datetime('24-Apr 07:47:22','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 60; % In minutes, 'duracion' before and after the current time
precision = 10 / 60; % Precision in minutes

[time, DTtime] = initTimes(inst, duracion, precision);

[ydata, epochd, epoch_og] = noiseGen(time,25544);

fun = @(x, xdata) paramsToDop(x, xdata);

x0 = [epoch_og];

res = lsqcurvefit(@paramsToDop,x0,time,ydata);
fprintf('Epoch diff = %s, noisy ->%10.10f, fitted ->%10.10f\n', (epochToUTC(epochd)-epochToUTC(res)),epochd,res)

fitted = paramsToDop(res,time);
og = paramsToDop(NaN,time);

figure(1)
plot(DTtime, ydata)
hold on
plot(DTtime,fitted)
plot(DTtime,og)
hold off
grid on
legend('Noisy','Fitted','Original')


