clc
clear
close all

inst = datetime('now', 'TimeZone', 'Local'); %Time origin
inst = datetime('3-Jun 16:43:26','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 60; % In minutes
precision = 30 / 60; % Precision in minutes
propagateB4andafter = true; % False propaga hacia delante desde inst

[time, DTtime] = initTimes(inst, duracion, precision,propagateB4andafter);

ID = 25544;

updateEOP(36); % Update EOP 
updateTLE(1, ID); % Update TLE 

n = 100;

RX1 = [40.45206046037957, -3.726407299669201, 670];
% RX2 = [61.45206046037957, -3.726407299669201, 670];

seed = zeros(1,n);

Bist1 = zeros(6,n);
Bist2 = zeros(6,n);
Bist3 = zeros(6,n);
Bist4 = zeros(6,n);

tic;
for lat = 0:5:20
    RX2 = [RX1(1)+lat RX1(2:end)];
    for i = 1:n 
    
        % Caso bistatico receptor1
        fprintf('Bistatic1=%i------------------------------------------------Lat=%i\n',i,lat)
    
        seed(1,i) = randi(2^32,1);
    
        [ydata, epochd,epoch_og,f_dev] = noiseGen(time,ID,RX1,seed(1,i)); 
        x0 = [epoch_og,0];
        fun = @(optimize,time) paramsToDop(optimize(1),time,RX1)+optimize(2);
        res = lsqcurvefit(fun,x0,time,ydata);
        t = (epochToUTC(epochd)-epochToUTC(res(1)));
        t.Format = 'hh:mm:ss.SSS';
        Bist1(1,i) = epochd; %Ground truth
        Bist1(2,i) = res(1); %Fitted epoch
        Bist1(3,i) = f_dev; %SDR freq deviation
        Bist1(4,i) = res(2); %Fitted freq deviation
        Bist1(5,i) = seconds(t); %Diferencia
        Bist1(6,i) = res(2)-f_dev;
    
        % Caso bistatico receptor2
        fprintf('Bistatic2=%i------------------------------------------------Lat=%i\n',i,lat)
        [ydata, epochd,epoch_og,f_dev] = noiseGen(time,ID,RX2,seed(1,i)); 
        x0 = [epoch_og,0];
        fun = @(optimize,time) paramsToDop(optimize(1),time,RX2)+optimize(2);
        res = lsqcurvefit(fun,x0,time,ydata);
        t = (epochToUTC(epochd)-epochToUTC(res(1)));
        t.Format = 'hh:mm:ss.SSS';
        Bist2(1,i) = epochd; %Ground truth
        Bist2(2,i) = res(1); %Fitted epoch
        Bist2(3,i) = f_dev; %SDR freq deviation
        Bist2(4,i) = res(2); %Fitted freq deviation
        Bist2(5,i) = seconds(t); %Epoch diff
        Bist2(6,i) = res(2)-f_dev; %Freq deviation diff
    
        %Caso multiestatico
        fprintf('Multistatic=%i------------------------------------------------Lat=%i\n',i,lat)
        [ydataMult, epochdMult,epoch_ogMult,f_devMult] = noiseGen(time,ID,RX1,seed(1,i),RX2);
        x0Mult = [epoch_ogMult,0];
        funMult = @(optimize,time) paramsToDop(optimize(1),time,RX1,RX2)+optimize(2);
        resMult = lsqcurvefit(funMult,x0Mult,time,ydataMult);
        tMult = (epochToUTC(epochdMult)-epochToUTC(resMult(1)));
        tMult.Format = 'hh:mm:ss.SSS';
        Bist4(1,i) = epochdMult; %Ground truth
        Bist4(2,i) = resMult(1); %Fitted epoch
        Bist4(3,i) = f_devMult; %SDR freq deviation
        Bist4(4,i) = resMult(2); %Fitted freq deviation
        Bist4(5,i) = seconds(tMult);
        Bist4(6,i) = resMult(2)-f_devMult;
    
    end
    
    %Media entre ambos
    
    Bist3(1,:) = Bist2(1,:);
    Bist3(2,:) = (Bist1(2,:)+Bist2(2,:))/2;
    Bist3(3,:) = Bist2(3,:);
    Bist3(4,:) = (Bist1(4,:)+Bist2(4,:))/2;
    Bist3(5,:) = seconds(epochToUTC(Bist3(1,:))-epochToUTC(Bist3(2,:)));
    % Bist3(5,:) = (Bist1(5,:)+Bist2(5,:))./2;
    Bist3(6,:) = Bist3(4,:)-Bist3(3,:);
    
    %------------------
    
    dur = toc;
    
    fprintf('Total duration = %fs\n',dur)
    
    avgBist1EpochDiffLat20 = sqrt(mean(Bist1(5,:).^2,2));%Solo RX1
    avgBist2EpochDiffLat20 = sqrt(mean(Bist2(5,:).^2,2));%Solo RX2
    avgBist3EpochDiffLat20 = sqrt(mean(Bist3(5,:).^2,2));%Media entre ambos
    avgBist4EpochDiffLat20 = sqrt(mean(Bist4(5,:).^2,2));%Ajuste ambos
    
    avgBist1FreqDiffLat20 = sqrt(mean(Bist1(6,:).^2,2));
    avgBist2FreqDiffLat20 = sqrt(mean(Bist2(6,:).^2,2));
    avgBist3FreqDiffLat20 = sqrt(mean(Bist3(6,:).^2,2));
    avgBist4FreqDiffLat20 = sqrt(mean(Bist4(6,:).^2,2));

save(['ResultsMontecarlo_' char(datetime('now','Format','ddMMyyyyHHmm')) 'Lat_',num2str(lat), '.mat'],'Bist1','Bist2','Bist3','Bist4')
end
% figure
% histogram(Bist1(5,:),20)
% figure
% histogram(Bist2(5,:),20)
% figure
% histogram(Bist3(5,:),20)
% figure
% histogram(Bist4(5,:),20)

% bists = 1:4;
% lats = 0:5:20;
% 
% figure; hold on;
% 
% for b = bists
%     y = [];
%     for lat = lats
%         varname = sprintf('avgBist%dFreqDiffLat%d', b, lat);
%         if evalin('base', sprintf('exist(''%s'', ''var'')', varname))
%             val = evalin('base', varname);
%             y(end+1) = mean(val);  % Change to suit your data (mean if vector)
%         else
%             warning('Variable %s does not exist.', varname);
%             y(end+1) = NaN;
%         end
%     end
%     plot(lats, y, '-o', 'DisplayName', sprintf('Bist %d', b));
% end
% 
% xlabel('Latitude');
% ylabel('Value');
% legend show;
% title('avgBist Epoch Diff vs Latitude');
% grid on;
