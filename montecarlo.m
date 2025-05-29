clc
clear
close all

inst = datetime('now', 'TimeZone', 'Local'); %Time origin
inst = datetime('28-May 13:28:26','InputFormat','dd-MMM HH:mm:ss',TimeZone='Local');
duracion = 60; % In minutes
precision = 30 / 60; % Precision in minutes
propagateB4andafter = true; % False propaga hacia delante desde inst

[time, DTtime] = initTimes(inst, duracion, precision,propagateB4andafter);

ID = 25544;

updateEOP(36); % Update EOP 
updateTLE(1, ID); % Update TLE 

n = 100;

RX1 = [40.45206046037957, -3.726407299669201, 670];
RX2 = [61.45206046037957, -3.726407299669201, 670];

seed = zeros(1,n);

Bist1 = zeros(6,n);
Bist2 = zeros(6,n);
Bist3 = zeros(6,n);
Bist4 = zeros(6,n);

tic;
for i = 1:n % Caso bistatico receptor1

    fprintf('Bistatic1=%i------------------------------------------------\n',i)

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

end

for i = 1:n % Caso bistatico receptor2

    fprintf('Bistatic2=%i------------------------------------------------\n',i)
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
    
end

%Media entre ambos

Bist3(1,:) = Bist2(1,:);
Bist3(2,:) = (Bist1(2,:)+Bist2(2,:))/2;
Bist3(3,:) = Bist2(3,:);
Bist3(4,:) = (Bist1(4,:)+Bist2(4,:))/2;
Bist3(5,:) = seconds(epochToUTC(Bist3(1,:))-epochToUTC(Bist3(2,:)));
Bist3(6,:) = Bist3(4,:)-Bist3(3,:);

%------------------

for i = 1:n

    fprintf('Multistatic=%i------------------------------------------------\n',i)
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

dur = toc;

fprintf('Total duration = %fs\n',dur)

avgBist1EpochDiff = mean(Bist1(5,:),2);%Solo RX1
avgBist2EpochDiff = mean(Bist2(5,:),2);%Solo RX2
avgBist3EpochDiff = mean(Bist3(5,:),2);%Media entre ambos
avgBist4EpochDiff = mean(Bist4(5,:),2);%Ajuste ambos

avgBist1FreqDiff = mean(Bist1(6,:),2);
avgBist2FreqDiff = mean(Bist2(6,:),2);
avgBist3FreqDiff = mean(Bist3(6,:),2);
avgBist4FreqDiff = mean(Bist4(6,:),2);

% save('ResultsN100','avgBist1EpochDiff','avgBist2EpochDiff','avgBist3EpochDiff','avgBist4EpochDiff','avgBist1FreqDiff','avgBist2FreqDiff','avgBist3FreqDiff','avgBist4FreqDiff')

figure
histogram(Bist1(5,:),20)
figure
histogram(Bist2(5,:),20)
figure
histogram(Bist3(5,:),20)
figure
histogram(Bist4(5,:),20)