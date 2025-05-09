clear
clc
close all
load('cropped.mat')
figure(1)
imagesc(20.*log10(abs(cropped)));
clim([77 85])


orden = 2;

w = 0.0002;

T = 15;
cte = 65;

pesos = [w (1-w).*exp(-0.3.*(0:orden-2))/(sum(exp(-(0:orden-2))));];
filtered = zeros(size(cropped));
for i = 1:size(cropped,2)

    if i<orden
        filtered(:,i) = 0;
    else
        filtered(:,i) = (cropped(:,i).*pesos(1));
        for j = 2:orden
            filtered(:,i) = filtered(:,i)+(filtered(:,i-j+1).*pesos(j));
        end
    end

end

figure(2)
imagesc(20.*log10(abs(filtered)));
% clim([77 85])

figure(3)
plot(20.*log10(abs(cropped(26,:))))
grid on
hold on
plot(20.*log10(abs(filtered(26,:)))+cte)
hold off