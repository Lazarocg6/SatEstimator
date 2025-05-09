close all
clc
clear
load("040520250129.mat")
f = 143.05e6;
t = 0:1/(4*f):1e4/f;
faxis = linspace(-2*f,2*f,length(t));
seno = sin(2*pi*f*t);
f_dom = fft(seno);
% figure
% plot(faxis,abs(f_dom))
fIQ = fft(IQ')';
conv = fIQ.*f_dom;
figure
plot(conv)
