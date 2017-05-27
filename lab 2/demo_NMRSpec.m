%% ECE 311 class demonstration, NMR spectrum from a single point
% Chao Ma
%
clear all; 
close all;
clc;

fignum = 600;
load NMRSpec.mat;

%% 128 point data
st1 = st(1:128);
% plot time domain signal
figure(fignum);
plot((0:length(st1)-1)*dt*1e3,abs(st1));
xlabel('Time (ms)');
ylabel('NMR signal');
title(['Time domain NMR signal: ',num2str(length(st1)), ' points']);
fignum = fignum + 1;

% DFT spectrum
sf1 = fftshift(fft(st1));
f = linspace(-BW/2,BW/2,length(sf1));
f = f(:);
figure(fignum);
plot(f,abs(sf1));
xlabel('Frequency (Hz)');
ylabel('NMR spectrum (magnitude)');
title(['Frequency domain NMR signal: ',num2str(length(sf1)), ' points spectrum']);
fignum = fignum + 1;

% DFT spectrum after zero-padding
sf1ZP = fftshift(fft(st1,1024));
f = linspace(-BW/2,BW/2,length(sf1ZP));
f = f(:);
figure(fignum);
plot(f,abs(sf1ZP));
xlabel('Frequency (Hz)');
ylabel('NMR spectrum (magnitude)');
title(['Frequency domain NMR signal: ',num2str(length(sf1ZP)), ' points spectrum']);
fignum = fignum + 1;
