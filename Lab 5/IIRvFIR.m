
close all
clear all
% FIR filter
rp = 1;           % Passband ripple
rs = 80;          % Stopband ripple
fs = 2000;        % Sampling frequency
f = [400 600];    % Cutoff frequencies
a = [1 0];        % Desired amplitudes
% Compute deviations
dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)]; 
% Estimate the filter order
[n,fo,ao,w] = firpmord(f,a,dev,fs);
% design the filter
b1 = firpm(n,fo,ao,w);
a1 = 1;

% IIR filter
Wp = 0.4;
Ws = 0.6;
Rp = 1;
Rs = 80;
[n,Wp]=ellipord(Wp,Ws,Rp,Rs);
[b2,a2] = ellip(n,Rp,Rs,Wp);

% input signal
n = 0:127;
x1 = cos(n*0.2*pi);
x2 = cos(n*0.4*pi);

y11 = filter(b1,1,x1);
y12 = filter(b2,a2,x1);

y21 = filter(b1,1,x2);
y22 = filter(b2,a2,x2);

figure; plot(n, x1,'--r','LineWidth',2); xlabel('n'); title('');
hold all; plot(n, y11,'b','LineWidth',2); 
xlabel('n','FontSize',14);
title('x_1[n] and y_{11}[n]','FontSize',14);
axis([0,128,-1.5,1.5]);
set(gca,'FontSize',14);

figure; plot(n, x1,'--r','LineWidth',2); xlabel('n'); title('');
hold all; plot(n, y12,'k','LineWidth',2); 
xlabel('n','FontSize',14);
title('x_1[n] and y_{12}[n]','FontSize',14);
axis([0,128,-1.5,1.5]);
set(gca,'FontSize',14);

figure; plot(n, x2,'--r','LineWidth',2); xlabel('n'); title('');
hold all; plot(n, y21,'b','LineWidth',2); 
xlabel('n','FontSize',14);
title('x_2[n] and y_{21}[n]','FontSize',14);
axis([0,128,-1.5,1.5]);
set(gca,'FontSize',14);

figure; plot(n, x2,'--r','LineWidth',2); xlabel('n'); title('');
hold all; plot(n, y22,'k','LineWidth',2); 
xlabel('n','FontSize',14);
title('x_2[n] and y_{22}[n]','FontSize',14);
axis([0,128,-1.5,1.5]);
set(gca,'FontSize',14);