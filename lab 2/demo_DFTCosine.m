%% ECE 311: Lab 2, demo 1: DFT spectra analysis
% Chao Ma, 09-12-2011
%

%% create a discrete cosine signal
W1 = 2*pi*3*50/16;        % frequency
T = 0.02;           % sampling period
   
N = 32;             % length of the signal

% create the cosine signal
n = 0:N-1;          
x = cos(W1*T*n);

% plot the signal
figure;
stem(n,x);
axis([-1,N,-1.1,1.1]);
xlabel('n','FontSize',18);
ylabel('x','FontSize',18);
title(['Cosine signal of length N = ',num2str(N)],'FontSize',18);
set(gca,'FontSize',14);

%% calculate the DFT of the signal
X = fft(x,128);
n = 0:128-1;   
% plot the magnitude of the DFT of the signal
figure;
plot(n,abs(X));
xlabel('k','FontSize',18);
ylabel('|x|','FontSize',18);
title(['The DFT of a cosine signal of length ', num2str(N)],'FontSize',18);
set(gca,'FontSize',14);



%% create a discrete cosine signal
W1 = 2*pi*3*50/16;        % frequency
T = 0.02;           % sampling period
   
N = 1024;             % length of the signal

% create the cosine signal
n = 0:N-1;          
x = cos(W1*T*n);

% plot the signal
figure;
plot(n,x);
axis([-1,N,-1.1,1.1]);
xlabel('n','FontSize',18);
ylabel('x','FontSize',18);
title(['Cosine signal of length N = ',num2str(N)],'FontSize',18);
set(gca,'FontSize',14);

% calculate the DFT of the signal
X = fft(x,128);
n = 0:128-1;   
% plot the magnitude of the DFT of the signal
figure;
plot(n,abs(X));
xlabel('k','FontSize',18);
ylabel('|x|','FontSize',18);
title(['The DFT of a cosine signal of length ', num2str(N)],'FontSize',18);
set(gca,'FontSize',14);


%% create a discrete cosine signal
W1 = 2*pi*3*50/16;        % frequency
T = 0.02;           % sampling period
   
N = 32;             % length of the signal

% create the cosine signal
n = 0:N-1;          
x = cos(W1*T*n).*hamming(N).';

% plot the signal
figure;
stem(n,x);
axis([-1,N,-1.1,1.1]);
xlabel('n','FontSize',18);
ylabel('x','FontSize',18);
title(['Cosine signal of length N = ',num2str(N)],'FontSize',18);
set(gca,'FontSize',14);

% calculate the DFT of the signal
X = fft(x,128);
n = 0:128-1;   
% plot the magnitude of the DFT of the signal
figure;
plot(n,abs(X));
xlabel('k','FontSize',18);
ylabel('|x|','FontSize',18);
title({['The DFT of a cosine signal of length ', num2str(N)]; 'with Hamming Window'},'FontSize',18);
set(gca,'FontSize',14);