%% load data
[x,Fs] = wavread('Overture01.wav');
soundsc(x,Fs);
pause
%% pass the data through an LSI system, i.e. low-pass filter
% display the time-domain signal
figure;
plot(x);
title('Time-domain signal, x[n]');
axis tight;

% calculate and display the DFT spectrum
Xk = fftshift(fft(x));
f = linspace(-Fs/2,Fs/2,length(x));
figure;
plot(f,abs(Xk));
axis([-3000,3000,0,120]);
xlabel('Frequency (Hz)');
title('DFT spectrum, X[k]')

% design a low-pass filter
d = fdesign.highpass('Fst,Fp,Ast,Ap',300,350,40,1,Fs);
hd = design(d,'equiripple');
b = hd.Numerator;
a = 1;
xFiltered = filter(b,a,x);
soundsc(xFiltered,Fs);

% frequency response
freqz(b,a,1024);

% calculate and display the DFT spectrum of the filtered signal
XkFiltered = fftshift(fft(xFiltered));
f = linspace(-Fs/2,Fs/2,length(xFiltered));
figure;
plot(f,abs(XkFiltered));
axis([-3000,3000,0,120]);
xlabel('Frequency (Hz)');
title('DFT spectrum, X[k]')

%% pass the data through a non-LSI system
w0 = 2*pi*700;
t = (0:length(x)-1)/Fs;
temp = cos(w0*t);
mx=zeros(size(temp));
mx(temp>0) = 1;
mx = mx(:);
xModulated = x.*mx;
soundsc(xModulated,Fs);

% display the time-domain signal
figure;
plot(xModulated);
title('Modulated signal xModulated[n]');
axis tight;

% calculate and display the DFT spectrum
XkModulated = fftshift(fft(xModulated));
f = linspace(-Fs/2,Fs/2,length(xModulated));
figure;
plot(f,abs(XkModulated));
xlabel('Frequency (Hz)');
title('DFT spectrum, X[k]')