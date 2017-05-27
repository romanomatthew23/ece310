%coDe for lab 3 
fignum = 1;
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 1                                    %
%                                                                        %
%------------------------------------------------------------------------%
% 1.	Tradeoffs between T and N
% HINT: This problem is similar to problem 3 of HW 3 from ECE 310
% Let xa(t) = cos(2?f1t) + 1.25cos(2?f2t) where f1 = 4Hz and f2 = 7Hz. In 
% this exercise, we investigate how our choices of the signal length N and 
% sampling period T affect the analysis of the spectrum. We will use a 
% 256-point DFT to get a fine sampling of Xd(?) for each choice of N and T.

% (a)	Compute and plot the magnitude of the 256-point DFT   for each of 
% the following cases:
% (i) N = 64, T = 1/30,
N = 64;
T = 1/30;
n = [0:N-1];
x = cos(2*pi*4*n*T) + 1.25*cos(2*pi*7*n*T);
X = fft(x,256);
h = [0:255];
figure(fignum);
stem(h,abs(X));
title('N = 64, T = 1/30');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

%(ii) N = 64, T = 1/6, 
N = 64;
T = 1/6;
n = [0:N-1];
x = cos(2*pi*4*n*T) + 1.25*cos(2*pi*7*n*T);
X = fft(x,256);
h = [0:255];
figure(fignum);
stem(h,abs(X));
title('N = 64, T = 1/6');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;


% (iii) N = 32, T = 1/30. 
N = 32;
T = 1/30;
n = [0:N-1];
x = cos(2*pi*4*n*T) + 1.25*cos(2*pi*7*n*T);
X = fft(x,256);
h = [0:255];
figure(fignum);
stem(h,abs(X));
title('N = 32, T = 1/30');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;


% In one case, aliasing is present - which is it?
%they are all 256 point DFTs so to obtain the frequency from the index k,
%we say wk = (2*pi*k/N)
% and w = wk/T
% so w = (2*pi*k/N/T)
%just look at the sampling period and the highest frequency.
% the highest frequency is 7Hz. Therefore, the nyquist rate is
% 14Hz. Theerefore, T needs to be at most 1/7.
% for part b T=1/6 therefore there will be aliasing



% (b)	For each of the cases in part (a), determine the analog frequencies
% corresponding to k = 61 and k = 161.
% HINT: Remember that the DFT is computed over the [0,2?] interval of the 
% DTFT!
%no MATLAB code for this

% (c)	Given that only N = 128 samples of xa(t) are to be acquired, how 
% would you choose T to best resolve the sinusoidal components? Find the 
% value of T which causes the main lobes to overlap by 1/2. Plot the
% result.
T = 1/384;
N = 128;
n = [0:N-1];
x = cos(2*pi*4*n*T) + 1.25*cos(2*pi*7*n*T);
X = fft(x,256);
h = [0:255];
figure(fignum);
stem(h,abs(X));
%plot(h,abs(X));
title('N = 128, T = 1/384');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

% (d)	Given that T = 1/30, find the N which causes the main lobes to 
% overlap by 1/2. Also find the minimum Nmin < N signal length at which 
% you can resolve these two sinusoids, via plotting the DFT in Matlab. 
% Plot both results.
T = 1/30;
N = 10;
n = [0:N-1];
x = cos(2*pi*4*n*T) + 1.25*cos(2*pi*7*n*T);
X = fft(x,256);
h = [0:255];
figure(fignum);
%stem(h,abs(X));
plot(h,abs(X));
title('N = 10, T = 1/30');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

T = 1/30;
N = 9;
n = [0:N-1];
x = cos(2*pi*4*n*T) + 1.25*cos(2*pi*7*n*T);
X = fft(x,256);
h = [0:255];
figure(fignum);
%stem(h,abs(X));
plot(h,abs(X));
title('N = 9, T = 1/30');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;


%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 2                                    %
%                                                                        %
%------------------------------------------------------------------------%

% 2.	Study sampling effects
% (a)	Consider a signal x1(t) = sinc2(f0t),f0 = 32Hz. What is the Nyquist
% sampling rate for x1(t)? Sample the signal x1(t) at 
% i) the Nyquist sampling rate, 
%N = 1024;
%n = [0:N-1];
%T = 1/64;
%f = 8;
%y = sinc(f*n*T)
%Y = fft(y);
%plot(n,abs(Y));

% twice the largest frequency right?
% so 64Hz?

%im having doubts whether or not the DFT i am displaying with the indices
% ranging from -256 to +256 or whatever the indices are
% it should go from 0 to some positive number and that should
% represent samples of the DTFT from [0,2pi]
% I was thinking I might need to use the command fftshift
% which shifts the values from the range [-pi,pi] to [0,2pi]

%if i sample the signal in the range t=[0,10] seconds
% and then do the DFT like I know how then compare the results
% i can figure out if im going wrong
% because the DFT will be right, its just that I'm plotting it against
% and arbitrary set of n's that I used to sample in the time domain
% i predict I was doing it wrong :) lets see




f0 = 32;
T=1/64;
N=640;
n=[1-N:N-1];
%n=[0:N-1];
x1 = (sinc(f0*n*T)).^2;
X1 = fft(x1);
% 
figure(fignum);
%plot(n(610:670),x1(610:670));
plot(n,x1);
title('x1 Sampled at Nyquist Rate');
ylabel('Amplitude');
xlabel('Index (n)');
fignum = fignum + 1;

k = [0:2*N-2];
figure(fignum);
plot(k,abs(X1));
title('DFT of x1 Sampled at Nyquist Rate');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

% ii) 3 times the Nyquist sampling rate, and 
T=1/192;
N=1920;
n=[1-N:N-1];
x1 = (sinc(f0*n*T)).^2;
X1 = fft(x1);
figure(fignum);
plot(n,x1);
title('x1 Sampled at 3 Times Nyquist Rate');
ylabel('Amplitude');
xlabel('Index (n)');
fignum = fignum + 1;

k = [0:2*N-2];
figure(fignum);
plot(k,abs(X1));
title('DFT of x1 Sampled at 3 TimesNyquist Rate');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

% iii) 0.8 times the Nyquist sampling rate, for |t| < 10. Plot the sampled signals. Plot 
% and compare the resulting DFT spectra.
T=1/(64*.8);
N=(64*.8)*10;
n=[1-N:N-1];
x1 = (sinc(f0*n*T)).^2;
X1 = fft(x1);
figure(fignum);
plot(n,x1);
title('x1 Sampled at 0.8 Times Nyquist Rate');
ylabel('Amplitude');
xlabel('Index (n)');
fignum = fignum + 1;

k = [0:2*N-2];
figure(fignum);
plot(k,abs(X1));
title('DFT of x1 Sampled at 0.8 Times Nyquist Rate');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

% HINT: Use the ‘sinc’ function as defined by Matlab, sinc(x) = sin(?x)/?x 
% for x 6= 0, and sinc(x) = 1 for x = 0. The time domain signal sinc(f0t) 
% = sin(?f0t)/(?f0t) will give you a frequency domain ‘rect’ function of 
% total width f0 (cutoff frequency f0/2). Recall that the fourier transform
% of the sinc2 function is a convolution of two rectangles.



% (b)	Modulate the signal by a complex exponential function:
% x2(t) = e?i2?f1tsinc2(f0t),f1 = 8Hz. What is the Nyquist sampling rate 
% for x2(t)? Sample the signal x2(t) at the Nyquist sampling rate of x2(t)
% and at the Nyquist sampling rate of x1(t), for |t| < 10. Plot the sampled
% signals. Plot and compare the resulting DFT spectra.
T = 1/80; %(32+8)*2
N = 800;
n = [1-N:N-1];
f1=8;
x2 = exp(-1*i*2*pi*f1*n*T).*(sinc(f0*n*T)).^2;
figure(fignum);
plot(n,abs(x2));
title('x2 Sampled at 80Hz');
ylabel('Amplitude');
xlabel('Index (n)');
fignum = fignum + 1;
%
X2 = fft(x2);
k = [0:2*N-2];
figure(fignum);
plot(k,abs(X2));
title('DFT of x2 Sampled at 80Hz');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;
%
T = 1/64; %(32)*2
N = 800;
n = [1-N:N-1];
f1=8;
x2 = exp(-1*i*2*pi*f1*n*T).*(sinc(f0*n*T)).^2;
figure(fignum);
plot(n,abs(x2));
title('x2 Sampled at 64Hz');
ylabel('Amplitude');
xlabel('Index (n)');
fignum = fignum + 1;
%
X2 = fft(x2);
k = [0:2*N-2];
figure(fignum);
plot(k,abs(X2));
title('DFT of x2 Sampled at 64Hz');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;



% (c)	Load the “rubberducky.jpg” image. Perform a 2D DFT to get the 
% spectrum of the image. Undersample the spectrum by a factor of 2 along 
% the i) rows, ii) columns, and iii) both rows and columns. Perform an 
% inverse 2D DFT on the undersampled spectra to reconstruct images. Display
% the reconstructed images. What do you observe and why?

L = imread('rubberducky.jpg');
Lk = fft2(L);
Lkxy = Lk(1:end, 1:2:end); 
kxy = ifft2(Lkxy);
figure(fignum); 
imagesc(abs(kxy)); 
title('Undersampled Rubber Ducky Along Rows');
colormap('gray'); 
axis image off;
fignum = fignum + 1;

L = imread('rubberducky.jpg');
Lk = fft2(L);
Lkxy = Lk(1:2:end, 1:end); 
kxy = ifft2(Lkxy);
figure(fignum); 
imagesc(abs(kxy)); 
title('Undersampled Rubber Ducky Along Columns');
colormap('gray'); 
axis image off;
fignum = fignum + 1;

L = imread('rubberducky.jpg');
Lk = fft2(L);
Lkxy = Lk(1:2:end, 1:2:end); 
kxy = ifft2(Lkxy);
figure(fignum); 
imagesc(abs(kxy)); 
title('Undersampled Rubber Ducky Along Rows and Columns');
colormap('gray'); 
axis image off;
fignum = fignum + 1;



%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 3                                    %
%                                                                        %
%------------------------------------------------------------------------%

% 3.	Window functions 
% (a)	In the time domain, plot the following
% windows of 32-point length: 
% (i) Rectangular window, 
t = [0:31];
w_r = ones(1,32);
figure(fignum)
stem(t,w_r);
title('Rectangular Window');
xlabel('Index (n)');
ylabel('Amplitude');
fignum = fignum + 1;

%(ii) Triangular window, 
t = [0:31];
w_t = zeros(1,32);
%0 to 16
for l = 0:16
    w_t(l+1) = 2*l/32;
end
%17 to 31
for l = 17:31
    w_t(l+1) = w_t(33-l)
end

figure(fignum)
stem(t,w_t);
title('Triangular Window');
xlabel('Index (n)');
ylabel('Amplitude');
fignum = fignum + 1;


% (iii) Hamming window, 
t = [0:31];
w_Hamm = zeros(1,32);
%0 to 16
for l=0:31
    w_Hamm(l+1) = .54 - .46*cos(2*pi*l/32);
end

figure(fignum)
stem(t,w_Hamm);
title('Hamming Window');
xlabel('Index (n)');
ylabel('Amplitude');
fignum = fignum + 1;


% (iv) Hanning (Hann) window, 
t = [0:31];
w_Hann = zeros(1,32);
%0 to 16
for l=0:31
    w_Hann(l+1) = .5 -.5*cos(2*pi*l/32);
end

figure(fignum)
stem(t,w_Hann);
title('Hanning Window');
xlabel('Index (n)');
ylabel('Amplitude');
fignum = fignum + 1;


% (v) Kaiser window. 
t = [0:31];
w_k = kaiser(32,3);

figure(fignum)
stem(t,w_k);
title('Kaiser Window');
xlabel('Index (n)');
ylabel('Amplitude');
fignum = fignum + 1;


% (b)	Plot the 256-point DFT of the above 32-point windows. More
% specifically, plot the magnitude of the DFT spectra of the windows in
% decibels (dB). HINT: First normalize the DFT magnitude spectrum to 1,
% then display the spectrum in dB using 20*log10(X). Also note that
% 20*log10(0) = -Inf. Typically Matlab will ignore infinite values when
% plotting. You can add a very small non-zero value to the magnitude
% spectrum to avoid this, but make sure it does not significantly affect
% the side lobe height (try eps which matlab recognizes as the number
% 2?52). 
W_r = fft(w_r,256);
W_r = abs(W_r)./max(W_r);            %normalize to 1
W_r = 20*log10(W_r + eps);            %dB
k = [0:255];
figure(fignum);
plot(W_r);
title('DFT of Rectangular Window');
xlabel('Index (k)');
ylabel('Magnitude (dB)');
fignum = fignum + 1;

W_t = fft(w_t,256);
%W_t = fftshift(W_t);
W_t = abs(W_t)./max(W_t);            %normalize to 1
W_t = 20*log10(W_t + eps);            %dB
k = [0:255];
figure(fignum);
plot(W_t);
title('DFT of Triangular Window');
xlabel('Index (k)');
ylabel('Magnitude (dB)');
fignum = fignum + 1;

W_Hamm = fft(w_Hamm,256);
W_Hamm = abs(W_Hamm)./max(W_Hamm);            %normalize to 1
W_Hamm = 20*log10(W_Hamm + eps);            %dB
k = [0:255];
figure(fignum);
plot(W_Hamm);
title('DFT of Hamming Window');
xlabel('Index (k)');
ylabel('Magnitude  (dB)');
fignum = fignum + 1;


W_Hann = fft(w_Hann,256);
W_Hann = abs(W_Hann)./max(W_Hann);            %normalize to 1
W_Hann = 20*log10(W_Hann + eps);            %dB
k = [0:255];
figure(fignum);
plot(W_Hann);
title('DFT of Hanning Window');
xlabel('Index (k)');
ylabel('Magnitude (dB)');
fignum = fignum + 1;

W_k = fft(w_k,256);
W_k = abs(W_k)./max(W_k);            %normalize to 1
%W_k = 20*log10(W_k);            %dB
W_k = 20*log10(W_k + eps);            %dB
k = [0:255];
figure(fignum);
plot(W_k);
title('DFT of Kaiser Window');
xlabel('Index (k)');
ylabel('Magnitude (dB)');
fignum = fignum + 1;

% (c)	Use the plot trace tool to determine the main lobe width
% and the magnitude of the highest side lobe for each of the above windows.

%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 4                                    %
%                                                                        %
%------------------------------------------------------------------------%
% 
% 4.	Study windowing effects: 1D case 
% (a)	Download the signal.mat
% file from the course website, where x is an unknown signal of length 32.
% It is known that the sampling period is 0.02 seconds, which satisfies the
% Nyquist criterion, and the signal contains 5 frequency components, two of
% which have much weaker magnitude compared to the other three components.
% Plot the magnitude DFT spectrum using each of the 5 windows from the
% previous problem to compare windowing effects. Find the five frequencies
% present in the signal. 
load signal.mat;
x_r = x.*w_r';
x_t = x.*w_t';
x_Hamm = x.*w_Hamm';
x_Hann = x.*w_Hann';
x_k = x.*w_k;

N = 2^10;
X_r = fft(x_r,N);
X_t = fft(x_t,N);
X_Hamm = fft(x_Hamm,N);
X_Hann = fft(x_Hann,N);
X_k = fft(x_k,N);

k = [0:N-1];
figure(fignum);
plot(k,abs(X_r));
title('DFT of x using Rectangular Window');
xlabel('Index (k)');
ylabel('Magnitude');
fignum = fignum + 1;

figure(fignum);
plot(k,abs(X_t));
title('DFT of x using Triangular Window');
xlabel('Index (k)');
ylabel('Magnitude');
fignum = fignum + 1;

figure(fignum);
plot(k,abs(X_Hamm));
title('DFT of x using Hamming Window');
xlabel('Index (k)');
ylabel('Magnitude');
fignum = fignum + 1;

figure(fignum);
plot(k,abs(X_Hann));
title('DFT of x using Hanning Window');
xlabel('Index (k)');
ylabel('Magnitude');
fignum = fignum + 1;

figure(fignum);
plot(k,abs(X_k));
title('DFT of x using Kaiser Window');
xlabel('Index (k)');
ylabel('Magnitude');
fignum = fignum + 1;



% (b)	Thus far in ECE 311, we have typically
% considered signals with frequency components that do not change over time
% (even in the telephone problem, each dialed number had two constant
% tones). Consider a signal with two time-dependent, logarithmic chirp
% frequency components, given below. 
t = 0:0.001:6; Fs = 1000; % 6 seconds @ 1kHz sample rate 
% make two logarithmic chirps offset by 20 Hz 
fo = 45;
f1 = 450; offset = 20;
ya = chirp(t,fo,6,f1,'logarithmic'); 
yb = chirp(t,fo+offset,6,f1+offset,'logarithmic'); y = ya + yb; 
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024
figure(fignum);
spectrogram(y,L,128,nfft,Fs,'yaxis');% view chirps 
soundsc(y,Fs) % listen to chirps 
fignum = fignum + 1;
% The spectrogram windows the time signal with a length L=256
% hamming window starting every 128 samples. For each windowed signal, a
% 1024-point DFT is taken, and the magnitude spectrum (one-sided) is shown
% as a heat map; each DFT is given as a vertical slice of the spectrogram.

% Show the spectrogram result in your lab report. Then, consider one short
% segment of the time-domain signal beginning at t=4 seconds. Plot the
% magnitude 1024-point DFT spectrum for multiple window lengths, including
% L=256, L=512 and L=1024. What happens to the magnitude spectrum as you
% increase the window length? Why? This is why the ‘short-time’ fourier
% transform is often critical in signal processing!

% y has 6001 points, starting at 4 seconds means starting at index 4000.

y_4 = y(4000:5023);
%256 length window
L = 256;
k = [0:1023];
y_4w = y_4(1:L)'.*hamming(L);
Y_4w = fft(y_4w,1024);
figure(fignum);
plot(k,abs(Y_4w));
title('Short Time FFT Using 256 Length Window');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

L = 512;
k = [0:1023];
y_4w = y_4(1:L)'.*hamming(L);
Y_4w = fft(y_4w,1024);
figure(fignum);
plot(k,abs(Y_4w));
title('Short Time FFT Using 512 Length Window');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;

L = 1024;
k = [0:1023];
y_4w = y_4(1:L)'.*hamming(L);
Y_4w = fft(y_4w,1024);
figure(fignum);
plot(k,abs(Y_4w));
title('Short Time FFT Using 1024 Length Window');
ylabel('Magnitude');
xlabel('Index (k)');
fignum = fignum + 1;
%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 5                                    %
%                                                                        %
%------------------------------------------------------------------------%

% 5. Study windowing effects: 2D case 
% (a)	Generate a 256×256 image using
% the command I = phantom(256). Display the image in gray scale. 
I = phantom(256);
figure(fignum);
colormap(gray(256));
imagesc(I);
title('Gray Scale Phantom');
fignum = fignum + 1;

% (b)
% Perform 2D DFT to get the spectrum of the image using the command s =
% fftshift(fft2(I)). Truncate the spectrum s by taking only the central or
% lower frequency parts of the spectrum (64 ×64 out of 256×256), and
% zero-padding the rest. Reconstruct and display the resulting image from
% the truncated spectrum. What is your observation and why? HINT: Use
% ifft2(ifftshift(st)), where st is the truncated spectrum, to reconstruct
% the image. 
s = fftshift(fft2(I));

%truncate the magnitude spectrum
MagSpec = abs(s);
NewMagSpec = zeros(256,256);
NewMagSpec(1:64,1:64) = MagSpec(1:64,1:64);
%reconstruct the image
PhSpec = angle(s);
ImRecon = ifft2(ifftshift(NewMagSpec.*exp(1i*PhSpec)));
figure(fignum);
imagesc(abs(ImRecon));
title('Reconstructed Image After Truncation Using Rectangular Window');
colormap(gray(256));
axis image off;
fignum = fignum + 1;


% (c)	Apply a Hamming window of length 64 along 
%i) row,
HammRowMagSpec = zeros(256,256);
for l = 1:64;
    HammRowMagSpec(l,1:64) = hamming(64).*NewMagSpec(l,1:64)';
end

%ii) column, and 
HammColMagSpec = zeros(256,256);
for m = 1:64;
    HammColMagSpec(1:64,m) = hamming(64).*NewMagSpec(1:64,m);
end

% iii) both row and column of the truncated spectrum in b.
HammBothMagSpec = zeros(256,256);
for m = 1:64;
    HammBothMagSpec(1:64,m) = hamming(64).*NewMagSpec(1:64,m);
end
for l = 1:64;
    HammBothMagSpec(l,1:64) = hamming(64).*HammBothMagSpec(l,1:64)';
end

%reconstruct the images
PhSpec = angle(s);
ImRecon = ifft2(ifftshift(HammRowMagSpec.*exp(1i*PhSpec)));
figure(fignum);
imagesc(abs(ImRecon));
title('Reconstructed Image After Truncation By Row Hamming Window');
colormap(gray(256));
axis image off;
fignum = fignum + 1;

PhSpec = angle(s);
ImRecon = ifft2(ifftshift(HammColMagSpec.*exp(1i*PhSpec)));
figure(fignum);
imagesc(abs(ImRecon));
title('Reconstructed Image After Truncation By Column Hamming Window');
colormap(gray(256));
axis image off;
fignum = fignum + 1;

PhSpec = angle(s);
ImRecon = ifft2(ifftshift(HammBothMagSpec.*exp(1i*PhSpec)));
figure(fignum);
imagesc(abs(ImRecon));
title('Reconstructed Image After Truncation By Row and Column Hamming Window');
colormap(gray(256));
axis image off;
fignum = fignum + 1;


% Reconstruct and display the resulting images from the spectra. How does
% this differ from using the rectangular window in part (b), and why?

