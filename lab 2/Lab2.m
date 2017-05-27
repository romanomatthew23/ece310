%Lab 2!!!!!

%------------------------------------------------------------------------%
%                                                                        %
%                              Part 1                                    %
%                                                                        %
%------------------------------------------------------------------------%
%Write a function called funcMyDFT(x,M) to compute the DFT of the sequence
%x. Note M is the DFT length which may be different from the sequence 
%length N. If M > N the function should zero pad the input before computing
%the DFT. If M < N then the input must be truncated. Test your function by 
%comparing the output of your function with that obtained by fft.m. You can
%use the randn.m function to generate the test vectors.
test = randn(4,1)
X = fft(test,5)
[X_func] = funcMyDFT(test, 5)
%it works!! :)

%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 2                                    %
%                                                                        %
%------------------------------------------------------------------------%
%2. Consider the following signal
% x[n] = [?3,8,8,12,?3,12,8,8]
% (a)	Compute the 8-point DFT of x[n]. Plot the magnitude and the phase 
% using the stem function (also for the other questions of this problem). Is the DFT real? Why?

x = [-3,8,8,12,-3,12,8,8];
[X] = funcMyDFT(x,8);
%X2 = fft(x,8)
k=0:7;

figure(1)
subplot(2,1,1)
stem(k,abs(X))
title('Magnitude of the DFT (a)');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,phase(X));
title('Phase of the DFT (a)');
ylabel('Phase (rad)');
xlabel('Index');


%%




% For parts (b)-(e), plot the magnitude and the phase of the DFT described,
% using the stem function. Compared with part (a), what do you observe? 
% Explain your observation in terms of the DFT and its properties.
% (b)	Shift (circular shift) the signal by 3, compute the 8-point DFT of the shifted signal.
%x[n] = [12,?3,12,8,8,-3,8,8]
x = [12,-3,12,8,8,-3,8,8]

[X] = funcMyDFT(x,8)
%X_DFT = fft(x,8)

k=0:7;
%X2 = fft(x,8)
figure(2)
subplot(2,1,1)
stem(k,abs(X))
title('Magnitude of the DFT (b)');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,phase(X));
title('Phase of the DFT (b)');
ylabel('Phase (rad)');
xlabel('Index');

%%
% (c)	Shift (circular shift) the signal by -3, compute the 8-point DFT of the shifted signal.
x = [12,8,8,-3,8,8,12,-3];
X = funcMyDFT(x,8);
%X = fft(x,8)
k=0:7;
figure(3)
subplot(2,1,1)
stem(k,abs(X))
title('Magnitude of the DFT (c)');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,angle(X));
title('Phase of the DFT (c)');
ylabel('Phase (rad)');
xlabel('Index');

%%
% (d)	Compute the 16-point DFT of the signal using zero padding.
x = [-3,8,8,12,-3,12,8,8];
[X] = funcMyDFT(x,16);
k=0:15;
figure(4)
subplot(2,1,1)
stem(k,abs(X))
title('Magnitude of the DFT (d)');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,angle(X));
title('Phase of the DFT (d)');
ylabel('Phase (rad)');
xlabel('Index');

%%
% (e)	Compute the 32-point DFT of the zero interpolated signal

% x[n] = [?3,0,0,8,0,0,8,0,0,12,0,0,?3,0,0,12,0,0,8,0,0,8,0,0].
x = [-3,0,0,8,0,0,8,0,0,12,0,0,-3,0,0,12,0,0,8,0,0,8,0,0];
[X] = funcMyDFT(x,32);
k=0:31;
figure(5)
subplot(2,1,1)
stem(k,abs(X))
title('Magnitude of the DFT (e)');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,phase(X));
title('Phase of the DFT (e)');
ylabel('Phase (rad)');
xlabel('Index');

%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 3                                    %
%                                                                        %
%------------------------------------------------------------------------%
% 3.	Let x[n] be a discrete time sequence:
% ((0.8)n,	if 0 ? n ? 8 x[n] =
% 	0,	else
% (a)	Determine the analytical expression for the DTFT of x[n] and 
% plot the magnitude and phase of the DTFT.
%after doing it on paper...
w = [0:2*pi/250:2*pi]
X_d = (1 - ((0.8).*exp(-1.*j.*w)).^9)./(1-0.8.*exp(-1.*j.*w))
figure(6);
subplot(2,1,1)
plot(w,abs(X_d))
title('Magnitude of DTFT');
xlabel('w (radians)');
ylabel('Magnitude');

subplot(2,1,2)
plot(w,phase(X_d))
title('Phase of DTFT');
xlabel('w (radians)');
ylabel('Phase');

%%
% (b)	In MATLAB, compute the 9-point DFT of x[n],0 ? n ? 8 using the 
% DFT function you wrote in Problem 1. Plot the magnitude and phase using the stem function.
x = zeros(9,1);
for p = 1:9
    x(p) = (0.8)^(p-1)
end
k=0:8;

[X_9] = funcMyDFT(x,9);
figure(7)
subplot(2,1,1)
stem(k,abs(X_9))
title('Magnitude of the DFT');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,phase(X_9));
title('Phase of the DFT');
ylabel('Phase (rad)');
xlabel('Index');



% (c)	Compute the 16-point DFT of x[n],0 ? n ? 15 using the DFT function 
% you wrote in Problem 1. Plot the magnitude and phase using the stem function. 
% Comment on the effect of zeropadding the signal on its DFT.
x = zeros(9,1);
for p = 1:9
    x(p) = (0.8)^(p-1)
end
k=0:15;
[X_16] = funcMyDFT(x,16);
figure(8)
subplot(2,1,1)
stem(k,abs(X_16))
title('Magnitude of the 16pt DFT');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
stem(k,phase(X_16));
title('Phase of the 16pt DFT');
ylabel('Phase (rad)');
xlabel('Index');


% (d)	Compute the 128-point DFT of x[n],0 ? n ? 127 using the DFT 
% function you wrote in Problem 1. Plot the magnitude and phase using the 
% plot function instead of stem when many dense points exist, to avoid the 
% appearance of a solid blue area.
x = zeros(9,1);
for p = 1:9
    x(p) = (0.8)^(p-1)
end
k=0:127;

[X_128] = funcMyDFT(x,128);
figure(9)
subplot(2,1,1)
plot(k,abs(X_128))
title('Magnitude of the 128pt DFT');
ylabel('Magnitude');
xlabel('Index');
subplot(2,1,2);
plot(k,phase(X_128));
title('Phase of the 128pt DFT');
ylabel('Phase (rad)');
xlabel('Index');
% (e)	Compare the results from part (d) to the plots of part (a). How 
% does this relate to the relationship between the digital frequency ? and the DFT index k?

%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 4                                    %
%                                                                        %
%------------------------------------------------------------------------%
% 4.	Download the file tones.mat from the course web page (this should 
% be included in the zipped assignment file for lab 2). The file contains 
% a signal that has multiple tones in it. Load the signal using the command
% load tones;. The variable x should now contain the signal.
% (a)	Truncate the length of x to be 16. Compute the 16-point DFT of x. 
% Plot the magnitude of the DFT using the plot command. How many distinct 
% tone frequencies do you see?
load tones
y = x(1:16);
X = funcMyDFT(y,16);
figure(8)
plot(abs(X));
title('Magnitude of 16 pt DFT of Truncated Signal');
ylabel('Magnitude');
xlabel('Index (K)');
% (b)	Try to improve the resolution of DFT spectrum by using zero padding
% on the length-16 signal. Can you find more tone frequencies? How many 
% tones can you distinguish and what are their values?
X1 = funcMyDFT(y,32);
figure(9);
plot(abs(X1));
title('Magnitude of 32 pt DFT of ZPed Truncated Signal');
ylabel('Magnitude');
xlabel('Index (K)');
%%
% (c)	Compute the DFT of x without truncation of x (M = N = 36). Try to 
% improve the resolution of the DFT spectrum by using zero padding. Can you
% find more tone frequencies? How many tones can you distinguish and what 
% are their values?
load tones;
X = funcMyDFT(x,36);
figure(10);
plot(abs(X));
title('Magnitude of 36 pt DFT (no ZP)');
ylabel('Magnitude');
xlabel('Index (K)');

X1 = funcMyDFT(x,128);
figure(11);
plot(abs(X1));
title('Magnitude of 128 pt DFT (ZP)');
ylabel('Magnitude');
xlabel('Index (K)');

X2 = funcMyDFT(x,512);
figure(12);
plot(abs(X2));
title('Magnitude of 512 pt DFT (ZP)');
ylabel('Magnitude');
xlabel('Index (K)');

%%
% (d)	When we truncate the signal x, we are using a rectangular, or 
% ‘boxcar’ window. Try at least two other windowing methods, such as the 
% Hamming and Hanning (Hann) windows. These windows may be applied using 
% the hamming.m and hann.m functions as follows
% L = length(x); x windowed = hamming(L).*x;
% Compute the DFT of xwindowed with and without zero padding. How many 
% tones can you distinguish and what are their values? Comment on the 
% differences between results obtained using various windowing methods, 
% and explain them in terms of the DTFT of the windows. What are the 
% tradeoffs in using different windows?
load tones
L = length(x);
%hamming(L)
x_Hamm = hamming(L).*x'
%no ZP
X_Hamm = funcMyDFT(x_Hamm,L);
figure(13)
plot(abs(X_Hamm));
title('Magnitude Plot Hamming Window (no ZP) 36pt DFT');
ylabel('Magnitude');
xlabel('Index (K)');

X_Hamm_ZP = funcMyDFT(x_Hamm,128);
figure(14)
plot(abs(X_Hamm_ZP));
title('Magnitude Plot Hamming Window (ZP) 128pt DFT');
ylabel('Magnitude');
xlabel('Index (K)');

X_Hamm_ZP_512 = funcMyDFT(x_Hamm,512);
figure(15)
plot(abs(X_Hamm_ZP_512));
title('Magnitude Plot Hamming Window (ZP) 512pt DFT');
ylabel('Magnitude');
xlabel('Index (K)');

%Hanning Window
x_Hann = hann(L).*x'
X_Hann = funcMyDFT(x_Hann,L);
figure(16)
plot(abs(X_Hann));
title('Magnitude Plot Hanning Window (no ZP) 36pt DFT');
ylabel('Magnitude');
xlabel('Index (K)');

X_Hann_ZP = funcMyDFT(x_Hann,128);
figure(17)
plot(abs(X_Hann_ZP));
title('Magnitude Plot Hanning Window (ZP) 128pt DFT');
ylabel('Magnitude');
xlabel('Index (K)');

X_Hann_ZP_512 = funcMyDFT(x_Hann,512);
figure(18)
plot(abs(X_Hann_ZP_512));
title('Magnitude Plot Hanning Window (ZP) 512pt DFT');
ylabel('Magnitude');
xlabel('Index (K)');


%%
%------------------------------------------------------------------------%
%                                                                        %
%                              Part 5                                    %
%                                                                        %
%------------------------------------------------------------------------%
% 5.	Finish the following tasks related to the lab demo:
% (a)	Load the file DialNumber.mat. The variable dialTone records a piece
% of 10-digit dialing signal, which is sampled at Fs = 1/T = 8000 samples/s
% . Calculate the DFT spectrum of every 100-ms segment of the signal. What 
% is the dialed telephone number? Calculate and display the short-time 
% Fourier transform of the signal using the following commands:
% figure; windowLength = 128; noverlap = 100;
% spectrogram(dialTone,windowLength,noverlap,512,8000);
load Dialnumber
%how long is each 100-ms segment index wise? Well we know the sampling rate
% 8000 samples per sec. There are 32000 samples therefore 32/8=4seconds
% passed by and there are 40 100ms segments. Right?
windowLength = 128;
noverlap = 100;
spectrogram(dialTone,windowLength,noverlap,512,8000);
%%
load Dialnumber
%x1_Hamm = hamming(100).*dialTone(
X1 = funcMyDFT(dialTone(1:100),1024);
figure(19)
plot(abs(X1));
title('DFT Time 0-100ms');
ylabel('Magnitude');
xlabel('Index (K)');
%%
X2 = funcMyDFT(dialTone(101:200),1024);
figure(20)
plot(abs(X2));
title('DFT Time 101-200ms');
ylabel('Magnitude');
xlabel('Index (K)');

X3 = funcMyDFT(dialTone(201:300),1024);
figure(21)
plot(abs(X3));
title('DFT Time 201-300ms');
ylabel('Magnitude');
xlabel('Index (K)');

X4 = funcMyDFT(dialTone(301:400),1024);
figure(22)
plot(abs(X4));
title('DFT Time 301-400ms');
ylabel('Magnitude');
xlabel('Index (K)');

X5 = funcMyDFT(dialTone(401:500),1024);
figure(23)
plot(abs(X5));
title('DFT Time 401-500ms');
ylabel('Magnitude');
xlabel('Index (K)');


%%
% (b)	Load the file NMRSpec.mat. The variable st contains an NMR signal. 
% The data-sampling rate is Fs = 1/T = 2000 samples/s. Calculate the DFT 
% spectrum of the signal. Plot the real and the imaginary parts of the 
% spectrum. Plot the magnitude and the phase of the spectrum. Calculate 
% the 32-point DFT spectrum of the first 32 samples of the signal. Can you 
% distinguish the peaks corresponding to creatine (around 209 Hz) and 
% choline (around 185 Hz)? Will zero-padding help you out? Why?
load NMRSpec
X_st = funcMyDFT(st,1024)
%Real and Imag
figure(24)
subplot(2,1,1)
plot(real(X_st));
title('Real Part of NMR DFT');
subplot(2,1,2)
plot(imag(X_st));
title('Imaginary Part of NMR DFT');

%Mag and Phase
figure(25)
subplot(2,1,1)
plot(abs(X_st));
title('Magnitude of NMR DFT');
subplot(2,1,2)
plot(phase(X_st));
title('Phase of NMR DFT');

%%32 Point DFT
k = 0:31;
X_st_32 = funcMyDFT(st(1:32),32);
figure(26)
stem(k,abs(X_st_32));
title('Magnitude of 32 point DFT of NMR');
xlabel('Index (K)')
ylabel('Magnitude');
%%
%remember wK = 2piK/N sooo
k = 1:32
w = 2*pi*k/32
f = w./(2*pi)

%%
% (c)	Load the image cameraman.tif by calling 
% Im = imread('cameraman.tif');. You may view the image using imagesc(Im). 
% Calculate the magnitude spectrum and the phase spectrum of the image.
%% load and dispaly an image
fignum = 27;
Im = imread('cameraman.tif');
figure(fignum);
imagesc(Im);
colormap(gray(256));
axis image off;
fignum = fignum + 1;

%% Calculate the fft of the image
s = fftshift(fft2(Im));

%% Display the magntiude spectrum
MagSpec = abs(s);
MagSpecLog = log10(abs(s)+eps)*20;
figure(fignum);
imagesc(MagSpecLog);
colormap(jet(256));
axis image off;
fignum = fignum + 1;
%% Display the phase spectrum
PhSpec = angle(s);
figure(fignum);
imagesc(PhSpec);
colormap(jet(256));
axis image off;
fignum = fignum + 1;

% (i)	Add a linear phase along each row to the original phase spectrum 
% (add a phase term changing linearly with the column index). Reconstruct 
% the image using the magnitude and the new phase spectrum. What has 
% changed about the image? Can you explain this change in terms of the DFT?
%% Distort the phase
%%want to add a constant value to the previous entry by row
%i.e. PhSpecL(n,m) = PhSpec(n,m) + (constant)*m
PhSpecL = zeros(size(s))
constant = 1.7;
for n=1:256
    for m=1:256
        PhSpecL(n,m) = PhSpec(n,m) + constant*m;
    end
end



ImPhRand = ifft2(ifftshift(abs(s).*exp(1i*PhSpecL)));

figure(fignum);
imagesc(abs(ImPhRand));
colormap(gray(256));
axis image off;
fignum = fignum + 1;
%%
% (ii)	Add a linear phase along each column to the original phase spectrum
% (add a phase term changing linearly with the row index). Reconstruct the 
% image using the magnitude and the new phase spectrum. What has changed 
% about the image? Can you explain this change in terms of the DFT?

%%want to add a constant value to the previous entry by column
%i.e. PhSpecL(n,m) = PhSpec(n,m) + (constant)*m
PhSpecL = zeros(size(s))
constant = 1.7;
for n=1:256
    for m=1:256
        PhSpecL(n,m) = PhSpec(n,m) + constant*n;
    end
end



ImPhRand = ifft2(ifftshift(abs(s).*exp(1i*PhSpecL)));

figure(fignum);
imagesc(abs(ImPhRand));
colormap(gray(256));
axis image off;
fignum = fignum + 1;
% You can find useful code in the MATLAB file demoMagnPhSpec.m. Add 
%enough linear phase to the original phase spectrum that you can see a distinct effect.

%%

h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['figure' num2str(length(h)+1-i)], 'jpeg'); 
end
