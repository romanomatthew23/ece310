

[x,Fs] = wavread('Flute.wav');
soundsc(x,Fs);

figure();
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
spectrogram(x,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content
%set(gca,'ylim',[0 2500]); % get a closer look at the frequency content
return

%% pass the data through an LSI system, i.e. low-pass filter
% design a low-pass filter
d = fdesign.lowpass('N,Fc',100,440*2,Fs);
hd = design(d);
b = hd.Numerator; a = 1;
% plot frequency response of filter
figure(); freqz(b,a,1024);
% filter signal
xFiltered = filter(b,a,x);
soundsc(xFiltered,Fs);

% show new spectrogram
figure();
spectrogram(xFiltered,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content


%% pass the data through a non-LSI (Shift-varying) system
w0 = 2*pi*2000;
t = (0:length(x)-1)/Fs;
temp = cos(w0*t);
xModulated = x.*temp.'; % modulation is a shift-varying operation
soundsc(xModulated,Fs);

figure();
spectrogram(xModulated,L,128,nfft,Fs,'yaxis'); 

% show modulation of the time-domain signal
figure; 
subplot(211); hold on
plot(t,x); xlabel('t');
closeup = (80000:84000); plot(t(closeup),x(closeup),'r'); 
title('Time Signal'); axis tight
subplot(212)
plot(t(closeup),x(closeup),'b'); hold on
plot(t(closeup),xModulated(closeup),'g'); legend('x','xModulated')
title('{\color{red} Segment} of xModulated[n] vs. orginal'); xlabel('t');
axis tight;

