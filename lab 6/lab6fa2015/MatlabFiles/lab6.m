%lab 6
%1a LPF wc = pi/3 rect window
N=25;
M = (N-1)/2;
n=0:N-1;
w=zeros(1,100);
w=1;
h=(1/3).*sinc((1/3).*(n-M)).*w;

figure, stem(0:N-1,h), title('h[n] Rectangular Window'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w)');
figure, stem(n,impz(h)), title('Impulse Response Rectangular Window'); xlabel('n');
%%
%1b LPF wc = pi/3 hamming window
N=25;
M = (N-1)/2;
n=0:N-1;
w=(0.54-0.46*cos(pi*n/12));
h=(1/3)*sinc((1/3)*(n-M)).*w;
figure, stem(0:N-1,h), title('h[n] Hamming Window'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w)');
figure, stem(n,impz(h)), title('Impulse Response Hamming Window'); xlabel('n');
%%
%1c LPF wc = pi/3 kaiser window
N=25;
M = (N-1)/2;
n=0:N-1;
beta = 2.05  %play around with beta until filter has stopband attenutation of 30dB
w=kaiser(N,beta);
h=(1/3)*sinc((1/3)*(n-M)).*w';
figure, stem(0:N-1,h), title('h[n] Kaiser Window'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w)');
figure, stem(n,impz(h)), title('Impulse Response Kaiser Window'); xlabel('n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%
%1d  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1d
%1a LPF wc = pi/3 rect window
N=50;
M = (N-1)/2;
n=0:N-1;
w=1;
h=(1/3).*sinc((1/3).*(n-M)).*w;
%stem(n,impz(h))
figure, stem(0:N-1,h), title('h[n] Rectangular Window N=50'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w) N=50');
figure, stem(n,impz(h)), title('Impulse Response Rectangular Window N=50'); xlabel('n');
%%
%1b LPF wc = pi/3 hamming window
N=50;
M = (N-1)/2;
n=0:N-1;
w=(0.54-0.46.*cos(pi*n/12));
h=(1/3).*sinc((1/3)*(n-M)).*w;
figure, stem(0:N-1,h), title('h[n] Hamming Window N=50'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w) N=50');
figure, stem(n,impz(h)), title('Impulse Response Hamming Window N=50'); xlabel('n');
%%
%1c LPF wc = pi/3 kaiser window
N=50;
M = (N-1)/2;
n=0:N-1;
beta = 2.05  %play around with beta until filter has stopband attenutation of 30dB
w=kaiser(N,beta);
h=(1/3)*sinc((1/3)*(n-M)).*w';
figure, stem(0:N-1,h), title('h[n] Kaiser Window N=50'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w) N=50');
figure, stem(n,impz(h)), title('Impulse Response Kaiser Window N=50'); xlabel('n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%
%%
%1e
% design HPF with wc = 2pi/3
%1b is below
N=25;
M = (N-1)/2;
n=0:N-1;
w=(0.54-0.46*cos(pi*n/12));
h=(1/3)*sinc((1/3).*(n-M)).*w;
%however need to modulate digitally by pi so
h=h.*(2*cos(pi*n))
figure, stem(0:N-1,h), title('h[n] Modulated to HPF'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w) Modulated to HPF');
figure, stem(n,impz(h)), title('Impulse Response Modulated to HPF'); xlabel('n');

%%
% part f now N=50
%1e
% design HPF with wc = 2pi/3
%1b is below
N=50;
M = (N-1)/2;
n=0:N-1;
w=(0.54-0.46*cos(pi*n/12));
h=(1/3)*sinc((1/3).*(n-M)).*w;
%however need to modulate digitally by pi so
h=h.*(2*cos(pi*n))
figure, stem(0:N-1,h), title('h[n] Modulated to HPF N=50'); xlabel('n'); axis tight;
figure, freqz(h,1), title('H_d(w) Modulated to HPF N=50');
figure, stem(n,impz(h)), title('Impulse Response Modulated to HPF N=50'); xlabel('n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                                                                        %
%                             Part 2                                     %   
%                                                                        %
%                                                                        %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%2b
N=25;
n = 0:N-1;
Gd2 = [ones(1,5) zeros(1,16) ones(1,4)];

g2 = ifft(Gd2);

figure, stem(n,abs(g2)), title('|h[n]| 2b'); xlabel('n'); axis tight;
figure, freqz(g2), title('H_d(w) 2b');
figure, stem(n,impz(g2)), title('Impulse Response 2b'); xlabel('n');

%%
%2d
N=25;
M=(N-1)/2;
n=0:N-1;
k1 = 0:4;
k2 = 21:24;
Gd2 = [exp(-1i*M*2*pi.*k1./N) zeros(1,16) exp(-1i*M*2*pi.*k2./N)];
g2 = ifft(Gd2);
figure, stem(n, abs(g2)), title('|h[n]| 2d');
figure, freqz(g2), title('Hd(w) 2d');
figure, stem(n,impz(g2)), title('Impulse Response 2d'); xlabel('n');
%%
a=.5; %keep trying different values for a
Gd2(5+1) = a*exp(-1i*M*2*pi.*6./N);
Gd2(20+1) = a*exp(-1i*M*2*pi.*21./N);

g2 = ifft(Gd2);
figure, stem(n, abs(g2)), title('|h[n]| 2e');
figure, freqz(g2), title('Hd(w) 2e');
figure, stem(n,impz(g2)), title('Impulse Response 2e'); xlabel('n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                                                                        %
%                             Part 3                                     %   
%                                                                        %
%                                                                        %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%a
%designed using fdatool wp=0.3pi, ws=0.36pi
N=25;
n=0:N-1;
b = [0.0126   -0.0642   -0.0301    0.0016    0.0326    0.0350   -0.0031   -0.0547   -0.0656    0.0033    0.1386  0.2736    0.3300    0.2736    0.1386    0.0033   -0.0656   -0.0547   -0.0031    0.0350    0.0326    0.0016 -0.0301   -0.0642    0.0126];
figure, stem(n, abs(b)), title('|h[n]| 3a');
figure, freqz(b), title('Hd(w) 3a');
figure, stem(n,impz(b)), title('Impulse Response 3a'); xlabel('n');
%%
%b same as a but attenuation at 30dB at cost of transition band (.3 and
%.42 now)
N=25;
n=0:N-1;
b = [  0.0139   -0.0075   -0.0210   -0.0180    0.0086    0.0339    0.0208   -0.0330   -0.0722   -0.0256    0.1198  0.2866    0.3602    0.2866    0.1198   -0.0256   -0.0722   -0.0330    0.0208    0.0339    0.0086   -0.0180 -0.0210   -0.0075    0.0139];


figure, stem(n, abs(b)), title('|h[n]| 3b');
figure, freqz(b), title('Hd(w) 3b');
figure, stem(n,impz(b)), title('Impulse Response 3b'); xlabel('n');
%%
%c same as b but adjusting weights instead to get to 30dB had to do 10 on
%Wstop and 1 on Wpass
N=25;
n=0:N-1;
b = [  -0.0239   -0.0100    0.0161    0.0538    0.0762    0.0590    0.0038   -0.0532   -0.0606    0.0117    0.1440 0.2719    0.3246    0.2719    0.1440    0.0117   -0.0606   -0.0532    0.0038    0.0590    0.0762    0.0538  0.0161   -0.0100   -0.0239];

figure, stem(n, abs(b)), title('|h[n]| 3c');
figure, freqz(b), title('Hd(w) 3c');
figure, stem(n,impz(b)), title('Impulse Response 3c'); xlabel('n');
%%
%d changed a to N=50
N=50;
n=0:N-1;
b=[ 0.0059   -0.0128   -0.0098   -0.0034    0.0054    0.0093    0.0035   -0.0078   -0.0132   -0.0053    0.0104  0.0184    0.0080   -0.0140   -0.0261   -0.0123    0.0196    0.0394    0.0204   -0.0305   -0.0689   -0.0419  0.0660    0.2116    0.3153    0.3153    0.2116    0.0660   -0.0419   -0.0689   -0.0305    0.0204    0.0394  0.0196   -0.0123   -0.0261   -0.0140    0.0080    0.0184    0.0104   -0.0053   -0.0132   -0.0078    0.0035 0.0093    0.0054   -0.0034   -0.0098   -0.0128    0.0059];
figure, stem(n, abs(b)), title('|h[n]| 3d');
figure, freqz(b), title('Hd(w) 3d');
figure, stem(n,impz(b)), title('Impulse Response 3d'); xlabel('n');

%%
%e design multiband 
%your design currently works but theres no reason to think its minimal
%N=60;
%n=0:N;
%f = [0 .1  .15   .3   .4  .7   .75   1];
%a = [1  1   0    0    1    1    0    0];
%w = [1 1.5 .5 5];
%b = firpm(N,f,a,w);
%figure, stem(n, abs(b)), title('|h[n]|');
%figure, freqz(b), title('Hd(w)');
%figure, stem(n,impz(b)), title('Impulse Response'); xlabel('n');
%%
%should us fircls(n,f,a,up,lo)
N=250;  %nothing actually works, they just stop when they get close enough :(
n=0:N;
% -30dB is .0316227766
% -40dB is .01
% 2dB is 1.258925412
% -2dB is .794328347
% .5dB is 1.059253725
% -.5dB is .9440608763
%    PASSBAND       STOPBAND        PASSBAND        STOPBAND
f = [0     .1      .15    .3       .4     .7      .75     1];
a = [  1        1       0      0       1       1        0  ];
up =[ 1.259    100    .03163   100  1.059254   100     .01]; 
lo =[ 0.7944  -100      0     -100  .944061   -100      0 ];
b = fircls(N,f,a,up,lo);
figure, stem(n, abs(b)), title('|h[n]| 3e');
figure, freqz(b), title('Hd(w) 3e');
figure, stem(n,impz(b)), title('Impulse Response 3e'); xlabel('n');

%%
%Problem 4
%a
%find a reasonable sampling freq between 6000 and 44100Hz
load('speechsig.mat');
Fs = 9250;  %10000 slightly too high
%8500 slightly too low
%anywhere between 8500 to 10000 seems possible
% i would guess 9250
soundsc(x,Fs);
%soundsc(x)
%soundsc(xnoise);

%%
%b calculate and plot mag of DFT for both signals


figure();
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
spectrogram(x,L,128,nfft,Fs,'yaxis'); 
title('Line up! (Without Noise)');
soundsc(x,Fs);
%%
figure();
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
spectrogram(xnoise,L,128,nfft,Fs,'yaxis');
title('Line Up! (With Noise)');
soundsc(xnoise,Fs);
%%
% part c
%desinging a LPF
%fdatool
b = [-0.0217411346989065,-0.00145848753357376,0.0155739633020352,0.0153206544667661,-0.00289783134587723,-0.00845130491296203,0.0100737085738482,0.0204253596031618,-0.000478343950468720,-0.0203417848172653,-0.00105524575249998,0.0290217529581835,0.0131386103847575,-0.0312150383043354,-0.0248406764628302,0.0371702830069087,0.0494524623769139,-0.0389279004751593,-0.0962737197514837,0.0420427603567735,0.315646460435112,0.458275598878023,0.315646460435112,0.0420427603567735,-0.0962737197514837,-0.0389279004751593,0.0494524623769139,0.0371702830069087,-0.0248406764628302,-0.0312150383043354,0.0131386103847575,0.0290217529581835,-0.00105524575249998,-0.0203417848172653,-0.000478343950468720,0.0204253596031618,0.0100737085738482,-0.00845130491296203,-0.00289783134587723,0.0153206544667661,0.0155739633020352,-0.00145848753357376,-0.0217411346989065];
n=0:42;
figure, stem(n, abs(b)), title('|h[n]| 4c');
figure, freqz(b), title('Hd(w) 4c');
figure, stem(n,impz(b)), title('Impulse Response 4c'); xlabel('n');

%%
%part d
xfiltered = filter(b,1,xnoise);
figure();
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
spectrogram(xfiltered,L,128,nfft,Fs,'yaxis');
title('Line Up! (Filtered Noise)');
soundsc(xfiltered,Fs);
%%
soundsc(x,Fs);
%%
soundsc(xnoise,Fs);
%%
%part e mean square thing
N=43;
a=1;
xmean = filter(b,a,xnoise);
b = (1/N).*ones(1,N);
figure();
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
spectrogram(xmean,L,128,nfft,Fs,'yaxis');
title('Line Up! (Filtered Noise with Mean Filter)');
soundsc(xmean,Fs);