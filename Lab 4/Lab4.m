%% This is The Code for Lab 4
% Problem 1
% a) write a func y=syseqn(x,yn1,a) 

%b)
a=1;
yn1=0;
%compute response for delta and unit step
%delta
x1= zeros(31,1);
x1(1) = 1;
%
y1 = syseqn(x1,yn1,a);

%unit step
x2=ones(31,1);
y2=syseqn(x2,yn1,a);

n=0:30;
figure();
stem(n,y1);
title('Impulse Response of System');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(n,y2);
title('Unit Step Response of System');
xlabel('index (n)');
ylabel('Amplitude');

%%
% c)
a = 1;
yn1 = -1;
%unit step
x1=ones(31,1);
y1=syseqn(x1,yn1,a);

%3*unit step
x2=3*ones(31,1);
y2=syseqn(x2,yn1,a);

n=0:30;
figure();
stem(n,y1);
title('Response of System (y_1[n])');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(n,y2);
title('Response of System (y_2[n])');
xlabel('index (n)');
ylabel('Amplitude');


%3*y1[n] - y2[n]
y3 = 3.*y1 - y2;
figure();
stem(n,y3);
title('Response of System (3y_1[n]-y_2[n])');
xlabel('index (n)');
ylabel('Amplitude');


%%
% d)
%when is it BIBO Stable?
a=1/5;
%unit step
x1=ones(31,1);
yn1 = 0;
y1=syseqn(x1,yn1,a);

yn1=3;
y2=syseqn(x1,yn1,a);



n=0:30;
figure();
stem(n,y1);
title('Response of System y[-1] =0');
xlabel('index (n)');
ylabel('Amplitude');


figure();
stem(n,y2);
title('Response of System y[-1]=3');
xlabel('index (n)');
ylabel('Amplitude');



%%
%p Problem 2
a = 2/5;
yn1 =0;
%delta
x1= zeros(20,1);
x1(1) = 1;
%calculate the impulse response of both systems
h1=syseqn(x1,yn1,a);
h2=syseqn2(x1,yn1,a);

n=0:19;
figure();
stem(n,h1);
title('h_1[n] (Impulse Response of System 1)');
xlabel('index (n)');
ylabel('Amplitude');


figure();
stem(n,h2);
title('h_2[n] (Impulse Response of System 2)');
xlabel('index (n)');
ylabel('Amplitude');


% part b) same thing for unit steps
a = 2/5;
yn1 =0;
%unit step
x1=ones(20,1);
%calculate the impulse response of both systems
s1=syseqn(x1,yn1,a);
s2=syseqn2(x1,yn1,a);

n=0:19;
figure();
stem(n,s1);
title('s_1[n] (Unit Step Response of System 1)');
xlabel('index (n)');
ylabel('Amplitude');


figure();
stem(n,s2);
title('s_2[n] (Unit Step Response of System 2)');
xlabel('index (n)');
ylabel('Amplitude');

%%
% convolve to obtain z1, z2
unit = ones(1,20);

z1 = conv(h1,unit);
z2 = conv(h2,unit);

z1=z1(1:20);
z2=z2(1:20);

n=0:19;
figure();
stem(n,z1);
title('z_1[n]');
xlabel('index (n)');
ylabel('Amplitude');


figure();
stem(n,z2);
title('z_2[n]');
xlabel('index (n)');
ylabel('Amplitude');

%part d
%plot s1 and z1 on same axes
figure();
subplot(2,1,1)
stem(n,s1);
title('s_1[n]');
xlabel('index (n)');
ylabel('Amplitude');
subplot(2,1,2);
stem(n,z1);
title('z_1[n]');
xlabel('index (n)');
ylabel('Amplitude');


%plot s2 and z2 on same axes
figure();
subplot(2,1,1)
stem(n,s2);
title('s_2[n]');
xlabel('index (n)');
ylabel('Amplitude');
subplot(2,1,2);
stem(n,z2);
title('z_2[n]');
xlabel('index (n)');
ylabel('Amplitude');

%%
% Problem 3

bm=[-0.0872, 
    0.0370, 
    0.0985, 
    0.1755, 
    0.2388, 
    0.2632, 
    0.2388, 
    0.1755, 
    0.0985, 
    0.0370,
    -0.0872]; 

% y[n] = summation m=0 to M bm*x[n-m]
% M is really min(10,n) but in matlab indices...
n = 0:31;
x1 = cos(pi*n/10);
x2 = cos(3*pi*n/10);

y1 = zeros(1,32);
for h=1:32
    for m=0:min(10,h-1)
        y1(h)=y1(h) + bm(m+1)*x1(h-m);
    end
end


y2 = zeros(1,32);
for h=1:32
    for m=0:min(10,h-1)
        y2(h)=y2(h) + bm(m+1)*x2(h-m);
    end
end

n=0:31;
figure();
stem(n,y1);
title('y_1[n]');
xlabel('index (n)');
ylabel('Amplitude');

n=0:31;
figure();
stem(n,y2);
title('y_2[n]');
xlabel('index (n)');
ylabel('Amplitude');

%compute the DFT of x1,x2,y1,y2
N = 128;
X1 = fft(x1,N);
X2 = fft(x2,N);
Y1 = fft(y1,N);
Y2 = fft(y2,N);

figure();
stem(abs(X1));
title('Magnitude of DFT of x1');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(X1));
title('Phase of DFT of x1');
xlabel('index (n)');
ylabel('Phase (radians)');

figure();
stem(abs(X2));
title('Magnitude of DFT of x2');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(X2));
title('Phase of DFT of x2');
xlabel('index (n)');
ylabel('Phase (radians)');

figure();
stem(abs(Y1));
title('Magnitude of DFT of y1');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(Y1));
title('Phase of DFT of y1');
xlabel('index (n)');
ylabel('Phase (radians)');

figure();
stem(abs(Y2));
title('Magnitude of DFT of y2');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(Y2));
title('Phase of DFT of y2');
xlabel('index (n)');
ylabel('Phase (radians)');

%%
% part b
n = 0:31;
x1 = cos(pi*n/10);
x2 = cos(3*pi*n/10);
y1 = x1.*cos(3*pi*n/5);
y2 = x2.*cos(3*pi*n/5);

%compute the DFT of x1,x2,y1,y2
N = 128;
X1 = fft(x1,N);
X2 = fft(x2,N);
Y1 = fft(y1,N);
Y2 = fft(y2,N);

figure();
stem(abs(X1));
title('Magnitude of DFT of x1');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(X1));
title('Phase of DFT of x1');
xlabel('index (n)');
ylabel('Phase (radians)');

figure();
stem(abs(X2));
title('Magnitude of DFT of x2');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(X2));
title('Phase of DFT of x2');
xlabel('index (n)');
ylabel('Phase (radians)');

figure();
stem(abs(Y1));
title('Magnitude of DFT of y1');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(Y1));
title('Phase of DFT of y1');
xlabel('index (n)');
ylabel('Phase (radians)');

figure();
stem(abs(Y2));
title('Magnitude of DFT of y2');
xlabel('index (n)');
ylabel('Amplitude');

figure();
stem(angle(Y2));
title('Phase of DFT of y2');
xlabel('index (n)');
ylabel('Phase (radians)');

%%
%problem 4
z = tf('z',-1);
H1 = (z-.5)*(z+5)/(z^2+1/16)/(z-0.3);
H2 = z/(z^3 + z^2 -4.75*z + 3.75);
H3 = (z^2 + 3*z + 1)/(z + exp(-1*j*pi/3))/(z + exp(1*j*pi/3))/(z + exp(-2*j*pi/3))/(z + exp(2*j*pi/3))
% plot the pole zero diagram
figure();
axis auto;
pzplot(H1);

figure();
axis auto;
pzplot(H2);

figure();
axis auto;
pzplot(H3);

%calulcate and plot unit impulse response
h1 = impulse(H1);
h2 = impulse(H2);
h3 = impulse(H3);
%%
n = 0:63;
h1(22:64) = 0;
figure();
stem(n,h1);
title('Unit Impulse Response System 1');
xlabel('index (n)');
ylabel('Amplitude');

h2(58:64) = 0;
figure();
stem(n,h2);
title('Unit Impulse Response System 2');
xlabel('index (n)');
ylabel('Amplitude');


figure();
stem(n,h3(1:64));
title('Unit Impulse Response System 3');
xlabel('index (n)');
ylabel('Amplitude');
%%
%sum the magnitude of values
magh1 = 0;
for r=1:64
    magh1 = abs(h1(r))+ magh1;
end
   
magh2 = 0;
for r=1:64
    magh2 = abs(h2(r))+ magh2;
end

magh3 = 0;
for r=1:64
    magh3 = abs(h3(r))+ magh3;
end

%%
%find input to make 3 go unstable
n=0:63;
x_explode = cos(pi*n/3);
[num,den] = tfdata(H3,'v');
y_explode = filter(num, den, x_explode);
%
n=0:63;
figure();
plot(abs(y_explode));
title('System 3 x[n]=cos(\pin/3)');

%%
% problem 5
[x,Fs] = wavread('Flute.wav');
soundsc(x,Fs);

% show spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(x,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content

%
a = .99;
w0 = 2*pi*880/Fs;

%create the transfer functionn for the filter
z = tf('z',-1);
H = (1 - 2*a*cos(w0) + a^2)/(z - a*exp(-j*w0))/(z - a*exp(j*w0));
[num,den] = tfdata(H,'v');
xFiltered = filter(num, den, x);
soundsc(xFiltered,Fs);

% show new spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(xFiltered,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content

%
% now calculate and apply the inverse filter
Hinv = 1/((1 - 2*a*cos(w0) + a^2)/(z - a*exp(-j*w0))/(z - a*exp(j*w0)))

[num,den] = tfdata(Hinv,'v');
denc = den(end);            %thing to get it to work
xRecovered = filter(num, denc, xFiltered);
soundsc(xRecovered,Fs);

% show recovered spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(xRecovered,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content

%%
%part b
[x,Fs] = wavread('Flute.wav');
soundsc(x,Fs);
%%
f0 = 440; %f0 = [440,494,554,587,659,880] Hz
a = 1;
w0 = 2*pi*f0/Fs;
%create the transfer functionn for the filter
z = tf('z',-1);
H = (1 - 2*a*cos(w0) + a^2)/(z - a*exp(-j*w0))/(z - a*exp(j*w0));
[num,den] = tfdata(H,'v');
xFiltered440 = filter(num, den, x);
soundsc(xFiltered440,Fs);

% show recovered spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(xFiltered440,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content
title('f_0 = 440Hz');
%%
f0 = 659; %f0 = [440,494,554,587,659,880] Hz
a = 1;
w0 = 2*pi*f0/Fs;
%create the transfer functionn for the filter
z = tf('z',-1);
H = (1 - 2*a*cos(w0) + a^2)/(z - a*exp(-j*w0))/(z - a*exp(j*w0));
[num,den] = tfdata(H,'v');
xFiltered659 = filter(num, den, x);
soundsc(xFiltered659,Fs);


% show recovered spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(xFiltered659,L,128,nfft,Fs,'yaxis'); 
set(gca,'ylim',[0 5000]); % get a closer look at the frequency content
title('f_0 = 659Hz');
