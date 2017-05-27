close all
%% System A
% Create the system by specifying the cofficients of an LCCDE equation
b = [0.1199,0.3801,0.3801,0.1199];
a = 1;
sysA = tf(b,a,-1);

% Plot the zero-pole digram
figure; axis auto; pzplot(sysA);



% Create the input signal
n = 0:19;
x = cos(pi*n/10) + cos(8*pi*n/10);

% Calculate the output signal
y = filter(b, a, x);

% Plot the output signal
figure;stem(y);title('y[n]');

% Calculate and plot the DFT spectra
num = 1024;
Xk = abs(fft(x,num));
figure; plot(Xk); title('X[k]');
axis tight;
Yk = abs(fft(y,num));
figure; plot(Yk); title('Y[k]');
axis tight;

%% System B
% Create the system by specifying the cofficients of an LCCDE equation
a = poly([exp(-1i*8*pi/10) exp(1i*8*pi/10)]);
b = poly(1);
sysB = tf(b,a,-1);

% Plot the zero-pole digram
figure; axis auto; pzplot(sysB);

% Create the first input signal
n = 0:19;
x = cos(pi*n/10);

% Calculate the output signal
y = filter(b, a, x);

% Plot the output signal
figure;stem(y);

% Create the second input signal
n = 0:19;
x = cos(8*pi*n/10);

% Calculate the output signal
y = filter(b, a, x);

% Plot the output signal
figure;stem(y);