close all

rp = 1;           % Passband ripple (in dB)
rs = 80;          % Stopband ripple (in dB)
f = [.4 .6];      % Normalized cutoff frequencies
a = [1 0];        % Desired amplitudes

% Compute deviations (not in dB)
dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)]; 

% Estimate the filter order
[n,fo,ao,w] = firpmord(f,a,dev); % It turns out n = 25

% Design the filter
b = firpm(n,fo,ao,w);

% Display the filter
fvtool(b,1);