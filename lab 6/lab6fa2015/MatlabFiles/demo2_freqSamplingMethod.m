%% original (sampled boxcar)

k1 = 0:2;
k2 = 19:20;
Gd2 = [exp(-1i*10*2*pi.*k1./21) zeros(1,16) exp(-1i*10*2*pi.*k2./21)];

g2 = ifft(Gd2);
figure, stem(0:20, abs(g2)), title('|h[n]|');
figure, freqz(g2), title('Hd(w)');

%% with transition band
Gd3 = [exp(-1i*10*2*pi.*k1./21) zeros(1,16) exp(-1i*10*2*pi.*k2./21)];
a = 0.3;
Gd3(4) = a*exp(-1i*10*2*pi*3/21); % k = 3
Gd3(19) = a*exp(-1i*10*2*pi*18/21); % k = 18
g3 = ifft(Gd3);
figure, stem(0:20, abs(g3)), title('h[n]');
figure, freqz(g3), title('Hd(w)');

