%% FIR filter designed using Parks-McClellan

N = 21;
n = 0:N-1;
g = (1/4).*sinc(1/4.*(n-10));

%%designed using fdatool, limited to order 20, wp = 0.25, ws = 7/20, Wp = 2, Ws=1
b = [   -0.0275    0.0385    0.0360    0.0162   -0.0233   -0.0569   -0.0465    0.0280    0.1477    0.2589    0.3042    0.2589    0.1477    0.0280   -0.0465   -0.0569   -0.0233    0.0162    0.0360    0.0385   -0.0275];

stem(0:20, b); xlabel('n'); title('h[n]');
freqz(b,1);

%% FIR filter designed using windowing method - Kaiser Window
kk = kaiser(1,N);

hkais = g.*kk;
figure, stem(0:N-1,hkais), title('h[n] Kaiser Window'); xlabel('n'); axis tight;
figure, freqz(hrect,1), title('H_d(w)');