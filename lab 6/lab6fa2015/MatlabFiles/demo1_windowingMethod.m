%% windowing method - rectangular window
N = 21;
n = 0:N-1;
g = (1/4).*sinc(1/4.*(n-10));

rect = ones(1,N);
hrect = rect.*g;
figure, stem(0:N-1,hrect), title('h[n] Rectangular Window'); xlabel('n'); axis tight;
figure, freqz(hrect,1), title('H_d(w)');

%% windowing method - hamming window
ham = hamming(N);
hham = g.*ham';
figure, stem(0:N-1,hham), title('h[n] Hamming Window'); xlabel('n'); axis tight;
figure, freqz(hham,1), title('H_d(w)');


%% Windowing method - Hamming window, various lengths
N1 = 21;
n1 = 0:N1-1;
N2 = 41;
n2 = 0:N2-1;
N3 = 61;
n3 = 0:N3-1;
g1 = (1/4).*sinc(1/4.*(n1-10));
g2 = (1/4).*sinc(1/4.*(n2-20));
g3 = (1/4).*sinc(1/4.*(n3-30));

ham1 = hamming(N1);
ham2 = hamming(N2);
ham3 = hamming(N3);
hham1 = g1.*ham1';
hham2 = g2.*ham2';
hham3 = g3.*ham3';
figure, stem(0:N1-1,hham1), title('h[n] Hamming Window'); xlabel('n'); axis tight;
figure, freqz(hham1,1), title('H_d(w)');
figure, stem(0:N2-1,hham2), title('h[n] Hamming Window'); xlabel('n'); axis tight;
figure, freqz(hham2,1), title('H_d(w)');
figure, stem(0:N3-1,hham3), title('h[n] Hamming Window'); xlabel('n'); axis tight;
figure, freqz(hham3,1), title('H_d(w)');