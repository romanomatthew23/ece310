%% Lab 5 Code
%1a
Fs = 6000;
d = fdesign.highpass('Fst,Fp,Ast,Ap',300,350,40,1,Fs);
hd = design(d,'kaiserwin');
b = hd.Numerator; a = 1;
fvtool(b,a)

%%
%1b
Fs = 6000;
d = fdesign.highpass('Fst,Fp,Ast,Ap',300,350,40,1,Fs);
hd = design(d,'equiripple');
b1b = hd.Numerator; a1b = 1;
%%
fvtool(b1b,a1b)

%%
%1d
Fs = 6000;
d = fdesign.lowpass('Fp,Fst,Ap,Ast',900,950,1,40,Fs);
hd = design(d,'equiripple');
b1d = hd.Numerator; a1d = 1;
%%
fvtool(b1d,a1d);
%%
%1e
Fs = 6000;
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',500,550,650,700,40,1,40,Fs);
hd = design(d,'equiripple');
b1e = hd.Numerator; a1e = 1;
%%
fvtool(b1e,a1e);
%%
%1f
%load flute2.wav and filter it with b,d, and e
[x,Fs] = wavread('flute2.wav');
soundsc(x,Fs);
% show spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(x,L,128,nfft,Fs,'yaxis'); 
title('Original Spectrogram');
%%
%1b filter
x1b = filter(b1b, a1b, x);
soundsc(x1b,Fs);
% show recovered spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(x1b,L,128,nfft,Fs,'yaxis'); 
title('Flute2 Filter 1b');
%%
%1d filter
x1d = filter(b1d, a1d, x);
soundsc(x1d,Fs);
% show recovered spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(x1d,L,128,nfft,Fs,'yaxis');
title('Flute2 Filter 1d');
%%
%1e filter
x1e = filter(b1e, a1e, x);
soundsc(x1e,Fs);
% show recovered spectrogram
L = 256; nfft = 1024; %Hamming (default) window length is 256, DFT length is 1024 for spectrogram
figure();
spectrogram(x1e,L,128,nfft,Fs,'yaxis'); 
title('Flute2 Filter 1e');


%%
%part 2
%%2a
Fs = 6000;
d = fdesign.highpass('Fst,Fp,Ast,Ap',300,350,40,1,Fs);
hd = design(d,'ellip');
hd = convert(hd,'df1');
b = hd.coeffs.Numerator; a = hd.coeffs.Denominator;
fvtool(b,a);
%%
%2c
Fs = 6000;
d = fdesign.highpass('Fst,Fp,Ast,Ap',300,350,40,1,Fs);
hd = design(d,'cheby2');
hd = convert(hd,'df1');
b = hd.coeffs.Numerator; a = hd.coeffs.Denominator;
fvtool(b,a);

%%
%2e
Fs = 6000;
d = fdesign.lowpass('Fp,Fst,Ap,Ast',900,950,1,40,Fs);
hd = design(d,'ellip');
hd = convert(hd,'df1');
b = hd.coeffs.Numerator; a = hd.coeffs.Denominator;
fvtool(b,a);
%%
%2f
Fs = 6000;
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',500,550,650,700,40,1,40,Fs);
hd = design(d,'ellip');
hd = convert(hd,'df1');
b = hd.coeffs.Numerator; a = hd.coeffs.Denominator;
fvtool(b,a);





