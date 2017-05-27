%% ECE 311 class demonstration, dual-tone multifrequency signal detection
% Chao Ma,
%

clear all;
close all;
clc;

fignum = 500;

%% Generate the twelve frequency pairs
symbol = ['1';'2';'3';'4';'5';'6';'7';'8';'9';'*';'0';'#'];
lfg = [697 770 852 941]; % Low frequency group
hfg = [1209 1336 1477];  % High frequency group
f  = [];
for c=1:4,
    for r=1:3,
        f = [ f [lfg(c);hfg(r)] ];
    end
end

%% Generate and visualize the DTMF tones
Fs  = 8000;       % Sampling frequency 8 kHz
N = 800;          % Tones of 100 ms
t   = (0:N-1)/Fs; % 800 samples at Fs

tones = zeros(N,size(f,2));
figure(fignum);
for toneChoice=1:12,
    % Generate tone
    tones(:,toneChoice) = sum(sin(f(:,toneChoice)*2*pi*t))';
    % Plot tone
    subplot(4,3,toneChoice),plot(t*1e3,tones(:,toneChoice));
    title(['Symbol "', symbol(toneChoice),'": [',num2str(f(1,toneChoice)),',',num2str(f(2,toneChoice)),']'])
    set(gca, 'Xlim', [0 25]);
    ylabel('Amplitude');
    if toneChoice>9, xlabel('Time (ms)'); end
end
set(gcf,'Color',[1,1,1]);
fignum = fignum + 1;

%% Spectra of the DTMF tone
Nfft = 1024;
figure(fignum);
faxis = Fs*(0:Nfft-1)/Nfft;
indf = find(faxis>600&faxis<1600);
for toneChoice=1:12,
    % Select tone
    tone=tones(:,toneChoice);
    % Estimate DFT using Goertzel
    ydft(:,toneChoice) = fft(tone,Nfft); % Goertzel use 1-based indexing
    % Plot magnitude of the DFT
    subplot(4,3,toneChoice),plot(faxis(indf),abs(ydft((indf),toneChoice)));
    title(['Symbol "', symbol(toneChoice),'": [',num2str(f(1,toneChoice)),',',num2str(f(2,toneChoice)),']'])
%     set(gca, 'XTick', estim_f, 'XTickLabel', estim_f, 'Xlim', [650 1550]);
    ylabel('DFT Magnitude');
    if toneChoice>9, xlabel('Frequency (Hz)'); end
end
set(gcf,'Color',[1,1,1]);
fignum = fignum + 1;
