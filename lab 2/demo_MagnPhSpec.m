%% ECE 311 class demonstration, magntidue spectrum, phase spectrum
% Chao Ma
%

clear all; 
close all;
clc
fignum = 500;

%% load and dispaly an image
Im = imread('cameraman.tif');
figure(fignum);
imagesc(Im);
colormap(gray(256));
axis image off;
fignum = fignum + 1;

%% Calculate the fft of the image
s = fftshift(fft2(Im));

%% Display the magntiude spectrum
MagSpec = abs(s);
MagSpecLog = log10(abs(s)+eps)*20;
figure(fignum);
imagesc(MagSpecLog);
colormap(jet(256));
axis image off;
fignum = fignum + 1;

%% Display the phase spectrum
PhSpec = angle(s);
figure(fignum);
imagesc(PhSpec);
colormap(gray(256));
axis image off;
fignum = fignum + 1;

%% Reconstruct the image
ImRecon = ifft2(ifftshift(abs(s).*exp(1i*PhSpec)));
figure(fignum);
imagesc(abs(ImRecon));
colormap(gray(256));
axis image off;
fignum = fignum + 1;

%% Distort the phase
PhSpecRand = PhSpec + 0.7*randn(size(s));
ImPhRand = ifft2(ifftshift(abs(s).*exp(1i*PhSpecRand)));

figure(fignum);
imagesc(abs(ImPhRand));
colormap(gray(256));
axis image off;
fignum = fignum + 1;
