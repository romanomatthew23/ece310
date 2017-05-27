%% image blurring
% ECE 311
% Chao
%

clear all;
close all;
clc;

fignum = 500;

%% load and display an image
[temp,map] = imread('car_plate.jpg');
% crop the image
temp2 = temp(1:1024,1:1536,:);
Im = imresize(temp2,[512,768]);
figure(fignum);
imagesc(Im);
colormap(map);
axis image off;
fignum = fignum +1;

%% Create a Gaussian window
w = linspace(-pi,pi,1023);
H = exp(-w.^2/0.4).*exp(-1i*(length(w)-1)/2*w);
temp = ifft(ifftshift(H))*sqrt(1023);
h = temp(511-20:511+20);
h = real(h);

figure(fignum);
plot(w,abs(H));
fignum = fignum + 1;

figure(fignum);
stem(real(h));
fignum = fignum + 1;

%% blur the image
Nx = size(Im,1);
Ny = size(Im,2);

temp = [];
for ind = 1:3
    for indx = 1:Nx
        temp(indx,:,ind) = cconv(h,Im(indx,:,ind));
    end
end
for ind = 1:3
    for indy = 1:Ny
        ImBlurred(:,indy,ind) = cconv(h(:),temp(:,indy,ind));
    end
end

figure(fignum);
imagesc(abs(ImBlurred)./max(ImBlurred(:)));
% colormap(map);
axis image off;
fignum = fignum + 1;

%% Design the deblurring filter
G = conj(H)./(abs(H).^2+0.01);
figure(fignum);
plot(w,abs(G));
fignum = fignum + 1;
% design the deblurring filter
temp = ifft(ifftshift(G))*sqrt(1023);
g = temp(511-20:511+20);
figure(fignum);
stem(real(g));
fignum = fignum + 1;

%% Deblurred the image

Nx = size(ImBlurred,1);
Ny = size(ImBlurred,2);

temp = [];
for ind = 1:3
    for indx = 1:Nx
        temp(indx,:,ind) = cconv(g,ImBlurred(indx,:,ind));
    end
end
for ind = 1:3
    for indy = 1:Ny
        ImDeBlurred(:,indy,ind) = cconv(g(:),temp(:,indy,ind));
    end
end

%%
figure(fignum);
imagesc(abs(ImDeBlurred)./max(abs(ImDeBlurred(:))));
colormap(map);
axis image off;
fignum = fignum + 1;