%Lab 1 pre-assignment stuff
%MATLAB is funny, % are for comments

%%
%this double % seperates code

%%
%exercise 1. Find functions on MATLAB that generate random numbers
%rand works and can take 2 arguments to make a rand matrix
A = rand(3,3)
%this creates a random 3x3 matrix

%%
%exercise 2. Find functions that calculate the conjugate, absolute value,
%real part, and imaginary part of a complex number
%abs, conj, imag, real are the functions


%calculate the conjugate, absolute value, real part, and imaginary part of
% 3 + 4i
x = 3 + 4*i;
x_conj = conj(x);
x_r = real(x);
x_im = imag(x);

%%
%exercise 3
%generate the signal shown on the lab
n = [0:20];
x = exp(i*3*pi*n/10);
%display the real and imaginary parts
x_r = real(x);
x_im = imag(x);

subplot(2,1,1)
stem(n,x_r, 'r')
title('Real Part of x[n]');
xlabel('Index (n)');
ylabel('Amplitude');

subplot(2,1,2)
stem(n,x_im, 'b')
title('Imaginary Part of x[n]');
xlabel('Index (n)');
ylabel('Amplitude');