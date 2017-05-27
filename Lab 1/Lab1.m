% Number 1)
% a)
%on a single plot (e.g. if x is complex plot(real(x),imag(x)) representing
%the complex plane, graph
%     the complex number z = 1 -3j (mark the location with an asterisk (*)
z = 1 - 3*1j;

%     (also give the phase and magnitude of this number
z_real = real(z);
z_imag = imag(z);

figure(1);
subplot(2,1,1)
plot(z,'*');
%plot(real(z),imag(z),'*');
%
hold on;
% The complex number exp(1-.2j) (marked with an 'o')
y = exp(1 - .2*j);
plot(y,'o')

theta = [0:2*pi/1000:2*pi];
unit_circle = exp(j*theta);
plot(unit_circle)
axis equal
title('Plot of 2 Complex Numbers With Unit Circle Reference');
xlabel('real');
ylabel('imag');
hold off;



subplot(2,1,2);
%
%part b on same figure
n = [0:99];
w = exp(.03*(1+3*j)*n);
plot(w,'x');
title('Part b');
xlabel('Real Part');
ylabel('Imaginary Part');
%
%  part c)
w_real = real(w);
w_im = imag(w);
figure(2);
%real part
subplot(2,1,1);
stem(n,w_real,'b');
title('Real');
xlabel('Samples');
ylabel('Real Amplitude');
%imag part
subplot(2,1,2);
stem(n,w_im,'r');
title('Imaginary');
xlabel('Samples');
ylabel('Imaginary Amplitude');
%
w_mag = abs(w);
w_phase = phase(w);
%vs unwrap and angle
w_angle = angle(w);

figure(3);
subplot(3,1,1);
stem(n,w_mag);
title('Magnitude Plot');
xlabel('Samples');
ylabel('Magnitude');

subplot(3,1,2);
stem(n,w_phase);
title('Phase Plot using "Phase"');
xlabel('Samples');
ylabel('Phase');

subplot(3,1,3);
stem(n,w_angle);
title('Phase Plot using "Angle"');
xlabel('Samples');
ylabel('Phase');

%%
%   Number 2
%M = 2,5,7,12 
N = 14;
n = [0:1:2*N-1]

figure(4);
subplot(2,2,1);
%M=2
x_m_2 = sin(2*pi*2*n/N);
stem(n,x_m_2,'b');
title('M=2');
xlabel('Samples (n)');
ylabel('Amplitude');

subplot(2,2,2);
%M=5;
x_m_5 = sin(2*pi*5*n/N);
stem(n,x_m_5,'r');
title('M=5');
xlabel('Samples (n)');
ylabel('Amplitude');

subplot(2,2,3);
%M=7;
x_m_7 = sin(2*pi*7*n/N);
stem(n,x_m_7,'y');
title('M=7');
xlabel('Samples (n)');
ylabel('Amplitude');

subplot(2,2,4);
%M=12;
x_m_12 = sin(2*pi*12*n/N);
stem(n,x_m_12,'g');
title('M=12');
xlabel('Samples (n)');
ylabel('Amplitude');

%%
% Number 3
t = [-10:.1:9.9];
b = sinc(t);
zero = 0*t;


figure(5);
plot(t,b);
title('Plot of Sinc Function');
xlabel('Time');
ylabel('Amplitude');
hold on
plot(t,zero,'r');
for i = 1:10
    plot(i,0,'o')
end
for i = -10:-1
    plot(i,0,'o')
end

%%
%b plot the dirichlet function
N = 7
t = [-4*pi:.1:4*pi]
d_7 = diric(t,N)
N = 4;
d_4 = diric(t,N);
figure(6);
subplot(2,1,1);
plot(t,d_7);
title('Diric Function n=7');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(t,d_4);
title('Diric Function n=4');
xlabel('Time');
ylabel('Amplitude');

%%
%Number 4 

%z-.35/(z-.8*exp(-j*pi/3));
%create a 100x100 matrix for the entries f(z) for Re(z)={0,.01,.02...} and
%similarly for Im(z)
%i = 0:.01:.99;
%j = 0:.01:.99;
A = zeros(100,100);
%matlab matrix indices start at 1

a = 0:.01:.99;
b = 0:.01:.99;

for l = 1:100
    for m = 1:100
        z = a(l) +j*b(m);
        A(l,m)= abs((z-.35/(z-.8*exp(-j*pi/3))));
    end
end
figure(7)
surf(a,b,A);
title('Magnitude Plot of Complex Function');
xlabel('Imaginary');
ylabel('Real');
zlabel('Magnitude');

%%
pha = phantom(256);
figure(8);
imagesc(pha);
title('Phantom');
n = 1:1:256;
%plotting the 150th row
figure(9);
plot(n, pha(150,n));
title('150th Row of Phantom');
xlabel('column index');
ylabel('Amplitude');
%%
load clown;
figure(10);
imagesc(X);
title('Clown');

%%
figure(11);
MAP = colormap;
NEWMAP = rgb2gray(MAP);
colormap(NEWMAP);
imagesc(X);
title('Gray Scale Clown');
%%
m = 1:1:320;
figure(12);
plot(m, X(150,m));
title('150th Row of clown');
xlabel('column index');
ylabel('Amplitude');

%what do I have to do for colormap???

%%
%part c write a function to flip an image
%upside down, and right to left
load clown;

clown_updown = Flip_Image(X,0);
clown_lr = Flip_Image(X,1);
clown_udlr = Flip_Image(clown_updown,1);
figure(13);
imagesc(X);
title('Regular Clown');

figure(14);
imagesc(clown_updown);
title('Upside Down Clown');

figure(15);
imagesc(clown_lr);
title('Flipped Right to Left Clown');

figure(16);
imagesc(clown_udlr);
title('Upside Down Flipped Right to Left Clown');

%%
rando = rand(200,320);
new_clown = rando.*X;
figure(17);
imagesc(new_clown);
title('Random Matrix Multiplied with Clown');

R = normrnd(.5,.1,200,320);
gauss_clown = R.*X;
figure(18);
imagesc(gauss_clown);
title('Gauss Clown');


%random matrix with 0s or 1s only
Q = rand(200,320);
for e = 1:1:200
    for p = 1:1:320
        Q(e,p) = round(Q(e,p));
    end
end
bin_clown = Q.*X;
figure(19);
imagesc(bin_clown);
title('Binary Clown');

%%
% Number 5
nx = [-3:1:7];
x = nx*0;
x(4) = 2;
x(6) = -1;
x(7) = 1;
x(8) = 3;
figure(20);
stem(nx,x);
title('Unshifted Signal');
xlabel('Samples');
ylabel('Amplitude');



y_1 = x;
y_2 = x;
y_3 = x;
y_4 = x;
ny1= nx + 2;
ny2= nx - 1;
ny3= nx*-1;
ny4= nx*-1 +1;

ii = 4;
jj = 1;
figure(21);
subplot(ii,jj,1);
stem(ny1,y_1,'b');
title('y_1 = x[n-2]');
xlabel('Samples');
ylabel('Amplitude');

subplot(ii,jj,2);
stem(ny2,y_2,'r');
title('y_2 = x[n+1]');
xlabel('Samples');
ylabel('Amplitude');

subplot(ii,jj,3);
stem(ny3,y_3,'g');
title('y_3 = x[-n]');
xlabel('Samples');
ylabel('Amplitude');

subplot(ii,jj,4);
stem(ny4,y_4,'y');
title('y_4 = x[-n+1]');
xlabel('Samples');
ylabel('Amplitude');