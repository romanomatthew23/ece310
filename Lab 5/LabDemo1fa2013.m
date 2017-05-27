clc, clear all, close all
ts = 14;
b = [1 1];
a = [1 -0.25838];
[h,w] = freqz(b,a,1024);

figure;
subplot(2,1,1);
plot(w,mag2db(abs(h)),'b','linewidth',2);
grid on;
set(gca,'GridLineStyle','-');
axis tight;
xlabel('\omega','fontsize',ts);
ylabel('|H_d(\omega)|','fontsize',ts);
set(gca,'xtick',[0 pi/4 pi/2 3*pi/4 pi]); set(gca,'xlim',[0,pi*1.01]);
set(gca,'xticklabel','0| p/4 | p/2| 3p/4 | p','fontname','symbol','fontsize',ts)

subplot(2,1,2);
plot(w,angle(h),'r','linewidth',2);
grid on;
set(gca,'GridLineStyle','-');
axis tight;
xlabel('\omega','fontsize',ts);
ylabel('\angle H_d(\omega)','fontsize',ts);
set(gca,'xtick',[0 pi/4 pi/2 3*pi/4 pi]); set(gca,'xlim',[0,pi*1.01]);
set(gca,'xticklabel','0| p/4 | p/2| 3p/4 | p','fontname','symbol','fontsize',ts)

n = 0:40-1;
x1n = cos(pi*n/10);
x2n = cos(8*pi*n/10);

y1n = filter(b,a,x1n);
y2n = filter(b,a,x2n);

figure;
subplot(2,2,1);
stem(n,x1n);
ylabel('x_1[n]');
xlabel('n');
set(gca,'ylim',[-3,3])

subplot(2,2,2);
stem(n,y1n);
ylabel('y_1[n]');
xlabel('n'); set(gca,'ylim',[-3,3])

subplot(2,2,3);
stem(n,x2n);
ylabel('x_2[n]');
xlabel('n'); set(gca,'ylim',[-3,3])

subplot(2,2,4);
stem(n,y2n);
ylabel('y_2[n]');
xlabel('n'); set(gca,'ylim',[-3,3])

%%
X1 = fftshift(fft(x1n,2^(5+nextpow2(length(x1n)))));
X2 = fftshift(fft(x2n,2^(5+nextpow2(length(x2n)))));
Y1 = fftshift(fft(y1n,2^(5+nextpow2(length(y1n)))));
Y2 = fftshift(fft(y2n,2^(5+nextpow2(length(y2n)))));

figure;
subplot(2,2,1);
plot(2*pi/length(X1)*(0:length(X1)-1)-pi,abs(X1));
axis tight;
xlabel('\omega');
ylabel('X_1(\omega)');
set(gca,'xtick',[-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]); set(gca,'xlim',[-pi*1.01,pi*1.01]);
set(gca,'xticklabel','-p || -p/2 || 0 || p/2 || p','fontname','symbol','fontsize',ts)

subplot(2,2,2);
plot(2*pi/length(X2)*(0:length(X2)-1)-pi,abs(X2));
axis tight;
xlabel('\omega');
ylabel('X_2(\omega)');
set(gca,'xtick',[-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]); set(gca,'xlim',[-pi*1.01,pi*1.01]);
set(gca,'xticklabel','-p || -p/2 || 0 || p/2 || p','fontname','symbol','fontsize',ts)

subplot(2,2,3);
plot(2*pi/length(Y1)*(0:length(Y1)-1)-pi,abs(Y1));
axis tight;
xlabel('\omega');
ylabel('Y_1(\omega)');
set(gca,'xtick',[-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]); set(gca,'xlim',[-pi*1.01,pi*1.01]);
set(gca,'xticklabel','-p || -p/2 || 0 || p/2 || p','fontname','symbol','fontsize',ts)

subplot(2,2,4);
plot(2*pi/length(Y2)*(0:length(Y2)-1)-pi,abs(Y2));
axis tight;
xlabel('\omega');
ylabel('Y_2(\omega)');
set(gca,'xtick',[-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]); set(gca,'xlim',[-pi*1.01,pi*1.01]);
set(gca,'xticklabel','-p || -p/2 || 0 || p/2 || p','fontname','symbol','fontsize',ts)