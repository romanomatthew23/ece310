% ECE 310 HW #3 Plotting Help
w = [0:pi/4:2*pi];

y_d = 1 - ((.5*exp(-j*w)).^4);
d = 1-.5*exp(-j*w)
%Y_d = (1 - (.5*exp(-j*w)).^4)/(1-.5*exp(-j*w))
Y_d = y_d./d;

%%
figure(1)
%subplot(2,1,1);
plot(w,abs(Y_d));
%subplot(2,1,2);
%plot(w,phase(Y_d));
Mag = abs(Y_d)

%%
w = [-pi:2*pi/100:pi];
v1 = 14*sinc(14*(w+(pi/4))/pi);
v2 = 14*sinc(14*(w-(pi/4))/pi);
v = v1 + v2;
figure(1);
subplot(3,1,1);
plot(w,v1);

subplot(3,1,2);
plot(w,v2);
             
subplot(3,1,3);
plot(w,v);

%%
%next part
%plotting the dirichlet function (not actually diirhahc)

w = 3*pi/28;
w1 = w - pi/4;
w2 = w + pi/4;
ve1 = 14*sin(14*w1)/sin(w1/2)
ve2 = 14*sin(14*w2)/sin(w2/2)
vee = ve1 + ve2

