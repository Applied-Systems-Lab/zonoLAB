clear
clc
close all

figure

a = 10;

Z1 = zono([a/2 a/2; a/2 0], [0;a/2]);
Z2 = zono([a/2 a/2; a/2 0], [-a;-a/2]);
Z3 = zono([a/2 a/2; a/2 0], [a;a/2]);

plot(Z1,'b',0.2)
hold on
plot(Z2,'r',0.2)
plot(Z3,'r',0.2)

plot([0 -10],[0 0],'m','LineWidth',3)
plot([0 10],[0 10],'m','LineWidth',3)

axis([-20 20 ...
    -20 20])
grid on

saveas(gcf,'relunn.png')

%%

clear
clc
close all
a = 5;
x = linspace(-a,a,50)

f = figure
f.Position = [704 818 1296 420]
subplot(1,3,1)
y2 = (exp(x)-exp(-x))./(exp(x)+exp(-x));
plot(x,y2,'k',LineWidth=3)
axis([-a a -a a])
grid on
title('Tanh')

subplot(1,3,2)
y3 = 1./(1+exp(-x));
plot(x,y3,'k',LineWidth=3)
axis([-a a -a a])
grid on
title('Sigmoid')

subplot(1,3,3)
y1 = x./(1+exp(-x));
plot(x,y1,'k',LineWidth=3)
axis([-a a -a a])
grid on
title('Sigmoid Linear Unit')


saveas(gcf,'activationFunctions.png')