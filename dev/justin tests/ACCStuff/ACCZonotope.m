clear
clc
close all

figure

Z = zono([1 2; 2 1], [1;1]);

plot(Z,'r',0.8)
axis([-5 5 ...
    -5 5])
yline(1)
xline(1)
%%
f = figure
f.Position = [345 248 560 420];

plot(1,1,'ok');
annotation(f,'textbox',[0.59 0.3 0.4 0.3], 'String','c','EdgeColor','none','Interpreter','latex','FontSize',12)
axis([-5 5 ...
    -5 5])
F(1) = getframe(gcf);
%%

delete(findall(gcf,'type','annotation'))
plot(1,1,'ok')
hold on
plot([1,2],[1,3],'-r','LineWidth',3)
plot([1,3],[1,2],'-b','LineWidth',3)
annotation(f,'textbox',[0.59 0.3 0.4 0.3], 'String','c','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.65 0.43 0.4 0.3], 'String','$G_r$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.67 0.35 0.4 0.3], 'String','$G_b$','EdgeColor','none','Interpreter','latex','FontSize',12)
hold off
axis([-5 5 ...
    -5 5])

F(2) = getframe(gcf);
%%

delete(findall(gcf,'type','annotation'))
plot(1,1,'ok')
hold on
plot([1,2],[1,3],'-r','LineWidth',3)
plot([1,3],[1,2],'-b','LineWidth',3)

plot([1,0],[1,-1],'--r','LineWidth',3)
plot([1,-1],[1,0],'--b','LineWidth',3)
hold off
axis([-5 5 ...
    -5 5])
annotation(f,'textbox',[0.59 0.3 0.4 0.3], 'String','c','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.65 0.43 0.4 0.3], 'String','$G_r$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.67 0.35 0.4 0.3], 'String','$G_b$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.55 0.23 0.4 0.3], 'String','$-G_r$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.47 0.33 0.4 0.3], 'String','$-G_b$','EdgeColor','none','Interpreter','latex','FontSize',12)

F(3) = getframe(gcf);

%%

delete(findall(gcf,'type','annotation'))
plot(1,1,'ok')
hold on
plot([0,2],[-1,3],'-r','LineWidth',3)
plot([-1,3],[0,2],'-b','LineWidth',3)
hold off
axis([-5 5 ...
    -5 5])

annotation(f,'textbox',[0.59 0.3 0.4 0.3], 'String','c','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.65 0.43 0.4 0.3], 'String','$|G_r|$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.67 0.35 0.4 0.3], 'String','$|G_b|$','EdgeColor','none','Interpreter','latex','FontSize',12)

F(4) = getframe(gcf);

%%

delete(findall(gcf,'type','annotation'))
plot(1,1,'ok')
hold on
plot([0,2],[-1,3],'-r','LineWidth',3)
plot([-2,2],[-2,0],'-b','LineWidth',3)
plot([0,4],[2,4],'-b','LineWidth',3)
hold off
axis([-5 5 ...
    -5 5])
annotation(f,'textbox',[0.59 0.3 0.4 0.3], 'String','c','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.65 0.43 0.4 0.3], 'String','$|G_r|$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.6 0.52 0.4 0.3], 'String','$|G_b|$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.5 0.14 0.4 0.3], 'String','$|G_b|$','EdgeColor','none','Interpreter','latex','FontSize',12)

F(5) = getframe(gcf);

%%

delete(findall(gcf,'type','annotation'))
plot(1,1,'ok')
hold on
plot([-2,0],[-2,2],'-r','LineWidth',3)
plot([2,4],[0,4],'-r','LineWidth',3)
plot([-2,2],[-2,0],'-b','LineWidth',3)
plot([0,4],[2,4],'-b','LineWidth',3)
hold off
axis([-5 5 ...
    -5 5])
annotation(f,'textbox',[0.59 0.3 0.4 0.3], 'String','c','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.75 0.40 0.4 0.3], 'String','$|G_r|$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.38 0.27 0.4 0.3], 'String','$|G_r|$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.6 0.52 0.4 0.3], 'String','$|G_b|$','EdgeColor','none','Interpreter','latex','FontSize',12)
annotation(f,'textbox',[0.5 0.14 0.4 0.3], 'String','$|G_b|$','EdgeColor','none','Interpreter','latex','FontSize',12)

F(6) = getframe(gcf);

%%
hold on
delete(findall(gcf,'type','annotation'))
plot(Z,'g',0.3)
annotation(f,'textbox',[0.5 0.14 0.4 0.3], 'String','$Z$','EdgeColor','none','Interpreter','latex','FontSize',12)
F(7) = getframe(gcf);


%%
f2 = figure
f2.Position = [345 248 560 420];
for i = 1: length(F)
    for j = 1: 10
        imshow(F(i).cdata)
        exportgraphics(gcf,'zonoAni.gif',Append=true)
    end
end