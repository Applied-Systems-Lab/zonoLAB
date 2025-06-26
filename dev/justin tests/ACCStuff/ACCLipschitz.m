clear
clc
close all

figure

x = linspace(-5, 5, 100);
y = tanh(x);

f = figure
slope = (2./(exp(x)+exp(-x))).^2;
annotation(f,'textbox',[0.5 0.14 0.4 0.3], 'String','$tanh(x)$','EdgeColor','none','Interpreter','latex','FontSize',12)
L = 0;

for i = 2: size(x,2)
    plot(x,y,'b','LineWidth',3)
    hold on
    y2 = slope(i)*(x);
    plot(x(i),y(i),'or',LineWidth=2)
    % plot(x+x(i),-y2+y(i),'--r','LineWidth',2)
    plot(x+x(i),y2+y(i),'--r','LineWidth',2)
    if(slope(i) == max(slope))
        L = slope(i)*(x);
        xL = x(i);
        yL = y(i);
    end
    if (L ~= 0)
        plot(xL,yL,'ok',LineWidth=2)
        plot(x+xL,-L+yL,'-k','LineWidth',2)
        plot(x+xL,L+yL,'-k','LineWidth',2)
    end
    hold off
    axis([-5 5 -1.5 1.5])
    % pause(0.25)
    % exportgraphics(gcf,'lipschitzani2.gif',Append=true)
end

%%
annotation(f,'textbox',[0.5 0.14 0.4 0.3], 'String','$tanh(x)$','EdgeColor','none','Interpreter','latex','FontSize',12)
for i = 2: size(x,2)
    y2 = slope(i)*(x);
    plot(x,y,'b','LineWidth',3)
    hold on
    plot(x(i),y(i),'ok',LineWidth=2)
    plot(x+x(i),-L+y(i),'-k','LineWidth',2)
    plot(x+x(i),L+y(i),'-k','LineWidth',2)
    plot(x+x(i),y2+y(i),'--r','LineWidth',2)
    hold off
    axis([-5 5 -1.5 1.5])
    exportgraphics(gcf,'lipschitzani3.gif',Append=true)
end

% f2 = figure
% f2.Position = [345 248 560 420];
% for i = 1: length(F)
%     imshow(F(i).cdata)
%     exportgraphics(gcf,'lipschitzani.gif',Append=true)
% end


%
% for a = 0:5
%         plot([-5+a 5-a],[-5 5],'b','LineWidth',3)
%         plot([-5+a 5-a],[-5 5],'b','LineWidth',3)
%         title(sprintf('Slope: %.2f', 10/(10-(2*a))))
%         axis([-5 5 ...
%             -5 5])
%         xlabel('input')
%         ylabel('output')
%         grid on
%         pause(1)
%         F(a+1) = getframe(gcf);
% end
%
