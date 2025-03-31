rng(42)
data = 2*rand(2,1000)-1;
f = @(x) 15.*(x-.5).*(x+1).*x.*(x+.5).*sin(1.2*(x-.25))-.35;

class1 = []; % 0s
class2 = [];  % 1s
for i = 1:size(data, 2)
    this = data(:,i);
    if this(2) >= f(this(1))+.5*randn
        class1 = [class1, this];
    else
        class2 = [class2, this];
    end
end

figure
hold on, grid on, grid minor
plot(class1(1,:), class1(2,:), 'b.', 'MarkerSize', 10)
plot(class2(1,:), class2(2,:), 'r.','MarkerSize',10)
fplot(f, [-1 1])
axis([-1 1 -1 1])