clear
clc
close all

f = figure;
f.Position = [345 248 560 420];

cZ = conZono([eye(3)],[ones(3,1)]);
plot(cZ, 'r', 0.8)

axis([-3 3 -3 3 -3 3])