clear
clc
close all

%%
load("mnist.mat")

training.images = imresize(training.images(9:22,8:21,:),0.33);
test.images = imresize(test.images(9:22,8:21,:),0.33);

save('5x5mnist.mat')