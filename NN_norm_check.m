clc
clear all
close all

X = [0 0 1 1
    0 1 0 1];
Y = [0 1 1 0];

theta_nn0(2)= -4.481706269421458e+01;
theta_nn0(3)= -2.311066575047726e+01;
theta_nn0(4)= -4.584524260712102e+01;
theta_nn0(5)=   4.522449892595682e+01;
theta_nn0(6)= -2.293930732337328e+01;

z_norm = 0;
theta_nn0(1) = 4.510866097865148e+01;
for x = 1:.01:50
    theta_nn0(1) = x;

    Theta = [theta_nn0(1:3)' theta_nn0(4:6)'];
    for kk = 1:4
        Yhat_f(1,kk) = [1 1]*neural_layer(X(:,kk),Theta);
    end
    z_norm(end+1) = norm(Y-Yhat_f);
end

plot(z_norm)