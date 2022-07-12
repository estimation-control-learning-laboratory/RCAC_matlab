clc
clear all
close all
format longe
% randn('state',1)

theta_0 = randn(1,6);


% myfun(theta)
options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-15)
theta = fsolve(@calc_error,theta_0,options);
theta'
% calc_error(theta)

Theta = [theta(1:3)' theta(4:6)'];

X = [0 0 1 1
     0 1 0 1];
Y = [0 1 1 0];
for ii = 1:4
    Yhat(1,ii) = [1 1]*neural_layer(X(:,ii),Theta); 
end
Yhat



