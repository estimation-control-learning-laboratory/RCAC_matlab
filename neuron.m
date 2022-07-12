function [y] = neuron(x,theta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


y = sigmoid([x' 1]*theta);

end



function y = sigmoid(x)
y = 1/(1+exp(-x));

end