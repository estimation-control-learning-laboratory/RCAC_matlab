function [z] = calc_error(theta)

theta = theta(:);

Theta = [theta(1:3) theta(4:6)];

X = [0 0 1 1
     0 1 0 1];
Y = [0 1 1 0];
Yhat = 0*Y;
for ii = 1:4
    Yhat(1,ii) = [1 1]*neural_layer(X(:,ii),Theta); 
end

z = Y-Yhat;


end





% function [x_ip1] = neural_layer(x_i,Theta)
% 
% 
% for jj = 1:size(Theta,2)
%     x_ip1(jj,1) = neuron(x_i,Theta(:,jj));
% end
% 
% end
% 
% function [y] = neuron(x,theta)
% %UNTITLED4 Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% y = sigmoid([x' 1]*theta);
% 
% end
% 
% 
% 
% function y = sigmoid(x)
% y = 1/(1+exp(-x));
% 
% end