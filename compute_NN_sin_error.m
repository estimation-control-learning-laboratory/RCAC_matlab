function [z] = compute_NN_sin_error(Theta)

jj = randi([1,4]);
x = rand*pi-pi/2;

%         for kk = 1:4
%             Yhat_f(1,kk) = [1 1]*neural_layer(X(:,kk),Theta);
%         end

Yhat = [1 1]*neural_layer(x,Theta);
z     = sin(x) - Yhat;

% keyboard

end