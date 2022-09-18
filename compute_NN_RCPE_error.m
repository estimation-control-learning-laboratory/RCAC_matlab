function [z] = compute_NN_RCPE_error(Theta)

X = [0 0 1 1
    0 1 0 1];
Y = [0 1 1 0];


jj = randi([1,4]);


%         for kk = 1:4
%             Yhat_f(1,kk) = [1 1]*neural_layer(X(:,kk),Theta);
%         end

Yhat = [1 1]*neural_layer(X(:,jj),Theta);
z     = Y(:,jj) - Yhat;

if 1
    z=0;
    for jj=1:4
        Yhat = [1 1]*neural_layer(X(:,jj),Theta);
        z(jj)     = Y(:,jj) - Yhat;
    end

%     z = z(randi([1 4])); 
    z = norm(z);
%     z = sum(z);
end


end