function [x_ip1] = neural_layer(x_i,Theta)


for jj = 1:size(Theta,2)
    x_ip1(jj,1) = neuron(x_i,Theta(:,jj));
end

end