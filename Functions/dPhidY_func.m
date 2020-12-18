function [dPhidY] = dPhidY_func(y, lure)

    % unpack parameters
    W1 = lure.W1; 
    b1 = lure.b1;
    W2 = lure.W2;
    I2 = ones(size(lure.W1));

    dPhidY = zeros(length(y),1);
    for i = 1:length(y)
        % dphi/dy
        M           = diag(I2 - tanh(W1*y(i) + b1).^2);
        dPhidY(i,1) = W2.'*M*W1;
    end
end

