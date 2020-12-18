function [dPhidTheta] = dPhidTheta_func(y, lure)

    % unpack parameters
    W1 = lure.W1; 
    b1 = lure.b1;
    W2 = lure.W2;
    I2 = ones(size(lure.W1));

    dPhidTheta = zeros(length(y),44);                                      % !GET RID OF THIS 14 LATER ON!
    for i = 1:length(y)
        % dphi/dW1
        dphidW1 = W2.'*diag(I2 - tanh(W1*y(i) + b1).^2)*y(i);

        % dphi/db1
        dphidb1 = W2.'*diag(tanh(b1).^2 - tanh(W1*y(i)+b1).^2);

        %dphi/dW2
        dphidW2 = ( tanh(W1*y(i)+b1) - tanh(b1) ).';

        dPhidTheta(i,:) = [zeros(1,35), dphidW1, dphidb1, dphidW2];
    end
end

