function [lure] = thetaToLure(theta)
    
    global scaleswitch
    
    % LTI Block
    lure.A  = reshape(theta(1:16), 4,4);
    lure.B  = reshape(theta(17:20),4,1);      
    lure.L  = reshape(theta(21:24),4,1);                               
    lure.C  = reshape(theta(25:28),1,4);                                        
    lure.D  = theta(29); 
    lure.F  = reshape(theta(30:33),1,4);      
    lure.G  = theta(34);
    lure.H  = theta(35);   

    % Nonlinearity
    lure.W1     = reshape(theta(36:38),3,1);  
    lure.b1     = reshape(theta(39:41),3,1);  
    lure.W2     = reshape(theta(42:44),3,1);
    lure.b2     = -lure.W2.'*tanh(lure.b1);
    lure.gamma  = lure.W2.'*lure.W1;
    
    if scaleswitch
        %%% Scale to unit norm of W2
        scaling       = norm(lure.W2);
        lure.W2       = lure.W2./scaling;
        lure.b2       = lure.b2./scaling;
        lure.B        = lure.B.*scaling;
        lure.G        = lure.G.*scaling;
        lure.gamma    = lure.gamma./(scaling);
        lure.scaling  = scaling;
        
        %%% Scale to unit norm of W1
        scaling1      = norm(lure.W1);
        lure.W1       = lure.W1./scaling1;
        lure.C        = lure.C.*scaling1;
        lure.D        = lure.D.*scaling1;
        lure.gamma    = lure.gamma./(scaling1);
        lure.scaling1 = scaling1;
    end

    % Define the nonlinearity and its partial derivatives
    lure.nonlin     = @(y) phi_func(y, lure);
    lure.dphidtheta = @(y) dPhidTheta_func(y, lure);
    lure.dphidy     = @(y) dPhidY_func(y, lure);
end

