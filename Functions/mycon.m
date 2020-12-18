function [c, ceq, gradc, gradceq] = mycon(x)

    % Transform decision variables to a lure-type model
    lure    = thetaToLure(x);

    % Nonlinear inequality constraints
    c  = [];
        
    % Nonlinear equality constraints
    ceq(1,1) = norm(lure.W1) - 1;
    ceq(2,1) = norm(lure.W2) - 1;

    if nargout > 2
        % Nonlinear inequality constraint gradient
        gradc        = [];
        
        % Nonlinear equality constraint gradient
        Jac_W1        = zeros(3, length(x));
        Jac_W2        = Jac_W1;
        Jac_W1(1,36)  = 1;
        Jac_W1(2,37)  = 1;
        Jac_W1(3,38)  = 1;
        Jac_W2(1,42)  = 1;
        Jac_W2(2,43)  = 1;
        Jac_W2(3,44)  = 1;
        gradceq(1,:)  = Jac_W1.'*lure.W1/norm(lure.W1);
        gradceq(2,:)  = Jac_W2.'*lure.W2/norm(lure.W2);
        
        gradceq = gradceq.';                                               %Fmincon requires Jacobian.', not sure why.
    end  
end