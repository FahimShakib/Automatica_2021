function [theta] = lureToTheta(lure)

    theta = [lure.A(:);
             lure.B;
             lure.L;
             lure.C(:);  
             lure.D;      
             lure.F(:);
             lure.G;
             lure.H;
             lure.W1; 
             lure.b1; 
             lure.W2];
end

