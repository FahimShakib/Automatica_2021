function stop = outfun(in,optimValues,state)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Output function for the nonlinear optimization method   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    persistent history
    stop = false;
 
     switch state
         case 'init'
             hold on
             history.x = [];
             history.f = [];
             history.t = [];  

         case 'iter'
            history.f = [history.f; optimValues.fval];
            history.x = [history.x; in(:)'];                               % Ensure in is a row vector
            %history.t = [history.t; toc];
            %tic;

         case 'done'
             hold off
             assignin('base','history',history);
         otherwise
     end
end

