function thetainit = controlledRandomSearch(crs_pars, J)

    % Unpack parameters 
    kmax    = crs_pars.kmax;
    tolCRS  = crs_pars.tolCRS;
    alpha   = crs_pars.alpha;
    tolMinJ = crs_pars.tolMinJ;
    
    % Initialize population
    disp('CRS: Started initialization')
    [R, JJ]          = initializePopulation(J, crs_pars);
    [maxJ, maxIndex] = max(JJ);                                            % Find the worst and best agents in the population
    [minJ, ~]        = min(JJ);   
    disp('CRS: Finished initialization')

    % Repetitively update the population until stopping conditions
    k  = 0; nf = 0; ns = 0; dontStop = true;
    
    while dontStop
        % Print status in the command window
        disp(['CRS: k = ',num2str(k), ', Conv = ',num2str(((maxJ - minJ)/maxJ),3), ', minJ = ', num2str(minJ,3)])

        % Draw a convergent primary point
        [rPrim,agts] = drawPrimaryPoint(J, R,crs_pars);   
        JPrim = J(rPrim);
        k     = k  + 1; 
        nf    = nf + 1;

        % Replace the worst agent, if the primary point is an improvement
        if JPrim < maxJ           
            R(:, maxIndex)   = rPrim;                                      % Replace the worst agent by the primary point
            JJ(maxIndex)     = JPrim;
            [maxJ, maxIndex] = max(JJ);                                    % Update the worst and best agents in the population
        
        % Calculate the secondary point from the selection
        elseif ( ns/nf < alpha )
            pDraw = R(:,agts);
            gBar  = mean(pDraw(:,1:end-1),2);                              % Compute the centroid of the first pn agents 
            rSec  = 0.5*(gBar + pDraw(:,end));                             % Find the secondary point
            JSec  = J(rSec);                                               % Find the secondary point costs

            % Check if the secondary point is non-convergent
            if JSec == 10000
                disp('CRS: Non-convergent secondary point')
            end
            
            % Replace the worst agent, if the secondary point is an improvement
            if JSec < maxJ 
                R(:, maxIndex)   = rSec;                                   % Replace the worst agent by the secondary point
                JJ(maxIndex)     = JSec;
                [maxJ, maxIndex] = max(JJ);                                % Update the worst and best agents in the population
                [minJ, ~] = min(JJ);
                ns = ns + 1;
                nf = nf + 1;
            end
        end     

        % Evaluate stopping conditions
        dontStop_1 = (maxJ - minJ)/maxJ > tolCRS; 
        dontStop_2 = k < kmax; 
        dontStop_3 = minJ > tolMinJ;
        dontStop   = all([dontStop_1;dontStop_2;dontStop_3]); 
    end

    % Define the best agend in the population as initial parameter vector
    [~, minIndex] = min(JJ);
    thetainit     = R(:,minIndex).';
end

function [rSub, agts] = drawPrimaryPoint(J, R, crs_pars)  

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%   Draw primary points that yield convergent Lur'e-type models.   %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   jPrim = 10000;
   
   % Draw again in case the primary point is non-convergent.
   while jPrim > 9999
        % Draw np+1 random agents from the population
        agts    = randi(crs_pars.N, [crs_pars.np+1,1]);
        pDraw   = R(:,agts);
        
        % Calculate primary point from this selection
        gBar    = mean(pDraw(:,1:end-1),2);                                % Compute the centoid of the first pn agents 
        rSub    = 2*gBar - pDraw(:,end);                                   % Substract the last agent

        % Constrain primary points to the feasible domain
        rSub = max(rSub,crs_pars.pLow);
        rSub = min(rSub,crs_pars.pHigh);    
        
        %  Calculate primary point costs to check their convergence
        jPrim = J(rSub);
   end
end

function [R, JJ] = initializePopulation(J, crs_pars)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Initialize a population R of N trial points that are     %%%
    %%%    1) uniformly distributed over search space X.           %%%
    %%%    2) correspond to a convergent model                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Unpack parameters
    N     = crs_pars.N;
    pLow  = crs_pars.pLow;
    pHigh = crs_pars.pHigh;
    np    = crs_pars.np;
    
    % Initialize random agents in the search space
    R   = pLow*ones(1,N) + ((pHigh - pLow)*ones(1,N)).*rand(np, N); 
        
    % Re-initialize non-convergent agents
    JJ = zeros(N,1);
    k = 1;
    for i = 1:N
        JJ(i) = J(R(:,i));
        while JJ(i) > 9999
            R(:,i) = pLow + (pHigh - pLow).*rand(np,1); 
            JJ(i)  = J(R(:,i));
        end
        % Display initialization progress
        if i/N*100 >= k*10
            disp(['CRS: Initialized ',num2str(k*10), ' %'])
            k  = k+1;
        end
    end
end