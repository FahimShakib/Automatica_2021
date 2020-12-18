function stop = plotFunction(x, optimValues, state, thetatrue)

    persistent plotdata
    
   %% Iteration logs initialization
   if strcmp(state,'init')
        plotdata.iteration      = [];
        plotdata.fval           = [];
        plotdata.firstorderopt	= [];
        plotdata.stepsize   	= [];
        
        plotdata.alpha0	= [];
        plotdata.Rlung	= [];
        plotdata.Rleak	= [];
        plotdata.Rlin	= [];
        plotdata.Clung	= [];
   end
   
   %% Iteration logs extension
   if strcmp(state,'iter') || strcmp(state,'done') 
       
       % Extend iteration logs
       plotdata.iteration       = [plotdata.iteration,      optimValues.iteration       ];
       plotdata.fval            = [plotdata.fval,           optimValues.fval            ];
       plotdata.firstorderopt   = [plotdata.firstorderopt,  optimValues.firstorderopt   ];
       
       if optimValues.iteration == 0
           plotdata.stepsize    = [plotdata.stepsize,       NaN                         ];
       else
          plotdata.stepsize     = [plotdata.stepsize,       optimValues.stepsize        ];
       end
       
       % Plot iteration nonlinearity
       lure     = thetaToLure(x);
       luretrue = thetaToLure(thetatrue);
%        
%         subplot(12,2,[1 3 5 7 9 11 13 15 17 19 21 23])
%             bar([lure.A ,luretrue.A;
%                  lure.B ,luretrue.B;
%                  lure.C ,luretrue.C;
%                  lure.D ,luretrue.D;
%                  lure.F ,luretrue.F;
%                  lure.G ,luretrue.G;
%                  lure.H ,luretrue.H;
%                  lure.L ,luretrue.L;
%                  lure.W1,luretrue.W1;
%                  lure.b1 ,luretrue.b1;
%                  lure.W2 ,luretrue.W2])
%             legend('\theta','\theta_{true}')
%             xticklabels({'A','B','C','D','F','G','H','L','W11','W12','b11','b12','W21','W22'})
        subplot(12,2,[2 4 6 8])
            semilogy(plotdata.iteration ,plotdata.fval,'kd','MarkerFaceColor',[1 0 1]);
            title(['Objective function value: ', num2str(optimValues.fval,4)])
            set(gca, 'YScale', 'log')            
            s = get(gca, 'Position');
            set(gca, 'Position', [s(1), s(2), s(3), s(4) * 0.6])
        subplot(12,2,[10 12 14 16])
            semilogy(plotdata.iteration, plotdata.firstorderopt,'kd','MarkerFaceColor',[1 0 1]);
            title(['First order optimality: ', num2str(optimValues.firstorderopt,4)])
            set(gca, 'YScale', 'log')
            s = get(gca, 'Position');
            set(gca, 'Position', [s(1), s(2), s(3), s(4) * 0.6])       
        subplot(12,2,[18 20 22 24])
            semilogy(plotdata.iteration,plotdata.stepsize,'kd','MarkerFaceColor',[1 0 1]);
            xlabel('Iteration number')
            title(['Step size: ', num2str(optimValues.stepsize,4)])
            s = get(gca, 'Position');
            set(gca, 'Position', [s(1), s(2), s(3), s(4) * 0.6])        
   end
   
   %% Iteration logs ending
   if strcmp(state,'done')
       stop = 1;
       assignin('base', 'plotdata', plotdata);
   else
       stop = 0;
   end

end
