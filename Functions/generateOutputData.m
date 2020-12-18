function IOData = generateOutputData(lure, IData, mtf_pars, nfi_pars, method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Periodic steady-state solution simulation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isstruct(lure) == 0
    lure = thetaToLure(lure);
end

% Solve the periodic steady-state solution
switch method
    case 'nfi'
        IOData = NFI_Solution(lure, IData, nfi_pars);
    case 'mtf'
        IOData = MTF_Solution(lure, IData, mtf_pars);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Compare NFI to MTF for a single solution   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compare_nfi_mtf = false;

if compare_nfi_mtf
    tic
    IOData_NFI = NFI_Solution(lure, IData, nfi_pars);
    toc
    tic
    IOData_MTF = MTF_Solution(lure, IData, mtf_pars);
    toc
    tic
    IOData_NFISim = NFISim_Solution(lure, IData, nfi_pars, mtf_pars);
    toc

    figure
    subplot(5,1,1)
        plot(IOData_MTF.t,IOData_MTF.w)  
        hold on
        plot(IOData_MTF.t,IOData_NFI.w,'--')
        plot(IOData_NFISim.t,IOData_NFISim.w,'-.')  
        legend('MTF','NFI','NFI\_Simulink')
        ylabel('w')
    subplot(5,1,2)
        plot(IOData_MTF.t,IOData_MTF.z)  
        hold on
        plot(IOData_MTF.t,IOData_NFI.z,'--')
        plot(IOData_NFISim.t,IOData_NFISim.z,'-.')  
        ylabel('z')
    subplot(5,1,3)
        plot(IOData_MTF.t,IOData_MTF.y)  
        hold on
        plot(IOData_MTF.t,IOData_NFI.y,'--')
        plot(IOData_NFISim.t,IOData_NFISim.y,'-.')  
        ylabel('y')
    subplot(5,1,4)
        plot(IOData_MTF.t,IOData_MTF.x)  
        hold on
        plot(IOData_MTF.t,IOData_NFI.x,'--')
        plot(IOData_NFISim.t,IOData_NFISim.x,'-.')  
        ylabel('x')
    subplot(5,1,5)
        plot(IOData_MTF.t,IOData_MTF.u)  
        hold on
        plot(IOData_MTF.t,IOData_NFI.u,'--')
        plot(IOData_NFISim.t,IOData_NFISim.u,'-.')     
        ylabel('u')       
        xlabel('t [s]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Compare NFI approaches   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compare_runtimes = false;

if compare_runtimes
    
    % NFI variable output stepper with zero order hold
    nfi_pars.hold    = 'zero order';
    nfi_pars.stepper = 'variable';
    tic
    IOData_NFI_zoh_v = NFI_Lure1(lure, IData, nfi_pars);
    t_run.NFI_zoh_v  = toc;
    
    % NFI dense output stepper with zero order hold
    nfi_pars.hold    = 'zero order';
    nfi_pars.stepper = 'dense';
    tic
    IOData_NFI_zoh_d = NFI_Lure1(lure, IData, nfi_pars);
    t_run.NFI_zoh_d  = toc;
 
    % NFI variable output stepper with first order hold    
    nfi_pars.hold    = 'first order';
    nfi_pars.stepper = 'variable';
    tic
    IOData_NFI_foh_v = NFI_Lure1(lure, IData, nfi_pars);
    t_run.NFI_foh_v  = toc;

    % NFI dense output stepper with first order hold    
    nfi_pars.hold    = 'first order';
    nfi_pars.stepper = 'dense';
    tic
    IOData_NFI_foh_d = NFI_Lure1(lure, IData, nfi_pars);
    t_run.NFI_foh_d  = toc;
    
    % MTF
    tic
    IOData_MTF = MTF_Lure1(lure, IData, mtf_pars);
    t_run.MTF  = toc;
    
    % compare nfi to MTF solutions
    figure(12)
    subplot(2,1,1)
        plot(IOData_NFI_zoh_v.t, IOData_NFI_zoh_v.w, '--b.')
        hold on
        plot(IOData_NFI_foh_v.t, IOData_NFI_foh_v.w, '--r.')
    subplot(2,1,2)
        scatter(IOData_NFI_zoh_v.t, IOData_NFI_zoh_v.z, 'b.')
        plot(IOData_NFI_zoh_d.t, IOData_NFI_zoh_d.z, 'b--')
        hold on
        scatter(IOData_NFI_foh_v.t, IOData_NFI_foh_v.z, 'r.')
        plot(IOData_NFI_foh_d.t, IOData_NFI_foh_d.z, 'r--')
        plot(IOData_NFI_foh_d.t, IOData_MTF.z,       'g--')                % ! a bit sloppy to use the NFI Time here for plotting!
        legend('ZOH: V','ZOH: D','FOH: V','FOH: D', 'MTF')
        
        figure(13)
        subplot(3,1,1)
            plot(IOData_NFI_zoh_d.t, IOData_NFI_zoh_d.w)
            hold on
            stairs(IOData_NFI_foh_d.t, IOData_NFI_zoh_d.w)
            ylabel('$\bar{w}$','Interpreter','Latex')
        subplot(3,1,2)
            plot(IOData_NFI_zoh_d.t, IOData_NFI_zoh_d.z)
            hold on
            plot(IOData_NFI_foh_d.t, IOData_NFI_foh_d.z)
            plot(IOData_NFI_foh_d.t, IOData_MTF.z, '--')                   % ! a bit sloppy to use the NFI Time here for plotting!
            legend('NFI: ZOH','NFI: FOH','MTF')
            ylabel('$\bar{z}$','Interpreter','Latex')
        subplot(3,1,3)
            plot(IOData_NFI_zoh_d.t, IOData_MTF.z - IOData_NFI_zoh_d.z)
            hold on
            plot(IOData_NFI_zoh_d.t, IOData_MTF.z - IOData_NFI_foh_d.z)
            legend('ZOH','FOH')    
            ylabel('$\bar{z}_{MTF} - \bar{z}_{NFI}$ ','Interpreter','Latex')     
end

end

function [IOData] = NFI_Solution(luresys, IData, nfi_pars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Periodic steady-state solution  %%%
%%%             solved by             %%%
%%%    forward numeric integration    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implicit settings
plot_fig  = true;
IData_sol = IData;
Period    = length(IData.t);

% Run the algorithm
yerror   = inf;
k        = 1;

while (yerror > nfi_pars.convtol)
    % Solve the first k periods of IData (it may be better to solve entire sequence at once)
    IData_sol.w = IData.w_tot(1:k*Period);
    IData_sol.t = IData.t_tot(1:k*Period);
    
    % Solve ODE by numeric forward integration
    options	= odeset('RelTol', nfi_pars.RelTol,'AbsTol', nfi_pars.AbsTol);
    ode     = @(t,x) lure_ex(t, x, luresys, IData_sol, nfi_pars);
    
    switch nfi_pars.stepper
        case 'dense'
            [t, x]	= ode45(ode, IData_sol.t, luresys.x0, options);
            w       = IData_sol.w.';
            
        case 'variable'
            [t, x]	= ode45(ode, [0 IData_sol.t(end)], luresys.x0, options);
            
            w = zeros(size(x)).';
            for i = 1:length(t)
                w(i) = intersample_behavior(IData_sol, t(i), nfi_pars);
            end
    end
    
    % reconstruct outputs
    x   = x.';
    y   = luresys.C*x + luresys.D*w;
    phi = luresys.W2.'*tanh(luresys.W1*y + luresys.b1) + luresys.b2;
    u   = -phi;
    z   = luresys.F*x + luresys.G*u + luresys.H*w;
    
    % Check for convergence
    if strcmp(nfi_pars.stepper,'dense') == 1
        Ts   = IData_sol.t(2);
        z_ss = z(t >= (k-1)*Period*Ts);
        Z_ss = fft(z_ss);
        
        if k>1
            errordat{k-1}   = abs((z_ss - z_ss_old)./z_ss_old);
            yerror          = norm(Z_ss - Z_ss_old)/norm(Z_ss_old);
            nrm(k-1)        = yerror;
        end
        
        Z_ss_old = Z_ss;
        z_ss_old = z_ss;
    end
    k = k + 1;
    disp(['NFI: k = ',num2str(k),', yerror = ', num2str(yerror)])
end

% Define output struct for the periodic solution (last period).
Ts  = IData_sol.t(2);
idx = find(t >= (k-2)*Period*Ts);

IOData.t       = t(idx);
IOData.w       = w(idx).';
IOData.x       = x(:,idx).';
IOData.y       = y(idx).';
IOData.u       = u(idx).';
IOData.z       = z(idx).';
IOData.w_freqs = IData.w_freqs;

% Show how the simulation converges to the steady-state solution
if (plot_fig) && strcmp(nfi_pars.stepper,'dense') == 1
    figure
    ii = 1;
    for i = 2:k-1
        semilogy(IOData.t, errordat{i-1})
        hold on
        plot([IOData.t(1),IOData.t(end)],[nrm(i-1),nrm(i-1)])
        legends{ii} = ['$k = ', num2str(i),': \frac{|z_k-z_{k-1}|}{|z_{k-1}|}$'];
        legends{ii+1} = ['$k = ', num2str(i),':||\frac{z_k-z_{k-1}}{z_{k-1}}||_2 $'];
        legend(legends,'Interpreter','Latex')
        ii = ii + 2;
        xlim([IOData.t(1), IOData.t(end)])
    end
end

end

function xdot = lure_ex(t, x, lure, IData, nfi_pars)

    % unpack parameters
    A  = lure.A;
    B  = lure.B;
    C  = lure.C;
    D  = lure.D;
    L  = lure.L;
    W1 = lure.W1;
    b1 = lure.b1;
    W2 = lure.W2;
    b2 = lure.b2;

    % Intersample behavior
    w   = intersample_behavior(IData, t, nfi_pars);

    % Nonlinear SS equations
    y    = C*x + D*w;
    phi  = W2.'*tanh(W1*y + b1) + b2;
    u    = -phi;
    xdot = A*x + B*u + L*w;
end


function w = intersample_behavior(IData,t, nfi_pars)

    % Find leading index
    Ts  = IData.t(2);
    ind = ceil(t/Ts);

    switch nfi_pars.hold
        % Calculate excitation based on zero-order hold
        case 'zero order'
            if ind >= length(IData.w)
                w = IData.w(end);
            elseif ind <= 0
                w = IData.w(1);
            else
                w = IData.w(ind);
            end
            % Calculate excitation based on first-order hold
        case 'first order'
            if ind >= length(IData.w)
                w = IData.w(end);
            elseif ind <= 0
                w = IData.w(1);
            else
                w = interp1([IData.t(ind); IData.t(ind+1)],...
                    [IData.w(ind); IData.w(ind+1)],...
                    t);
            end
    end
end