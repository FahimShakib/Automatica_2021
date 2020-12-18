function [IO] = MTF_Solution(lure, IData, mtf_pars)

%% System definition

% Loop-transform the lure-type system
lure.A_tf      = lure.A - lure.gamma/2*lure.B*lure.C;
lure.F_tf      = lure.F - lure.gamma/2*lure.G*lure.C;
lure.L_tf      = lure.L - lure.gamma/2*lure.B*lure.D;
lure.H_tf      = lure.H - lure.gamma/2*lure.G*lure.D;
lure.nonlin_tf = @(y) lure.nonlin(y) - lure.gamma/2*y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Lur'e-type system:     %%%
%%%                             %%%
%%%   x_dot = A*x + B*u + L*w   %%%
%%%   y     = C*x       + D*w   %%%
%%%   z     = F*x + G*u + H*w   %%%
%%%   u     = -phi(y)           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define linear transfers
xdim = size(lure.A,2);
G_yw = ss(lure.A_tf, lure.L_tf,    lure.C,        lure.D);
G_yu = ss(lure.A_tf,    lure.B,    lure.C,             0);
G_zw = ss(lure.A_tf, lure.L_tf, lure.F_tf,     lure.H_tf);
G_zu = ss(lure.A_tf,    lure.B, lure.F_tf,        lure.G);
G_xw = ss(lure.A_tf, lure.L_tf, eye(xdim), zeros(xdim,1));
G_xu = ss(lure.A_tf,    lure.B, eye(xdim), zeros(xdim,1));

%% MTF algorithm setup

% Extract parameters
n = mtf_pars.n;
T = mtf_pars.T;
f = (0 : 1 : n-1) * (1/T);                                                 % Corresponding frequency vector with frequency steps 1/T

% Evaluate linear transfers
G_yu_vec = squeeze(freqresp(G_yu,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));       % Notes: 1) These spectrum repeat themselves
G_yw_vec = squeeze(freqresp(G_yw,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));       %        2) The negative frequencies are stored in fft order behind the positive ones)
G_zu_vec = squeeze(freqresp(G_zu,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));       %        3) Use squeeze to make it a vector
G_zw_vec = squeeze(freqresp(G_zw,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));       
G_xu_vec = squeeze(freqresp(G_xu,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));
G_xw_vec = squeeze(freqresp(G_xw,[f(1:1:n/2+1),-f(n/2:-1:2)]*2*pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%   Calculate initial response of linear dynamics     %%%
%%%     1) Subject to excitation w(t)                   %%%
%%%     2) Assume u(t) = 0,                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W  = fft(IData.w);                                                         % Fourier coefficients W of excitation w(t). 
Y0 = G_yw_vec.*W;
Z0 = G_zw_vec.*W;

% MTF algorithm execution
y0      = zeros(n,1);                                                      % initial guess for steady state response for e
y       = y0;                                                              % output signal (error e in this case)
Y       = zeros(n,1);                                                      % fft of error signal
Y_old   = zeros(n,1);                                                      % used for determining convergence in MTF algorithm
y_old   = zeros(n,1);                                                      % used for determining convergence in MTF algorithm
U       = zeros(n,1);                                                      % fft of output of the nonlinearity

yerror  = 1;                                                               % set initial error > tol
k       = 0;
while k < mtf_pars.max_iter && yerror > mtf_pars.convtol
    if any(isnan(Y))
        disp('Response diverges')
        break;
    end
    u     = -lure.nonlin_tf(y.').';
    U     = fft(u);                                                        % fft of output of nonlinearity
    Y     = Y0 + G_yu_vec.*U;                                              % linear dynamics in frequency domain , superposition ?!
    y     = real(ifft(Y));                                                 % use real value because there is a very small imaginary part present (1e-21)
    Z     = Z0 + G_zu_vec.*U;
    z     = real(ifft(Z));
    
    % Error based on Fourier coefficients to check convergence
    yerror   = norm(Y-Y_old)/norm(Y_old);
    %yerror   = norm((y-y_old)./y_old);
    nrm(k+1) = yerror;
    if k < mtf_pars.max_iter - 1
    Y_old    = Y;
    y_old    = y;
    end
    
    k     = k + 1;                                                         % update the iteration index
end

% calculate periodic steady-state state solutions
for k = 1:size(lure.A,2)
    X(k,:) = G_xw_vec(k,:).'.*W + G_xu_vec(k,:).'.*U;
    x(k,:) = real(ifft(X(k,:)));
end

% Formulate steady state solution
IO.w        = IData.w;
IO.w_freqs  = IData.w_freqs;
IO.t        = IData.t;
IO.y        = y;
IO.x        = x.';
IO.u        = -lure.nonlin(IO.y);
IO.z        = z;

%% Optional: Analysis on MTF in case it didn't converge
show_figure = false;
if k == mtf_pars.max_iter
    disp(['MTF Sol: Max Iter Reached: yerror = ', num2str(yerror)])
    if show_figure
        % Show MTF-history and analyze last vs second last iteration
        figure
        subplot(4,1,1)
            plot(IData.t,real(ifft(Y)))
            hold on
            plot(IData.t,real(ifft(Y_old)),'--')
            legend('y(t)','y_{old}(t)')
            xlabel('t [s]')
        subplot(4,1,2)
            title('')
            plot(abs(Y-Y_old))
            legend('|Y-Y_{old}|')
            xlabel('samples [-]')
        subplot(4,1,3)
            plot(Y,'.')
            hold on
            plot(Y_old,'o')
            xlabel('Re [-]')
            ylabel('Im [-]')
            legend('Y','Y_{old}')
        subplot(4,1,4)
            plot(nrm)
            title('norm history')
            xlabel('k')

        % Analyze third circle-criterion LIKE condition
        fgr = logspace(-6,2,1e3);
        G_yu_check = frd(G_yu,fgr*2*pi);
        diff       = abs(G_yu_check.ResponseData(1)) - 1/(lure.gamma/2);     

        figure
            bodemag_custom({G_yu_check},{'G_{yu}'})
            hold on
            plot([fgr(1),fgr(end)],mag2db([1/(lure.gamma/2),1/(lure.gamma/2)]));
            ax2 = gca;
            ax2.Legend.String{2} = '1/\gamma';
            xlim([fgr(1) fgr(end)])
            title(['Difference: $| ||G_{uy}||_{\infty} - \frac{1}{\bar{\gamma}} | = $', num2str(diff)], 'Interpreter','Latex')
    end
end
end