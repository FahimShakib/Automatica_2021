function [SSGradData] = MTF_Sensitivity(lure, IData_sens, mtf_pars)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Lur'e-type system:     %%%     Loop-transformed Lur'e-type system:   %%%
%%%                              %%%                                           %%%
%%%   x_dot = A*x + B*Ũ + Wbar   %%%   x_dot = A_tf*x + B*Û + I*Wbar +B_tf*Ȳ   %%%
%%%   y     = C*x       + Ȳ      %%%   y     =    C*x                  + I*Ȳ   %%%
%%%   z     = F*x + G*Ũ + Zbar   %%%   z     = F_tf*x + G*Û + I*Zbar +G_tf*Ȳ   %%%     
%%%   Ũ     = -Psi*y             %%%   Û     = - (Psi + gamma/2)*y             %%%
%%%                              %%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loop-transform the Lure-type system
lure.A_tf      = lure.A - lure.gamma/2*lure.B*lure.C;
lure.F_tf      = lure.F - lure.gamma/2*lure.G*lure.C;
lure.B_tf      =        - lure.gamma/2*lure.B;
lure.G_tf      =        - lure.gamma/2*lure.G;
lure.nonlin    = @(y_theta) IData_sens.Psi.*y_theta;
lure.nonlin_tf = @(y_theta) lure.nonlin(y_theta) - lure.gamma/2*y_theta;

%% Define LTI-block transfer functions
% Linear transfers to y
G_yu    = ss(lure.A_tf,	     lure.B, lure.C,          0);
G_ywbar = ss(lure.A_tf,      eye(4), lure.C, zeros(1,4));                  % eye or ones?   
G_yybar = ss(lure.A_tf,   lure.B_tf, lure.C,          1);
G_yzbar = ss(lure.A_tf,  zeros(4,1), lure.C,          0);

% Linear transfers to z
G_zu    = ss(lure.A_tf,     lure.B, lure.F_tf,     lure.G);
G_zwbar = ss(lure.A_tf,     eye(4), lure.F_tf, zeros(1,4));                % eye or ones?
G_zybar = ss(lure.A_tf,  lure.B_tf, lure.F_tf,  lure.G_tf);
G_zzbar = ss(lure.A_tf, zeros(4,1), lure.F_tf,          1);

% Linear transfers to x
I2      = eye(4);

G_xu    = ss(lure.A_tf,    lure.B, I2, zeros(4,1));
G_xwbar = ss(lure.A_tf,   eye(4), I2,  zeros(4,4));                         % eye or ones?
G_xybar = ss(lure.A_tf, lure.B_tf, I2, zeros(4,1));

%% Evaluate LTI-block transfer functions

% Define frequency vector
n = mtf_pars.n;
f = (0 : 1 : mtf_pars.n-1) * (1/mtf_pars.T);                              
                                                                           
% Calculate Fourier coefficients of G_eu_bar and G_er_bar 
G_yu_vec    = squeeze(freqresp(G_yu,   [f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_ywbar_vec = squeeze(freqresp(G_ywbar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_yybar_vec = squeeze(freqresp(G_yybar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_yzbar_vec = squeeze(freqresp(G_yzbar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));       

G_zu_vec    = squeeze(freqresp(G_zu,   [f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_zybar_vec = squeeze(freqresp(G_zybar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_zzbar_vec = squeeze(freqresp(G_zzbar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_zwbar_vec = squeeze(freqresp(G_zwbar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));

G_xu_vec    = squeeze(freqresp(G_xu,   [f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_xwbar_vec = squeeze(freqresp(G_xwbar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));
G_xybar_vec = squeeze(freqresp(G_xybar,[f(1:1:n/2+1), -f(n/2:-1:2)]*2*pi));

if size(G_xu_vec,1) ~= n
    G_xu_vec    = G_xu_vec.';
    G_xybar_vec = G_xybar_vec.';
end
                                                                           
%% Execute the MTF Algorithm

% Fourier transform the excitation
Wi_fft  = fft(IData_sens.Wi);            
Yi_fft  = fft(IData_sens.Yi);
Zi_fft  = fft(IData_sens.Zi);

% Calculate initial response assuming zero nonlinearity output (u = 0)
Y_theta0 = sum(G_ywbar_vec.'.*Wi_fft, 2) + G_yybar_vec.*Yi_fft;                          % ASK FAHIM!
Z_theta0 = sum(G_zwbar_vec.'.*Wi_fft, 2) + G_zybar_vec.*Yi_fft + G_zzbar_vec.*Zi_fft;    % ASK FAHIM!

% Pre-allocate memory
y_theta0 = zeros(n,1);                                                     % initial guess for steady state response for e
Y_theta  = zeros(n,1);                                                     % fft of error signal
Y_old    = zeros(n,1);                                                     % used for determining convergence in MTF algorithm

% Run the MTF algorithm
yerror   = 1;                                                              % set initial error > tol
k        = 0;                                                              % set iteration index of while loop
y_theta	 = y_theta0;                                                       % output signal
while k < mtf_pars.max_iter && yerror > mtf_pars.convtol
    if any(isnan(Y_theta))
        disp('Response diverges')
        break;
    end
    utilde  = -lure.nonlin_tf(y_theta);
    UTilde  = fft(utilde);                                                 % fft of output of nonlinearity
    Y_theta = Y_theta0 + G_yu_vec.*UTilde;                                 % linear dynamics in frequency domain , superposition ?!
    y_theta	= real(ifft(Y_theta));                                         % use real value because there is a very small imaginary part present (1e-21)
    Z_theta = Z_theta0 + G_zu_vec.*UTilde;
    z_theta = real(ifft(Z_theta));
    
    % error based on Fourier coefficients to check convergence
    yerror      = norm(Y_theta - Y_old)/norm(Y_old);
    nrm(k+1)    = yerror;
    Y_old       = Y_theta;
    
    % update the iteration index
    k     = k + 1;
end

% Simulate steady-state state trajectory
for k = 1:length(lure.A)
    X_theta(k,:) = sum(squeeze(G_xwbar_vec(:,k,:)).'.*Wi_fft, 2) + G_xu_vec(:,k).*UTilde + G_xybar_vec(:,k).*Yi_fft;
    x_theta(k,:) = real(ifft(X_theta(k,:)));
end

% Formulate output data
SSGradData.x = x_theta.';
SSGradData.y = y_theta;
SSGradData.z = z_theta;
SSGradData.t = 0:mtf_pars.T/mtf_pars.n:(mtf_pars.T-mtf_pars.T/mtf_pars.n);

%% Optional: Analysis on MTF in case it didn't converge
show_figure = false;
if k == mtf_pars.max_iter
    disp(['MTF Sens: Max Iter Reached: yerror = ', num2str(yerror)])
    if show_figure
        % Show MTF-history and analyze last vs second last iteration
        figure
        subplot(4,1,1)
            plot(IData.t,real(ifft(Y_theta)))
            hold on
            plot(IData.t,real(ifft(Y_old)),'--')
            legend('y(t)','y_{old}(t)')
            xlabel('t [s]')
        subplot(4,1,2)
            title('')
            plot(abs(Y_theta-Y_old))
            legend('|Y-Y_{old}|')
            xlabel('samples [-]')
        subplot(4,1,3)
            plot(Y_theta,'.')
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