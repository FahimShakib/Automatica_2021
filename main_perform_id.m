clc; clear all; close all; addpath('Functions');addpath('HelperFunctions')

global t_logs scaleswitch

% Simulation method: 1 for MTF, 0 for NFI
sim_method = 1;

% Define noise scenario: 0 for no noise, 1 for SNR of 10 dB, 2 for SNR of
% 20 dB, 3 for SNR of 40 dB and 4 for SNR of 60 dB
noise_scenario = 4;

% Compare timing between MTF and NFI (takes a lot of time)
comptime = 0;

%            %
%            %
% Save Data? %
%            %
%            %
save_data = 0;


%% Define the ground-truth Lur'e-type model that is to be identified

% Define model parameters
m1 = 0.15;
m2 = 0.45;
d1 = 0.4;
d2 = 0.4;
k1 = 1000;
k2 = 2200;
W1 = [0.5;    1;  1];
W2 = [  1;  0.5;  1];
b1 = [  3;   -3;  0];

% Define true lure-type system
luretrue.A       = [           0,           1,      0,      0;
                     -(k1+k2)/m1, -(d1+d2)/m1,  k2/m1,  d2/m1;
                               0,           0,      0,      1;
                           k2/m2,       d2/m2, -k2/m2, -d2/m2];
luretrue.B       = [0; 1; 0; 0];
luretrue.L       = [0; 0; 0; 1];
luretrue.C       = [20, 0, 0, 0];
luretrue.D       = 0;
luretrue.F       = [0, 0, 1, 0];
luretrue.G       = 0;
luretrue.H       = 0;
luretrue.W1      = W1./norm(W1);  
luretrue.W2      = W2./norm(W2);    
luretrue.b1      = b1;                       
luretrue.b2      = -luretrue.W2.'*tanh(luretrue.b1);
luretrue.gamma   = luretrue.W2.'*luretrue.W1;
luretrue.x0      = [1;0;1.5;0];

luretrue.nonlin     = @(y) phi_func(y, luretrue);
luretrue.dphidtheta = @(y) dPhidTheta_func(y, luretrue);
luretrue.dphidy     = @(y) dPhidY_func(y, luretrue); 

% Define true parameter vector
thetatrue  = lureToTheta(luretrue);

%% Identification algorithm user settings

% Multisine excitation signal settings
input.Ts        = 2e-3;                                                    % Sampling time       [s]
input.Period    = 1000;                                                    % Signal length       [samples]
input.NumPeriod = 20;                                                      % Number of periods   [-]
input.Range     = [-2000 2000];                                            % Amplitude           [-]
input.NumTrials = 50;                                                      % Random-phase trials [-]
input.df        = 0.5;                                                     % Delta frequency     [Hz]
input.fMin      = 0.5;                                                     % Slowest frequency   [Hz]
input.fMax      = 200;                                                     % Fastest frequency   [Hz]

% Controlled random search parameter initialization settings
crs_pars.np       = 44;                                                    % Number of parameters to be estimated
crs_pars.N        = 25;                                                    % Number of agents in the population
crs_pars.alpha    = 0.5;                                                   % Success-rate threshold
crs_pars.tolCRS   = 1e-3;                                                  % Maximum. difference between the best and worst agent
% crs_pars.tolMinJ  = 1e-4;                                                  % Cost threshold for initial parameter vector
crs_pars.kmax     = 500;                                                   % Maximum number of iterations
crs_pars.convLow  = [-inf*ones(1,35),zeros(1,3),-inf(1,3),zeros(1,3)].';   % Hard lower-bound on the parameters
crs_pars.convHigh = inf(1,44).';                                           % Hard upper-bound on the parameters
crs_pars.pLow     = max(crs_pars.convLow   + 1e-3, thetatrue - 0.1);       % Actual lower-bound on the parameters
crs_pars.pHigh    = min(crs_pars.convHigh  - 1e-3, thetatrue + 0.1);       % Actual upper-bound on the parameters

% MTF algorithm settings
mtf_pars.max_iter = 1000;                                                  % Max number of iterations
mtf_pars.convtol  = 1e-16;                                                 % Convergence accuracy
mtf_pars.n        = input.Period;                                          % Number of points used to approximate response (>n, more accurate)
mtf_pars.T        = input.Period*input.Ts;                                 % Period time of input signal

% Numeric forward integration settings
nfi_pars.hold     = 'first order';                                         % Intersample behavior of the excitation signal
nfi_pars.stepper  = 'dense';                                               % Dense or variable output stepper?
nfi_pars.P        = 20;                                                    % Number of simulated periods
nfi_pars.convtol  = 1e-6;                                                  % Convergence accuracy
nfi_pars.RelTol   = 1e-12;                                                 % Relative error stepper tolerance
nfi_pars.AbsTol   = 1e-12;                                                 % Absolute error stepper tolerance        
nfi_pars.comptime = comptime;                                              % Run alternating method in parallel to compare computation time
nfi_pars.n        = mtf_pars.n;
nfi_pars.T        = mtf_pars.T;

% Initialize timing-logs if required
if nfi_pars.comptime
        t_logs.lure1.mtf = [];
        t_logs.lure1.nfi = [];
        t_logs.lure2.mtf = [];
        t_logs.lure2.nfi = [];
end

%% Algorithm set-up

% Define the cost-function by simulating IO-behavior
IData  = generateInputData(input);
SSData = generateOutputData(luretrue, IData, mtf_pars, nfi_pars, 'mtf');

% Add noise to the output
switch noise_scenario
    case 0 % No noise (do nothing)
        S = -inf;
    case 1 % 10 dB SNR 
        S = 10;
    case 2
        S = 20;
    case 3
        S = 40;
    case 4
        S = 60;
end
vare        = sqrt(var(SSData.z-mean(SSData.z))/(10^(S/10)));
if noise_scenario
    e           = randn(length(IData.w),1)*vare;
else 
    e           = zeros(length(IData.w),1);
end
SSData.z0   = SSData.z;
SSData.z    = SSData.z0 + e;

% Adapt the search criterion of the CRS based on the amount of noise
if noise_scenario
    crs_pars.tolMinJ  = vare^2*10;                                         % Cost threshold for initial parameter vector
else
    crs_pars.tolMinJ  = 1e-4;                                              % Cost threshold for initial parameter vector
end

% Cost-function definition
if sim_method 
    J      = @(theta) costFunction(theta, mtf_pars, nfi_pars, SSData, 'mtf');
else
    J      = @(theta) costFunction(theta, mtf_pars, nfi_pars, SSData, 'NFI');
end
    
display(['Cost at true model parameters ' num2str(J(thetatrue))])

% Parameter-vector initialization
init_switch = 'crs';
switch init_switch
    case 'true'
        thetainit = thetatrue;
    case 'close'    
        Jc = inf;
        while (Jc > 9999)
            dtheta    = -1 + 2*rand(size(thetatrue));
            thetainit = thetatrue + 0.01*dtheta;
            Jc = J(thetainit);
        end 
    case 'crs'
        scaleswitch = true;                                                % Scale parameters online such that ||W_1|| = ||W_2|| = 1
        thetainit_unscaled = controlledRandomSearch(crs_pars, J);
        lureinit           = thetaToLure(thetainit_unscaled);
        thetainit          = lureToTheta(lureinit);                    
end

if save_data
    save('SNR20_CRS')
end

%% Gradient-based search using fmincon
scaleswitch = false;                                                       % Scaling already incorporated in fmincon constraint function

% Optimization options
options                           = optimoptions(@fmincon);
options.Display                   = 'final-detailed';
options.GradObj                   = 'on';
options.GradConstr                = 'on';
options.DerivativeCheck           = 'off';
options.MaxIter                   = 500;
options.MaxFunctionEvaluations    = 1e6;
options.TolFun                    = 1e-14;
options.FiniteDifferenceType      = 'central';
options.FiniteDifferenceStepSize  = 1e-7;
options.OptimalityTolerance       = 1e-14;
options.ConstraintTolerance       = 1e-6;
options.ObjectiveLimit            = J(thetatrue)*0.99;
options.PlotFcn                   = @(x, optimValues, state) plotFunction(x, optimValues, state, thetatrue);
options.OutputFcn                 = @outfun;

% Timing logs initialization for running NFI parallel to MTF
nfi_pars.comptime = comptime;

% Cost-function updated definition
J = @(theta) costFunction(theta, mtf_pars, nfi_pars, SSData, 'mtf');

tStart = tic;
[thetafinal,~,~,~,~,~,~] = fmincon(J,thetainit,[],[],[],[],crs_pars.pLow,crs_pars.pHigh,@mycon,options);
tel = toc(tStart);

% Define output as a lure-type model
lurefinal = thetaToLure(thetafinal);

%% Extract similarity transformation from identified LTI block
dummy = true;
if dummy
    T = eye(4);  
else
    T = luretrue.B*pinv(lurefinal.B);
end
Tinv = inv(T);

lurefinal_trans   = lurefinal;
lurefinal_trans.A = T*lurefinal.A*Tinv;
lurefinal_trans.B = T*lurefinal.B;
lurefinal_trans.L = T*lurefinal.L;
lurefinal_trans.C = lurefinal.C*Tinv;
lurefinal_trans.F = lurefinal.F*Tinv;

%% Visualize results

% % Visualize timing results
% if nfi_pars.comptime
%     its = size(t_logs.lure1(1,:),2);
%     figure
%     semilogy(t_logs.lure1(1,:))
%     hold on
%     plot([1 its], [mean(t_logs.lure1(1,:)) mean(t_logs.lure1(1,:))],'--')
%     semilogy(t_logs.lure1(2,:))
%     hold on
%     plot([1 its], [mean(t_logs.lure1(2,:)) mean(t_logs.lure1(2,:))],'--')
%     xlabel('Iterations [-]')
%     ylabel('Time [s]')
%     legend('MTF','MTF-mean','NFI','NFI-mean')
% end

% SS-solutions
lureinit  = thetaToLure(thetainit);
SSLInit   = generateOutputData(lureinit,  IData, mtf_pars, nfi_pars, 'mtf');
SSLTrue   = generateOutputData(luretrue,  IData, mtf_pars, nfi_pars, 'mtf');
SSLFinal  = generateOutputData(lurefinal, IData, mtf_pars, nfi_pars, 'mtf'); 
SSLFinalT = generateOutputData(lurefinal_trans, IData, mtf_pars, nfi_pars, 'mtf'); 

% Plot SS-solutions
figure
subplot(5,1,1)
    plot(SSLInit.t,SSLInit.w)
    hold on
    plot(SSLFinal.t,SSLFinal.w)   
    plot(SSLFinalT.t,SSLFinalT.w)   
    plot(SSLTrue.t,SSLTrue.w,'--')
    legend('Init','Final','Final-trans','True')
    ylabel('w')
subplot(5,1,2)
    plot(SSLInit.t,SSLInit.z)
    hold on
    plot(SSLFinal.t,SSLFinal.z)   
    plot(SSLFinalT.t,SSLFinalT.z)   
    plot(SSLTrue.t,SSLTrue.z,'--')
    ylabel('z')
subplot(5,1,3)
    plot(SSLInit.t,SSLInit.y)
    hold on
    plot(SSLFinal.t,SSLFinal.y)   
    plot(SSLFinalT.t,SSLFinalT.y)   
    plot(SSLTrue.t,SSLTrue.y,'--')
    ylabel('y')
subplot(5,1,4)
    plot(SSLInit.t,SSLInit.x)
    hold on
    plot(SSLFinal.t,SSLFinal.x)   
    plot(SSLFinalT.t,SSLFinalT.x)   
    plot(SSLTrue.t,SSLTrue.x,'--')  
    ylabel('x')
subplot(5,1,5)
    plot(SSLInit.t,SSLInit.u)
    hold on
    plot(SSLFinal.t,SSLFinal.u)   
    plot(SSLFinalT.t,SSLFinalT.u)   
    plot(SSLTrue.t,SSLTrue.u,'--')  
    ylabel('u')       
    xlabel('t [s]')
        
% Plot nonlinearity
yrange = -15:0.001:15;
figure
    plot(yrange,lureinit.nonlin(yrange))
    hold on
    plot(yrange,lurefinal.nonlin(yrange))
    plot(yrange,luretrue.nonlin(yrange),'--')
    legend('$\theta_{init}$','$\hat{\theta}_N$','$\theta_0$')
    ylabel('$\varphi$','Interpreter','Latex')
    xlabel('$y$','Interpreter','Latex')
    grid minor
    set(findobj(gcf,'Type','Line'),'linewidth',2)

% Define LTI transfer functions
fgrid           = 2*pi*(0.01:0.01:100); 
tfs_init        = eval_tfs(lureinit, fgrid);
tfs_final       = eval_tfs(lurefinal, fgrid);
tfs_final_trans = eval_tfs(lurefinal_trans, fgrid);
tfs_true        = eval_tfs(luretrue, fgrid);

% Define linestyle
ls_par = {'-','-','--'};
figure()
subplot(2,2,1)
    bodemag_custom({tfs_init.G_yu_frd,                   ...
                    tfs_final.G_yu_frd,                  ...
                    tfs_true.G_yu_frd},                  ...
                   {'Init','Final','True'},ls_par  ...
                  );
    title('$G_{yu}$')    
subplot(2,2,2)
    bodemag_custom({tfs_init.G_yw_frd,                   ...
                    tfs_final.G_yw_frd,                  ...
                    tfs_true.G_yw_frd},                  ...
                   {'Init','Final','True'},ls_par  ...
                  );
    title('$G_{yw}$')    
lgd = legend('$\theta_{init}$','$\hat{\theta}_N$','$\theta_0$',...
    'orientation','vertical','location','SW');
% lgd.Position = [0.32 0.93 0.4 0.05];
subplot(2,2,3)
    bodemag_custom({tfs_init.G_zu_frd,                   ...
                    tfs_final.G_zu_frd,                  ...
                    tfs_true.G_zu_frd},                  ...
                   {'Init','Final','True'},ls_par  ...
                  );
    title('$G_{zu}$')
subplot(2,2,4)
    bodemag_custom({tfs_init.G_zw_frd,                   ...
                    tfs_final.G_zw_frd,                  ...
                    tfs_true.G_zw_frd},                  ...
                   {'Init','Final','True'},ls_par  ...
                  );
    title('$G_{zw}$')
set(findobj(gcf,'Type','Line'),'linewidth',2)

if save_data
    save('SNR20_GBS')
end