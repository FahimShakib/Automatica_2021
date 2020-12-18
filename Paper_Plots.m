%% Table 1 in paper (SNR20,40,60 and SNRminf_GBS)
display('SNR20')
load Paper_Data/SNR20
display(['sigma_e^2: ' num2str(vare^2)])
display(['Cost init: ' num2str(J(thetainit))])
display(['Cost final: ' num2str(J(thetafinal))])

display('SNR40')
load Paper_Data/SNR40
display(['sigma_e^2: ' num2str(vare^2)])
display(['Cost init: ' num2str(J(thetainit))])
display(['Cost final: ' num2str(J(thetafinal))])

display('SNR60')
load Paper_Data/SNR60
display(['sigma_e^2: ' num2str(vare^2)])
display(['Cost init: ' num2str(J(thetainit))])
display(['Cost final: ' num2str(J(thetafinal))])

display('SNRminf')
load Paper_Data/SNRminf_GBS
display(['sigma_e^2: ' num2str(vare^2)])
display(['Cost init: ' num2str(J(thetainit))])
display(['Cost final: ' num2str(J(thetafinal))])

return

%% Table 2 in paper (SNRminf_CRS and SNRminf_GBS)
load Paper_Data/SNRminf_CRS
display(['Nr of simulations: ' num2str(length(t_logs.lure1.nfi))])
display(['MTF time CRS: ' num2str(sum(t_logs.lure1.mtf))])
display(['NFI time CRS: ' num2str(sum(t_logs.lure1.nfi))])

load Paper_Data/SNRminf_GBS
display(['Nr of simulations: ' ...
    num2str(length(t_logs.lure1.nfi)+ numel(t_logs.lure2.nfi))])
display(['MTF time GBS: ' ...
    num2str(sum(t_logs.lure1.mtf)+sum(sum(t_logs.lure2.mtf)))])
display(['NFI time GBS: ' ...
    num2str(sum(t_logs.lure1.nfi)+sum(sum(t_logs.lure2.nfi)))])

return

%% Figure 3 in paper: Bode plots + NLity (of the 60dB noise case)
addpath('Functions');addpath('HelperFunctions')
load Paper_Data/SNR60

saveplot = 0;

% Define linestyle
ls_par = {'-','-','--'};
h = figure();
h.Position  = [680   200   1020   350];
subplot(2,3,1)
    bodemag_custom({tfs_true.G_yu_frd,                   ...
                    tfs_init.G_yu_frd,                   ...
                    tfs_final.G_yu_frd},                 ...
                   {''},ls_par  ...
                  );
    title('$G_{yu}$')    
xlabel('')
subplot(2,3,2)
    bodemag_custom({tfs_true.G_yw_frd,                    ...
                    tfs_init.G_yw_frd,                    ...
                    tfs_final.G_yw_frd},                  ...
                   {},ls_par  ...
                  );    
    title('$G_{yw}$')    
xlabel('')
ylabel('')
ylim([-100 0])
subplot(2,3,4)
    bodemag_custom({tfs_true.G_zu_frd,                   ...
                    tfs_init.G_zu_frd,                  ...
                    tfs_final.G_zu_frd},                  ...
                   {},ls_par  ...
                  );    
    title('$G_{zu}$')
ylim([-140 -20])
subplot(2,3,5)
    bodemag_custom({tfs_true.G_zw_frd,                   ...
                    tfs_init.G_zw_frd,                  ...
                    tfs_final.G_zw_frd},                  ...
                   {},ls_par  ...
                  );    
    title('$G_{zw}$')
ylabel('')
ylim([-140 -20])

subplot(2,3,[3 6])
plot(yrange,luretrue.nonlin(yrange),'color',[0 0.4470 0.7410])
hold on
plot(yrange,lureinit.nonlin(yrange),'k')
plot(yrange,lurefinal.nonlin(yrange),'--','color',[0.85000,0.32500,0.09800])

legend('$\theta_{0}$','$\theta_{init}$','$\hat{\theta}_N$','location','NW')
ylabel('$\varphi$','Interpreter','Latex')
xlabel('$y$','Interpreter','Latex')
xlim([-15 10])
% grid minor
set(findobj(gcf,'Type','Line'),'linewidth',2)


if saveplot
    cleanfigure('minimumPointsDistance', 0.1)
    matlab2tikz(['../../Automatica - Second Submission/Figures/EX_Sim_Bode_NL.tex'],'width','\fwidth','height','\fheight');
end

return

%% Time response plot + validation (of the 20dB noise case)
clear all; clc
addpath('Functions');addpath('HelperFunctions')

saveplot    = 0;

load Paper_Data/SNR20

figure
subplot(211)
plot(SSLTrue.t,SSLTrue.z,'color',[0 0.4470 0.7410])
hold all
plot(SSLInit.t,SSLTrue.z-SSLInit.z,'k')
hold on
plot(SSLFinal.t,SSLTrue.z-SSLFinal.z,'color',[0.85000,0.32500,0.09800])   
set(findobj(gcf,'Type','Line'),'linewidth',2)
legend('$\bar{\zz}(t)$','$\bar \epsilon(t,\theta_{init})$',...
    '$\bar \epsilon(t,\hat{\theta}_N)$','orientation','horizontal')
ylabel('Position [m]')
ylim([-1 1.3])
title('Identification data')

% Perform validation test
input_val       = input;
input_val.Range = input_val.Range*10;
IData_val       = generateInputData(input_val);
SSData_val      = generateOutputData(luretrue, IData_val, mtf_pars, nfi_pars, 'mtf');
SSLInit_val     = generateOutputData(lureinit, IData_val, mtf_pars, nfi_pars, 'mtf');
SSLfinal_val    = generateOutputData(lurefinal, IData_val, mtf_pars, nfi_pars, 'mtf');

subplot(212)
plot(SSData_val.t,SSData_val.z,'color',[0 0.4470 0.7410])
hold all
plot(SSLInit_val.t,SSData_val.z-SSLInit_val.z,'k')
hold on
plot(SSLfinal_val.t,SSData_val.z-SSLfinal_val.z,'color',[0.85000,0.32500,0.09800])   
set(findobj(gcf,'Type','Line'),'linewidth',2)
xlabel('Time [sec]')
ylabel('Position [m]')
title('Validation data')

if saveplot
    cleanfigure('minimumPointsDistance', 0.1)
    matlab2tikz(['../../Automatica - Second Submission/Figures/Ex_Sim_TimeResp.tex'],'width','\fwidth','height','\fheight');
end

return