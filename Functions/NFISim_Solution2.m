function [IOData] = NFISim_Solution2(luresys, IData, nfi_pars)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Approximate gradient ss-solution via a Simulink NFI implementation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define simulation length
    P     = nfi_pars.P;                                                    % Number of simulated periods  [-]
    n     = nfi_pars.n;                                                    % Number of samples per period [-]
    T     = nfi_pars.T;                                                    % Period time [s]

    % Define simulation inputs
    tvec  = 0:T/n:T*P-T/n;
    Win   = [tvec;kron(ones(1,P),IData.Wi.')]';                            % Input signal for simulink
    Yin   = [tvec;kron(ones(1,P),IData.Yi.')]';                            % Input signal for simulink
    Zin   = [tvec;kron(ones(1,P),IData.Zi.')]';                            % Input signal for simulink
    Psiin = [tvec;kron(ones(1,P),IData.Psi.')]';
    
    % Simulation
    options = simset('SrcWorkspace','current');
    sim('Lure_model_gradient_simulation',[],options)
    
    % Select the last period as steady-state solution
    idx      = find(NFISim_t >= (P-1)*T);
    IOData.t = NFISim_t(idx)-(P-1)*T;
    IOData.y = NFISim_y(idx);
    IOData.x = NFISim_x(idx);
    IOData.u = NFISim_u(idx);
    IOData.z = NFISim_z(idx);
end