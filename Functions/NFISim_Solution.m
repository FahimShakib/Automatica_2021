function [IOData] = NFISim_Solution(luresys, IData, nfi_pars)

    % NFI via Simulink 
    P     = nfi_pars.P;                                                    % Simulate P periods
    n     = nfi_pars.n;                                                    % Number of samples per period [-]
    T     = nfi_pars.T;                                                    % Period time [s]

    % Simulink input signal definition
    win	= [0:T/n:T*P-T/n;kron(ones(1,P),IData.w.')]';                      
    
    % Solve the Simulink model by numeric forward integration
    options = simset('SrcWorkspace','current');
    sim('Lure_model_simulation',[],options)
    
    % Select the last period as steady-state solution
    idx = find(NFISim_t>(P-1)*T);
    IOData.w = NFISim_w(idx);
    IOData.t = NFISim_t(idx)-(P-1)*T;                                      % Start at t = 0
    IOData.y = NFISim_y(idx);
    IOData.x = NFISim_x(idx,:);
    IOData.u = NFISim_u(idx);
    IOData.z = NFISim_z(idx);
end