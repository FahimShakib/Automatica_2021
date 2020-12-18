function [JJ, dJdtheta] = costFunction(theta, mtf_pars, nfi_pars, SSData_meas, method)

    global t_logs
    
    % Check convergence for the selected parameter vector
    lure     = thetaToLure(theta);

    A_tf   	= lure.A - lure.gamma/2*lure.B*lure.C;
    Gyu  	= ss(A_tf, lure.B, lure.C, 0);
    infnorm	= norm(Gyu, Inf);
    
    conv_1     = all([lure.W1.', lure.W2.'] > 0);                          % < or <= for the W's?? !STILL CHECK THIS!
    conv_2     = all(real(eig(lure.A)) < 0 );
    conv_3     = infnorm < 1/(lure.gamma/2) - 1e-4;                        % Small tolerance results in less MTF-convergence issues
    convergent = all([conv_1,conv_2,conv_3]);

    if convergent

        % Evaluate the model's steady-state solution
        switch method
            case 'mtf'
                t_var      = tic;
                SSData_sim = MTF_Solution(lure, SSData_meas, mtf_pars);
                t_mtf      = toc(t_var);
            case 'nfi'
                t_var      = tic;
                SSData_sim = NFISim_Solution(lure, SSData_meas, nfi_pars);
                t_nfi      = toc(t_var);
        end

        % Run alternating method in parallel if desired
        if nfi_pars.comptime
            switch method
                case 'mtf'
                    t_var = tic;
                    [~]   = NFISim_Solution(lure, SSData_meas, nfi_pars);
                    t_nfi = toc(t_var);
                case 'nfi'
                    t_var = tic;
                    [~]   = MTF_Solution(lure, SSData_meas, mtf_pars);
                    t_mtf = toc(t_var);
            end
            
            % Update timing logs
            t_logs.lure1.mtf = [t_logs.lure1.mtf, t_mtf];
            t_logs.lure1.nfi = [t_logs.lure1.nfi, t_nfi];
        end

        % Evaluate cost-function
        epsilon    = SSData_sim.z - SSData_meas.z;
        JJ         = 1/(mtf_pars.n)*(epsilon.'*epsilon);

        if nargout > 1   
            % Partial derivatives of LTI block w.r.t. theta
            Ajac = [            eye(16),  zeros(16,28)];
            Bjac = [zeros(4,16), eye(4),  zeros( 4,24)];
            Ljac = [zeros(4,20), eye(4),  zeros( 4,20)];
            Cjac = [zeros(4,24), eye(4),  zeros( 4,16)];
            Djac = [zeros(1,28), eye(1),  zeros( 1,15)];
            Fjac = [zeros(4,29), eye(4),  zeros( 4,11)];
            Gjac = [zeros(1,33), eye(1),  zeros( 1,10)];
            Hjac = [zeros(1,34), eye(1),  zeros( 1, 9)];

            A_theta = reshape(Ajac,4,4,44);
            B_theta = reshape(Bjac,4,1,44);
            L_theta = reshape(Ljac,4,1,44);
            C_theta = reshape(Cjac,1,4,44);
            D_theta = reshape(Djac,1,1,44);
            F_theta = reshape(Fjac,1,4,44);
            G_theta = reshape(Gjac,1,1,44);
            H_theta = reshape(Hjac,1,1,44);

            % partial derivative of the nonlinearity w.r.t. 
            % y and theta, evaluated at steady-state y.
            dphidTheta     = lure.dphidtheta(SSData_sim.y);   
            IData_sens.Psi = lure.dphidy(SSData_sim.y);                        

            % Jacobian of the cost function evaluated at x
            J = zeros(mtf_pars.n, length(theta));
            for i=1:length(theta)
                % SS-inputs for gradient computation theta_i
                IData_sens.Wi = A_theta(:,:,i)*SSData_sim.x.' + B_theta(:,:,i)*SSData_sim.u.' + L_theta(:,:,i)*SSData_sim.w.' - lure.B*dphidTheta(:,i).';
                IData_sens.Yi = C_theta(:,:,i)*SSData_sim.x.'                                 + D_theta(:,:,i)*SSData_sim.w.';
                IData_sens.Zi = F_theta(:,:,i)*SSData_sim.x.' + G_theta(:,:,i)*SSData_sim.u.' + H_theta(:,:,i)*SSData_sim.w.' - lure.G*dphidTheta(:,i).';

                IData_sens.Wi = IData_sens.Wi';
                IData_sens.Yi = IData_sens.Yi';
                IData_sens.Zi = IData_sens.Zi';

                % Evaluate the model's steady-state solution
                switch method
                    case 'mtf'
                        t_var      = tic;
                        SSGradData = MTF_Sensitivity(lure,  IData_sens, mtf_pars);
                        t_mtf(i,1) = toc(t_var);
                    case 'nfi'
                        t_var      = tic;
                        SSGradData = NFISim_Solution2(lure, IData_sens, nfi_pars);
                        t_nfi(i,1) = toc(t_var);
                end

                % Run alternating method in parallel if desired
                if nfi_pars.comptime
                    switch method
                        case 'mtf'
                            t_var       = tic;
                            SSGradData2 = NFISim_Solution2(lure, IData_sens, nfi_pars);
                            t_nfi(i,1)  = toc(t_var);
                        case 'nfi'
                            t_var = tic;
                            SSGradData2  = MTF_Sensitivity(lure, IData_sens, mtf_pars);
                            t_mtf(i,1)   = toc(t_var);
                    end

                    % Compare NFI- to MTF-solution for gradient checks
                    compfig = false;
                    if compfig
                            figure
                            plot(SSGradData.t, SSGradData.z)
                            hold on
                            plot(SSGradData2.t, SSGradData2.z,'--');
                            xlabel('Time [s]')
                            ylabel(['$\nabla z_{\theta_',num2str(i),'}$'],'Interpreter','Latex')
                            legend('NFI-Simulink','MTF')  
                    end
                end

                J(:,i) = SSGradData.z;
            end

            % Update timing logs
            if nfi_pars.comptime
                t_logs.lure2.mtf  = [t_logs.lure2.mtf, t_mtf];
                t_logs.lure2.nfi  = [t_logs.lure2.nfi, t_nfi];
            end

            % Evaluate cost function gradient
            dJdtheta = 2/(mtf_pars.n)*epsilon.'*J;
        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Large costs for non-convergent models   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        JJ       = 10000;
        dJdtheta = zeros(1,44);  
    end
end