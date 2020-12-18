function [IData] = generateInputData(input)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Periodic excitation design   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate some input design settings
    Fs           = 1/input.Ts;                                             % [Hz]
    fN           = Fs/2;                                                   % [Hz]
    Band         = [input.fMin, input.fMax]/fN;                            % Fraction of the nyquist frequency
    NumSinusoids = round((input.fMax - input.fMin)/input.df + 1);          % [-]
    dFgrid       = 1/(input.Ts*input.Period);                              % [Hz]
    GridSkip     = round(input.df/dFgrid);                                 % [samples]

    % validate which set of frequencies is to be selected
    test_tol      = 1e-6;
    fGrid_Init    = dFgrid*(0:GridSkip:fix(input.Period/2)-1);
    fGrid         = fGrid_Init(input.fMin - fGrid_Init <=  test_tol & ...
                               input.fMax - fGrid_Init >= -test_tol);

    % generate the multisine type input sequence
    [w, ~]    = idinput_custom(                      ...
                    [input.Period 1 input.NumPeriod],   ...
                    'sine',                             ...
                    Band,                               ...
                    input.Range,                        ...
                    [NumSinusoids, input.NumTrials, GridSkip]);   
       
    % Define output struct for this function 
    IData.t_tot   = input.Ts*(0:input.NumPeriod*input.Period-1).';
    IData.w_tot   = w;
    IData.t       = IData.t_tot(1:input.Period);
    IData.w       = IData.w_tot(1:input.Period);
    IData.w_freqs = fGrid;                                                 %freq/(pi*dFgrid); % Original freq vector assumes unit sampling --> to Hz! 
    
    % analyze the excitation signal design
    fig_inputanalysis = false;
    if fig_inputanalysis

        % 1) Analyze the frequency selection of idinput
        figure
        plot(fGrid,'-*')
        hold on
        plot([1 length(fGrid)], [IData.w_freqs; IData.w_freqs].', '--*')
        plot([1 length(fGrid)], [input.fMin input.fMin]     , '-k')
        plot([1 length(fGrid)], [input.fMax input.fMax]     , '-k')
        ylim([input.fMin, input.fMax])
        
        % analyze the full excitation signal
        wfft            = fft(IData.w_tot);
        L               = length(IData.w_tot);
        magTS           = abs(wfft/L);
        magSS           = magTS(1:L/2+1);
        magSS(2:end-1)  = 2*magSS(2:end-1);
        f_fft           = (0:(L/2))/(L*input.Ts);

        % analyze the first period of excitation signal
        wfft_per            = fft(IData.w);
        L_per               = input.Period;
        magTS_per           = abs(wfft_per/L_per);
        magSS_per           = magTS_per(1:L_per/2+1);
        magSS_per(2:end-1)  = 2*magSS_per(2:end-1);
        f_fft_per           = (0:(L_per/2))/(L_per*input.Ts);

        % plot analysis for full and periodic signal
        figure
            subplot(2,1,1)
                plot(IData.t,IData.w)                     
            subplot(2,1,2)
                stem(f_fft_per, magSS_per) 
                title('Single-Sided Amplitude Spectrum of w(t) period')
                xlabel('Frequency [Hz]')
                ylabel('Amplitude')    
                xticks(IData.w_freqs)
                xlim([IData.w_freqs(1) IData.w_freqs(end)])
    end
end
