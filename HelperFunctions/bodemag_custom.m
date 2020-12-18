function [h] = bodemag_custom(FRD_OL,legends,ls)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         Custom bode-plot for FRD-objects         %%%%%%%%%%%%%%
%%%%%%%%%%%         Input:  1) Cell array of FRD-Objects     %%%%%%%%%%%%%%
%%%%%%%%%%%                 2) Cell array of legend entries  %%%%%%%%%%%%%%
%%%%%%%%%%%                 3) linestyle                     %%%%%%%%%%%%%%
%%%%%%%%%%%         Output: No outputs, Figure was built     %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot data
c(1,:) = [0 0.4470 0.7410];
c(2,:) = [0 0 0];
c(3,:) = [0.2660    1    0.2880];
c(3,:) = [0.85000   0.32500 0.09800];

for i = 1:length(FRD_OL)
    complex_nrs = squeeze(FRD_OL{i}.ResponseData);
    freqs       = FRD_OL{i}.Frequency/(2*pi);
    h = semilogx(freqs,mag2db(abs(complex_nrs)),ls{i},'color',c(i,:));
    hold on
end

% Figure Layout
grid minor
ylabel('Magnitude [dB]')
xlabel('Frequency [Hz]')
% legend(legends)
axis tight
xlim([freqs(2) freqs(end)])
end