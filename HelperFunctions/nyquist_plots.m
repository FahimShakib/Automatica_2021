function [ax] = nyquist_plots(FRD_objs,Leg)

% Settings
small_axis = 'no';

% Set figure axis
k   = 1;
fig = gcf;
ax1 = gca;

%% Plot in big axis ax1

% Nyquist diagram default plots
fig.CurrentAxes = ax1;
colors = get(gca,'colororder');
hold(ax1,'on');
[~] = circles(0,0,1,'facecolor','none');
hold on
[~] = circles(-1,0,0.5,'FaceColor',[0.9 0.9 0.9]);
hold on
scatter(ax1,-1,0,40,'k.','HandleVisibility','off')

for i = 1:length(FRD_objs)
    complex_nrs = squeeze(FRD_objs{i}.ResponseData);
    hold(ax1,'on')
    h(k) = plot(real(complex_nrs),imag(complex_nrs),'Color',colors(i,:));
    legends{k} = [Leg{i},'~$:~f>0$'];
    k = k+1;
    h(k) = plot(real(complex_nrs),-imag(complex_nrs),'--','Color',colors(i,:));
    legends{k} = [Leg{i},'~$:~f<0$'];
    k = k+1;
end

%% Plot in small axis x2
if strcmp(small_axis,'yes') == 1
    ax2 = axes('Position',[.5 .2 .27 .27]);
    fig.CurrentAxes = ax2;
    box on;
    for i = 1:length(FRD_objs)
        complex_nrs = squeeze(FRD_objs{i}.ResponseData);
        hold(ax2,'on')
        plot(real(complex_nrs),imag(complex_nrs),'Color',colors(i,:));
        plot(real(complex_nrs),-imag(complex_nrs),'--','Color',colors(i,:));
    end
end

%% Figure layout
fig.CurrentAxes     = ax1;
ax1.XAxisLocation   = 'origin';
ax1.YAxisLocation   = 'origin';
grid(ax1,'minor')
ylim(ax1,[-3 3])
xlim(ax1,[-3 3])
xlabel(ax1,'Real axis [-]')
ylabel(ax1,'Imaginary axis [-]')
axis(ax1,'square')
legend(h,legends,'Interpreter','Latex')
ax{1} = ax1;    

if strcmp(small_axis,'yes') == 1
    fig.CurrentAxes = ax2;
    xlabel(ax2,'Re [-]')
    ylabel(ax2,'Im [-]')
    axis(ax2,'square')
    grid(ax2,'minor')
    ax2.XAxisLocation = 'origin';
    ax2.YAxisLocation = 'origin';
    ax{2} = ax2;
end


end

