%% Plot polar Threshold
% Written by Tobias Jüterbock
close all
clear all

savePlot = true;
plotpath = 'Plots (full grid)';
plotfile = 'polarThreshold_full';


%% Load polar threshold

filepath = './Results/Threshold';
filename = 'polarThreshold_all_VPs_absoluteLatAngle_lat90added_a396ae7cf8afe4c00387625fa804ebe0';

load (fullfile(filepath,filename));

%% Plot

AKf(12,12)

x = [-lat_ref_abs(end:-1:1),0, lat_ref_abs];
y = pol_ref;
C = [-polarThreshold(:,end:-1:1,2),zeros(numel(pol_ref),1), polarThreshold(:,:,1)];

i = imagesc(x,y,C);
hold on
xline(0,'k','linewidth',6,'alpha',1)
set(gca,'YDir','normal')
        AKpSetColormapAndColorbar_gca('Reds', 'EastOutside', 20, [0 180], 'polar Threshold')
c = colorbar;
set(c,'Ticklabelinterpreter','latex')
ylabel(c,'polar Threshold [$^\circ$]','interpreter','latex')

text(-45,230,{'\textbf{Negative}','\textbf{opening}','\textbf{angles}'},...
     'color','w','FontSize',12, ...
     'HorizontalAlignment','center', ...
     'EdgeColor','w',...
     'interpreter','latex')
text( 45,230,{'\textbf{Positive}','\textbf{opening}','\textbf{angles}'},...
     'color','w','FontSize',12, ...
     'HorizontalAlignment','center', ...
     'EdgeColor','w',...
     'interpreter','latex')

xt = -80:20:80;
% xticks(xt)
% xticklabels(abs(xt))
% xticks(linspace(-90,90,39))
xticks([-80.53,-61.58,-42.36, -23.68, -4.737,4.747,23.68,42.36, 61.58,80.53])
xticklabels([80 , 60, 40,20,0,0,20,40,60,80])
set(gca,'Ticklabelinterpreter','latex')


xlabel('$|\phi_\mathrm{ref}|$ [$^\circ$]','interpreter','latex')
ylabel('$\theta_\mathrm{ref}$ [$^\circ$]','interpreter','latex')

set(gca, 'Position', [0.103,0.09,.736,0.89])

if savePlot
    savefig(fullfile(plotpath,plotfile))
    saveas(gcf, fullfile(plotpath,[plotfile,'.pdf']))
    print('-depsc', fullfile(plotpath,[plotfile,'.eps']))
end

