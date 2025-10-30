% Publication-ready MATLAB plot of p(z|o) vs orientations z (peak at 45 deg)
clear; close all; clc;

% --- Parameters ---
mu_deg   = 45;                  % peak orientation in degrees
kappa    = 18;                  % concentration (higher -> narrower peak)
z_deg    = linspace(0,90,1000); % orientations from 0–90 deg
deg2rad  = pi/180;
mu       = mu_deg*deg2rad;
z_rad    = z_deg*deg2rad;

% -- von Mises pdf implementation (circular, normalized)
I0 = besseli(0,kappa);
p = exp(kappa*cos(z_rad - mu)) ./ (2*pi*I0);

% Normalize over 0–90 so integral = 1 in this range
p = p ./ trapz(z_deg,p);

% --- Figure styling (publication ready) ---
f = figure('Units','inches','Position',[1 1 6 6], 'Color','w'); 
ax = axes('Parent',f); hold(ax,'on');

% Plot line (thicker)
plot(z_deg, p, 'LineWidth',3, 'Color',[0 0.2 0.6]);

% Axes / labels (bold + large)
xlabel('\textbf{Orientation $\mathbf{z}$}', ...
    'Interpreter','latex', 'FontSize',40);
ylabel('\textbf{Posterior} $\mathbf{p(z \mid o)}$', ...
    'Interpreter','latex', 'FontSize',40);

% Ticks, limits
xlim([0 90]);
ylim([0 max(p)*1.15]);
xticks(0:45:90);

% Remove y-ticks but keep axis + label
yticks([]);

% Make axes and ticks thicker
set(ax, 'LineWidth',2, ...           % axis lines
        'TickLength',[0.03 0.03], ...
        'Box','off', ...
        'TickDir','out', ...
        'FontName','Times New Roman', ...
        'FontSize',30);  % tick label size

% Adjust axis position (smaller plotting box -> more space for labels)
ax.Position = [0.25 0.25 0.65 0.65];  % shrink and center

% Square aspect ratio
axis square;

% Save files
print(f, 'pz_given_o_plot.png', '-dpng', '-r600');  
print(f, 'pz_given_o_plot.pdf', '-dpdf');
