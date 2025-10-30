close all
clear all

% Publication-ready schematic for posterior representation

% Define posterior (Gaussian)
mu = 1; 
sigma = 0.8;
x = linspace(-5,5,500);
y = normpdf(x,mu,sigma);

% Generate samples
numSamples = 1000;
samples = mu + sigma*randn(numSamples,1);

figure('Color','w','Position',[100 100 900 260]);

%% --- 1. Posterior distribution ---
subplot(1,3,1);
plot(x,y,'k','LineWidth',2);
xlim([-4 4]);
ylim([0 max(y)*1.1]);
xlabel('Stimulus ( s )','FontSize',16);
ylabel('p(s|o2)','FontSize',16);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.2,'FontSize',18);
% title("Posterior")
box off;

% Assuming sigma is defined
sigma = 0.8;

hold on

% % Draw horizontal line at y = height corresponding to sigma^2
% % For illustration, let's show y at p(mu+sigma)
% y_sigma = normpdf(mu + 1.5*sigma, mu, sigma);
% h = plot([-1.2 1.2], [y_sigma y_sigma], 'r--', 'LineWidth', 2);
% 
% % Optional: add text to label variance
% text(0, y_sigma + 0.05, '\sigma^2', 'FontSize',18, 'Color','r', 'HorizontalAlignment','center');
% 
% hold off

% Set sparse x-ticks
xticks([-4 -4+8/3 4-8/3 4]);                     % you can add intermediate if needed
xticklabels({'0','45','90', '135'});      % map -4 -> 0, 0 -> 67.5, 4 -> 135


% % Add axis arrows
% ax = gca;
% axesPos = ax.Position;
% annotation('arrow',[axesPos(1) axesPos(1)+axesPos(3)],...
%                    [axesPos(2) axesPos(2)]);  % x-axis
% annotation('arrow',[axesPos(1) axesPos(1)],...
%                    [axesPos(2) axesPos(2)+axesPos(4)]); % y-axis


%% --- 2. Samples from posterior ---
idx = randsample(length(samples), 30);
sel_samples = samples(idx);

subplot(1,3,2);
scatter(sel_samples, zeros(size(sel_samples)), ...
    100, 'o', ...
    'MarkerFaceColor',[0 0 0], ...
    'MarkerEdgeColor',[0 0 0], ...
    'MarkerFaceAlpha',0.1, ...
    'LineWidth',1.2);
xlim([-4 4]);
ylim([-0.1 0.1]);
xlabel('Stimulus (s)','FontSize',16);
ylabel('Samples','FontSize',16);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.2,'FontSize',18);
box off;

% Set sparse x-ticks
xticks([-4 -4+8/3 4-8/3 4]);                     % you can add intermediate if needed
xticklabels({'0','45','90', '135'});      % map -4 -> 0, 0 -> 67.5, 4 -> 135


% ax = gca;
% axesPos = ax.Position;
% annotation('arrow',[axesPos(1) axesPos(1)+axesPos(3)],...
%                    [axesPos(2) axesPos(2)]);
% annotation('arrow',[axesPos(1) axesPos(1)],...
%                    [axesPos(2) axesPos(2)+axesPos(4)]);

%% --- 3. Histogram of samples ---
subplot(1,3,3);
histogram(samples, 'Normalization','pdf','FaceColor',[0.2 0.2 0.2], ...
    'EdgeColor','none','BinWidth',0.4);
xlim([-4 4]);
xlabel('Stimulus (s)','FontSize',16);
ylabel('Frequency','FontSize',16);
% title("Histogram of Samples")
set(gca,'XTick',[],'YTick',[],'LineWidth',1.2,'FontSize',18);
box off;

% Set sparse x-ticks
xticks([-4 -4+8/3 4-8/3 4]);                     % you can add intermediate if needed
xticklabels({'0','45','90', '135'});      % map -4 -> 0, 0 -> 67.5, 4 -> 135


% ax = gca;
% axesPos = ax.Position;
% annotation('arrow',[axesPos(1) axesPos(1)+axesPos(3)],...
%                    [axesPos(2) axesPos(2)]);
% annotation('arrow',[axesPos(1) axesPos(1)],...
%                    [axesPos(2) axesPos(2)+axesPos(4)]);