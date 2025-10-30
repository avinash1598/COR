close all; clear; clc;

%% Parameters
mu = 0; 
sigma = 1; 
numBins = 50;      % number of time bins
binWidth = 0.1;    % seconds per bin
maxRate = 20;      % spikes/s

%% Generate posterior samples (firing rate per bin)
samples = mu + sigma*randn(numBins,1);
rates = maxRate * exp(-0.5*((samples-mu)/sigma).^2); % Gaussian tuning

%% Generate spike times (inhomogeneous Poisson)
spikeTimes = [];
for b = 1:numBins
    nSpikes = poissrnd(rates(b)*binWidth);
    spikeTimes = [spikeTimes, b*binWidth + rand(1,nSpikes)*binWidth];
end

%% Figure
figure('Color','w','Position',[100 100 1000 300]);

%% 1) Spike raster
subplot(1,3,1); hold on;
plot(spikeTimes, ones(size(spikeTimes)), 'k.', 'MarkerSize',12);
xlabel('Time (s)','FontSize',16);
ylabel('Spikes','FontSize',16);
xlim([0 numBins*binWidth]);
ylim([0 1.5]);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.2,'FontSize',14);
box off;
axis square;

%% Add axis arrows for raster
% quiver(0,0, numBins*binWidth,0,0,'k','LineWidth',1.5,'MaxHeadSize',0.2);
% quiver(0,0,0,1.2,0,'k','LineWidth',1.5,'MaxHeadSize',0.2);

title('Single Neuron spikes','FontSize',16);

%% 2) Firing rate / samples
subplot(1,3,2); hold on;
bar((1:numBins)*binWidth, samples, 'FaceColor',[0 0 0], ...
    'FaceAlpha',0.2, 'EdgeColor','k', 'LineWidth',1.2);
xlabel('Time (s)','FontSize',16);
ylabel('Firing rate','FontSize',16);
xlim([0 numBins*binWidth]);
ylim([min(samples)-0.5 max(samples)+0.5]);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.2,'FontSize',14);
box off;
axis square;

% Axis arrows
% quiver(0,0, numBins*binWidth,0,0,'k','LineWidth',1.5,'MaxHeadSize',0.2);
% quiver(0,min(samples)-0.5,0,max(samples)+1,0,'k','LineWidth',1.5,'MaxHeadSize',0.2);

title('Firing rate = sample','FontSize',16);

%% 3) Histogram → posterior
subplot(1,3,3); hold on;
histogram(samples,'Normalization','pdf','FaceColor',[0.2 0.2 0.2], ...
    'EdgeColor','none','FaceAlpha',0.5,'BinWidth',0.4);
xplot = linspace(min(samples)-1,max(samples)+1,200);
% plot(xplot,normpdf(xplot,mu,sigma),'r-','LineWidth',2);
xlabel('Firing rate','FontSize',16);
ylabel('Frequency','FontSize',16);
set(gca,'XTick',[],'YTick',[],'LineWidth',1.2,'FontSize',14);
box off;
axis square;
% quiver(min(samples)-1,0,max(samples)-min(samples)+2,0,0,'k','LineWidth',1.5,'MaxHeadSize',0.2);
% quiver(0,0,0,max(normpdf(xplot,mu,sigma))+0.1,0,'k','LineWidth',1.5,'MaxHeadSize',0.2);

% title('Firing rate distribution → Posterior','FontSize',16);
