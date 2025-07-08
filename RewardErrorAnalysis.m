close all
clear all

% data = load('COR01_Ranjan random ori.mat');
% data = load('COR01_600_trials.mat');
data = load('COR01.mat');

stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

absErr = abs(stimOri - reportedOri);
absOriError = min(absErr, 180 - absErr);  % ensures error is in [0, 90]
rawError = stimOri - reportedOri;
rawOriError = mod(rawError + 90, 180) - 90;

data.dat.absOriError = absOriError;
data.dat.rawOriError = rawOriError;

%% Analytical reward function
maxTolerableError = 20; % In degrees
sigmaHC = sqrt(30);     % HC reward function std deviation
valLC   = 0.3;          % LC constant reward
x = linspace(-20, 20, 100);

mu = 0;
sigma = sigmaHC;
y1 = exp( - (abs(x)).^2 / (2 * sigma^2) );
y2 = ones(1, numel(x))*valLC;

x1 = sqrt( - log(0.3)* 2 * sigma^2 );
x2 = -x1;

figure
hold on
plot(x, y1, DisplayName="HC", LineWidth=2)
plot(x, y2, DisplayName="LC", LineWidth=2)
xlabel("Perceptual error")
ylabel("Reward")
legend

xline(x1, LineWidth=2, LineStyle="--", HandleVisibility="off")
xline(x2, LineWidth=2, LineStyle="--", HandleVisibility="off")

hold off

%% Histogram
theoreticalPropHC = zeros(1, numel(unique(data.dat.stimContrast)));
actualPropHC = zeros(1, numel(unique(data.dat.stimContrast)));

tableSummary = groupsummary(data.dat, {'stimContrast', 'reportedConf'}, ...
    {'mean', 'std', 'numel'}, 'rawOriError');
tableSummary = tableSummary(tableSummary.reportedConf == 1, : );
actualPropHC(:) = tableSummary.GroupCount / 120;

[G, contrastLevels] = findgroups(data.dat.stimContrast);
grpIdxes = unique(G);

% Sort conditions by uncertainty difference
[B_, sidx] = sort(contrastLevels);

% sSpread   = spreadLevels(sidx);
% sContrast = contrastLevels(sidx);

figure
for i=1:numel(sidx)
    gidx_ = sidx(i);

    subplot(2, 2, i)
    hold on

    grpOriErr = rawOriError(G == gidx_);
    histogram(grpOriErr, -90:3:90, "Normalization", "pdf")
    
    % fltOriErr = grpOriErr( zscore(grpOriErr) < 3 ); % Data within 3 SD
    fltOriErr = grpOriErr( abs(grpOriErr) < 50 );
    pd = fitdist(fltOriErr, 'Normal');
    x_values = -90:1:90;
    y = pdf(pd, x_values);
    % Compute from actual data as well - HC and LC
    
    plot(x_values, y, LineWidth=1.5)
    
    xline(x1, LineWidth=2, LineStyle="--", HandleVisibility="off")
    xline(x2, LineWidth=2, LineStyle="--", HandleVisibility="off")
    
    mu = pd.mu;
    sigma = pd.sigma;
    y = cdf('Normal', x1, mu, sigma) - cdf('Normal', x2, mu, sigma);

    theoreticalPropHC(i) = y;

    % Prop of error above the threshold (without considering confidence report)
    aboveThreshodl = sum(abs(grpOriErr) > x1) / numel(grpOriErr);  % Theoretical low confidence
    belowThreshodl = sum(abs(grpOriErr) <= x1) / numel(grpOriErr); % Theoretical high confidence

    xlabel("Orientation (deg)")
    ylabel("count")
    title(sprintf("Contrast: %.3f, \n Prop HC: (%.3f), Prop LC: (%.3f) " + ...
        "\n Below Thr: %.2f, Above Thr: %.2f " + ...
        "\n sigma = %.3f", ...
        contrastLevels(gidx_), y, 1-y, ...
        belowThreshodl, aboveThreshodl, sigma))
    % xlim([-50, 50])
    hold off
%     disp(pd)

    
end

%% Histogram by confidence

[G, contrastLevels, confReports] = findgroups(data.dat.stimContrast, data.dat.reportedConf);
grpIdxes = unique(G);

figure
for i=1:numel(unique(contrastLevels))
    contrLevel_ = contrastLevels(2*i);
    % confLevel_ = confReports(i);
    fltIdxLC = (data.dat.stimContrast == contrLevel_) & (data.dat.reportedConf == 0);
    fltIdxHC = (data.dat.stimContrast == contrLevel_) & (data.dat.reportedConf == 1);
    rawErrLC = rawOriError(fltIdxLC);
    rawErrHC = rawOriError(fltIdxHC);
    rewardLC = data.dat.reward(fltIdxLC);
    rewardHC = data.dat.reward(fltIdxHC);
    
    subplot(2, 4, i)
    hold on

    histogram(abs(rawErrLC), 0:3:90, FaceColor="red", EdgeColor="red", DisplayName="LC")
    histogram(abs(rawErrHC), 0:3:90, FaceColor="green", EdgeColor="green", DisplayName="HC", ...
        FaceAlpha=0.7, EdgeAlpha=0.5)

    xline(x1, LineWidth=2, LineStyle="--", HandleVisibility="off")
    % xline(x2, LineWidth=2, LineStyle="--", HandleVisibility="off")
    
    xlabel("Orientation (deg)")
    ylabel("count")
    title(sprintf("Contrast: %.3f, Conf: %.3f", contrastLevels(i), confReports(i)))
    legend
    % xlim([-50, 50])
    hold off
    
    % Avg earned reward in each bin
    binIdxesHC = discretize(abs(rawErrHC), 0:3:90);
    binIdxesLC = discretize(abs(rawErrLC), 0:3:90);
    earnedRewardHC = zeros(1, numel(0:3:90) - 1);
    earnedRewardLC = zeros(1, numel(0:3:90) - 1);
    
    for k=1:numel(earnedRewardHC)
        fltidx_ = binIdxesLC == k;
        avgearnedReward = sum( rewardLC(fltidx_) );
        earnedRewardLC(k) = avgearnedReward;
        
        fltidx_ = binIdxesHC == k;
        avgearnedReward = sum( rewardHC(fltidx_) );
        earnedRewardHC(k) = avgearnedReward;
    end
    
    subplot(2, 4, 4 + i)
    edges = 0:3:90;
    binCenters = ( edges(1:end-1) + edges(2:end) ) / 2;
    hold on
    plot(binCenters, earnedRewardLC, LineWidth=2, DisplayName="LC")
    plot(binCenters, earnedRewardHC, LineWidth=2, DisplayName="HC")
    xlabel("Perceptual error (abs)")
    ylabel("Total earned reward")
    legend
    hold off
end

%% Actual HC report vs theoreticala report
figure
scatter(theoreticalPropHC, actualPropHC)
xlabel("Theoretical prop of HC reports")
ylabel("Actual prop of HC reports")
xlim([0 1])
ylim([0 1])

% Add text labels near each point
for i = 1:length(theoreticalPropHC)
    text(theoreticalPropHC(i) + 0.01, actualPropHC(i), ...
         sprintf(' Contrast %.3f', tableSummary.stimContrast(i)));
end



