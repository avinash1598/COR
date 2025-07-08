% close all
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

% datOK = data.dat(data.dat.trialStatus == 1, :);

%% Sample spread and mean distribution
spread = unique(data.dat.stimSpread);
idx = (data.dat.stimSpread == spread(1));
actualOris = data.dat.stimOri(idx);
actualSpreads = data.dat.stimSpread(idx);
sampleMeans = data.dat.stimSampleMeanOri(idx);
sampleSpreads = data.dat.stimSampleSpread(idx);

figure

% Mean ori distribution
rawError1 = actualOris - sampleMeans;
rawOriError1 = mod(rawError1 + 90, 180) - 90; % Error between -90 and 90
subplot(1, 2, 1)
histogram(rawOriError1)
std(rawOriError1(~isnan(rawOriError1)))

% Sample Spread distribution
spreadErrs = actualSpreads - sampleSpreads;
subplot(1, 2, 2)
histogram(spreadErrs)
std(~isnan(spreadErrs))

%% Figure 1 - Raw errors by group - Estimation error (Human subject)
% Get elements in the group
[G, contrastLevels] = findgroups(data.dat.stimContrast);
grpIdxes = unique(G);

% Sort conditions by uncertainty difference
[B_, sidx] = sort(contrastLevels);

sContrast = contrastLevels(sidx);

% Generate x-axis labels as combined strings like "0.2|5"
xLabels = strcat(string(sContrast));

% Define visual properties
colors = [0.85 0.33 0.10; 0 0.45 0.74; 0.85 0.33 0.10; 0 0.45 0.74]; % colorblind-friendly red/blue
alphas = [1, 1, 0.3, 0.3];

% Create figure
figure('Color','w');
pbaspect([3 5 1]) 
hold on

% Plot each group
for i = 1:numel(grpIdxes)
    grpOriErr = rawOriError(G == sidx(i));
    meanErr = mean(grpOriErr(~isnan(grpOriErr)));
    xPts = i + 0.4*(rand(1, numel(grpOriErr)) - 0.5);

    scatter(xPts, grpOriErr, 60, ...
        'MarkerEdgeColor', 'w', ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerFaceAlpha', alphas(i), ...
        'MarkerEdgeAlpha', alphas(i), ...
        'LineWidth', 0.5);

    plot([i - 0.35, i + 0.35], [meanErr, meanErr], LineStyle="-", LineWidth=2, Color='black')
end

% Labels and styling
ylabel("Estimation error (°)", 'FontSize', 14)
xlabel("Uncertainty level", 'FontSize', 14)
xticks(1:length(xLabels));
xticklabels(xLabels);
yticks([-90, 0, 90])
ylim([-90, 90])
xlim([0.3, 5])

% Improve figure aesthetics
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
hold off

exportgraphics(gcf, 'Human1.eps', 'ContentType', 'vector')
% print(gcf, 'Human1', '-depsc2', '-r300')  % -r300 = 300 DPI

%% Figure 2 - Error by uncertainty - take MSE error instead (maybe)
% === Figure 2 ===
% Clean data
errHC = rawOriError(data.dat.reportedConf == 1);
errLC = rawOriError(data.dat.reportedConf == 0);
errHC = errHC(~isnan(errHC));
errLC = errLC(~isnan(errLC));

abErrHC = abs(errHC);
abErrLC = abs(errLC);

% Prepare values
means = [mean(abErrLC), mean(abErrHC)];
sems  = [std(abErrLC)/sqrt(numel(abErrLC)), std(abErrHC)/sqrt(numel(abErrHC))]; % No bias in this plot needs to be shown
colors = [0.85 0.33 0.10; 0.47 0.67 0.19]; 

% Create figure
figure('Color','w', 'Units','inches', 'Position',[1, 1, 4.5, 5]);
pbaspect([3 5 1]) 
hold on

labels = ["Uncertain", "Certain"];
% Bar plot with custom colors
for i = 1:2
    b = bar(i+0.5, means(i), 0.6, DisplayName=labels(i)); % 0.6 = bar width
    b.FaceColor = colors(i,:);
    b.EdgeColor = 'none';
end

% Add error bars
errorbar([1.5 2.5], means, sems, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10, HandleVisibility='off');

% Styling
ylabel("|Estimation Error| (°)", 'FontSize', 14)
% yticks([0, round( max(means + sems) * 1.2 ) ])
yticks([0, 30])
xticks([])
% xticklabels({'Low Conf.', 'High Conf.'})
lgd = legend('\color[rgb]{0.85, 0.33, 0.10} Uncertain','\color[rgb]{0.47, 0.67, 0.19} Certain');
lgd.Box = "off";
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
% ylim([0, max(means + sems) * 1.2])  % Leave headroom
ylim([0, 30]) 

% Grid and finish
hold off

exportgraphics(gcf, 'Human2.eps', 'ContentType', 'vector')
% print(gcf, 'Human2', '-depsc', '-r300')  % -r300 = 300 DPI

%% Orientation bias
stimOri = data.dat.stimOri;
bins = linspace(0, 180, 18);
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2;
binIdx = discretize(stimOri, bins);

mseErrs = zeros(1, numel(unique(binIdx)));
avgErrs = zeros(1, numel(unique(binIdx)));
stdErrs = zeros(1, numel(unique(binIdx)));
stdDevs = zeros(1, numel(unique(binIdx)));

for i=1:numel(unique(binIdx))
    errors_ = rawOriError(binIdx == i);
    errors_ = errors_(~isnan(errors_));
    mseErrs(i) = sum(errors_.^2)/numel(errors_);
    avgErrs(i) = mean(errors_);
    stdErrs(i) = std(errors_) / sqrt(numel(errors_));
    stdDevs(i) = std(errors_);
end

figure 
subplot(2, 1, 1)
errorbar(binCenters, avgErrs, stdErrs, ...
    'k', 'LineStyle', '-', 'LineWidth', 1.5, 'CapSize', 10)
xlabel("Orientation (deg)")
ylabel("Error")
title("Orientation bias")

subplot(2, 1, 2)
x = 1:numel(unique(binIdx));
width = 0.25;
hold on
bar(x - width/2, stdDevs, width, DisplayName="Std Dev")
bar(x + width/2, mseErrs, width, DisplayName="MSE Err")
xlabel("Orientation (deg)")
xticks(x)
xticklabels( round(binCenters))
legend
hold off

%% Histogram

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
    grpOriErr = rawOriError(G == gidx_);
    histogram(grpOriErr, -90:3:90)
    
    xlabel("Orientation (deg)")
    ylabel("count")
    title(sprintf("Contrast: %.2f", contrastLevels(gidx_)))
    xlim([-50, 50])
end

title("Distribution of errors")


%% MSE error by contrast levels and confidence
datOK = data.dat(data.dat.trialStatus == 1, :);

summaryTable = groupsummary(datOK, {'stimContrast', 'reportedConf'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'mean', 'std', 'numel'}, 'rawOriError');

disp(summaryTable)

summaryTable = groupsummary(datOK, {'stimContrast'}, ...
                            {'mean', 'std', 'numel'}, 'absOriError');


sContrasts = summaryTable.stimContrast;
[B_, sidx] = sort(sContrasts);

sContrast = sContrasts(sidx);
meanError = summaryTable.mean_absOriError(sidx);
stdError = summaryTable.std_absOriError(sidx);

% Generate x-axis labels as combined strings like "0.2|5"
xLabels = strcat(string(sContrast));

% Plot with error bars
figure;
subplot(2, 2, 1)
errorbar(meanError, stdError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('Absolute Error by Contrast and Spread');
grid on;


% By confidence
summaryTable = groupsummary(datOK, {'stimContrast'}, ...
                            {'mean', 'std', 'numel'}, 'absOriError');

sContrasts = summaryTable.stimContrast;
[B_, sidx] = sort(sContrasts);
sContrast = sContrasts(sidx);

summaryTable = groupsummary(datOK, {'stimContrast', 'reportedConf'}, ...
                            {'mean', 'std', 'numel'}, 'absOriError');

mErrHC = zeros(1, numel(sContrast));
stdErrHC = zeros(1, numel(sContrast));
mErrLC = zeros(1, numel(sContrast));
stdErrLC = zeros(1, numel(sContrast));

for i = 1:numel(sContrast)
    rowHC = summaryTable(summaryTable.stimContrast == sContrast(i) & ...
        summaryTable.reportedConf == 1, :);
    mErrHC(i) = rowHC.mean_absOriError;
    stdErrHC(i) = rowHC.std_absOriError;

    rowLC = summaryTable(summaryTable.stimContrast == sContrast(i) & ...
        summaryTable.reportedConf == 0, :);
    if ~isempty(rowLC)
        mErrLC(i) = rowLC.mean_absOriError;
        stdErrLC(i) = rowLC.std_absOriError;
    else
        mErrLC(i) = NaN;
        stdErrLC(i) = NaN;
    end

end

subplot(2, 2, 2)

% Generate x-axis labels as combined strings like "0.2|5"
xLabels = strcat(string(sContrast));

% Plot with error bars
subplot(2, 2, 2)
hold on
errorbar(mErrHC, stdErrHC, 'o-', LineWidth=2, DisplayName='HC');
errorbar(mErrLC, stdErrLC, 'o-', LineWidth=2, DisplayName='LC');

xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('Absolute Error by Contrast and Spread');
grid on;
legend
hold off


% MSE error
summaryTable = groupsummary(datOK, {'stimContrast'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'numel'}, 'rawOriError');

meanError = summaryTable.fun1_rawOriError(sidx);

subplot(2, 2, 3)
plot(1:numel(meanError), meanError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimContrast" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('MSE error');
grid on;


% MSE err split by confidence
summaryTable = groupsummary(datOK, {'stimContrast', 'reportedConf'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'numel'}, 'rawOriError');

mseHC = zeros(1, numel(sContrast));
mseLC = zeros(1, numel(sContrast));

for i = 1:numel(sContrast)
    rowHC = summaryTable(summaryTable.stimContrast == sContrast(i) & ...
        summaryTable.reportedConf == 1, :);
    mseHC(i) = rowHC.fun1_rawOriError;

    rowLC = summaryTable(summaryTable.stimContrast == sContrast(i) & ...
        summaryTable.reportedConf == 0, :);
    if ~isempty(rowLC)
        mseLC(i) = rowLC.fun1_rawOriError;
    else
        mseLC(i) = NaN;
    end

end

subplot(2, 2, 4)
hold on
plot(1:numel(mseHC), mseHC, 'o-', LineWidth=2, DisplayName="HC");
plot(1:numel(mseLC), mseLC, 'o-', LineWidth=2, DisplayName="LC");
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('MSE error');
grid on;
hold off
