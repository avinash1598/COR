close all
clear all

% data = load('COR14.mat');
% data = load('COR15.mat');
data = load('COR16.mat');

stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

absErr = abs(stimOri - reportedOri);
absOriError = min(absErr, 180 - absErr);  % ensures error is in [0, 90]
rawError = stimOri - reportedOri;
rawOriError = mod(rawError + 90, 180) - 90;

data.dat.rawOriErrorOrg = rawOriError;

% Filter data in appropriate range
median_val = median(rawOriError);
mad_val = median(abs(rawOriError - median_val));
mad_val = mad_val / 0.6745;
lower_bound = median_val - 3*mad_val;
upper_bound = median_val + 3*mad_val;
rawOriError(rawOriError < lower_bound | rawOriError > upper_bound) = NaN;

data.dat.absOriError = absOriError;
data.dat.rawOriError = rawOriError;

% Err wrt sample mean
stimOri_S = data.dat.stimSampleMeanOri;
reportedOri_S = data.dat.reportedOri;

absErr_S = abs(stimOri_S - reportedOri_S);
absOriError_S = min(absErr_S, 180 - absErr_S);  % ensures error is in [0, 90]
rawError_S = stimOri_S - reportedOri_S;
rawOriError_S = mod(rawError_S + 90, 180) - 90;

% Filter data in appropriate range
median_val = median(rawOriError_S);
mad_val = median(abs(rawOriError_S - median_val));
mad_val = mad_val / 0.6745;
lower_bound = median_val - 3*mad_val;
upper_bound = median_val + 3*mad_val;
rawOriError_S(rawOriError_S < lower_bound | rawOriError_S > upper_bound) = NaN;

data.dat.absOriError_S = absOriError_S;
data.dat.rawOriError_S = rawOriError_S;


% datOK = data.dat(data.dat.trialStatus == 1, :);

%% Sample spread and mean distribution
% !!!
% IMP
% Note: small spread at all angles will be reflected in
% huge differences in actual spread vs sample spread even for low contrast.
% !!!
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
title("Sample mean error dsitribution")

% Sample Spread distribution
spreadErrs = actualSpreads - sampleSpreads;
subplot(1, 2, 2)
histogram(spreadErrs)
std(~isnan(spreadErrs))
title("Sample spread error distribution")

%% Common to all
% TODO: compute group level metrics here
% Get elements in the group
[G, contrastLevels, spreadLevels, stimDur] = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
grpIdxes = unique(G);
fitStds = zeros(1, numel(grpIdxes));

grpRawErrs_Mean = zeros(1, numel(grpIdxes));
grpRawErrs_Std  = zeros(1, numel(grpIdxes));
grpAbsErrs_Mean = zeros(1, numel(grpIdxes));
grpAbsErrs_Std  = zeros(1, numel(grpIdxes));

grpRawErrs_Mean_by_conf = zeros(2, numel(grpIdxes));
grpRawErrs_Std_by_conf  = zeros(2, numel(grpIdxes));
grpAbsErrs_Mean_by_conf = zeros(2, numel(grpIdxes));
grpAbsErrs_Std_by_conf  = zeros(2, numel(grpIdxes));

grpRawErrs_Mean_S = zeros(1, numel(grpIdxes));
grpRawErrs_Std_S  = zeros(1, numel(grpIdxes));
grpAbsErrs_Mean_S = zeros(1, numel(grpIdxes));
grpAbsErrs_Std_S  = zeros(1, numel(grpIdxes));

grpRawErrs_Mean_by_conf_S = zeros(2, numel(grpIdxes));
grpRawErrs_Std_by_conf_S  = zeros(2, numel(grpIdxes));
grpAbsErrs_Mean_by_conf_S = zeros(2, numel(grpIdxes));
grpAbsErrs_Std_by_conf_S  = zeros(2, numel(grpIdxes));

% TODO: plot RMSE error as well to account for bias

for i=1:numel(grpIdxes)
    gidx_ = grpIdxes(i);
    grpOriErr = rawOriError(G == gidx_);
    grpOriErr_Abs = absOriError(G == gidx_);

    fltOriErr = grpOriErr;
    pd = fitdist(fltOriErr(~isnan(fltOriErr)), 'Normal');
    
    mu = pd.mu;
    sigma = pd.sigma;
    fitStds(i) = sigma;

    grpRawErrs_Mean(i) = mu;
    grpRawErrs_Std(i)  = sigma;
    grpAbsErrs_Mean(i) = mean(grpOriErr_Abs);
    grpAbsErrs_Std(i)  = std(grpOriErr_Abs) / sqrt(numel(grpOriErr_Abs));

    % Split by confidence
    reportedConf = data.dat.reportedConf(G == gidx_);
    idxHC = reportedConf == 1;
    idxLC = reportedConf == 0;

    rawErrHC = grpOriErr(idxHC);
    rawErrLC = grpOriErr(idxLC);
    absErrHC = grpOriErr_Abs(idxHC);
    absErrLC = grpOriErr_Abs(idxLC);
    
    grpAbsErrs_Mean_by_conf(1, i) = mean(absErrLC);
    grpAbsErrs_Mean_by_conf(2, i) = mean(absErrHC);
    grpAbsErrs_Std_by_conf(1, i)  = std(absErrLC) / sqrt(numel(absErrLC));
    grpAbsErrs_Std_by_conf(2, i)  = std(absErrHC) / sqrt(numel(absErrHC));

    pdHC = fitdist(rawErrHC(~isnan(rawErrHC)), 'Normal');
    pdLC = fitdist(rawErrLC(~isnan(rawErrLC)), 'Normal');
    
    grpRawErrs_Std_by_conf(1, i)  = pdLC.sigma;
    grpRawErrs_Std_by_conf(2, i)  = pdHC.sigma;
    grpRawErrs_Mean_by_conf(1, i) = pdLC.mu;
    grpRawErrs_Mean_by_conf(2, i) = pdHC.mu;

    % Wrt sample mean
    grpOriErr_S = rawOriError_S(G == gidx_);
    grpOriErr_Abs_S = absOriError_S(G == gidx_);
    
    fltOriErr_S = grpOriErr_S;
    pd = fitdist(fltOriErr_S(~isnan(fltOriErr_S)), 'Normal');
    
    mu = pd.mu;
    sigma = pd.sigma;

    grpRawErrs_Mean_S(i) = mu;
    grpRawErrs_Std_S(i)  = sigma;
    grpAbsErrs_Mean_S(i) = mean(grpOriErr_Abs_S);
    grpAbsErrs_Std_S(i)  = std(grpOriErr_Abs_S) / sqrt(numel(grpOriErr_Abs_S));

    % Split by confidence
    reportedConf = data.dat.reportedConf(G == gidx_);
    idxHC = reportedConf == 1;
    idxLC = reportedConf == 0;

    rawErrHC = grpOriErr_S(idxHC);
    rawErrLC = grpOriErr_S(idxLC);
    absErrHC = grpOriErr_Abs_S(idxHC);
    absErrLC = grpOriErr_Abs_S(idxLC);
    
    grpAbsErrs_Mean_by_conf_S(1, i) = mean(absErrLC);
    grpAbsErrs_Mean_by_conf_S(2, i) = mean(absErrHC);
    grpAbsErrs_Std_by_conf_S(1, i)  = std(absErrLC) / sqrt(numel(absErrLC));
    grpAbsErrs_Std_by_conf_S(2, i)  = std(absErrHC) / sqrt(numel(absErrHC));
    
    pdHC = fitdist(rawErrHC(~isnan(rawErrHC)), 'Normal');
    pdLC = fitdist(rawErrLC(~isnan(rawErrLC)), 'Normal');
    
    grpRawErrs_Std_by_conf_S(1, i)  = pdLC.sigma;
    grpRawErrs_Std_by_conf_S(2, i)  = pdHC.sigma;
    grpRawErrs_Mean_by_conf_S(1, i) = pdLC.mu;
    grpRawErrs_Mean_by_conf_S(2, i) = pdHC.mu;


end

% Sort conditions by uncertainty difference
[B_, idx] = sort(fitStds);
sidx = grpIdxes(idx);

sSpread   = spreadLevels(sidx);
sContrast = contrastLevels(sidx);
sDur      = stimDur(sidx);

grpRawErrs_Mean = grpRawErrs_Mean(sidx);
grpRawErrs_Std  = grpRawErrs_Std(sidx);
grpAbsErrs_Mean = grpAbsErrs_Mean(sidx);
grpAbsErrs_Std  = grpAbsErrs_Std(sidx);

grpAbsErrs_Mean_by_conf(1, :) = grpAbsErrs_Mean_by_conf(1, sidx);
grpAbsErrs_Mean_by_conf(2, :) = grpAbsErrs_Mean_by_conf(2, sidx);
grpAbsErrs_Std_by_conf(1, :)  = grpAbsErrs_Std_by_conf(1, sidx);
grpAbsErrs_Std_by_conf(2, :)  = grpAbsErrs_Std_by_conf(2, sidx);

grpRawErrs_Mean_by_conf(1, :) = grpRawErrs_Mean_by_conf(1, sidx);
grpRawErrs_Mean_by_conf(2, :) = grpRawErrs_Mean_by_conf(2, sidx);
grpRawErrs_Std_by_conf(1, :)  = grpRawErrs_Std_by_conf(1, sidx);
grpRawErrs_Std_by_conf(2, :)  = grpRawErrs_Std_by_conf(2, sidx);

% WRT sample mean
grpRawErrs_Mean_S = grpRawErrs_Mean_S(sidx);
grpRawErrs_Std_S  = grpRawErrs_Std_S(sidx);
grpAbsErrs_Mean_S = grpAbsErrs_Mean_S(sidx);
grpAbsErrs_Std_S  = grpAbsErrs_Std_S(sidx);

grpAbsErrs_Mean_by_conf_S(1, :) = grpAbsErrs_Mean_by_conf_S(1, sidx);
grpAbsErrs_Mean_by_conf_S(2, :) = grpAbsErrs_Mean_by_conf_S(2, sidx);
grpAbsErrs_Std_by_conf_S(1, :)  = grpAbsErrs_Std_by_conf_S(1, sidx);
grpAbsErrs_Std_by_conf_S(2, :)  = grpAbsErrs_Std_by_conf_S(2, sidx);

grpRawErrs_Mean_by_conf_S(1, :) = grpRawErrs_Mean_by_conf_S(1, sidx);
grpRawErrs_Mean_by_conf_S(2, :) = grpRawErrs_Mean_by_conf_S(2, sidx);
grpRawErrs_Std_by_conf_S(1, :)  = grpRawErrs_Std_by_conf_S(1, sidx);
grpRawErrs_Std_by_conf_S(2, :)  = grpRawErrs_Std_by_conf_S(2, sidx);

% Generate x-axis labels as combined strings like "0.2|5"
xLabels = strcat(string(sDur), ', ', string(sContrast), ', ', string(sSpread));

%% Figure 1 - Raw errors by group - Estimation error (Human subject)

% Define visual properties
colors = [0.85 0.33 0.10; 0 0.45 0.74; 0.85 0.33 0.10; 0 0.45 0.74; 0 0.45 0.74; 0 0.45 0.74]; % colorblind-friendly red/blue
alphas = [1, 1, 1, 0.3, 0.3, 0.3, 0.3];

% Create figure
figure('Color','w');
pbaspect([3 5 1]) 
hold on

% Plot each group
for i = 1:numel(grpIdxes)
    grpOriErr = data.dat.rawOriErrorOrg(G == sidx(i));
    meanErr = mean(grpOriErr(~isnan(grpOriErr)));
    xPts = i + 0.4*(rand(1, numel(grpOriErr)) - 0.5);

    scatter(xPts, grpOriErr, 60, ...
        'LineWidth', 0.5);
%         'MarkerEdgeColor', 'w', ...
%         'MarkerFaceColor', colors(i,:), ...
%         'MarkerFaceAlpha', alphas(i), ...
%         'MarkerEdgeAlpha', alphas(i), ...
        

    plot([i - 0.35, i + 0.35], [meanErr, meanErr], LineStyle="-", LineWidth=2, Color='black')
end

% Labels and styling
ylabel("Estimation error (°)", 'FontSize', 14)
xlabel("Uncertainty level", 'FontSize', 14)
xticks(1:length(xLabels));
xticklabels(xLabels);
yticks([-90, 0, 90])
ylim([-90, 90])
% xlim([0.3, 7])

% Improve figure aesthetics
set(gca, 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out', 'Box', 'off')
hold off

% exportgraphics(gcf, 'Human1.eps', 'ContentType', 'vector')
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

% exportgraphics(gcf, 'Human2.eps', 'ContentType', 'vector')
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
hold on
errorbar(binCenters, avgErrs, stdErrs, ...
    'k', 'LineStyle', '-', 'LineWidth', 1.5, 'CapSize', 10)
xlabel("Orientation (deg)")
ylabel("Error")
yline(0, LineStyle="--")
title("Orientation bias")
hold off

subplot(2, 1, 2)
x = 1:numel(unique(binIdx));
width = 0.25;
hold on
bar(x - width/2, stdDevs, width, DisplayName="Std Dev")
bar(x + width/2, sqrt(mseErrs), width, DisplayName="RMSE Err")
xlabel("Orientation (deg)")
xticks(x)
xticklabels( round(binCenters))
legend
title("Orientation dependent std dev")
hold off

%% Histogram

% [G, contrastLevels, spreadLevels, stimDur] = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
% grpIdxes = unique(G);

% grpIdxes = zeros(1, numel(sidx));
% fitStds = zeros(1, numel(sidx));

figure
% for i=1:numel(grpIdxes)
for i=1:numel(sidx)
    gidx_ = sidx(i);
%     gidx_ = grpIdxes(i);

    subplot(2, numel(sidx)/2, i)
    grpOriErr = rawOriError(G == gidx_);
    histogram(data.dat.rawOriErrorOrg(G == gidx_), -90:3:90, Normalization="pdf")
    
    pd = fitdist(grpOriErr(~isnan(grpOriErr)), 'Normal');
    x_values = -90:1:90;
    y = pdf(pd, x_values);

    mu = pd.mu;
    sigma = pd.sigma;
    % fitStds(i) = sigma;

    hold on
    plot(x_values, y, LineWidth=1.5)

    val = std(grpOriErr);

    xlabel("Orientation (deg)")
    ylabel("count")
    title(sprintf("Durtion: %.2f, Contrast: %.2f, \n Spread: %.2f \n fit std %.2f, std: %.2f", ...
        stimDur(gidx_), contrastLevels(gidx_), spreadLevels(gidx_), sigma, val))

    hold off
end

% title("Distribution of errors")


%% Histogram by confidence
figure
for confVal = [0, 1]
    for i=1:numel(sidx)
        gidx_ = sidx(i);
    
        subplot(2, numel(sidx), (numel(sidx))*confVal + i)
        grpOriErr = rawOriError(G == gidx_);
        reportedConf = data.dat.reportedConf(G == gidx_);
        grpOriErr = grpOriErr(reportedConf == confVal);
        % histogram(grpOriErr, -90:3:90, Normalization="pdf")
        histogram(data.dat.rawOriErrorOrg((G == gidx_) & (data.dat.reportedConf == confVal)), -90:3:90, Normalization="pdf")
        
        pd = fitdist(grpOriErr(~isnan(grpOriErr)), 'Normal');
        x_values = -90:1:90;
        y = pdf(pd, x_values);
    
        mu = pd.mu;
        sigma = pd.sigma;
        % fitStds(i) = sigma;
    
        hold on
        plot(x_values, y, LineWidth=1.5)
    
        val = std(grpOriErr);
    
        xlabel("Orientation (deg)")
        ylabel("count")
        title(sprintf("fit std %.2f, std: %.2f", ...
             sigma, val))
%         title(sprintf("Durtion: %.2f, Contrast: %.2f, \n Spread: %.2f \n fit std %.2f, std: %.2f", ...
%             stimDur(gidx_), contrastLevels(gidx_), spreadLevels(gidx_), sigma, val))
    
        hold off
    end
end


%% Histogram wrt sample mean

figure
% for i=1:numel(grpIdxes)
for i=1:numel(sidx)
    gidx_ = sidx(i);

    subplot(2, numel(sidx)/2, i)
    grpOriErr = rawOriError_S(G == gidx_);
    histogram(grpOriErr, -90:3:90, Normalization="pdf")
    
%     grpOriErr_rad = deg2rad(grpOriErr);
%     [mu_hat, kappa_hat] = circ_vmpar(grpOriErr_rad);
%     % grpOriErr should be in radians
%     
%     x = linspace(-pi/2, pi/2, 101);
%     y = exp(kappa_hat .* cos(x - mu_hat)) ./ (2*pi*besseli(0, kappa_hat));
    
%     gm = fitgmdist(grpOriErr, 2, 'Options', statset('MaxIter',1000));
%     mu_gm = gm.mu;       % 2 means
%     sigma_gm = squeeze(sqrt(gm.Sigma));   % 2 variances
%     p_comp = gm.ComponentProportion;
%     
%     sigma_net = sqrt( sum( squeeze(gm.Sigma).*p_comp' ) );

%     x = linspace(-90, 90, 101)';  % column vector
%     y = pdf(gm, x);

    % Filter based on MAD estimate (https://www.sciencedirect.com/science/article/pii/S0022103113000668?ref=pdf_download&fr=RR-2&rr=96dd498e3f2970e0)
    % MAD-based estimate
    median_val = median(grpOriErr);
    mad_val = median(abs(grpOriErr - median_val));
    sigma_mad = mad_val / 0.6745;
    mad_val = sigma_mad;
    
    % fltOriErr = grpOriErr( zscore(grpOriErr) < 3 ); % Data within 3 SD
    % fltOriErr = grpOriErr(abs(grpOriErr) < 50);
    % fltOriErr = grpOriErr( abs(grpOriErr) < 50 );
    % fltOriErr = grpOriErr( grpOriErr > (median_val - 3*mad_val) & grpOriErr < (median_val + 3*mad_val) );
    pd = fitdist(grpOriErr(~isnan(grpOriErr)), 'Normal');
    x = -90:1:90;
    y = pdf(pd, x_values);

    hold on
    plot(x, y, LineWidth=1.5)

    val = std(grpOriErr);

    % IQR-based estimate
    Q1 = prctile(grpOriErr, 25);
    Q3 = prctile(grpOriErr, 75);
    iqr_val = Q3 - Q1;
    sigma_iqr = iqr_val / 1.349;
    
    % MAD-based estimate
    median_val = median(grpOriErr);
    mad_val = median(abs(grpOriErr - median_val));
    sigma_mad = mad_val / 0.6745;
    
    xlabel("Orientation (deg)")
    ylabel("count (wrt sample mean)")
    title(sprintf("Durtion: %.2f, Contrast: %.2f, \n Spread: %.2f \n IQR std %.2f, MAD std %.2f, fit std %.2f, std: %.2f", ...
        stimDur(gidx_), contrastLevels(gidx_), spreadLevels(gidx_), sigma_iqr, sigma_mad, pd.std, val))

    hold off
end

%% MSE error by contrast levels and confidence
datOK = data.dat(data.dat.trialStatus == 1, :);

summaryTable = groupsummary(datOK, {'stimDur', 'stimContrast', 'stimSpread'}, ...
                            {'mean', 'std', 'numel'}, 'absOriError');

disp(summaryTable)

meanError = grpAbsErrs_Mean;
stdError = grpAbsErrs_Std;

% Plot with error bars
figure;
subplot(2, 4, 1)
% plot(stdError)
errorbar(meanError, stdError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('Absolute Error');

grid on;

% By confidence
summaryTable = groupsummary(datOK, {'stimDur', 'stimContrast', 'stimSpread', 'reportedConf'}, ...
                            {'mean', 'std', 'numel'}, 'absOriError');

disp(summaryTable)

% Plot with error bars
subplot(2, 4, 2)
hold on
errorbar(grpAbsErrs_Mean_by_conf(1, :), grpAbsErrs_Std_by_conf(1, :), 'o-', LineWidth=2, DisplayName='LC');
errorbar(grpAbsErrs_Mean_by_conf(2, :), grpAbsErrs_Std_by_conf(2, :), 'o-', LineWidth=2, DisplayName='HC');

xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('Absolute Error');
grid on;
legend
hold off

% grpRawErrs_Mean_by_conf = zeros(2, numel(grpIdxes));
% grpRawErrs_Std_by_conf  = zeros(2, numel(grpIdxes));

% meanError = grpRawErrs_Mean;
stdError = grpRawErrs_Std;

% Plot with error bars
subplot(2, 4, 3)
plot(stdError, LineWidth=2)
% errorbar(meanError, stdError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Std dev (deg)');
title('Raw Error');

grid on;

% By confidence
subplot(2, 4, 4)
hold on
plot(grpRawErrs_Std_by_conf(1, :), LineWidth=2, DisplayName='LC')
plot(grpRawErrs_Std_by_conf(2, :), LineWidth=2, DisplayName='HC')

xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Std dev (deg)');
title('Raw Error');
grid on;
legend
hold off


% sample mean
% If subject is using sample mean to report orientation then variability
% wrt sample mean should ideally be less than wrt actual mean.
meanError = grpAbsErrs_Mean_S;
stdError = grpAbsErrs_Std_S;

% Plot with error bars
subplot(2, 4, 5)
% plot(stdError)
errorbar(meanError, stdError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('Absolute Error (sample mean)');

grid on;


% Plot with error bars
subplot(2, 4, 6)
hold on
errorbar(grpAbsErrs_Mean_by_conf_S(1, :), grpAbsErrs_Std_by_conf_S(1, :), 'o-', LineWidth=2, DisplayName='LC');
errorbar(grpAbsErrs_Mean_by_conf_S(2, :), grpAbsErrs_Std_by_conf_S(2, :), 'o-', LineWidth=2, DisplayName='HC');

xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Orientation Error (deg)');
title('Absolute Error');
grid on;
legend
hold off

% grpRawErrs_Mean_by_conf = zeros(2, numel(grpIdxes));
% grpRawErrs_Std_by_conf  = zeros(2, numel(grpIdxes));

% meanError = grpRawErrs_Mean;
stdError = grpRawErrs_Std_S;

% Plot with error bars
subplot(2, 4, 7)
plot(stdError, LineWidth=2)
% errorbar(meanError, stdError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Std dev (deg)');
title('Raw Error');

grid on;

% By confidence
subplot(2, 4, 8)
hold on
plot(grpRawErrs_Std_by_conf_S(1, :), LineWidth=2, DisplayName='LC')
plot(grpRawErrs_Std_by_conf_S(2, :), LineWidth=2, DisplayName='HC')

xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("StimDur, StimContrast, StimSpread" + newline + "(increasing uncertainty)")
ylabel('Std dev (deg)');
title('Raw Error');
grid on;
legend
hold off


%% RMSE error

summaryTable = groupsummary(datOK, {'stimContrast', 'stimSpread', 'stimDur'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'numel'}, 'rawOriError');

meanError = summaryTable.fun1_rawOriError(sidx);

figure
subplot(1, 2, 1)
plot(1:numel(meanError), meanError, 'o-', LineWidth=2);
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("increasing uncertainty")
ylabel('Orientation Error (deg)');
title('RMSE error');
grid on;


% MSE err split by confidence
summaryTable = groupsummary(datOK, {'stimContrast', 'stimSpread', 'stimDur', 'reportedConf'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'numel'}, 'rawOriError');

mseHC = zeros(1, numel(sContrast));
mseLC = zeros(1, numel(sContrast));

for i = 1:numel(sContrast)
    rowHC = summaryTable(summaryTable.stimContrast == sContrast(i) & ...
        summaryTable.stimSpread == sSpread(i) & ...
        summaryTable.stimDur == sDur(i) & ...
        summaryTable.reportedConf == 1, :);
    mseHC(i) = rowHC.fun1_rawOriError;

    rowLC = summaryTable(summaryTable.stimContrast == sContrast(i) & ...
        summaryTable.stimSpread == sSpread(i) & ...
        summaryTable.stimDur == sDur(i) & ...
        summaryTable.reportedConf == 0, :);
    if ~isempty(rowLC)
        mseLC(i) = rowLC.fun1_rawOriError;
    else
        mseLC(i) = NaN;
    end

end

subplot(1, 2, 2)
hold on
plot(1:numel(mseHC), mseHC, 'o-', LineWidth=2, DisplayName="HC");
plot(1:numel(mseLC), mseLC, 'o-', LineWidth=2, DisplayName="LC");
xticks(1:length(xLabels));
xticklabels(xLabels);
xlabel("increasing uncertainty")
ylabel('Orientation Error (deg)');
title('MSE error');
grid on;
hold off



%%

x = linspace(0, 180, 100);
y = 9 + (16 - 9)*abs(sind(2*x));

figure
plot(x, abs(y))
