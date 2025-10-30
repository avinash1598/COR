close all
clear all

data = load('/Users/avinashranjan/Desktop/UT Austin/Goris lab/Uncertainty/Stimuli/COR/Data/COR14.mat'); % dur, contrast, dispersion
% data = load('../Data/COR15.mat'); % Changed reward function - affects perceptual variance
% data = load('../Data/COR16.mat'); % orientation dependent reward
% data = load('../Data/COR17.mat'); % ori dependent reward - changed reward function c1 = 5
% COR18 - c1=1, c2=0.3, -3 (all HC)
% COR19 - c1=1, c2=0.3, -0.5 (all HC)

stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

absErr = abs(reportedOri - stimOri);
absOriError = min(absErr, 180 - absErr);  % ensures error is in [0, 90]
rawError = reportedOri - stimOri;
rawOriError = mod(rawError + 90, 180) - 90;

data.dat.rawOriError = rawOriError;

% This might be needed if the extreme values represent random guesses.
% Filter data in appropriate range (should i do this? Is this really an outlier. Not really probably!)
median_val = median(rawOriError);
mad_val = median(abs(rawOriError - median_val));
mad_val = mad_val / 0.6745;
lower_bound = median_val - 3*mad_val;
upper_bound = median_val + 3*mad_val;
fltIdx = rawOriError < lower_bound | rawOriError > upper_bound;
% rawOriError(fltIdx) = NaN;
% absOriError(fltIdx) = NaN; % Filter abs ori error as well

data.dat.absOriErrorFlt = absOriError;
data.dat.rawOriErrorFlt = rawOriError;

% Error wrt sample mean
stimOri_S = data.dat.stimSampleMeanOri;
reportedOri_S = data.dat.reportedOri;

absErr_S = abs(reportedOri_S - stimOri_S);
absOriError_S = min(absErr_S, 180 - absErr_S);  % ensures error is in [0, 90]
rawError_S = reportedOri_S - stimOri_S;
rawOriError_S = mod(rawError_S + 90, 180) - 90;

data.dat.rawOriError_S = rawOriError_S;

% This might be needed if the extreme values represent random guesses.
% Filter data in appropriate range (should i do this? Is this really an outlier. Not really probably!)
median_val = median(rawOriError_S);
mad_val = median(abs(rawOriError_S - median_val));
mad_val = mad_val / 0.6745; % Does this use gaussian assumption?
lower_bound = median_val - 3*mad_val;
upper_bound = median_val + 3*mad_val;
fltIdx = rawOriError_S < lower_bound | rawOriError_S > upper_bound;
% rawOriError_S(fltIdx) = NaN;
% absOriError_S(fltIdx) = NaN; % Filter absolute ori error as well

data.dat.absOriError_S_Flt = absOriError_S;
data.dat.rawOriError_S_Flt = rawOriError_S;

%% Flag to decide whether to use filtered data or not
useFilteredData = true;

if ~useFilteredData
    selectedRawOriError = data.dat.rawOriError;
else
    selectedRawOriError = data.dat.rawOriErrorFlt;
end

%% Prepare data structures
% TODO: Arrange the data in following format nlevels, norientations, ntrialPerOri
% Once I have this datastructure I can just use the existing analysis
% script

[G, contrastLevels, spreadLevels, stimDur] = findgroups(data.dat.stimContrast, data.dat.stimSpread, data.dat.stimDur);
grpIdxes           = unique(G);
uncertainty_levels = numel(grpIdxes);
nTheta             = 10;
% nTrials            = numel(data.dat.rawOriError(G == grpIdxes(1))); %rawOriError
nTrials            = numel(selectedRawOriError(G == grpIdxes(1)));

% Arrange groups in increasing uncertainty level
fitStds = zeros(1, uncertainty_levels);

for i=1:uncertainty_levels
    grpIdx          = grpIdxes(i);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    % grpRawErr_Flt   = data.dat.rawOriError(G == grpIdx); % rawOriErrorFlt
    
    fltOriErr = grpRawErr_Flt;
    pd = fitdist(fltOriErr(~isnan(fltOriErr)), 'Normal');
    
    mu = pd.mu;
    sigma = pd.sigma;
    fitStds(i) = sigma;
end

[B_, idx] = sort(fitStds);
sidx = grpIdxes(idx);

% Arrange in structured format - same as the one used for analysis
theta_true_all        = zeros(uncertainty_levels, nTrials);
theta_resp_all        = zeros(uncertainty_levels, nTrials); % Recorded theta based on user response
confidence_report_all = zeros(uncertainty_levels, nTrials);
resp_err_all          = zeros(uncertainty_levels, nTrials);

for i=1:uncertainty_levels
    grpIdx          = sidx(i);
    % grpRawErr_Flt   = data.dat.rawOriError(G == grpIdx); % rawOriError 
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    grpStimOri      = data.dat.stimOri(G == grpIdx);
    grpOriResp      = data.dat.reportedOri(G == grpIdx);
    grpReportedConf = data.dat.reportedConf(G == grpIdx);
    
    theta_true_all(i, :)          = grpStimOri;
    theta_resp_all(i, :)          = grpOriResp;
    confidence_report_all(i, :)   = grpReportedConf;
    resp_err_all(i, :)            = grpRawErr_Flt;
end

% resp_err_all = (theta_resp_all - theta_true_all);
% resp_err_all = mod(resp_err_all + 90, 180) - 90; % Find minimum acute angle error

%% Check distribution of stim orientations for each uncertianty level (Should be uniform)
figure

for i=1:uncertainty_levels
    subplot(2, 4, i)
    histogram(theta_true_all(i, :), 'BinEdges', 0:10:180, Normalization='pdf');
    xlabel("Orientation")
    ylabel("Count")
    title(sprintf("Uncetainty level %.2f \nMean %d", i, mean(theta_true_all(i, :))));
end

%% Orientation bias
stimOri = data.dat.stimOri;
bins = linspace(0, 180, 18);
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2;
binIdx = discretize(stimOri, bins);

meanStimOriErr = zeros(1, numel(unique(binIdx)));
stdStimOriErr  = zeros(1, numel(unique(binIdx)));
mseErrs        = zeros(1, numel(unique(binIdx)));
avgErrs        = zeros(1, numel(unique(binIdx)));
stdErrs        = zeros(1, numel(unique(binIdx)));
stdDevs        = zeros(1, numel(unique(binIdx)));

for i=1:numel(unique(binIdx))
    stimOrisInThisBin = data.dat.stimOri(binIdx == i);
    meanStimOriErr(i) = mean(stimOrisInThisBin) - binCenters(i);
    stdStimOriErr(i)  = std(stimOrisInThisBin);
    
    % Calculate mean and std error
    % errorsInThisBin = data.dat.rawOriError(binIdx == i); % rawOriError 
    errorsInThisBin = selectedRawOriError(binIdx == i);
    errorsInThisBin = errorsInThisBin(~isnan(errorsInThisBin));
    
    mseErrs(i) = sum(errorsInThisBin.^2)/numel(errorsInThisBin);
    avgErrs(i) = mean(errorsInThisBin);
    stdErrs(i) = std(errorsInThisBin) / sqrt(numel(errorsInThisBin)); % Plot this as bar as well
    stdDevs(i) = std(errorsInThisBin);
end

figure 

subplot(4, 1, 1)
hold on
errorbar(binCenters, meanStimOriErr, stdStimOriErr, ...
    'k', 'LineStyle', '-', 'LineWidth', 1.5, 'CapSize', 10)
xlabel("Orientation (deg)")
ylabel(sprintf("Mean StimOri Err \n(deg)"))
yline(0, LineStyle="--")
title("Stimulus orientation")
hold off

subplot(4, 1, 2)
hold on
errorbar(binCenters, avgErrs, stdErrs, ...
    'k', 'LineStyle', '-', 'LineWidth', 1.5, 'CapSize', 10)
xlabel("Orientation (deg)")
ylabel("Error")
yline(0, LineStyle="--")
title("Orientation bias")
hold off

subplot(4, 1, 3)
hold on
errorbar(binCenters, avgErrs - meanStimOriErr, stdErrs, ...
    'k', 'LineStyle', '-', 'LineWidth', 1.5, 'CapSize', 10)
xlabel("Orientation (deg)")
ylabel("Error")
yline(0, LineStyle="--")
title(sprintf("Orientation bias \n(Corrected for stimOri distribution)"))
hold off

subplot(4, 1, 4)
x = 1:numel(unique(binIdx));
width = 0.25;
hold on 
bar(x - width/2, stdDevs, width, DisplayName="Std Dev")
% bar(x - width/2, stdErrs, width, DisplayName="Std Dev")
bar(x + width/2, sqrt(mseErrs), width, DisplayName="RMSE Err")
xlabel("Orientation (deg)")
xticks(x)
xticklabels( round(binCenters))
legend
title(sprintf("Orientation dependent std err \n(perceptual err)"))
hold off


%% Orientation depenent error - split by confidence
% Mean error split by confidence

stimOri_all_flat            = theta_true_all(:);
resp_err_all_flat           = resp_err_all(:); % resp_err_all % something to do with NaN
confidence_report_all_flat  = confidence_report_all(:);

bins = linspace(0, 180, 18);
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2;
binIdx = discretize(stimOri_all_flat, bins);

stdDevHC       = zeros(1, numel(unique(binIdx)));
stdDevLC       = zeros(1, numel(unique(binIdx)));
stdErrHC       = zeros(1, numel(unique(binIdx)));
stdErrLC       = zeros(1, numel(unique(binIdx)));
meanErrHC      = zeros(1, numel(unique(binIdx)));
meanErrLC      = zeros(1, numel(unique(binIdx)));
meanStimOriErrHC = zeros(1, numel(unique(binIdx)));
meanStimOriErrLC = zeros(1, numel(unique(binIdx)));

for i=1:numel(unique(binIdx))
    % Get all the orientations within this bin
    idxes_ = (binIdx == i);
    binned_stimOri     = stimOri_all_flat(idxes_);
    binned_err         = resp_err_all_flat(idxes_); % something to do with NaN
    binned_conf_report = confidence_report_all_flat(idxes_);
    
    fltIDx             = ~isnan(binned_err);
    binned_err         = binned_err(fltIDx);
    binned_stimOri     = binned_stimOri(fltIDx);
    binned_conf_report = binned_conf_report(fltIDx);
    
    % Split the reports by HC and LC wihtin bin
    idxHC = binned_conf_report == 1;
    idxLC = binned_conf_report == 0;
    
    meanErrHC(i) = mean(binned_err(idxHC)); 
    meanErrLC(i) = mean(binned_err(idxLC));
    
    stdDevHC(i) = std(binned_err(idxHC)); 
    stdDevLC(i) = std(binned_err(idxLC)); 
    
    stdErrHC(i) = std(binned_err(idxHC))/sqrt(numel(binned_err(idxHC)));
    stdErrLC(i) = std(binned_err(idxLC))/sqrt(numel(binned_err(idxLC))); 
    
    meanStimOriErrHC(i) = mean(binned_stimOri) - binCenters(i); 
    meanStimOriErrLC(i) = mean(binned_stimOri) - binCenters(i); 
end

figure

subplot(3, 1, 1)
hold on
plot(binCenters, meanStimOriErrHC, 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'HC')
plot(binCenters, meanStimOriErrLC, 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'LC')
xlabel("Orientation (deg)")
ylabel("Mean Stim Error")
yline(0, LineStyle="--", HandleVisibility="off")
title(sprintf("Stimulus distribution"))
legend
ylim([-5, 5])
hold off

subplot(3, 1, 2)
hold on
plot(binCenters, meanErrHC, 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'HC')
plot(binCenters, meanErrLC, 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'LC')
xlabel("Orientation (deg)")
ylabel("Mean Error")
yline(0, LineStyle="--", HandleVisibility="off")
title(sprintf("Orientation bias"))
legend
ylim([-10, 10])
hold off

subplot(3, 1, 3)
hold on
sx = 1:numel(unique(binIdx));
% bar(x, stdDevHC)
% bar(x, stdDevLC)
plot(binCenters, stdDevHC, 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'HC')
plot(binCenters, stdDevLC, 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'LC')
xlabel("Orientation (deg)")
ylabel("Std dev (Error)")
yline(0, LineStyle="--", HandleVisibility="off")
title(sprintf("Orientation dependent std dev"))
% xticks(x)
% xticklabels( round(binCenters))
legend
hold off


%% Orientation dependent std split by uncertainty level
bins = linspace(0, 180, 10);
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2;
stdDevByUncertaintyLevel = zeros(uncertainty_levels, numel(binCenters));

for i=1:uncertainty_levels % [1, 8]
    grpIdx          = sidx(i);
    grpStimOri      = data.dat.stimOri(G == grpIdx);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    % grpRawErr_Flt   = data.dat.rawOriError(G == grpIdx); % rawOriErrorFlt 
    
    binIdx          = discretize(grpStimOri, bins);
    
    for j=1:numel(unique(binIdx))
        idxes_ = (binIdx == j);
        x1 = grpRawErr_Flt(idxes_); % Remove NaN if it has it
        stdDevByUncertaintyLevel(i, j) = std(x1, 'omitnan');
    end
end


% % Fit a and b
% coeff_arr_a = zeros(1, uncertainty_levels);
% coeff_arr_b = zeros(1, uncertainty_levels);
% 
% % opts = fitoptions('Method', 'NonlinearLeastSquares', ...
% %     'StartPoint', [5*rand, 5*rand], ...
% %     'Lower', [0, 0]); % ensures a ≥ 0, b ≥ 0
% % 
% % ft = fittype('b + a*abs(sin(2*x))', ...
% %     'independent', 'x', ...
% %     'coefficients', {'a', 'b'}, ...
% %     'options', opts);
% 
% for i=1:uncertainty_levels
%     % f = fit(binCenters', stdDevByUncertaintyLevel(i, :)', ft);
%     
%     x = binCenters;
%     y = stdDevByUncertaintyLevel(i, :);
% 
%     model = @(p, x) p(2) + p(1).*abs(sin(2.*x)); % Define model: y_hat = b + a*abs(sin(2*x))
%     fitness = @(p) sum( (y - model(p,x)).^2 );   % Objective for GA: sum of squared errors
%     
%     % Bounds: force positivity if desired (change as needed)
%     lb = [0, 0];        % a >= 0, b >= 0
%     ub = [Inf, Inf];    % no upper bound (or choose sensible finite upper bounds)
%     
%     % GA options (tweak PopulationSize/MaxGenerations for speed vs quality)
%     gaOpts = optimoptions('ga', ...
%         'Display', 'iter', ...
%         'PopulationSize', 100, ...
%         'MaxGenerations', 200, ...
%         'UseParallel', false);   % set true if you have Parallel Toolbox
%     
%     % Run the genetic algorithm (2 variables)
%     nvars = 2;
%     [bestP, bestFval, exitflag, output] = ga(fitness, nvars, [], [], [], [], lb, ub, [], gaOpts);
% 
% 
%     coeff_arr_a(i) = bestP(1);
%     coeff_arr_b(i) = bestP(2);
% end

% Plot result

figure
subplot(1, 2, 1)

hold on

for i=[1, 8] %1:uncertainty_levels % [1, 8]
    % Define sine model: a*sin(b*x + c) + d
    plot(binCenters, stdDevByUncertaintyLevel(i, :) , 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', sprintf("%d", i))
    % scatter(binCenters, stdDevByUncertaintyLevel(i, :), 'filled');
    % plot(binCenters, coeff_arr_b(i) + coeff_arr_a(i)*abs(sind(2*binCenters)));
end

xlabel("Orientation (deg)")
ylabel("Std dev (Error)")
title(sprintf("Orientation dependent std dev \n(split by uncertainty)"))
legend
hold off

subplot(1, 2, 2)
coeff_arr_a = [3.7, 6, 6.1, 7.9, 6, 7, 5.3, 6];
coeff_arr_b = [5, 7, 7.9, 6.9, 8, 12, 13.7, 11.34];
avg_b_by_a = mean(coeff_arr_b./coeff_arr_a);
scatter(coeff_arr_a, coeff_arr_b, 'filled')
ylabel("Coefficient a")
xlabel("Coefficient b")
title(sprintf("Average b/a: %.2f", avg_b_by_a))

%% Orientation dependent bias split by uncertainty level

bins = linspace(0, 180, 10);
binCenters = ( bins(1:end-1) + bins(2:end) ) / 2;
biasByUncertaintyLevel = zeros(uncertainty_levels, numel(binCenters));

for i=1:uncertainty_levels % [1, 8]
    grpIdx          = sidx(i);
    grpStimOri      = data.dat.stimOri(G == grpIdx);
    grpRawErr_Flt   = selectedRawOriError(G == grpIdx);
    % grpRawErr_Flt   = data.dat.rawOriError(G == grpIdx); % rawOriErrorFlt 
    
    binIdx          = discretize(grpStimOri, bins);
    
    for j=1:numel(unique(binIdx))
        idxes_ = (binIdx == j);
        x1 = grpRawErr_Flt(idxes_); % Remove NaN if it has it
        x1 = x1(~isnan(x1));
        biasByUncertaintyLevel(i, j) = mean(x1) ;
    end
end

% Plot result
figure
hold on

for i=[1, 5, 6, 8] %1:uncertainty_levels % [1, 8]
    plot(binCenters, biasByUncertaintyLevel(i, :) , 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', sprintf("%d", i))
end

xlabel("Orientation (deg)")
ylabel("Mean Error)")
title(sprintf("Orientation dependent bias \n(split by uncertainty)"))
yline(0, LineStyle="--")
legend

