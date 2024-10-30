%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
wav_path = fullfile(file_path, sprintf('../../figures/%s/wavs/', exp_name)); % get path to raw waveforms
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
matched_permutation_test_file = fullfile(file_path, sprintf('../../results/%s/matched_permutation_test.mat', exp_name));
load(matched_permutation_test_file)
indices = [6, 8, 10]; % corresponds to timepoints 1500, 2000, 2500 ms
test_indices = setdiff(1:10, indices);
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/timing_analysis.mat', exp_name)); % get path to timing analysis
load(timing_analysis_file)
elbow_points = est_params(:, 3);

% load exp1 timing analysis data for comparison
timing_analysis_file = fullfile(file_path, '../../results/exp7/timing_analysis_exp1_and_exp7.mat'); % get path to timing analysis
load(timing_analysis_file)
elbow_points_comparison = est_params(:, 3);

load(fullfile(file_path, '../../assets/stationarity_scores/audioset_stationarity_scores.mat'))
temp_table = readtable(fullfile(file_path, '../../assets/stationarity_scores/spectrum_matched_stationarity.csv'));
stationarity_scores = stationarity_scores;
smn_stationarity = temp_table.stationarity_score(1:160);
texture_stationarity = temp_table.stationarity_score(161:320);

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% set colors
exp1_color = [34, 113, 178] / 255; % "honolulu blue" see http://mkweb.bcgsc.ca/colorblind/palettes.mhtml (Martin Krzywinski — Designing for Color Blindness) 
main_color = [247, 72, 165] / 255; % "barbie pink" see http://mkweb.bcgsc.ca/colorblind/palettes.mhtml (Martin Krzywinski — Designing for Color Blindness) 
clrs = [main_color; exp1_color];

% settings for cochleagrams...note must match saved cochleagrams in wav_path
stims = [22 66];
labels = {'wind', 'typing'};
x_spacing = 550;
y_spacing = 25;
x_end = 0;
y_end = -y_spacing;
scale = 2000;

%% plot diagram

% initialize diagram figure
figure('units','inches','innerposition',[0 0 10 4])
hold on

% loop over trials
for n = 1:length(stims)
    % load background cochleagrams    
    background_file = sprintf('cochleagram_%d.mat', stims(n));
    load(fullfile(wav_path, background_file))

    background_file = sprintf('smn_cochleagram_%d.mat', stims(n));
    load(fullfile(wav_path, background_file))
    
    % set x position
    x_start = 0;
    x_end = x_start + size(cochleagram, 1) - 1;
    
    % set y position
    y_start = y_end + y_spacing;
    y_end = y_start + size(cochleagram, 2) - 1;

    % plot texture cochleagram
    imagesc(x_start, y_start, cochleagram')
    rectangle('Position',[x_start y_start x_end-x_start y_end-y_start],'EdgeColor', exp1_color, 'LineWidth', 3)
    
    % add text label
    x_pos = x_start - 0.5*x_spacing;
    y_pos = mean([y_start y_end]);
    text(x_pos, y_pos, labels{n}, ...
     'FontSize', 18, ...
     'FontWeight', 'bold', ...
     'Color',ones(3,1)*0, ...
     'HorizontalAlignment', 'center', ...
     'Rotation', 90)

    % set x position
    x_start = x_end + x_spacing;      
    x_end = x_start + size(cochleagram, 1) - 1;    

    % plot spectrally matched noise cochleagram
    imagesc(x_start, y_start, smn_cochleagram')
    rectangle('Position',[x_start y_start x_end-x_start y_end-y_start],'EdgeColor', clrs(1,:), 'LineWidth', 3)

    % plot mean excitations
    x_pos = x_end + x_spacing*0.25;  
    plot(x_pos + mean(cochleagram)*scale, y_start:y_end , 'LineWidth', 3, 'Color', exp1_color)
    plot(x_pos + mean(smn_cochleagram)*scale, y_start:y_end , 'LineWidth', 3, 'Color', clrs(1,:))        
end

% reset x position
x_start = 0;
x_end = x_start + size(cochleagram, 1) - 1;

% add time axis label
x_pos = mean([x_start x_end]);
y_pos = y_start - 0.3*y_spacing;
t = text(x_pos, y_pos, 'time', ...
     'FontSize', 16, ...
     'FontAngle','normal', ...
     'Color',ones(3,1)*0, ...
     'HorizontalAlignment', 'center');
plot([x_start, t.Extent(1) - 0.07*x_spacing], [y_pos y_pos], 'k', 'LineWidth', 1.5)
plot([t.Extent(1) + t.Extent(3) + 0.07*x_spacing x_end], [y_pos y_pos], 'k', 'LineWidth', 1.5)

% add texture background title
y_pos = y_end + y_spacing;
text(x_pos, y_pos, sprintf('Texture\nBackgrounds'), ...
     'FontSize', 18, ...
     'FontWeight', 'bold', ...
     'Color', exp1_color, ...
     'HorizontalAlignment', 'center', ...
     'Rotation', 0)

% add frequency axis label
x_pos = x_end + 0.4*x_spacing;
y_pos = mean([y_start y_end]);
t = text(x_pos, y_pos, 'frequency', ...
     'FontSize', 16, ...
     'FontAngle','normal', ...
     'Color',ones(3,1)*0, ...
     'HorizontalAlignment', 'center', ...
     'Rotation', -90);
plot([x_pos x_pos], [y_start, t.Extent(2) - 0.07*y_spacing] , 'k', 'LineWidth', 1.5)
plot([x_pos x_pos], [t.Extent(2) + t.Extent(4) + 0.07*y_spacing y_end], 'k', 'LineWidth', 1.5)

% adjust x position
x_start = x_end + x_spacing;      
x_end = x_start + size(cochleagram, 1) - 1; 

% add spectrally matched noise background title
x_pos = mean([x_start x_end]);
y_pos = y_end + y_spacing;
text(x_pos, y_pos, sprintf('Spectrally Matched\nNoise Backgrounds'), ...
     'FontSize', 18, ...
     'FontWeight', 'bold', ...
     'Color', clrs(1,:), ...
     'HorizontalAlignment', 'center', ...
     'Rotation', 0)

% add mean excitation label
text(x_end + 0.9*x_spacing, y_pos, sprintf('mean\nexcitation'), ...
     'FontSize', 16, ...
     'FontAngle','normal', ...
     'Color',ones(3,1)*0, ...
     'HorizontalAlignment', 'center', ...
     'Rotation', 0)

% adjust colormap
colormap(flipud(colormap(gray)))
caxis([0, 0.4])

% make look nicer
ylim([-y_spacing*0.5, y_end + 1.8*y_spacing])
xlim([-1*x_spacing, x_end + 3*x_spacing])
axis off

% save diagram figure
set(gca,'LooseInset',get(gca,'TightInset')); % remove space around edges
ax = gca;
ax.Position = [0 0 1 1];
save_name = fullfile(save_path, sprintf('%s_diagram.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% plot results figure

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
d_primes = (d_primes(1:10,:) + d_primes(11:20,:)) / 2;
n_subs = size(d_primes,2);

% load comparison data from exp 1
comparison_exp = 'exp1';
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp)); % get path to results
load(results_file)
d_primes_comparison = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_comparison(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end
d_primes_comparison1 = d_primes_comparison(:, logical(keep_subs));
d_primes_comparison1 = (d_primes_comparison1(1:10, :) + d_primes_comparison1(11:20, :)) / 2;

% load comparison data from exp 7
comparison_exp = 'exp7';
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp)); % get path to results
load(results_file)
d_primes_comparison = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_comparison(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end
d_primes_comparison7 = d_primes_comparison(:, logical(keep_subs));

d_primes_comparison = [d_primes_comparison1, d_primes_comparison7];
n_subs_comparison = size(d_primes_comparison,2);

% fit elbow function
x = linspace(250,2500,10);
x_test = linspace(250, 2500, 1000);

y = mean(d_primes,2)';
est_params = fit_elbow_func(x,y);
actual_elbow_point = est_params(3);
y_pred = elbow_func(x_test, est_params);

y_comparison = mean(d_primes_comparison,2,'omitnan')';
est_params_comparison = fit_elbow_func(x,y_comparison);
actual_elbow_point_comparison = est_params_comparison(3);
y_pred_comparison = elbow_func(x_test, est_params_comparison);

% initialize figure
figure('units','inches','innerposition',[0 0 9 11])
hold on
timevec = linspace(250,2500,10);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% plot elbow point stats 
med = median(elbow_points);
stdev_ep = std(elbow_points);
ypos = -0.08;
plot([med - stdev_ep, med + stdev_ep], [ypos, ypos], 'LineWidth', 3, 'Color', clrs(1,:))
scatter(actual_elbow_point, ypos, 150, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)

med = median(elbow_points_comparison);
stdev_ep = std(elbow_points_comparison);
ypos = -0.14;
plot([med - stdev_ep, med + stdev_ep], [ypos, ypos], 'LineWidth', 3, 'Color', clrs(2,:))
scatter(actual_elbow_point_comparison, ypos, 150, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot shaded errorbars (SEM)
m = mean(d_primes,2);
s = std(d_primes,0,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub', flip(lb')];
ptch = patch(x_patch, y_patch, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

m = mean(d_primes_comparison,2);
s = std(d_primes_comparison,0,2)/sqrt(n_subs_comparison);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot elbow function fit to data
plot(x_test, y_pred, '--', 'LineWidth', 3, 'Color', clrs(1,:))
plot(x_test, y_pred_comparison, '--', 'LineWidth', 3, 'Color', clrs(2,:));

% plot human results
scatter(timevec, mean(d_primes,2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes_comparison,2), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot dividing line
plot([0 2750], [0, 0], 'k', 'LineWidth', 2)

% plot deltas
delta_x = 2650;
plot([delta_x, delta_x], [y_pred(1) y_pred(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred(end) y_pred(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred(1) y_pred(1)], 'LineWidth', 2, 'Color', clrs(1,:))
delta_x = 2600;
plot([delta_x, delta_x], [y_pred_comparison(1) y_pred_comparison(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_comparison(end) y_pred_comparison(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_comparison(1) y_pred_comparison(1)], 'LineWidth', 2, 'Color', clrs(2,:))

% add number of subjects
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',clrs(1,:))
text(75, 2.75, sprintf('n = %d', n_subs_comparison), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',clrs(2,:))

% stylize figure and add labels
box off
ylabel('Detection Performance (d'')', 'Position', [-300, 1.5])
xlabel('Foreground Onset Time (ms)')
ylim([-0.25 3])
yticks(0:0.5:3)
ytickformat('%.1f')
xlim([0 2750])
xticks((timevec))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% improve look of legend
[l, objh]  = legend({'  More stationary (Exp. 9)', '  Less stationary (Exps. 1 & 7)'},'Location','northeast', 'FontSize', 25);
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.EdgeColor = 'none';
text(830, 2.9, 'Background Type', ...
     'FontSize', 26, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', ...
     'Color',ones(3,1)*0)

% reposition
set(gca,'Position', [0.145, 0.165, (1 - 0.165), (1 - 0.195)]) % remove space around edges
l.Position(1) = l.Position(1)*1;
l.Position(2) = l.Position(2)*0.96;

% save as svg
save_name = fullfile(save_path, sprintf('%s_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% plot matched performance inset

d_primes_matched = d_primes(:,matched_subs);
d_primes_matched_comparison = d_primes_comparison(:, matched_subs_comparison);

% fit elbow function
x = linspace(250,2500,10);
x = x(test_indices);
y = mean(d_primes_matched(test_indices,:),2,'omitnan')';
y_comparison = mean(d_primes_matched_comparison(test_indices,:),2,'omitnan')';
est_params = fit_elbow_func(x,y);
est_params_comparison = fit_elbow_func(x,y_comparison);
x_test = linspace(250, 2500, 1000);
y_pred = elbow_func(x_test, est_params);
y_pred_comparison = elbow_func(x_test, est_params_comparison);

% initialize figure
y_scale = 1.2 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 9 11*y_scale])
hold on
timevec = linspace(250,2500,10);

% plot shaded errorbars (SEM)
m = mean(d_primes_matched(test_indices,:),2,'omitnan');
s = std(d_primes_matched(test_indices,:),0,2,'omitnan')/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(test_indices), flip(timevec(test_indices))];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

m = mean(d_primes_matched_comparison(test_indices,:),2,'omitnan');
s = std(d_primes_matched_comparison(test_indices,:),0,2,'omitnan')/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(test_indices), flip(timevec(test_indices))];
y_patch = [ub(:,1)', flip(lb(:,1)')];
ptch = patch(x_patch, y_patch, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot elbow function fit to data
p1 = plot(x_test, y_pred, '--', 'LineWidth', 3, 'Color', clrs(1,:));
p2 = plot(x_test, y_pred_comparison, '--', 'LineWidth', 3, 'Color', clrs(2,:));

% plot human results
scatter(timevec(test_indices), mean(d_primes_matched(test_indices,:),2,'omitnan'), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec(test_indices), mean(d_primes_matched_comparison(test_indices,:),2,'omitnan'), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot deltas
delta_x = 2650;
plot([delta_x, delta_x], [y_pred(1) y_pred(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred(end) y_pred(end)], 'LineWidth', 2, 'Color', clrs(1,:))
plot([delta_x - 25, delta_x + 5], [y_pred(1) y_pred(1)], 'LineWidth', 2, 'Color', clrs(1,:))
delta_x = 2600;
plot([delta_x, delta_x], [y_pred_comparison(1) y_pred_comparison(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_comparison(end) y_pred_comparison(end)], 'LineWidth', 2, 'Color', clrs(2,:))
plot([delta_x - 25, delta_x + 5], [y_pred_comparison(1) y_pred_comparison(1)], 'LineWidth', 2, 'Color', clrs(2,:))

% add text
text(1750, 1.56, sprintf('Performance\nmatched'), ...
     'FontSize', 32, ...
     'FontAngle', 'italic', ...
     'HorizontalAlignment', 'center', ...
     'Color',ones(3,1)*0)

% stylize figure and add labels
box off
ylabel('d''', 'Position', [-300, 1.75])
xlabel('Foreground Onset Time (ms)')
ylim([1.25 2])
yticks([1.25, 1.5, 2])
yticklabels({'0', '1.5', '2.0'})
xlim([0 2750])
xticks((timevec(test_indices)))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% trim edges
set(gca,'Position', [(0.145), (0.165/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges

% add axis break
ax1 = gca;
ax2 = axes('Position', ax1.Position);
px=[-50 50];
py1=[1.275 1.325];
height=0.04;
py2=py1+height;
hold on
fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','w');
plot(px,py1,'k','LineWidth',2)
plot(px,py2,'k','LineWidth',2)
ax2.YLim = ax1.YLim;
ax2.XLim = ax1.XLim;
ax2.Position(1) = ax2.Position(1)*0.945;
axis off
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf, 'Color', 'w');

save_name = fullfile(save_path, sprintf('%s_matched_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% Plot stationarity histogram

figure('units','inches','innerposition',[0 0 18 4])
hold on

% plot stationarity scores from audioset
histogram(stationarity_scores, 'FaceColor', 0.6*ones(3,1), 'EdgeColor', 'None', 'Normalization','probability', 'FaceAlpha', 0.6);
histogram(stationarity_scores, 'Normalization', 'probability', 'EdgeColor', 0.6*ones(3,1), 'DisplayStyle', 'stairs', 'LineWidth', 3)

% plot gray lines connecting textures to spectrally matched noise stimuli
x1 = 2.25/100;
x2 = 1.75/100;
plot([texture_stationarity, smn_stationarity]', [x1 x2], 'Color', [0.90*ones(1,3),0.25], 'LineWidth', 2)

% plot texture stationarity scores
scatter(texture_stationarity,x1*ones(size(texture_stationarity)), 250, 'MarkerFaceColor', exp1_color, 'MarkerEdgeColor','w', 'MarkerFaceAlpha',0.2)
scatter(mean(texture_stationarity), x1, 250, 'MarkerFaceColor', exp1_color, 'MarkerEdgeColor','w')
plot([mean(texture_stationarity), mean(texture_stationarity)],[0, x1], '-', 'LineWidth', 4, 'Color', exp1_color)

% plot spectrally matched noise stationarity scores
scatter(smn_stationarity,x2*ones(size(smn_stationarity)), 250, 'MarkerFaceColor', clrs(1,:), 'MarkerEdgeColor','w', 'MarkerFaceAlpha',0.2)
scatter(mean(smn_stationarity),x2,  250, 'MarkerFaceColor', clrs(1,:), 'MarkerEdgeColor','w')
plot([mean(smn_stationarity), mean(smn_stationarity)],[0, x2], '-', 'LineWidth', 4, 'Color', clrs(1,:))

% add text
text(0, 0.4/100, 'AudioSet', 'HorizontalAlignment','center','FontSize', 28, 'Color', 1*ones(1,3), 'FontAngle', 'italic', 'Rotation', 0)

% clean figure up
xlim([-1, 1.5])
xtickformat('%.1f')
ylbl = ylabel('Percentage');
ylbl.Position(2) = 1/100;
ylim([0, 2.45/100])
yticks((0:2)/100)
yticklabels([0, 1, 2])
xlabel_string = sprintf('Normalized Standard Deviation of Texture Statistics');
xlabel(xlabel_string)
set(gca, 'FontSize', 28)
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca, 'Layer','top');
ax = gca;
ax.Position(1) = ax.Position(1)*0.45;
ax.Position(3) = ax.Position(3)*1.19;
annotation('rectangle',[0.051 0.815 .014 .115], 'FaceColor','w', 'EdgeColor', 'none'); 

save_name = fullfile(save_path, sprintf('%s_stationarity.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');