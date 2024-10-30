%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
wav_path = fullfile(file_path, sprintf('../../figures/%s/wavs/', exp_name)); % get path to raw waveforms
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/timing_analysis.mat', exp_name)); % get path to timing analysis
load(timing_analysis_file)
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% set colors
main_color = [53, 155, 115]; % "ocean green" see http://mkweb.bcgsc.ca/colorblind/palettes.mhtml (Martin Krzywinski — Designing for Color Blindness) 
alt_color = main_color + (255 - main_color) * 0.4;
clrs = [main_color; alt_color] / 255;

%% plot results figure

azimuth_errors = nan(length(results(1).mean_azimuth_error),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    azimuth_errors(:,i) = results(i).mean_azimuth_error;
    keep_subs(i) = results(i).pass;
end

azimuth_errors = azimuth_errors(:, logical(keep_subs));
n_subs = size(azimuth_errors,2);

elbow_points = est_params(:, 3);

% fit elbow function
x = linspace(250,2250,5);
y = mean(azimuth_errors,2)';
est_params = fit_elbow_func_localization(x,y);
actual_elbow_point = est_params(3);
x_test = linspace(250, 2250, 1000);
y_pred = elbow_func(x_test, est_params);

% initialize figure
x_scale = 2.5 / 2.75; % adjust for size difference in x dimension since fewer foreground onset times
figure('units','inches','innerposition',[0 0 9*x_scale 11])
hold on
timevec = linspace(250,2250,5);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), '--', 'LineWidth', 3, 'Color', clrs(1,:)); % style of elbow function line
plot(nan(1), nan(1), 'o-', 'LineWidth', 3, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of elbow point location

% plot shaded errorbars (SEM)
m = mean(azimuth_errors,2);
s = std(azimuth_errors,0,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch = [ub', flip(lb')];
ptch = patch(x_patch, y_patch, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot human results
scatter(timevec, mean(azimuth_errors,2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)

% plot elbow function fit to data
p1 = plot(x_test, y_pred, '--', 'LineWidth', 3, 'Color', clrs(1,:));

% plot elbow point location
med = median(elbow_points);
stdev_ep = std(elbow_points);
ypos = 36;
plot([med - stdev_ep, med + stdev_ep], [ypos, ypos], 'LineWidth', 3, 'Color', clrs(1,:))
scatter(actual_elbow_point, ypos, 150, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)

% plot dividing line
plot([0 2500], [35, 35], 'k', 'LineWidth', 2)

% add number of subjects
text(75, 6, sprintf('n = %d', n_subs), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',ones(3,1)*0)

% stylize figure and add labels
box off
ylabel(sprintf('Localization Error (%s)', char(176)), 'Position', [-300*x_scale, 20])
xlabel('Foreground Onset Time (ms)')
ylim([5, 37.5])
yticks(5:5:35)
xlim([0 2500])
xticks((timevec))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)
set(gca,'YDir','reverse')

% improve look of legend
[l, objh]  = legend({'    Avg. performance', '    Elbow function fit', '    Elbow point location'},'Location','southeast');
objh(4).XData = [objh(4).XData(1), objh(4).XData(2)*1.18];
objh(4).Color = [objh(4).Color 0.5];
objh(5).XData = mean(objh(4).XData);
objh(5).LineWidth = 1;
objh(6).XData = [objh(6).XData(1), objh(6).XData(2)*1.5];
objh(8).XData = [objh(8).XData(1), objh(8).XData(2)*1.18];
objh(9).XData = mean(objh(8).XData);
objh(9).LineWidth = 1;
l.Position(1) = l.Position(1)*1.15;
l.Position(2) = l.Position(2)*1.4;
l.EdgeColor = 'none';

% save figure as svg
set(gca,'Position', [(0.145/x_scale), (0.165), (1 - 0.165/x_scale), (1 - 0.195)]) % remove space around edges
save_name = fullfile(save_path, sprintf('%s_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');