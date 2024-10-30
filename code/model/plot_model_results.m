%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
model_name = '500ms_present_1000ms_past';
results_file = fullfile(file_path, sprintf('../../results/model/%s/model_results.mat', model_name)); % get path to results
load(results_file)
save_path = fullfile(file_path,'../../figures/model/'); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

% paths = [pathsep, path, pathsep];
% utils_path = fullfile(file_path, '../utils/');
% on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
% if ~ on_path; addpath(utils_path); end

% set colors
main_color = [0, 0, 0];
alt_color = main_color + (255 - main_color) * 0.5;
clrs = [main_color; alt_color] / 255;

%% plot results figure

% initialize figure
y_scale = 2 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 9*1.1 11*y_scale*1.1])
hold on
timevec = linspace(250,2500,10);
SNRs = [-2, -6];

% plot shaded errorbars (stdev)
m = reshape(mean_DPs,10,2);
s = reshape(stdev_DPs,10,2);
ub = m + s;
lb = m - s;
x_patch = [timevec, flip(timevec)];
y_patch1 = [ub(:,1)', flip(lb(:,1)')];
y_patch2 = [ub(:,2)', flip(lb(:,2)')];
ptch = patch(x_patch, y_patch1, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';
ptch = patch(x_patch, y_patch2, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot model results
scatter(timevec, mean_DPs(1:10), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean_DPs(11:20), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% add SNR labels
text(timevec(end), max(ub(9:10,1))+.075, sprintf('%d dB', SNRs(1)), ...
     'FontSize', 24, ...
     'FontAngle','italic', ...
     'Color',clrs(1,:), ...
     'HorizontalAlignment','right')
text(timevec(end), max(ub(9:10,2))+0.075, sprintf('%d dB', SNRs(2)), ...
     'FontSize', 24, ...
     'FontAngle','italic', ...
     'Color',clrs(2,:), ...
     'HorizontalAlignment','right')

% stylize figure and add labels
box off
ylabel('Detection Performance (d'')', 'Position', [-300, 1])
xlabel('Foreground Onset Time (ms)')
ylim([0 2])
yticks(0:0.5:3)
ytickformat('%.1f')
xlim([0 2750])
xticks((timevec))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% save figure as svg
set(gca,'Position', [(0.145), (0.165/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges
save_name = fullfile(save_path, 'model_results.pdf');
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% plot scatter plot with human results

comparison_exp_name = 'exp1';
comparison_results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp_name)); % get path to results
load(comparison_results_file)

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

mean_human_d_primes = mean(d_primes,2);

y_scale = 2 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 7 11*y_scale*1.1])
axis square
hold on
plot(mean_human_d_primes(1:10), mean_DPs(1:10), 'o', 'MarkerSize', 13, 'MarkerEdgeColor', ones(3,1), 'MarkerFaceColor', clrs(1,:))
plot(mean_human_d_primes(11:20), mean_DPs(11:20), 'o', 'MarkerSize', 13, 'MarkerEdgeColor', ones(3,1), 'MarkerFaceColor', clrs(2,:))

%stylize
box off
ylabel('Model Performance (d'')')
xlabel('Human Performance (d'')')
ylim([0 2])
yticks(0:0.5:2)
ytickformat('%.1f')
xlim([0 2])
xticks(0:0.5:2)
xtickformat('%.1f')
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'Position', [(0.145), (0.11/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges

% save figure as svg
save_name = fullfile(save_path, 'model_human_scatterplot.pdf');
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');