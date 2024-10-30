%% initialization and settings

[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

matched_permutation_test_file = fullfile(file_path, sprintf('../../results/%s/matched_permutation_test.mat', exp_name));
load(matched_permutation_test_file)
indices = [6, 8, 10]; % corresponds to timepoints 1500, 2000, 2500 ms
test_indices = setdiff(1:10, indices);

% set colors
main_color = [247, 72, 165]; 
alt_color = [34, 113, 178];
clrs = [main_color; alt_color] / 255;

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end
d_primes = d_primes(:, logical(keep_subs));
d_primes = (d_primes(1:10, :) + d_primes(11:20, :)) / 2;

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

d_primes_matched = d_primes(:,matched_subs);
d_primes_matched_comparison = d_primes_comparison(:, matched_subs_comparison);
n_subs = size(d_primes_matched,2);

%% plot results figure
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
y_scale = 3 / 3.25; % adjust for size difference in y dimension since not plotting elbow point CI
figure('units','inches','innerposition',[0 0 9 11*y_scale])
hold on
timevec = linspace(250,2500,10);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

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

% add subject numbers
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 25, ...
     'FontAngle','italic', ...
     'Color','k')

% stylize figure and add labels
box off
ylabel('Detection Performance (d'')', 'Position', [-300, 1.5])
xlabel('Foreground Onset Time (ms)')
ylim([0 3])
yticks(0:0.5:3)
ytickformat('%.1f')
xlim([0 2750])
xticks((timevec))
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',25)
set(gca,'XTickLabelRotation',45)

% improve look of legend
[l, objh]  = legend({'    Stationary noise', '    Natural textures'},'Location','southeast');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.Position(1) = l.Position(1)*1.075;
l.Position(2) = l.Position(2)*1.2;
l.EdgeColor = 'none';
text(2010, 0.5, 'Background Type', ...
     'FontSize', 22.5, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', ...
     'Color',ones(3,1)*0)

% save figure as svg
set(gca,'Position', [(0.145), (0.145/y_scale), (1 - 0.165), (1 - 0.195/y_scale)]) % remove space around edges
save_name = fullfile(save_path, sprintf('%s_supplement_results.svg', exp_name));
saveas(gcf, save_name)