%% initialization and settings

[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/spectrotemporal_split_timing_analysis.mat', exp_name)); % get path to timing analysis
load(timing_analysis_file)
elbow_points_high = est_params_high(:,3);
elbow_points_low = est_params_low(:,3);
pairings_file = fullfile(file_path, sprintf('../../results/%s/spectrotemporal_split_pairings.mat', exp_name)); % get path to timing analysis
load(pairings_file)
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% set colors
main_color = ones(1,3)*0.2; 
alt_color = ones(1,3)*0.7;
clrs = [main_color; alt_color];

%% plot figure

low_pairs(:,2) = goodinds(low_pairs(:,2));
low_pairs(:, 3:4) = [d_spec_low, d_spectemp_low];
high_pairs(:,2) = goodinds(high_pairs(:,2));
high_pairs(:, 3:4) = [d_spec_high, d_spectemp_high];

all_groups = [];
for i = 1:3
    for j = 1:10
        all_groups(end+1, :) = [i j];
    end
end

d_primes = nan(30, length(results));
distances = nan(4, length(results));
remove_subs = [];
for i = 1:length(results)
    if results(i).pass
        responses = results(i).responses;
        on_trials = find(responses(:,5)==2);
        low_trials = on_trials(ismember(responses(on_trials,1:2), low_pairs(:,1:2), 'rows'));
        responses(low_trials,3) = 3;
        high_trials = on_trials(ismember(responses(on_trials,1:2), high_pairs(:,1:2), 'rows'));
        responses(high_trials,3) = 2;
        responses_to_use = responses(:, 2:end);
        responses_to_use(sum(responses_to_use,2)==0,:) = []; % remove trials with no data                           
        [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
        group_info = [ID1 ID2];
        hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
        mean_false_alarm_rate = compute_false_alarm_rate(responses_to_use(:, 4), responses_to_use(:, 5)); % compute false alarms over all "conditions" since they only apply with foregrounds present
        inds = find(ismember(all_groups, [ID1 ID2], 'rows'));
        d_primes(inds, i) = norminv(hit_rates) - norminv(mean_false_alarm_rate);
        
        low_inds = find(ismember(low_pairs(:, 1:2), responses(low_trials,1:2), 'rows'));
        high_inds = find(ismember(high_pairs(:, 1:2), responses(high_trials,1:2), 'rows'));
        distances(1, i) = mean(low_pairs(low_inds, 3));        
        distances(2, i) = mean(low_pairs(low_inds, 4));
        distances(3, i) = mean(high_pairs(high_inds, 3));
        distances(4, i) = mean(high_pairs(high_inds, 4));
    else
        remove_subs(end+1) = i;
    end
end

d_primes = d_primes(11:30, :);
d_primes(:, remove_subs) = [];
n_subs = size(d_primes,2);

% initialize figure
figure('units','inches','innerposition',[0 0 9 11])
hold on
timevec = linspace(250,2500,10);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% fit elbow function
x = linspace(250,2500,10);
y = mean(d_primes,2,'omitnan')';
est_params1 = fit_elbow_func(x,y(1:10));
actual_elbow_point_high = est_params1(3);
est_params2 = fit_elbow_func(x,y(11:20));
actual_elbow_point_low = est_params2(3);
x_test = linspace(250, 2500, 1000);
y_pred1 = elbow_func(x_test, est_params1);
y_pred2 = elbow_func(x_test, est_params2);

% plot shaded errorbars (SEM)
m = reshape(mean(d_primes,2,'omitnan'),10,2);
s = reshape(std(d_primes,0,2,'omitnan'),10,2)/sqrt(n_subs);
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

% plot elbow function fit to data
p1 = plot(x_test, y_pred1, '--', 'LineWidth', 3, 'Color', clrs(1,:));
p2 = plot(x_test, y_pred2, '--', 'LineWidth', 3, 'Color', clrs(2,:));

% plot human results
scatter(timevec, mean(d_primes(1:10,:),2,'omitnan'), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes(11:20,:),2,'omitnan'), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% add number of subjects
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 25, ...
     'FontAngle','italic', ...
     'Color',ones(3,1)*0)

% plot elbow point locations
med = median(elbow_points_high);
stdev_ep = std(elbow_points_high);
ypos = -0.08;
plot([med - stdev_ep, med + stdev_ep], [ypos, ypos], 'LineWidth', 3, 'Color', clrs(1,:))
scatter(actual_elbow_point_high, ypos, 150, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)

med = median(elbow_points_low);
stdev_ep = std(elbow_points_low);
ypos = -0.14;
plot([med - stdev_ep, med + stdev_ep], [ypos, ypos], 'LineWidth', 3, 'Color', clrs(2,:))
scatter(actual_elbow_point_low, ypos, 150, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot dividing line
plot([0 2750], [0, 0], 'k', 'LineWidth', 2)

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
set(gca,'FontSize',25)
set(gca,'XTickLabelRotation',45)

% improve look of legend
[l, objh]  = legend({'    Large Spectrotemporal Distance', '    Small Spectrotemporal Distance'},'Location','southeast');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.Position(1) = l.Position(1)*1.1;
l.Position(2) = l.Position(2)*1.5;
l.EdgeColor = 'none';
% text(1730, 0.41, 'Pairing Type', ...
%      'FontSize', 22.5, ...
%      'FontWeight', 'bold', ...
%      'HorizontalAlignment', 'center', ...
%      'Color',ones(3,1)*0)
text(575, 0.41, sprintf('Foreground-Background Pairing Type'), ...
     'FontSize', 22.5, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', ...
     'Color',ones(3,1)*0)
set(gca,'Position', [0.145, 0.145, (1 - 0.165), (1 - 0.195)]) % remove space around edges

% initialize figure
figure('units','inches','innerposition',[0 0 9 11])
hold on

distances(:, remove_subs) = [];
m = mean(distances,2);
s = std(distances, 0, 2);

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 5, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 5, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% errorbar(1:2, m(1:2), s(1:2), '.-', 'LineWidth', 5, 'CapSize', 0, 'MarkerSize', 30, 'Color', clrs(2,:))
% errorbar((1:2)+0.02, m(3:4), s(3:4), '.-', 'LineWidth', 5, 'CapSize', 0, 'MarkerSize', 30, 'Color', clrs(1,:))

errorbar(1:2, m(1:2), s(1:2), '-', 'LineWidth', 5, 'CapSize', 0, 'Color', clrs(2,:))
scatter(1:2, m(1:2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)
errorbar((1:2)+0.02, m(3:4), s(3:4), '-', 'LineWidth', 5, 'CapSize', 0, 'Color', clrs(1,:))
scatter((1:2)+0.02, m(3:4), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)

xlim([0.8, 2.2])
% ylim([0.25,0.55])
ylim([0.20,0.55])
yticks(linspace(0.2,0.55,8))

ylabel('Normalized Distance')
xlabel('Distance Type')
xticks(1:2)
xticklabels({'Spectral', 'Spectrotemporal'})
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',25)

% improve look of legend
[l, objh]  = legend({'    Large Spectrotemporal Distance', '    Small Spectrotemporal Distance'},'Location','southwest');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.Position(1) = l.Position(1)*1.15;
l.Position(2) = l.Position(2)*1.15;
l.EdgeColor = 'none';

text(0.875, 0.25, sprintf('Foreground-Background Pairing Type'), ...
     'FontSize', 22.5, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', ...
     'Color',ones(3,1)*0)
