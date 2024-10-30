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
main_color = [0, 129, 105]; % "deep sea" see http://mkweb.bcgsc.ca/colorblind/palettes.mhtml (Martin Krzywinski — Designing for Color Blindness) 
alt_color = main_color + (255 - main_color) * 0.4;
clrs = [main_color; alt_color] / 255;

% settings for waveforms...note must match saved waveforms in wav_path
nstims = 5;
stims = [123 99 32 40 22]; % stim IDs (bg/fg pairing)
outliers = [1 2 2 1 2]; % presence of outlier
times = [3 1 9 5 7]; % foreground onset times
snrs = [2 1 2 2 1]; % SNRs
trials = [1 2 3 319 320]; % trial numbers

% waveform plot settings
clip_thresh = 0.25; % value at which to clip waveforms to make them look a bit nicer
y_positions = flip((1:(nstims+1)))*clip_thresh*2; % vector of y-positions of waveforms
skip_position = y_positions(end-2); % position that is skipped and replaced with ellipsis
y_positions(end-2) = []; % actually remove the skipped position

% stimulus settings (taken from stimulus creation script)
seglen = 500; % of the foreground (in MS)
bg_length = 3.25; % length of background in seconds
audio_sr = 20000; % sample rate of all the audio
SNRvec = linspace(-5,-8,2); % set of possible SNRs
onset_times = linspace(0.25,2.50,10); % set of possible positions of the outlier

%% plot diagram

% initialize diagram figure
x_scale = 92/77;
figure('units','inches','innerposition',[0 0 7*x_scale 8])
ylim([0.25, 4])
xlim([-0.5*audio_sr, bg_length*audio_sr + 0.1*audio_sr + 0.75*audio_sr])
axis off
hold on
% loop over trials
for n = 1:nstims
    % load background waveform    
    bacgkround_file = sprintf('background_%d.wav', stims(n));
    background = audioread(fullfile(wav_path, bacgkround_file));
    background(abs(background) > clip_thresh) = clip_thresh;

    % load foreground waveform
    foreground_file = sprintf('foreground_%d.wav', stims(n));
    foreground = audioread(fullfile(wav_path, foreground_file));

    % get proper onset_time
    t = times(n);
    ind1 = onset_times(t)*audio_sr; % convert the ms timing to a sample

    % set SNR of foreground relative to background
    snr_ind = snrs(n);
    clip1 = background(ind1:ind1+(audio_sr*seglen/1000)-1); % extract corresponding segment of background            
    setsnr = SNRvec(snr_ind); % get SNR level
    c = (rms(clip1)/rms(foreground))*10^(setsnr/20); % compute the scaling factor to get the proper SNR
    outlier_switch = outliers(n);
    if outlier_switch == 1 % then no outlier
        c = 0; % will kill the outlier
    else % otherwise keep the right SNR level
    end

    foreground_cue = foreground;
    c_cue = (rms(clip1)/rms(foreground))*10^(SNRvec(1)/20); % compute the scaling factor to get the proper SNR
    foreground_cue = c_cue*foreground;
    foreground_cue(abs(foreground_cue) > clip_thresh) = clip_thresh; % apply clip threshold

    foreground = c*foreground; % actually set proper SNR
    foreground(abs(foreground) > clip_thresh) = clip_thresh; % apply clip threshold

    % plot stimulus
    temp_timevec = 1:length(background);
    bg_timevec = temp_timevec + 0.75*audio_sr;
    plot(bg_timevec, background + y_positions(n), 'LineWidth', 1.5, 'Color', 0.6*ones(3,1))
    plot(1:length(foreground), foreground_cue + y_positions(n), 'LineWidth', 1.5, 'Color', 'k')
    if outlier_switch==2 % then also plot the foreground
        temp_timevec = 1:length(background);
        fg_timevec = temp_timevec(ind1:ind1+(audio_sr*seglen/1000)-1) + 0.75*audio_sr;
        plot(fg_timevec, foreground + y_positions(n), 'LineWidth', 1.5, 'Color', 'k')
    end
    x_offset = 0;
    if n <= 3; x_offset = 0.025*audio_sr; end
    % add trial number label
    text(-0.35*audio_sr + x_offset, y_positions(n), sprintf('%d', trials(n)), 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')

end

% add more text to figure
text(-0.35*audio_sr + 0.025*audio_sr, skip_position, '. . .', 'HorizontalAlignment','center', 'FontSize', 24, 'Rotation', 90, 'VerticalAlignment','baseline')
text((length(background) + 0.75*audio_sr)/2, skip_position, '. . .', 'HorizontalAlignment','center', 'FontSize', 24, 'Rotation', 90, 'VerticalAlignment','baseline')
text(-0.35*audio_sr + 0.025*audio_sr, 3.4, 'trial', 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')
text(0.5*audio_sr/2, 3.4, 'cue', 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')
text((length(background))/2 + 0.75*audio_sr, 3.4, 'stimulus', 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')
text(mean([-0.5*audio_sr, bg_length*audio_sr + 0.1*audio_sr + 0.75*audio_sr]), 3.75, sprintf('Did you hear the cued sound?'), ...
    'HorizontalAlignment','center', 'FontSize', 24, 'VerticalAlignment','middle', 'FontWeight','bold', 'FontAngle','italic', 'Color', clrs(1,:))

% save diagram figure
set(gca,'LooseInset',get(gca,'TightInset')); % remove space around edges
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
n_subs = size(d_primes,2);

elbow_points = est_params(:, 3);

% fit elbow function
x = linspace(250,2500,10);
y = mean(d_primes,2)';
est_params1 = fit_elbow_func(x,y(1:10));
est_params2 = fit_elbow_func(x,y(11:20));
x_test = linspace(250, 2500, 1000);
y_pred1 = elbow_func(x_test, est_params1);
y_pred2 = elbow_func(x_test, est_params2);

d_primes_avg = (d_primes(1:10,:) + d_primes(11:20,:)) / 2;
y = mean(d_primes_avg,2)';
est_params_avg = fit_elbow_func(x,y);
actual_elbow_point = est_params_avg(3);

% initialize figure
figure('units','inches','innerposition',[0 0 9 11])
hold on
timevec = linspace(250,2500,10);
SNRs = SNRvec;

% plot dummy variables for legend
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), '--', 'LineWidth', 3, 'Color', clrs(1,:)); % style of elbow function line
plot(nan(1), nan(1), 'o-', 'LineWidth', 3, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of elbow point location

% plot shaded errorbars (SEM)
m = reshape(mean(d_primes,2),10,2);
s = reshape(std(d_primes,0,2),10,2)/sqrt(n_subs);
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
scatter(timevec, mean(d_primes(1:10,:),2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes(11:20,:),2), 125, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% add SNR labels
text(timevec(end), max(ub(9:10,1))+.075, sprintf('%d dB', SNRs(1)), ...
     'FontSize', 20, ...
     'FontAngle','italic', ...
     'Color',clrs(1,:), ...
     'HorizontalAlignment','right')
text(timevec(end), max(ub(9:10,2))+0.075, sprintf('%d dB', SNRs(2)), ...
     'FontSize', 20, ...
     'FontAngle','italic', ...
     'Color',clrs(2,:), ...
     'HorizontalAlignment','right')

% plot elbow point location
med = median(elbow_points);
stdev_ep = std(elbow_points);
ypos = -0.11;
plot([med - stdev_ep, med + stdev_ep], [ypos, ypos], 'LineWidth', 3, 'Color', clrs(1,:))
scatter(actual_elbow_point, ypos, 150, 'MarkerEdgeColor','w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)

% plot dividing line
plot([0 2750], [0, 0], 'k', 'LineWidth', 2)

% add number of subjects
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 25, ...
     'FontAngle','italic', ...
     'Color',ones(3,1)*0)

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
l.Position(2) = l.Position(2)*1.5;
l.EdgeColor = 'none';

% save figure as svg
set(gca,'Position', [0.145, 0.145, (1 - 0.165), (1 - 0.195)]) % remove space around edges
save_name = fullfile(save_path, sprintf('%s_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');