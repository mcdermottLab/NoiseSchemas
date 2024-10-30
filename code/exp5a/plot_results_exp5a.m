%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
wav_path = fullfile(file_path, sprintf('../../figures/%s/wavs/', exp_name)); % get path to raw waveforms
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% set random number generator
rng(1)

% set colors
clrs = [219, 133, 55; 239, 206, 89]/255;

% settings for waveforms...note must match saved waveforms in wav_path
nstims = 5;
stims = [123 99 32 40 22]; % stim IDs (bg/fg pairing)
outliers = [1 2 2 1 2]; % presence of outlier
times = [3 1 8 5 4]; % foreground onset times
snrs = [1 1 1 1 1]; % SNRs
ints = [2 1 2 1 1]; % interrupter types
trials = [1 2 3 319 320]; % trial numbers

% waveform plot settings
wav_scale = 2; % scales amplitude of everything for nicer plotting
clip_thresh = 0.25; % value at which to clip waveforms to make them look a bit nicer
y_positions = flip((1:(nstims+1)))*clip_thresh*2; % vector of y-positions of waveforms
skip_position = y_positions(end-2); % position that is skipped and replaced with ellipsis
y_positions(end-2) = []; % actually remove the skipped position

% stimulus settings (taken from stimulus creation script)
audio_sr = 20000; % sample rate of all the audio
seglen = 500; % of the foreground (in MS)
bg_length = 4; % length of background in seconds
int_length = 0.5; % length of noise in (s)
int_start = ((bg_length - int_length)/2)*audio_sr + 1; % starting position of interrupter
int_end = int_start + int_length*audio_sr - 1; % ending position of interrupter
SNRvec = -2; % set of possible SNRs
noise_dB = 12; % level of noise relative to background
desired_rms = .05 / (10^(noise_dB/20)); % sets noise to 0.05 rms and rest of texture lower
onset_times = [linspace(0.25,1,4), linspace(2.5,3.25,4)]; % set of possible positions of the outlier

%% plot diagram

% initialize diagram figure
% x_scale = 92/77;
figure('units','inches','innerposition',[0 0 6 8])
ylim([0.25, 4])
xlim([-0.5*audio_sr, bg_length*audio_sr + 0.1*audio_sr])
axis off
hold on
% loop over trials
for n = 1:nstims
    % load background waveform    
    bacgkround_file = sprintf('background_%d.wav', stims(n));
    background = audioread(fullfile(wav_path, bacgkround_file));
    background = background * wav_scale;
    background(abs(background) > clip_thresh) = clip_thresh;

    % load foreground waveform
    foreground_file = sprintf('foreground_%d.wav', stims(n));
    foreground = audioread(fullfile(wav_path, foreground_file));
    
    % create interrupter
    interrupter = randn(int_length*audio_sr, 1);
    interrupter = desired_rms * interrupter / rms(interrupter); % normalize rms level to desired_rms
    if ints(n)==1
        int_scale = 0;
    else
        int_scale = 10^(noise_dB/20);
    end
    interrupter = int_scale * interrupter;
    interrupter = interrupter * wav_scale;

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
    
    foreground = c*foreground; % actually set proper SNR
    foreground(abs(foreground) > clip_thresh) = clip_thresh; % apply clip threshold
    
    % plot stimulus
    plot(background + y_positions(n), 'LineWidth', 1.5, 'Color', 0.6*ones(3,1))    
    if outlier_switch==2 % then also plot the foreground
        temp_timevec = 1:length(background);
        fg_timevec = temp_timevec(ind1:ind1+(audio_sr*seglen/1000)-1);
        plot(fg_timevec, foreground + y_positions(n), 'LineWidth', 1.5, 'Color', 'k')
    end
    if int_scale==0; lw = 2; else; lw = 1.5; end
    plot(int_start:int_end, interrupter + y_positions(n), 'LineWidth', lw, 'Color', clrs(ints(n),:))
    x_offset = 0;
    if n <= 3; x_offset = 0.025*audio_sr; end
    % add trial number label
    text(-0.35*audio_sr + x_offset, y_positions(n), sprintf('%d', trials(n)), 'HorizontalAlignment','center', 'FontSize', 24, 'VerticalAlignment','middle')

end

% add more text to figure
label_offset = 0.45; % offset for labels to account for interrupting noise height
text(-0.35*audio_sr + 0.025*audio_sr, skip_position-0.1, '. . .', 'HorizontalAlignment','center', 'FontSize', 24, 'Rotation', 90, 'VerticalAlignment','baseline')
text(length(background)/2, skip_position-0.1, '. . .', 'HorizontalAlignment','center', 'FontSize', 24, 'Rotation', 90, 'VerticalAlignment','baseline')
text(-0.35*audio_sr + 0.025*audio_sr, y_positions(1) + label_offset, 'trial', 'HorizontalAlignment','center', 'FontSize', 24, 'VerticalAlignment','middle')
text(length(background)/2, y_positions(1) + label_offset, 'stimulus', 'HorizontalAlignment','center', 'FontSize', 24, 'VerticalAlignment','middle')

% save diagram figure
set(gca,'LooseInset',get(gca,'TightInset')); % remove space around edges
ax = gca;
ax.Position(4) = 1;
save_name = fullfile(save_path, sprintf('%s_diagram.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');

%% plot results

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

% initialize figure
x_scale = 3.5 / 2.75;
y_scale = 3 / 3.25;
figure('units','inches','innerposition',[0 0 6.3*x_scale 11*y_scale])
hold on
timevec = [linspace(250,1000,4), linspace(2500,3250,4)]; % scale so plot matches previous figure scale
plotting_timevec = [timevec, timevec];
interrupter_times = [1750, 2250];

% plot dummy stuff for legend styling
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(1,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:)); % style of performance data
plot(nan(1), nan(1), 'o-', 'LineWidth', 15, 'Color', clrs(2,:), 'MarkerSize', 12.5, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:)); % style of performance data

% plot shaded error bars & average performance
inds = [1:4, 9:12];
m = reshape(mean(d_primes(inds,:),2),4,2);
s = reshape(std(d_primes(inds,:),0,2),4,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(1:4), flip(timevec(1:4))];
y_patch1 = [ub(:,1)', flip(lb(:,1)')];
y_patch2 = [ub(:,2)', flip(lb(:,2)')];
ptch = patch(x_patch, y_patch1, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';
ptch = patch(x_patch, y_patch2, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

inds = [5:8, 13:16];
m = reshape(mean(d_primes(inds,:),2),4,2);
s = reshape(std(d_primes(inds,:),0,2),4,2)/sqrt(n_subs);
ub = m + s;
lb = m - s;
x_patch = [timevec(5:8), flip(timevec(5:8))];
y_patch1 = [ub(:,1)', flip(lb(:,1)')];
y_patch2 = [ub(:,2)', flip(lb(:,2)')];
ptch = patch(x_patch, y_patch1, clrs(1,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';
ptch = patch(x_patch, y_patch2, clrs(2,:));
ptch.FaceAlpha = 0.5;
ptch.EdgeColor = 'none';

% plot human results
scatter(timevec, mean(d_primes(1:8,:),2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(1,:), 'LineWidth', 0.75)
scatter(timevec, mean(d_primes(9:16,:),2), 125, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', clrs(2,:), 'LineWidth', 0.75)

% plot interrupter as shaded region
patch([interrupter_times(1), interrupter_times(1), interrupter_times(2), interrupter_times(2)], [0, 3, 3, 0], ones(1,3)*0.85, ...
    'FaceColor', ones(1,3)*0.85, 'FaceAlpha', 0.25, 'EdgeColor', 'None')

% text(mean(interrupter_times), 1.5, '500 ms interrupter', ...
%      'FontSize', 30, ...
%      'Rotation', 90, ...
%      'Color',ones(3,1)*0.6, ...
%      'HorizontalAlignment','center')

% add subject numbers
text(75, 2.9, sprintf('n = %d', n_subs), ...
     'FontSize', 30, ...
     'FontAngle','italic', ...
     'Color',ones(3,1)*0)

% stylize figure and add labels
box off
ylim([0 3])
yticks(0:0.5:3)
ytickformat('%.1f')
ylabel('Detection Performance (d'')', 'Position', [-450, 1.5])
xlim([0 3500])
xticks(timevec)
xlabel('Foreground Onset Time (ms)')
set(gca,'linewidth',2)
set(gca, 'TickDir', 'both')
set(gca,'FontSize',30)
set(gca,'XTickLabelRotation',45)

% add legend
[l, objh]  = legend({'  Silence', '  White Noise'},'Location','southwest');
objh(3).XData = [objh(3).XData(1), objh(3).XData(2)*1.18];
objh(3).Color = [objh(3).Color 0.5];
objh(4).XData = mean(objh(3).XData);
objh(4).LineWidth = 1;
objh(5).XData = [objh(5).XData(1), objh(5).XData(2)*1.18];
objh(5).Color = [objh(5).Color 0.5];
objh(6).XData = mean(objh(5).XData);
objh(6).LineWidth = 1;
l.Position(1) = l.Position(1)*1.35;
l.Position(2) = l.Position(2)*1.05;
l.EdgeColor = 'none';

% legend title
text(100, 0.5, 'Interrupter Type', ...
     'FontSize', 28, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', ...
     'Color',ones(3,1)*0)

% save figure as svg
set(gca,'Position', [(0.22/x_scale), (0.165/y_scale), (1 - 0.23/x_scale), (1 - 0.195/y_scale)]) % remove space around edges
save_name = fullfile(save_path, sprintf('%s_results.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');